#include <chrono>
#include <stdexcept>

#include "scidf.h"
#include "spade.h"
#include "proto/hywall_interop.h"

#include "typedef.h"
#include "calc_u_bulk.h"
#include "misc.h"
#include "bc.h"


int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    const std::size_t dim = 3;

    std::string input_filename = "none";
    for (auto i: range(0, argc))
    {
        std::string line(argv[i]);
        if (ends_with(line, ".sdf"))
        {
            input_filename = line;
            if (group.isroot()) print("Reading", input_filename);
        }
    }
    if (input_filename == "none")
    {
        if (group.isroot()) print("E: No input file name provided!");
        return 1;
    }
    
    scidf::node_t input;
    scidf::clargs_t clargs(argc, argv);
    try { scidf::read(input_filename, input, clargs); }
    catch (const std::exception& e)
    {
        if (group.isroot()) print(std::string("INPUT FILE ERROR:\n") + e.what());
        group.sync();
        return 143;
    }

    spade::ctrs::array<int, dim> num_blocks      = input["Grid"]["num_blocks"];
    spade::ctrs::array<int, dim> cells_in_block  = input["Grid"]["num_cells"];
    spade::ctrs::array<int, dim> exchange_cells  = input["Grid"]["num_exchg"];
    spade::bound_box_t<real_t, dim> bounds       = input["Grid"]["dims"];

	int        nt_backup        = 10;
    real_t     targ_cfl         = input["Time"]["cfl"];
    int        nt_max           = input["Time"]["nt_max"];
    int        nt_skip          = input["Time"]["nt_skip"];
    int        checkpoint_skip  = input["Time"]["ck_skip"];
    bool       output_timing    = input["Time"]["output_timing"];
    
    real_t   Twall         = input["WallModel"]["wallTemperature"];
    real_t   prandtl       = input["WallModel"]["fluidPrandtl"];
    real_t   mu_wall       = input["Fluid"]["mu_wall"];
    real_t   tau_wall      = input["Fluid"]["tau_wall"];
    real_t   rho_b         = input["Fluid"]["rho_b"];

    enum perturb_type
    {
        perturb_none = 0,
        perturb_freq = 1,
        perturb_rand = 2
    } ptype = (perturb_type)(scidf::menu_t<std::string>(input["Fluid"]["perturb"]).selected_index());
    
    bool     wm_enable     = input["Fluid"]["wm_enable"];
    int      couple_dist   = input["Fluid"]["couple_dist"];
    bool     const_init    = input["Fluid"]["const_init"];
    real_t   rgas          = input["WallModel"]["gasConstant"];
    real_t   fluidCp       = input["WallModel"]["fluidCp"];
    
    
    real_t             eps_ducr  = input["Num"]["eps_ducr"];
    real_t                eps_p  = input["Num"]["eps_p"];
    real_t                eps_T  = input["Num"]["eps_T"];
	real_t                cw     = input["Num"]["wale_cw"];
    bool               smooth_wm = input["Num"]["smooth_wm"];

    std::string        init_file = input["IO"]["init_file"];
    std::string        out_dir   = input["IO"]["out_dir"];
	std::string      hist_file   = input["IO"]["hist_file"];
	std::string      visc_file   = input["IO"]["visc_file"];
    bool             output_ducr = input["IO"]["output_ducr"];
    bool             output_sgs  = input["IO"]["output_sgs"];
    
    
    spade::coords::identity<real_t> coords;
    
    std::filesystem::path primary(out_dir);
    if (!std::filesystem::is_directory(primary)) std::filesystem::create_directory(primary);
    
    std::filesystem::path ck_dir(primary / "checkpoint");
    if (!std::filesystem::is_directory(ck_dir)) std::filesystem::create_directory(ck_dir);

    const std::string data_out = primary / "viz";
    
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    real_t delta = 0.5*(bounds.size(1));
    
    prim_t fill1 = 0.0;
    flux_t fill2 = 0.0;
    real_t fill3 = 0.0;
    
    spade::grid::grid_array prim       (grid, fill1);
	spade::grid::grid_array prim_crash (grid, fill1);
    spade::grid::grid_array rhs        (grid, fill2);
    spade::grid::grid_array scalar     (grid, fill3);
    
    spade::fluid_state::ideal_gas_t<real_t> air;
    air.R = rgas;
    air.gamma = fluidCp/(fluidCp-rgas);

    const real_t delta_g = grid.get_dx(1);
    spade::viscous_laws::power_law_t lam_visc_law(mu_wall, Twall, 0.76, prandtl);
    spade::subgrid_scale::wale_t eddy_visc   (air, cw, delta_g, 0.9);
    spade::viscous_laws::sgs_visc_t visc_law(lam_visc_law, eddy_visc);
    
    const real_t Lx  = bounds.size(0);
    const real_t Lz  = bounds.size(2);

    real_t Tref = 7.0*Twall;
    real_t p0 = rho_b*air.get_R()*Tref;
    real_t aw = std::sqrt(air.get_gamma()*air.get_R()*Twall);
    real_t u0 = 800*std::sqrt(tau_wall*rho_b);

    const int nmode = 11;
    using point_type = decltype(grid)::coord_point_type;
    auto ini = [&](const point_type& x) -> prim_t
    {
        real_t yh = x[1]/delta;
        prim_t output;
        output.p() = p0;
        output.T() = Tref - (Tref - Twall)*yh*yh;
        output.u() = (3.0/2.0)*u0*(1.0-yh*yh);
        if (const_init)
        {
            output.T() = Tref;
            output.u() = u0;
        }
        output.v() = 0;
        output.w() = 0;
        
        return output;
    };

    auto freq_perturb_func = [&](const prim_t& val, const point_type& x)
    {
        real_t up = 0.0;
        real_t vp = 0.0;
        real_t wp = 0.0;
        real_t yh = x[1]/delta;
        int imin = 1;
        for (int ii = imin; ii < imin + nmode; ++ii)
        {
            real_t ampl   = 0.1*u0*(1.0-yh*yh)/(ii*ii);
            real_t freq_x = 2.0*spade::consts::pi*ii/(0.5*bounds.max(0));
            real_t freq_y = 2.0*spade::consts::pi*ii/(0.5*delta);
            real_t freq_z = 2.0*spade::consts::pi*ii/(0.5*bounds.max(2));
            real_t phase_x = std::sin(14*ii)*2.0*spade::consts::pi;
            real_t phase_y = std::sin(10*ii)*2.0*spade::consts::pi;
            real_t phase_z = std::sin(17*ii)*2.0*spade::consts::pi;
            up += ampl*std::sin(freq_x*x[0]-phase_x)*std::sin(freq_y*x[1]+phase_x)*std::sin(freq_z*x[2]-phase_x);
            vp += ampl*std::sin(freq_x*x[0]+phase_y)*std::sin(freq_y*x[1]-phase_y)*std::sin(freq_z*x[2]+phase_y);
            wp += ampl*std::sin(freq_x*x[0]-phase_z)*std::sin(freq_y*x[1]+phase_z)*std::sin(freq_z*x[2]-phase_z);
        }
        prim_t output = val;
        output.u() += up;
        output.v() += vp;
        output.w() += wp;
        return output;
    };

    std::vector<int> glob(grid.get_partition().get_num_local_blocks());
    for (auto i: range(0, glob.size())) glob[i] = grid.get_partition().get_global_block(i);
    spade::utils::random_seed(1234);
    std::vector<real_t> u_perturb(grid.get_partition().get_num_global_blocks());
    std::vector<real_t> v_perturb(grid.get_partition().get_num_global_blocks());
    std::vector<real_t> w_perturb(grid.get_partition().get_num_global_blocks());

    const real_t perturb_ampl = 0.15*u0;
    std::generate(u_perturb.begin(), u_perturb.end(), [&](){return perturb_ampl*(1.0-2.0*spade::utils::unitary_random());});
    std::generate(v_perturb.begin(), v_perturb.end(), [&](){return perturb_ampl*(1.0-2.0*spade::utils::unitary_random());});
    std::generate(w_perturb.begin(), w_perturb.end(), [&](){return perturb_ampl*(1.0-2.0*spade::utils::unitary_random());});
    auto rand_perturb_func = [&](const prim_t& val, const point_type& x, const spade::grid::cell_idx_t& ii)
    {
        prim_t output = val;
        real_t yh = x[1]/delta;
        const real_t shape = (1.0-yh*yh*yh*yh);
        const auto lbglob = glob[ii.lb()];
        output.u() += u_perturb[lbglob]*shape;
        output.v() += v_perturb[lbglob]*shape;
        output.w() += w_perturb[lbglob]*shape;
        return output;
    };


    spade::algs::fill_array(prim, ini);

    if (ptype == perturb_freq) spade::algs::transform_inplace(prim, freq_perturb_func);
    if (ptype == perturb_rand) spade::algs::transform_inplace(prim, rand_perturb_func);

	
	real_t ub, rhob;
    calc_u_bulk(prim, air, ub, rhob);
    real_t ratio = rho_b/rhob;

    auto rhob_correct = [&](const prim_t& val)
    {
        cons_t w;
        prim_t q;
        spade::fluid_state::convert_state(val, w, air);
        w.rho() = ratio*w.rho();
        spade::fluid_state::convert_state(w, q, air);
        return q;
    };
    spade::algs::transform_inplace(prim, rhob_correct);
    calc_u_bulk(prim, air, ub, rhob);
    if (group.isroot())
    {
        print("Corrected rho_b with ratio", ratio);
        print("Specified: ", rho_b);
        print("Calculated:", rhob);
    }
    
    int nt_min = 0;
    if (init_file != "none")
    {
        if (group.isroot()) print("Loading", init_file+"...");
        if (!std::filesystem::exists(init_file))
        {
            if (group.isroot()) print("Invalid init file!");
            return 140;
        }
        spade::io::binary_read(init_file, prim);
        if (group.isroot()) print("Init done.");
        grid.exchange_array(prim);
        set_channel_slip(prim, Twall, wm_enable);

        nt_min = ext_nt(init_file, 8);
        if (group.isroot()) print("Detected nt_min =", nt_min);
    }
    
    spade::state_sensor::ducros_t ducr(eps_ducr);
    spade::convective::totani_lr        tscheme(air);
    spade::convective::pressure_diss_lr dscheme(air, ducr, eps_p, eps_T);
    spade::viscous::visc_lr             visc_scheme(visc_law, air);
    

    auto get_u = [&](const prim_t& val){return std::sqrt(air.gamma*air.R*val.T()) + std::sqrt(val.u()*val.u() + val.v()*val.v() + val.w()*val.w());};

    spade::reduce_ops::reduce_max<real_t> max_op;
    
    
    
    const real_t dx       = spade::utils::min(grid.get_dx(0), grid.get_dx(1), grid.get_dx(2));
    const real_t umax_ini = spade::algs::transform_reduce(prim, get_u, max_op);
    const real_t dt       = targ_cfl*dx/umax_ini;
    
    const real_t force_term = tau_wall/(delta*rho_b);
    auto source = [&](const prim_t& val)
    {
        cons_t w;
        spade::fluid_state::convert_state(val, w, air);
        flux_t output;
        output.continuity() = 0.0;
        output.energy()     = w.rho()*force_term*val.u();
        output.x_momentum() = w.rho()*force_term;
        output.y_momentum() = 0.0;
        output.z_momentum() = 0.0;
        return output;
    };
    
    spade::bound_box_t<bool, grid.dim()> boundary = true;
    boundary.min(1) = false;
    boundary.max(1) = false;
    
    spade::proto::hywall_binding_t wall_model(prim, rhs, air, couple_dist);
    wall_model.read(input["WallModel"]);
    for (auto& b: boundary) b = !b;
    wall_model.init(prim, boundary);
    for (auto& b: boundary) b = !b;
    wall_model.set_dt(dt);

    auto boundary_cond = [&](auto& q, const auto& t)
    {
        grid.exchange_array(q);
        set_channel_slip(q, Twall, wm_enable);
    };

    auto calc_rhs = [&](auto& rhsin, const auto& qin, const auto& tin) -> void
    {
        rhsin = 0.0;
        if (wm_enable)
        {
            auto policy = spade::pde_algs::block_flux_all;
            spade::pde_algs::flux_div(qin, rhsin, tscheme);
            spade::pde_algs::flux_div(qin, rhsin, policy, boundary, visc_scheme, dscheme);
            wall_model.sample(qin, lam_visc_law);
            wall_model.solve();
            if (smooth_wm) wall_model.smooth();
            wall_model.apply_flux(rhsin);
        }
        else
        {
            spade::pde_algs::flux_div(qin, rhsin, tscheme, visc_scheme, dscheme);
        }
        spade::pde_algs::source_term(qin, rhsin, source);
    };
    
    cons_t transform_state;
    spade::fluid_state::state_transform_t trans(transform_state, air);

    real_t time0 = 0.0;
    spade::time_integration::time_axis_t       axis    (time0, dt);
    spade::time_integration::ssprk3_t          alg;
    spade::time_integration::integrator_data_t q       (prim, rhs, alg);
    spade::time_integration::integrator_t      time_int(axis, alg, q, calc_rhs, boundary_cond, trans);

    boundary_cond(time_int.solution(), time0);

    spade::utils::avg_t<real_t> perf;
    
    spade::timing::mtimer_t tmr("advance");
    std::ofstream myfile(hist_file);
    std::ofstream myfile2(visc_file);
    for (auto nt: range(nt_min, nt_min+nt_max+1))
    {
        const real_t umax   = spade::algs::transform_reduce(time_int.solution(), get_u, max_op);
        calc_u_bulk(time_int.solution(), air, ub, rhob);
        
        if (group.isroot())
        {
            const real_t cfl = umax*dt/dx;
            const int pn = 10;
            print(
                "nt: ", spade::utils::pad_str(nt,   pn),
                "cfl:", spade::utils::pad_str(cfl,  pn),
                "u+a:", spade::utils::pad_str(umax, pn),
                "mb: ", spade::utils::pad_str(ub/aw,pn),
                "rb: ", spade::utils::pad_str(rhob, pn),
                "dx: ", spade::utils::pad_str(dx,   pn),
                "dt: ", spade::utils::pad_str(dt,   pn)
            );
            myfile << nt << " " << cfl << " " << umax << " " << ub << " " << rhob << " " << dx << " " << dt << std::endl;
            myfile.flush();
        }
        auto tau_avg = wall_model.mean_tau();
        auto qw_avg  = wall_model.mean_qw ();
        if (group.isroot())
        {
            print("Mean shear/heat:", tau_avg, qw_avg);
            myfile2 << time_int.time() << " " << tau_avg << " " << qw_avg << "\n";
            myfile2.flush();
        }
        if (nt%nt_skip == 0)
        {
            if (group.isroot()) print("Output solution...");
            std::string nstr = spade::utils::zfill(nt, 8);
            std::string filename = "prims"+nstr;
            spade::io::output_vtk(data_out, filename, grid, time_int.solution());
            if (group.isroot()) print("Done.");
            auto output_scalar = [&](const auto& kern, const std::string& scalar_name)
            {
                spade::algs::transform_to(time_int.solution(), scalar, kern);
                std::string filename_scl = scalar_name + nstr;
                spade::io::output_vtk(data_out, filename_scl, scalar);
            };
            if (output_ducr)
            {
                if (group.isroot()) print("Output sensor...");
                using sens_t = decltype(ducr);
                const auto kern = spade::omni::make_kernel([](const sens_t& sens, const auto& info)
                {
                    return sens.get_sensor(info);
                }, ducr);
                output_scalar(kern, "sensor");
                if (group.isroot()) print("Done");
            }
            if (output_sgs)
            {
                if (group.isroot()) print("Output sensor...");
                using lam_t  = decltype(lam_visc_law);
                using turb_t = decltype(eddy_visc);
                const auto kern = spade::omni::make_kernel([&](const lam_t& lam, const turb_t& turb, const auto& info)
                {
                    auto mut = turb.get_mu_t(info);
                    auto mu  = lam.get_visc(info);
                    return mut/mu;
                }, lam_visc_law, eddy_visc);
                output_scalar(kern, "turbratio");
                if (group.isroot()) print("Done");
            }
        }
        if (nt%checkpoint_skip == 0)
        {
            if (group.isroot()) print("Output checkpoint...");
            std::string nstr = spade::utils::zfill(nt, 8);
            std::string filename = "check"+nstr;
            filename = std::string(ck_dir)+"/"+filename+".bin";
            spade::io::binary_write(filename, time_int.solution());
            if (group.isroot()) print("Done.");
        }
		if (nt%nt_backup == 0)
		{
		  if (group.isroot()) print("Create backup...");
		  prim_crash = time_int.solution();
		  if (group.isroot()) print("Done"); 
		}
        tmr.start("advance");
        time_int.advance();
        tmr.stop ("advance");
        if (group.isroot()) print(tmr);

        auto   dur             = tmr.duration("advance");
        int    num_points      = cells_in_block[0]*cells_in_block[1]*cells_in_block[2]*num_blocks[0]*num_blocks[1]*num_blocks[2];
        int    num_ranks       = group.size();
        real_t updates_per_rank_per_sec = num_points/(num_ranks*dur);
        if (group.isroot()) print("Updates/core/s:", updates_per_rank_per_sec);
        perf << updates_per_rank_per_sec;
        if (group.isroot()) print("Mean:", perf.value(), "   St. dev:", perf.stdev());

        if (std::isnan(umax) || std::isnan(rhob))
        {
            if (group.isroot())
            {
                print("A tragedy has occurred!");
				print("Saving crash file...");
            }
            group.sync();
            std::string crashbase = "crash" + spade::utils::zfill(nt, 8);
			spade::io::output_vtk(data_out, crashbase, prim_crash);
            std::string cbfl = std::string(ck_dir / crashbase) + ".bin";
            spade::io::binary_write(cbfl, prim_crash);
			if (group.isroot()) print("See ya later");
            return 155;
        }
    }
    return 0;
}
