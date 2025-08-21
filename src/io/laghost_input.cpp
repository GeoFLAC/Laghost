#include <algorithm>  // For std::is_sorted
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "laghost_parameters.hpp"
#include "laghost_input.hpp"

std::map<std::string, int> bc_unit_map = {
    {"cm/yr", 1},
    {"mm/yr", 2},
    {"cm/s", 3}
};

// overarching function handling input parameters to be called during the main sequence.
void read_and_assign_input_parameters(OptionsParser& args, Param& param, const int &myid)
{
    const char* input_parameter_file = "./defaults.cfg";
    args.AddOption(&input_parameter_file, "-i", "--input", "Input parameter file to use.");

    // Let some options be overwritten by command-line options.
    args.AddOption(&param.sim.dim, "-dim", "--dimension", "Dimension of the problem.");
    args.AddOption(&param.sim.t_final, "-tf", "--t-final", "Final time; start time is 0.");
    args.AddOption(&param.sim.max_tsteps, "-ms", "--max-steps", "Maximum number of steps (negative means no restriction).");
    args.AddOption(&param.sim.visualization, "-vis", "--visualization", "-no-vis", "--no-visualization", "Enable or disable GLVis visualization.");
    args.AddOption(&param.sim.vis_steps, "-vs", "--visualization-steps", "Visualize every n-th timestep.");
    args.AddOption(&param.sim.visit, "-visit", "--visit", "-no-visit", "--no-visit", "Enable or disable VisIt visualization.");
    args.AddOption(&param.sim.paraview, "-paraview", "--paraview-datafiles", "-no-paraview", "--no-paraview-datafiles", "Save data files for ParaView (paraview.org) visualization.");
    args.AddOption(&param.sim.gfprint, "-print", "--print", "-no-print", "--no-print", "Enable or disable result output (files in mfem format).");
    // args.AddOption(param.sim.basename.c_str(), "-k", "--outputfilename", "Name of the visit dump files");
    // args.AddOption(param.sim.device.c_str(), "-d", "--device", "Device configuration string, see Device::Configure().");
    args.AddOption(&param.sim.dev, "-dev", "--dev", "GPU device to use.");
    args.AddOption(&param.sim.check, "-chk", "--checks", "-no-chk", "--no-checks", "Enable 2D checks.");
    args.AddOption(&param.sim.mem_usage, "-mb", "--mem", "-no-mem", "--no-mem", "Enable memory usage.");
    args.AddOption(&param.sim.fom, "-f", "--fom", "-no-fom", "--no-fom", "Enable figure of merit output.");
    args.AddOption(&param.sim.gpu_aware_mpi, "-gam", "--gpu-aware-mpi", "-no-gam", "--no-gpu-aware-mpi", "Enable GPU aware MPI communications.");
    // args.AddOption(&param.mesh.mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&param.mesh.rs_levels, "-rs", "--refine-serial", "Number of times to refine the mesh uniformly in serial.");
    args.AddOption(&param.mesh.rp_levels, "-rp", "--refine-parallel", "Number of times to refine the mesh uniformly in parallel.");
    args.AddOption(&param.mesh.partition_type, "-pt", "--partition", "Customized x/y/z Cartesian MPI partitioning of the serial mesh.\n\t"
                                                                       "Here x,y,z are relative task ratios in each direction.\n\t"
                                                                       "Example: with 48 mpi tasks and -pt 321, one would get a Cartesian\n\t"
                                                                       "partition of the serial mesh by (6,4,2) MPI tasks in (x,y,z).\n\t"
                                                                       "NOTE: the serially refined mesh must have the appropriate number\n\t"
                                                                       "of zones in each direction, e.g., the number of zones in direction x\n\t"
                                                                       "must be divisible by the number of MPI tasks in direction x.\n\t"
                                                                       "Available options: 11, 21, 111, 211, 221, 311, 321, 322, 432.");
    args.AddOption(&param.mesh.order_v, "-ok", "--order-kinematic", "Order (degree) of the kinematic finite element space.");
    args.AddOption(&param.mesh.order_e, "-ot", "--order-thermo", "Order (degree) of the thermodynamic finite element space.");
    args.AddOption(&param.mesh.order_q, "-oq", "--order-intrule", "Order  of the integration rule.");
    args.AddOption(&param.solver.ode_solver_type, "-s", "--ode-solver", "ODE solver: 1 - Forward Euler,\n\t"
                                                                         "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6,\n\t"
                                                                         "            7 - RK2Avg.");
    args.AddOption(&param.solver.cfl, "-cfl", "--cfl", "CFL-condition number.");
    args.AddOption(&param.solver.cg_tol, "-cgt", "--cg-tol", "Relative CG tolerance (velocity linear solve).");
    args.AddOption(&param.solver.ftz_tol, "-ftz", "--ftz-tol", "Absolute flush-to-zero tolerance.");
    args.AddOption(&param.solver.cg_max_iter, "-cgm", "--cg-max-steps", "Maximum number of CG iterations (velocity linear solve).");
    args.AddOption(&param.solver.p_assembly, "-pa", "--partial-assembly", "-fa", "--full-assembly", "Activate 1D tensor-based assembly (partial assembly).");
    args.AddOption(&param.solver.impose_visc, "-iv", "--impose-viscosity", "-niv", "--no-impose-viscosity", "Use active viscosity terms even for smooth problems.");

    // TMOP
    args.AddOption(&param.tmop.tmop, "-TMOP", "--enable-TMOP", "-no-TMOP", "--disable-TMOP", "Target Mesh Optimization Paradigm.");
    args.AddOption(&param.tmop.amr, "-amr", "--enable-amr", "-no-amr", "--disable-amr", "Adaptive mesh refinement.");
    args.AddOption(&param.tmop.remesh_steps, "-rstep", "--remesh_steps", "remeshing frequency.");
    args.AddOption(&param.tmop.jitter, "-ji", "--jitter", "Random perturbation scaling factor.");
    args.AddOption(&param.tmop.metric_id, "-mid", "--metric-id", "Mesh optimization metric:\n\t"
                                                                   "T-metrics\n\t"
                                                                   "1  : |T|^2                          -- 2D no type\n\t"
                                                                   "2  : 0.5|T|^2/tau-1                 -- 2D shape (condition number)\n\t"
                                                                   "7  : |T-T^-t|^2                     -- 2D shape+size\n\t"
                                                                   "9  : tau*|T-T^-t|^2                 -- 2D shape+size\n\t" // not applicable
                                                                   "14 : |T-I|^2                        -- 2D shape+size+orientation\n\t" // not applicable 
                                                                   "22 : 0.5(|T|^2-2*tau)/(tau-tau_0)   -- 2D untangling\n\t" // not applicable
                                                                   "50 : 0.5|T^tT|^2/tau^2-1            -- 2D shape\n\t"
                                                                   "55 : (tau-1)^2                      -- 2D size\n\t"
                                                                   "56 : 0.5(sqrt(tau)-1/sqrt(tau))^2   -- 2D size\n\t"
                                                                   "58 : |T^tT|^2/(tau^2)-2*|T|^2/tau+2 -- 2D shape\n\t"
                                                                   "77 : 0.5(tau-1/tau)^2               -- 2D size\n\t"
                                                                   "80 : (1-gamma)mu_2 + gamma mu_77    -- 2D shape+size\n\t"
                                                                   "85 : |T-|T|/sqrt(2)I|^2             -- 2D shape+orientation\n\t"
                                                                   "90 : balanced combo mu_50 & mu_77   -- 2D shape+size\n\t"
                                                                   "94 : balanced combo mu_2 & mu_56    -- 2D shape+size\n\t"
                                                                   "98 : (1/tau)|T-I|^2                 -- 2D shape+size+orientation\n\t"
                                                                   // "211: (tau-1)^2-tau+sqrt(tau^2+eps)  -- 2D untangling\n\t"
                                                                   // "252: 0.5(tau-1)^2/(tau-tau_0)       -- 2D untangling\n\t"
                                                                   "301: (|T||T^-1|)/3-1              -- 3D shape\n\t"
                                                                   "302: (|T|^2|T^-1|^2)/9-1          -- 3D shape\n\t"
                                                                   "303: (|T|^2)/3/tau^(2/3)-1        -- 3D shape\n\t"
                                                                   "304: (|T|^3)/3^{3/2}/tau-1        -- 3D shape\n\t"
                                                                   // "311: (tau-1)^2-tau+sqrt(tau^2+eps)-- 3D untangling\n\t"
                                                                   "313: (|T|^2)(tau-tau0)^(-2/3)/3   -- 3D untangling\n\t"
                                                                   "315: (tau-1)^2                    -- 3D no type\n\t"
                                                                   "316: 0.5(sqrt(tau)-1/sqrt(tau))^2 -- 3D no type\n\t"
                                                                   "321: |T-T^-t|^2                   -- 3D shape+size\n\t"
                                                                   "322: |T-adjT^-t|^2                -- 3D shape+size\n\t"
                                                                   "323: |J|^3-3sqrt(3)ln(det(J))-3sqrt(3)  -- 3D shape+size\n\t"
                                                                   "328: balanced combo mu_301 & mu_316   -- 3D shape+size\n\t"
                                                                   "332: (1-gamma) mu_302 + gamma mu_315  -- 3D shape+size\n\t"
                                                                   "333: (1-gamma) mu_302 + gamma mu_316  -- 3D shape+size\n\t"
                                                                   "334: (1-gamma) mu_303 + gamma mu_316  -- 3D shape+size\n\t"
                                                                   "328: balanced combo mu_302 & mu_318   -- 3D shape+size\n\t"
                                                                   "347: (1-gamma) mu_304 + gamma mu_316  -- 3D shape+size\n\t"
                                                                   // "352: 0.5(tau-1)^2/(tau-tau_0)     -- 3D untangling\n\t"
                                                                   "360: (|T|^3)/3^{3/2}-tau              -- 3D shape\n\t"
                                                                   "A-metrics\n\t"
                                                                   "11 : (1/4*alpha)|A-(adjA)^T(W^TW)/omega|^2 -- 2D shape\n\t"
                                                                   "36 : (1/alpha)|A-W|^2                      -- 2D shape+size+orientation\n\t"
                                                                   "107: (1/2*alpha)|A-|A|/|W|W|^2             -- 2D shape+orientation\n\t"
                                                                   "126: (1-gamma)nu_11 + gamma*nu_14a         -- 2D shape+size\n\t");
    args.AddOption(&param.tmop.target_id, "-tid", "--target-id", "Target (ideal element) type:\n\t"
                                                                   "1: Ideal shape, unit size\n\t"
                                                                   "2: Ideal shape, equal size\n\t"
                                                                   "3: Ideal shape, initial size\n\t"
                                                                   "4: Given full analytic Jacobian (in physical space)\n\t"
                                                                   "5: Ideal shape, given size (in physical space)");
    args.AddOption(&param.tmop.lim_const, "-lc", "--limit-const", "Limiting constant.");
    args.AddOption(&param.tmop.adapt_lim_const, "-alc", "--adapt-limit-const", "Adaptive limiting coefficient constant.");
    args.AddOption(&param.tmop.quad_type, "-qt", "--quad-type", "Quadrature rule type:\n\t"
                                                                   "1: Gauss-Lobatto\n\t"
                                                                   "2: Gauss-Legendre\n\t"
                                                                   "3: Closed uniform points");
    args.AddOption(&param.tmop.quad_order, "-qo", "--quad_order", "Order of the quadrature rule.");
    args.AddOption(&param.tmop.solver_type, "-st", "--solver-type", " Type of solver: (default) 0: Newton, 1: LBFGS");
    args.AddOption(&param.tmop.solver_iter, "-ni", "--newton-iters", "Maximum number of Newton iterations.");
    args.AddOption(&param.tmop.solver_rtol, "-rtol", "--newton-rel-tolerance", "Relative tolerance for the Newton solver.");
    args.AddOption(&param.tmop.solver_art_type, "-art", "--adaptive-rel-tol", "Type of adaptive relative linear solver tolerance:\n\t"
                                                                               "0: None (default)\n\t"
                                                                               "1: Eisenstat-Walker type 1\n\t"
                                                                               "2: Eisenstat-Walker type 2");
    args.AddOption(&param.tmop.lin_solver, "-ls", "--lin-solver", "Linear solver:\n\t"
                                                                   "0: l1-Jacobi\n\t"
                                                                   "1: CG\n\t"
                                                                   "2: MINRES\n\t"
                                                                   "3: MINRES + Jacobi preconditioner\n\t"
                                                                   "4: MINRES + l1-Jacobi preconditioner");
    args.AddOption(&param.tmop.max_lin_iter, "-li", "--lin-iter", "Maximum number of iterations in the linear solve.");
    args.AddOption(&param.tmop.move_bnd, "-bnd", "--move-boundary", "-fix-bnd", "--fix-boundary", "Enable motion along horizontal and vertical boundaries.");
    args.AddOption(&param.tmop.combomet, "-cmb", "--combo-type", "Combination of metrics options:\n\t"
                                                                   "0: Use single metric\n\t"
                                                                   "1: Shape + space-dependent size given analytically\n\t"
                                                                   "2: Shape + adapted size given discretely; shared target");
    args.AddOption(&param.tmop.bal_expl_combo, "-bec", "--balance-explicit-combo", "-no-bec", "--balance-explicit-combo", "Automatic balancing of explicit combo metrics.");
    args.AddOption(&param.tmop.hradaptivity, "-hr", "--hr-adaptivity", "-no-hr", "--no-hr-adaptivity", "Enable hr-adaptivity.");
    args.AddOption(&param.tmop.h_metric_id, "-hmid", "--h-metric", "Same options as metric_id. Used to determine refinement"
                                                                   " type for each element if h-adaptivity is enabled.");
    args.AddOption(&param.tmop.normalization, "-nor", "--normalization", "-no-nor", "--no-normalization", "Make all terms in the optimization functional unitless.");
    args.AddOption(&param.tmop.fdscheme, "-fd", "--fd_approximation", "-no-fd", "--no-fd-approx", "Enable finite difference based derivative computations.");
    args.AddOption(&param.tmop.exactaction, "-ex", "--exact_action", "-no-ex", "--no-exact-action", "Enable exact action of TMOP_Integrator.");
    args.AddOption(&param.tmop.verbosity_level, "-vl", "--verbosity-level", "Verbosity level for the involved iterative solvers:\n\t"
                                                                               "0: no output\n\t"
                                                                               "1: Newton iterations\n\t"
                                                                               "2: Newton iterations + linear solver summaries\n\t"
                                                                               "3: newton iterations + linear solver iterations");
    args.AddOption(&param.tmop.adapt_eval, "-ae", "--adaptivity-evaluator", "0 - Advection based (DEFAULT), 1 - GSLIB.");
    args.AddOption(&param.tmop.n_hr_iter, "-nhr", "--n_hr_iter", "Number of hr-adaptivity iterations.");
    args.AddOption(&param.tmop.n_h_iter, "-nh", "--n_h_iter", "Number of h-adaptivity iterations per r-adaptivity iteration.");
    args.AddOption(&param.tmop.mesh_node_ordering, "-mno", "--mesh_node_ordering", "Ordering of mesh nodes."
                                                                                   "0 (default): byNodes, 1: byVDIM");
    args.AddOption(&param.tmop.barrier_type, "-btype", "--barrier-type", "0 - None,"
                                                                           "1 - Shifted Barrier,"
                                                                           "2 - Pseudo Barrier.");
    args.AddOption(&param.tmop.worst_case_type, "-wctype", "--worst-case-type", "0 - None,"
                                                                                 "1 - Beta,"
                                                                                 "2 - PMean.");

    args.Parse();
    
    // Check if help was requested - if so, exit gracefully after showing usage
    if (!args.Good())
    {
        if (myid == 0)
        {
            args.PrintUsage(std::cout);
        }
        // Don't treat help request as an error - exit gracefully
        MPI_Finalize();
        exit(0);
    }
    
    // Now read the configuration file (after command line parsing)
    get_input_parameters(input_parameter_file, param);
    
    if (myid == 0)
    {
        args.PrintOptions(std::cout);
    }

    // Further determine some parameters based on the input.
    param.bc.vel_unit = 1.0; // m/s to m/s. 
    // See laghost_parameters.hpp for the mapping.
    switch (bc_unit_map[param.bc.bc_unit]) {
        case 1: // cm/yr to m/s. 
            param.bc.vel_unit = 1.0e-2 / YEAR2SEC;
            break;
        case 2: // mm/yr to m/s
            param.bc.vel_unit = 1.0e-3 / YEAR2SEC;
            break;
        case 3: // cm/s to m/s
            param.bc.vel_unit = 1.0e-2;
            break;
        default: // already in m/s
            break;
    }

    param.tmop.mesh_poly_deg = param.mesh.order_v;
    param.tmop.quad_order = 2*param.mesh.order_v - 1; // integration order = 2p  - 1

    if(param.sim.max_tsteps > -1)
        param.sim.t_final = 1.0e38; // set a large number to avoid the final time restriction.
    if(param.sim.year) {
        param.sim.t_final = param.sim.t_final * YEAR2SEC;
        if ( myid == 0)
            std::cout << "Use years in output instead of seconds is true" << std::endl;
    }
    else {
        if ( myid == 0 )
            std::cout << "Use seconds in output instead of years is true" << std::endl;
    }
}

static void declare_parameters(po::options_description &cfg,
                               Param &p)
{
    /* To have a new input parameter declared as such in parameters.hpp,
     *
     *     struct SectionType {
     *         type name;
     *     }
     *     struct Param {
     *         SectionType section;
     *     }
     *
     * add this line in this function:
     *
     *     ("section.name", po::value<type>(&p.section.name), "help string")
     *
     */

    cfg.add_options()
        ("sim.problem", po::value<int>(&p.sim.problem)->default_value(1),
         "Problem type: 1 to ??")
        ("sim.dim", po::value<int>(&p.sim.dim)->default_value(3),
         "Physical dimension of the problem: 2 or 3")
        ("sim.t_final", po::value<double>(&p.sim.t_final)->default_value(1.0),
         "Final time for the simulation in seconds")
        ("sim.max_tsteps", po::value<int>(&p.sim.max_tsteps)->default_value(-1),"Final time; start time is 0.") 
        ("sim.year", po::value<bool>(&p.sim.year)->default_value(false),
         "Is time unit year?")
        ("sim.visualization", po::value<bool>(&p.sim.visualization)->default_value(false),
         "Start in-app visualization? true or false") 
        ("sim.vis_steps", po::value<int>(&p.sim.vis_steps)->default_value(1000),
         "Visuaization output at every N step: N in integer")
        ("sim.visit", po::value<bool>(&p.sim.visit)->default_value(false),
         "Visualization outputs for VisIt? true or false")
        ("sim.paraview", po::value<bool>(&p.sim.paraview)->default_value(true),
         "Visualization outputs for ParaView? true or false")
        ("sim.gfprint", po::value<bool>(&p.sim.gfprint)->default_value(false),
         "Enable or disable result output (files in mfem format)") 
        ("sim.basename", po::value<std::string>(&p.sim.basename)->default_value("results/Laghost"),
         "Prefix for the output files")
        ("sim.device", po::value<std::string>(&p.sim.device)->default_value("cpu"),
         "Choose device, cpu or gpu")
        ("sim.dev", po::value<int>(&p.sim.dev)->default_value(0),
         "Choose a gpu device, 0, 1, ...")
        ("sim.check", po::value<bool>(&p.sim.check)->default_value(false),
         "Enable 2D checks.")
        ("sim.mem_usage", po::value<bool>(&p.sim.mem_usage)->default_value(false),
         "Print memory usage")
        ("sim.fom", po::value<bool>(&p.sim.fom)->default_value(false),
         "Enable figure of merit output.")
        ("sim.gpu_aware_mpi", po::value<bool>(&p.sim.gpu_aware_mpi)->default_value(false),
         "Enable GPU aware MPI communications.")
        ;

    cfg.add_options()
        ("solver.ode_solver_type", po::value<>(&p.solver.ode_solver_type)->default_value(7)," ")
        ("solver.cfl", po::value<double>(&p.solver.cfl)->default_value(0.5), " ")
        ("solver.cg_tol", po::value<double>(&p.solver.cg_tol)->default_value(1.0e-10)," ")
        ("solver.ftz_tol", po::value<double>(&p.solver.ftz_tol)->default_value(0.0)," ")
        ("solver.cg_max_iter", po::value<int>(&p.solver.cg_max_iter)->default_value(300)," ")
        ("solver.p_assembly", po::value<bool>(&p.solver.p_assembly)->default_value(false)," ")
        ("solver.impose_visc", po::value<bool>(&p.solver.impose_visc)->default_value(true)," ")
        ;

    cfg.add_options()
        ("control.pseudo_transient", po::value<bool>(&p.control.pseudo_transient)->default_value(false)," ")
        ("control.transient_num", po::value<int>(&p.control.transient_num)->default_value(5)," ")
        
        // ("control.flat_rate", po::value<double>(&p.control.flat_rate)->default_value(1.0e-7), " ")
        ("control.lithostatic", po::value<bool>(&p.control.lithostatic)->default_value(true)," ")
        ("control.init_dt", po::value<double>(&p.control.init_dt)->default_value(1.0), " ")
        ("control.mscale", po::value<double>(&p.control.mscale)->default_value(5.0e5), " ")
        ("control.gravity", po::value<double>(&p.control.gravity)->default_value(10.0), " ")
        ("control.thickness", po::value<double>(&p.control.thickness)->default_value(10.0e3), " ")
        ("control.mass_bal", po::value<bool>(&p.control.mass_bal)->default_value(false)," ")
        ("control.dyn_damping", po::value<bool>(&p.control.dyn_damping)->default_value(true)," ")
        ("control.dyn_factor", po::value<double>(&p.control.dyn_factor)->default_value(0.8), " ")
        // lower limit of maximum velocity, 0.1 mm/yr, to prevent inf. massscaling.
        ("control.max_vbc_val", po::value<double>(&p.control.max_vbc_val)->default_value(3.1709791983764588e-12), " ") 
        ;

    cfg.add_options()
        ("mesh.mesh_file", po::value<std::string>(&p.mesh.mesh_file)->default_value("default")," ")
        ("mesh.rs_levels", po::value<int>(&p.mesh.rs_levels)->default_value(2)," ")
        ("mesh.rp_levels", po::value<int>(&p.mesh.rp_levels)->default_value(0)," ")
        ("mesh.partition_type", po::value<int>(&p.mesh.partition_type)->default_value(0),
         "Partition type as an integer")
        ("mesh.order_v", po::value<int>(&p.mesh.order_v)->default_value(2),"Order (degree) of the kinematic finite element space.")
        ("mesh.order_e", po::value<int>(&p.mesh.order_e)->default_value(1),"Order (degree) of the thermodynamic finite element space.")
        ("mesh.order_q", po::value<int>(&p.mesh.order_q)->default_value(-1),"Order  of the integration rule.")
        ("mesh.local_refinement", po::value<bool>(&p.mesh.local_refinement)->default_value(false), " ")
        ("mesh.l2_basis", po::value<>(&p.mesh.l2_basis)->default_value(1)," ")
        ;

    cfg.add_options()
        ("bc.bc_unit", po::value<std::string>(&p.bc.bc_unit)->default_value("cm/yr"),"Unit of Velocity")
        ("bc.bc_ids", po::value<std::string>(&p.bc.bc_ids)->default_value("[0]"),"Boundary indicators '[d0, d1, d2, ...]")
        ("bc.bc_vxs", po::value<std::string>(&p.bc.bc_vxs)->default_value("[0]"), "Boundary veloicty x '[d0, d1, d2, ...]")
        ("bc.bc_vys", po::value<std::string>(&p.bc.bc_vys)->default_value("[0]"), "Boundary velocity y '[d0, d1, d2, ...]")
        ("bc.bc_vzs", po::value<std::string>(&p.bc.bc_vzs)->default_value("[0]"), "Boundary velocity z '[d0, d1, d2, ...]")
        ("bc.winkler_foundation", po::value<bool>(&p.bc.winkler_foundation)->default_value(false)," ")
        ("bc.winkler_flat", po::value<bool>(&p.bc.winkler_flat)->default_value(false)," ")
        ("bc.winkler_rho", po::value<double>(&p.bc.winkler_rho)->default_value(2700.0), " ")
        ("bc.surf_proc", po::value<bool>(&p.bc.surf_proc)->default_value(true)," ")
        ("bc.surf_diff", po::value<double>(&p.bc.surf_diff)->default_value(1.0e-7), " ")
        ("bc.surf_alpha", po::value<double>(&p.bc.surf_alpha)->default_value(0.0), " ")
        ("bc.base_proc", po::value<bool>(&p.bc.base_proc)->default_value(true)," ")
        ("bc.base_diff", po::value<double>(&p.bc.base_diff)->default_value(1.0e-7), " ")
        ("bc.base_alpha", po::value<double>(&p.bc.base_alpha)->default_value(0.0), " ")
        ;

    cfg.add_options()
        ("mat.plastic", po::value<bool>(&p.mat.plastic)->default_value(true), " ")
        ("mat.viscoplastic", po::value<bool>(&p.mat.viscoplastic)->default_value(false), " ")
        ("mat.rho", po::value<std::string>(&p.mat.rho)->default_value("[2700.0]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.lambda", po::value<std::string>(&p.mat.lambda)->default_value("[3e10]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.mu", po::value<std::string>(&p.mat.mu)->default_value("[3e10]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.tension_cutoff", po::value<std::string>(&p.mat.tension_cutoff)->default_value("[0.0]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.cohesion0", po::value<std::string>(&p.mat.cohesion0)->default_value("[44.0e6]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.cohesion1", po::value<std::string>(&p.mat.cohesion1)->default_value("[44.0e6]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.friction_angle0", po::value<std::string>(&p.mat.friction_angle0)->default_value("[30.0]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.friction_angle1", po::value<std::string>(&p.mat.friction_angle1)->default_value("[30.0]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.dilation_angle0", po::value<std::string>(&p.mat.dilation_angle0)->default_value("[0.0]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.dilation_angle1", po::value<std::string>(&p.mat.dilation_angle1)->default_value("[0.0]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.pls0", po::value<std::string>(&p.mat.pls0)->default_value("[0.0]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.pls1", po::value<std::string>(&p.mat.pls1)->default_value("[0.5]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.plastic_viscosity", po::value<std::string>(&p.mat.plastic_viscosity)->default_value("[1.0]"),"Material indicators '[d0, d1, d2, ...]")
        ("mat.weak_rad", po::value<double>(&p.mat.weak_rad)->default_value(1.0e3), "circular weakzone")//
        ("mat.weak_x", po::value<double>(&p.mat.weak_x)->default_value(50.0e3), " x coord of circular")//
        ("mat.weak_y", po::value<double>(&p.mat.weak_y)->default_value(2.0e3), "y coord of circular")  //
        ("mat.weak_z", po::value<double>(&p.mat.weak_z)->default_value(0.0), "z coord of circular")    //
        ("mat.ini_pls", po::value<double>(&p.mat.ini_pls)->default_value(0.5), "initial plasticity")   //
        ;
    cfg.add_options()
        ("tmop.tmop", po::value<bool>(&p.tmop.tmop)->default_value(false), " ")
        ("tmop.amr", po::value<bool>(&p.tmop.amr)->default_value(false), " ")
        ("tmop.ale", po::value<double>(&p.tmop.ale)->default_value(0.5), " ")
        ("tmop.remesh_steps", po::value<int>(&p.tmop.remesh_steps)->default_value(50000), " ")
        ("tmop.mesh_poly_deg", po::value<int>(&p.tmop.mesh_poly_deg)->default_value(2), " ")
        ("tmop.jitter", po::value<double>(&p.tmop.jitter)->default_value(0.0), " ")
        ("tmop.metric_id", po::value<int>(&p.tmop.metric_id)->default_value(2), " ")
        ("tmop.target_id", po::value<int>(&p.tmop.target_id)->default_value(1), " ")
        ("tmop.lim_const", po::value<double>(&p.tmop.lim_const)->default_value(0.0), " ")
        ("tmop.adapt_lim_const", po::value<double>(&p.tmop.adapt_lim_const)->default_value(0.0), " ")
        ("tmop.quad_type", po::value<int>(&p.tmop.quad_type)->default_value(1), " ")
        ("tmop.quad_order", po::value<int>(&p.tmop.quad_order)->default_value(8), " ")
        ("tmop.solver_type", po::value<int>(&p.tmop.solver_type)->default_value(0), " ")
        ("tmop.solver_iter", po::value<int>(&p.tmop.solver_iter)->default_value(20), " ")
        ("tmop.solver_rtol", po::value<double>(&p.tmop.solver_rtol)->default_value(1e-10), " ")
        ("tmop.solver_art_type", po::value<int>(&p.tmop.solver_art_type)->default_value(0), " ")
        ("tmop.lin_solver", po::value<int>(&p.tmop.lin_solver)->default_value(2), " ")
        ("tmop.max_lin_iter", po::value<int>(&p.tmop.max_lin_iter)->default_value(100), " ")
        ("tmop.move_bnd", po::value<bool>(&p.tmop.move_bnd)->default_value(false), " ")
        ("tmop.combomet", po::value<int>(&p.tmop.combomet)->default_value(0), " ")
        ("tmop.bal_expl_combo", po::value<bool>(&p.tmop.bal_expl_combo)->default_value(false), " ")
        ("tmop.hradaptivity", po::value<bool>(&p.tmop.hradaptivity)->default_value(false), " ")
        ("tmop.h_metric_id", po::value<int>(&p.tmop.h_metric_id)->default_value(-1), " ")
        ("tmop.normalization", po::value<bool>(&p.tmop.normalization)->default_value(false), " ")
        ("tmop.verbosity_level", po::value<int>(&p.tmop.verbosity_level)->default_value(0), " ")
        ("tmop.fdscheme", po::value<bool>(&p.tmop.fdscheme)->default_value(false), " ")
        ("tmop.adapt_eval", po::value<int>(&p.tmop.adapt_eval)->default_value(0), " ")
        ("tmop.exactaction", po::value<bool>(&p.tmop.exactaction)->default_value(false), " ")
        ("tmop.n_hr_iter", po::value<int>(&p.tmop.n_hr_iter)->default_value(5), " ")
        ("tmop.n_h_iter", po::value<int>(&p.tmop.n_h_iter)->default_value(1), " ")
        ("tmop.mesh_node_ordering", po::value<int>(&p.tmop.mesh_node_ordering)->default_value(0), " ")
        ("tmop.barrier_type", po::value<int>(&p.tmop. barrier_type)->default_value(0), " ")
        ("tmop.worst_case_type", po::value<int>(&p.tmop.worst_case_type)->default_value(0), " ")
        ("tmop.tmop_cond_num", po::value<double>(&p.tmop.tmop_cond_num)->default_value(0.5), " ")
        ;
}


static void read_parameters_from_file
(const char* filename,
 const po::options_description cfg,
 po::variables_map &vm)
{
    try {
        po::store(po::parse_config_file<char>(filename, cfg), vm);
        po::notify(vm);
    }
    catch (const boost::program_options::multiple_occurrences& e) {
        std::cerr << e.what() << " from option: " << e.get_option_name() << '\n';
        std::exit(1);
    }
    catch (std::exception& e) {
        std::cerr << "Error reading config_file '" << filename << "'\n";
        std::cerr << e.what() << "\n";
        std::exit(1);
    }
}

static void get_input_parameters(const char* filename, Param& p)
{
    po::options_description cfg("Config file options");
    po::variables_map vm;

    declare_parameters(cfg, p);
    // print help message
    if (std::strncmp(filename, "-h", 3) == 0 ||
        std::strncmp(filename, "--help", 7) == 0) {
        std::cout << cfg;
        std::exit(0);
    }
    read_parameters_from_file(filename, cfg, vm);
    // validate_parameters(vm, p);
}



#if 0
template<class T>
static int read_numbers(const std::string &input, std::vector<T> &vec, int len)
{
    /* Read 'len' numbers from input.
     * The format of input must be '[n0, n1, n2]' or '[n0, n1, n2,]' (with a trailing ,),
     * for len=3.
     */

    std::istringstream stream(input);
    vec.resize(len);

    char sentinel;

    stream >> sentinel;
    if (sentinel != '[') return 1;

    for (int i=0; i<len; ++i) {
        stream >> vec[i];

        if (i == len-1) break;

        // consume ','
        char sep;
        stream >> sep;
        if (sep != ',') return 1;
    }

    stream >> sentinel;
    if (sentinel == ',') stream >> sentinel;
    if (sentinel != ']') return 1;

    if (! stream.good()) return 1;

    // success
    return 0;
}

template<class T>
static void get_numbers(const po::variables_map &vm, const char *name,
                        std::vector<T> &values, int len, int optional_size=0)
{
    if ( ! vm.count(name) ) {
        std::cerr << "Error: " << name << " is not provided.\n";
        std::exit(1);
    }

    std::string str = vm[name].as<std::string>();
    int err = read_numbers(str, values, len);
    if (err && optional_size) {
        err = read_numbers(str, values, optional_size);
    }

    if (err) {
        std::cerr << "Error: incorrect format for " << name << ",\n"
                  << "       must be '[d0, d1, d2, ...]'\n";
        std::exit(1);
    }
}

static void validate_parameters(const po::variables_map &vm, Param &p)
{
    std::cout << "Checking consistency of input parameters...\n";

    //
    // stopping condition and output interval are based on either model time or step
    //
    if ( ! (vm.count("sim.max_steps") || vm.count("sim.max_time_in_yr")) ) {
        std::cerr << "Must provide either sim.max_steps or sim.max_time_in_yr\n";
        std::exit(1);
    }
    if ( ! vm.count("sim.max_steps") )
        p.sim.max_steps = std::numeric_limits<int>::max();
    if ( ! vm.count("sim.max_time_in_yr") )
        p.sim.max_time_in_yr = std::numeric_limits<double>::max();

    if ( ! (vm.count("sim.output_step_interval") || vm.count("sim.output_time_interval_in_yr")) ) {
        std::cerr << "Must provide either sim.output_step_interval or sim.output_time_interval_in_yr\n";
        std::exit(1);
    }
    if ( ! vm.count("sim.output_step_interval") )
        p.sim.output_step_interval = std::numeric_limits<int>::max();
    if ( ! vm.count("sim.output_time_interval_in_yr") )
        p.sim.output_time_interval_in_yr = std::numeric_limits<double>::max();

    //
    // These parameters are required when restarting
    //
    if (p.sim.is_restarting) {
        if ( ! vm.count("sim.restarting_from_modelname") ) {
            std::cerr << "Must provide sim.restarting_from_modelname when restarting.\n";
            std::exit(1);
        }
        if ( ! vm.count("sim.restarting_from_frame") ) {
            std::cerr << "Must provide sim.restarting_from_frame when restarting.\n";
            std::exit(1);
        }
    }

    //
    // these parameters are required in mesh.meshing_option == 2
    //
    if (p.mesh.meshing_option == 2) {
        if ( ! vm.count("mesh.refined_zonex") ||
#ifdef THREED
             ! vm.count("mesh.refined_zoney") ||
#endif
             ! vm.count("mesh.refined_zonez") ) {
        std::cerr << "Must provide mesh.refined_zonex, "
#ifdef THREED
                  << "mesh.refined_zoney, "
#endif
                  << "mesh.refined_zonez.\n";
        std::exit(1);
        }

        /* get 2 numbers from the string */
        double_vec tmp;
        int err;
        std::string str;
        str = vm["mesh.refined_zonex"].as<std::string>();
        err = read_numbers(str, tmp, 2);
        if (err || tmp[0] < 0 || tmp[1] > 1 || tmp[0] > tmp[1]) {
            std::cerr << "Error: incorrect value for mesh.refine_zonex,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zonex.first = tmp[0];
        p.mesh.refined_zonex.second = tmp[1];
#ifdef THREED
        str = vm["mesh.refined_zoney"].as<std::string>();
        err = read_numbers(str, tmp, 2);
        if (err || tmp[0] < 0 || tmp[1] > 1 || tmp[0] > tmp[1]) {
            std::cerr << "Error: incorrect value for mesh.refine_zoney,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zoney.first = tmp[0];
        p.mesh.refined_zoney.second = tmp[1];
#endif
        str = vm["mesh.refined_zonez"].as<std::string>();
        err = read_numbers(str, tmp, 2);
        if (err || tmp[0] < 0 || tmp[1] > 1 || tmp[0] > tmp[1]) {
            std::cerr << "Error: incorrect value for mesh.refine_zonez,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zonez.first = tmp[0];
        p.mesh.refined_zonez.second = tmp[1];
    }

    if (p.mesh.smallest_size > p.mesh.largest_size) {
        std::cerr << "Error: mesh.smallest_size is greater than mesh.largest_size.\n";
        std::exit(1);
    }

#ifdef THREED
    if (p.mesh.remeshing_option == 2) {
        std::cerr << "Error: mesh.remeshing_option=2 is not available in 3D.\n";
        std::exit(1);
    }
#endif

    //
    // bc
    //
    {
        if ( p.bc.has_winkler_foundation && p.control.gravity == 0 ) {
            p.bc.has_winkler_foundation = 0;
            std::cerr << "Warning: no gravity, Winkler foundation is turned off.\n";
        }
        if ( p.bc.has_winkler_foundation && p.bc.vbc_z0 != 0 ) {
            p.bc.vbc_z0 = 0;
            std::cerr << "Warning: Winkler foundation is turned on, setting bc.vbc_z0 to 0.\n";
        }
        if ( p.bc.has_water_loading && p.control.gravity == 0 ) {
            p.bc.has_water_loading = 0;
            std::cerr << "Warning: no gravity, water loading is turned off.\n";
        }
        if ( p.bc.has_water_loading && p.bc.vbc_z1 != 0 ) {
            p.bc.vbc_z1 = 0;
            std::cerr << "Warning: water loading is turned on, setting bc.vbc_z1 to 0.\n";
        }

        if ( p.bc.vbc_z0 > 3) {
            std::cerr << "Error: bc.vbc_z0 is not 0, 1, 2, or 3.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_z1 > 3) {
            std::cerr << "Error: bc.vbc_z0 is not 0, 1, 2, or 3.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_n0 != 1 && p.bc.vbc_n0 != 3 && p.bc.vbc_n0 != 11 && p.bc.vbc_n0 != 13 ) {
            std::cerr << "Error: bc.vbc_n0 is not 1, 3, 11, or 13.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_n1 != 1 && p.bc.vbc_n1 != 3 && p.bc.vbc_n1 != 11 && p.bc.vbc_n1 != 13 ) {
            std::cerr << "Error: bc.vbc_n1 is not 1, 3, 11, or 13.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_n2 != 1 && p.bc.vbc_n2 != 3 && p.bc.vbc_n2 != 11 && p.bc.vbc_n2 != 13 ) {
            std::cerr << "Error: bc.vbc_n2 is not 1, 3, 11, or 13.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_n3 != 1 && p.bc.vbc_n3 != 3 && p.bc.vbc_n3 != 11 && p.bc.vbc_n3 != 13 ) {
            std::cerr << "Error: bc.vbc_n3 is not 1, 3, 11, or 13.\n";
            std::exit(1);
        }
    }

    //
    // control
    //
    {
        if ( p.control.dt_fraction < 0 || p.control.dt_fraction > 1 ) {
            std::cerr << "Error: control.dt_fraction must be between 0 and 1.\n";
            std::exit(1);
        }
        if ( p.control.damping_factor < 0 || p.control.damping_factor > 1 ) {
            std::cerr << "Error: control.damping_factor must be between 0 and 1.\n";
            std::exit(1);
        }

    }

    //
    // ic
    //
    {
        if ( p.ic.mattype_option == 1) {
            get_numbers(vm, "ic.layer_mattypes", p.ic.layer_mattypes, p.ic.num_mattype_layers);
            get_numbers(vm, "ic.mattype_layer_depths", p.ic.mattype_layer_depths, p.ic.num_mattype_layers-1);
            // mattype_layer_depths must be already sorted
            if (! std::is_sorted(p.ic.mattype_layer_depths.begin(), p.ic.mattype_layer_depths.end())) {
                std::cerr << "Error: the content of ic.mattype_layer_depths is not ordered from"
                    " small to big values.\n";
                std::exit(1);
            }
        }
    }

    //
    // marker
    //
    {

    }

    //
    // material properties
    //
    {
        std::string str = vm["mat.rheology_type"].as<std::string>();
        if (str == std::string("elastic"))
            p.mat.rheol_type = MatProps::rh_elastic;
        else if (str == std::string("viscous"))
            p.mat.rheol_type = MatProps::rh_viscous;
        else if (str == std::string("maxwell"))
            p.mat.rheol_type = MatProps::rh_maxwell;
        else if (str == std::string("elasto-plastic"))
            p.mat.rheol_type = MatProps::rh_ep;
        else if (str == std::string("elasto-visco-plastic"))
            p.mat.rheol_type = MatProps::rh_evp;
        else {
            std::cerr << "Error: unknown rheology: '" << str << "'\n";
            std::exit(1);
        }

#ifdef THREED
        if ( p.mat.is_plane_strain ) {
            p.mat.is_plane_strain = false;
            std::cerr << "Warning: mat.is_plane_strain is not avaiable in 3D.\n";
        }
#endif

        if (p.mat.phase_change_option != 0 && p.mat.nmat == 1) {
            std::cerr << "Error: mat.phase_change_option is chosen, but mat.num_materials is 1.\n";
            std::exit(1);
        }
        if (p.mat.phase_change_option == 1 && p.mat.nmat < 8) {
            std::cerr << "Error: mat.phase_change_option is 1, but mat.num_materials is less than 8.\n";
            std::exit(1);
        }

        if (p.mat.nmat < 1) {
            std::cerr << "Error: mat.num_materials must be greater than 0.\n";
            std::exit(1);
        }

        if (p.mat.nmat == 1 && p.control.ref_pressure_option != 0) {
            p.control.ref_pressure_option = 0;
            std::cerr << "Warning: mat.num_materials is 1, using simplest control.ref_pressure_option.\n";
        }
        if (p.mat.nmat == 1 && p.markers.replenishment_option != 1) {
            p.markers.replenishment_option = 1;
            std::cerr << "Warning: mat.num_materials is 1, using simplest markers.replenishment_option.\n";
        }

        get_numbers(vm, "mat.rho0", p.mat.rho0, p.mat.nmat, 1);
        get_numbers(vm, "mat.alpha", p.mat.alpha, p.mat.nmat, 1);

        get_numbers(vm, "mat.bulk_modulus", p.mat.bulk_modulus, p.mat.nmat, 1);
        get_numbers(vm, "mat.shear_modulus", p.mat.shear_modulus, p.mat.nmat, 1);

        get_numbers(vm, "mat.visc_exponent", p.mat.visc_exponent, p.mat.nmat, 1);
        get_numbers(vm, "mat.visc_coefficient", p.mat.visc_coefficient, p.mat.nmat, 1);
        get_numbers(vm, "mat.visc_activation_energy", p.mat.visc_activation_energy, p.mat.nmat, 1);

        get_numbers(vm, "mat.heat_capacity", p.mat.heat_capacity, p.mat.nmat, 1);
        get_numbers(vm, "mat.therm_cond", p.mat.therm_cond, p.mat.nmat, 1);

        get_numbers(vm, "mat.pls0", p.mat.pls0, p.mat.nmat, 1);
        get_numbers(vm, "mat.pls1", p.mat.pls1, p.mat.nmat, 1);
        get_numbers(vm, "mat.cohesion0", p.mat.cohesion0, p.mat.nmat, 1);
        get_numbers(vm, "mat.cohesion1", p.mat.cohesion1, p.mat.nmat, 1);
        get_numbers(vm, "mat.friction_angle0", p.mat.friction_angle0, p.mat.nmat, 1);
        get_numbers(vm, "mat.friction_angle1", p.mat.friction_angle1, p.mat.nmat, 1);
        get_numbers(vm, "mat.dilation_angle0", p.mat.dilation_angle0, p.mat.nmat, 1);
        get_numbers(vm, "mat.dilation_angle1", p.mat.dilation_angle1, p.mat.nmat, 1);
    }

}
#endif
