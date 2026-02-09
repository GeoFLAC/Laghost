// laghost_input_cli11.cpp
// CLI11-based configuration file parser for Laghost
// This file replaces Boost.Program_options with CLI11 for TOML config parsing

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include "extern/CLI11.hpp"
#include "laghost_input.hpp"
#include "laghost_parameters.hpp"

using namespace mfem;

std::map<std::string, int> bc_unit_map = {
    {"m/s", 0},
    {"cm/yr", 1},
    {"mm/yr", 2},
    {"cm/s", 3}
};

// Forward declarations
static void validate_parameters(Param &p);
static void parse_material_arrays(Param &p);

// Helper to parse "[v1, v2, v3]" string format into mfem::Vector
template<class T>
static int read_numbers(const std::string &input, std::vector<T> &vec, int len)
{
    std::istringstream stream(input);
    vec.resize(len);
    char sentinel;

    stream >> sentinel;
    if (sentinel != '[') return 1;

    for (int i = 0; i < len; ++i) {
        stream >> vec[i];
        if (i == len - 1) break;
        char sep;
        stream >> sep;
        if (sep != ',') return 1;
    }

    stream >> sentinel;
    if (sentinel == ',') stream >> sentinel;
    if (sentinel != ']') return 1;
    if (!stream.good() && !stream.eof()) return 1;

    return 0;
}

template<class T>
static void parse_array_string(const std::string &str, const char* name, 
                               Vector &values, int len)
{
    std::vector<T> temp_values;
    int err = read_numbers(str, temp_values, len);
    
    if (err) {
        std::cerr << "Error: incorrect format for " << name << ",\n"
                  << "       must be '[d0, d1, d2, ...]'\n";
        std::exit(1);
    }
    
    values.SetSize(temp_values.size());
    for (size_t i = 0; i < temp_values.size(); ++i) {
        values[i] = temp_values[i];
    }
}

// ============================================================================
// Main entry point for reading input parameters
// ============================================================================
void read_and_assign_input_parameters(OptionsParser& args, Param& param, const int &myid)
{
    const char* input_parameter_file = "./defaults.cfg";
    args.AddOption(&input_parameter_file, "-i", "--input", "Input parameter file to use.");

    // Command-line options override config file values (handled by MFEM's OptionsParser)
    args.AddOption(&param.sim.dim, "-dim", "--dimension", "Dimension of the problem.");
    args.AddOption(&param.sim.t_final, "-tf", "--t-final", "Final time; start time is 0.");
    args.AddOption(&param.sim.max_tsteps, "-ms", "--max-steps", "Maximum number of steps.");
    args.AddOption(&param.sim.visualization, "-vis", "--visualization", "-no-vis", "--no-visualization", "Enable GLVis visualization.");
    args.AddOption(&param.sim.vis_steps, "-vs", "--visualization-steps", "Visualize every n-th timestep.");
    args.AddOption(&param.sim.visit, "-visit", "--visit", "-no-visit", "--no-visit", "Enable VisIt visualization.");
    args.AddOption(&param.sim.paraview, "-paraview", "--paraview-datafiles", "-no-paraview", "--no-paraview-datafiles", "Save ParaView data files.");
    args.AddOption(&param.sim.gfprint, "-print", "--print", "-no-print", "--no-print", "Enable MFEM format output.");
    args.AddOption(&param.sim.dev, "-dev", "--dev", "GPU device to use.");
    args.AddOption(&param.sim.check, "-chk", "--checks", "-no-chk", "--no-checks", "Enable 2D checks.");
    args.AddOption(&param.sim.mem_usage, "-mb", "--mem", "-no-mem", "--no-mem", "Enable memory usage.");
    args.AddOption(&param.sim.fom, "-f", "--fom", "-no-fom", "--no-fom", "Enable figure of merit output.");
    args.AddOption(&param.sim.gpu_aware_mpi, "-gam", "--gpu-aware-mpi", "-no-gam", "--no-gpu-aware-mpi", "Enable GPU aware MPI.");
    args.AddOption(&param.mesh.rs_levels, "-rs", "--refine-serial", "Serial mesh refinement levels.");
    args.AddOption(&param.mesh.rp_levels, "-rp", "--refine-parallel", "Parallel mesh refinement levels.");
    args.AddOption(&param.mesh.partition_type, "-pt", "--partition", "MPI partition type.");
    args.AddOption(&param.mesh.order_v, "-ok", "--order-kinematic", "Kinematic FE order.");
    args.AddOption(&param.mesh.order_e, "-ot", "--order-thermo", "Thermodynamic FE order.");
    args.AddOption(&param.mesh.order_q, "-oq", "--order-intrule", "Integration rule order.");
    args.AddOption(&param.solver.ode_solver_type, "-s", "--ode-solver", "ODE solver type.");
    args.AddOption(&param.solver.cfl, "-cfl", "--cfl", "CFL number.");
    args.AddOption(&param.solver.cg_tol, "-cgt", "--cg-tol", "CG tolerance.");
    args.AddOption(&param.solver.ftz_tol, "-ftz", "--ftz-tol", "Flush-to-zero tolerance.");
    args.AddOption(&param.solver.cg_max_iter, "-cgm", "--cg-max-steps", "Maximum CG iterations.");
    args.AddOption(&param.solver.p_assembly, "-pa", "--partial-assembly", "-fa", "--full-assembly", "Partial assembly mode.");
    args.AddOption(&param.solver.impose_visc, "-iv", "--impose-viscosity", "-niv", "--no-impose-viscosity", "Impose viscosity.");

    // TMOP options
    args.AddOption(&param.tmop.tmop, "-TMOP", "--enable-TMOP", "-no-TMOP", "--disable-TMOP", "Enable TMOP.");
    args.AddOption(&param.tmop.amr, "-amr", "--enable-amr", "-no-amr", "--disable-amr", "Enable AMR.");
    args.AddOption(&param.tmop.remesh_steps, "-rstep", "--remesh_steps", "Remeshing frequency.");
    args.AddOption(&param.tmop.jitter, "-ji", "--jitter", "Random perturbation scaling.");
    args.AddOption(&param.tmop.metric_id, "-mid", "--metric-id", "Mesh optimization metric.");
    args.AddOption(&param.tmop.target_id, "-tid", "--target-id", "Target element type.");
    args.AddOption(&param.tmop.lim_const, "-lc", "--limit-const", "Limiting constant.");
    args.AddOption(&param.tmop.adapt_lim_const, "-alc", "--adapt-limit-const", "Adaptive limiting constant.");
    args.AddOption(&param.tmop.quad_type, "-qt", "--quad-type", "Quadrature rule type.");
    args.AddOption(&param.tmop.quad_order, "-qo", "--quad_order", "Quadrature order.");
    args.AddOption(&param.tmop.solver_type, "-st", "--solver-type", "Solver type.");
    args.AddOption(&param.tmop.solver_iter, "-ni", "--newton-iters", "Newton iterations.");
    args.AddOption(&param.tmop.solver_rtol, "-rtol", "--newton-rel-tolerance", "Newton relative tolerance.");
    args.AddOption(&param.tmop.solver_art_type, "-art", "--adaptive-rel-tol", "Adaptive relative tolerance type.");
    args.AddOption(&param.tmop.lin_solver, "-ls", "--lin-solver", "Linear solver type.");
    args.AddOption(&param.tmop.max_lin_iter, "-li", "--lin-iter", "Maximum linear iterations.");
    args.AddOption(&param.tmop.move_bnd, "-bnd", "--move-boundary", "-fix-bnd", "--fix-boundary", "Move boundary.");
    args.AddOption(&param.tmop.combomet, "-cmb", "--combo-type", "Combination metrics type.");
    args.AddOption(&param.tmop.bal_expl_combo, "-bec", "--balance-explicit-combo", "-no-bec", "--no-balance-explicit-combo", "Balance explicit combo.");
    args.AddOption(&param.tmop.hradaptivity, "-hr", "--hr-adaptivity", "-no-hr", "--no-hr-adaptivity", "HR adaptivity.");
    args.AddOption(&param.tmop.h_metric_id, "-hmid", "--h-metric", "H-adaptivity metric.");
    args.AddOption(&param.tmop.normalization, "-nor", "--normalization", "-no-nor", "--no-normalization", "Normalization.");
    args.AddOption(&param.tmop.fdscheme, "-fd", "--fd_approximation", "-no-fd", "--no-fd-approx", "Finite difference scheme.");
    args.AddOption(&param.tmop.exactaction, "-ex", "--exact_action", "-no-ex", "--no-exact-action", "Exact action.");
    args.AddOption(&param.tmop.verbosity_level, "-vl", "--verbosity-level", "Verbosity level.");
    args.AddOption(&param.tmop.adapt_eval, "-ae", "--adaptivity-evaluator", "Adaptivity evaluator.");
    args.AddOption(&param.tmop.n_hr_iter, "-nhr", "--n_hr_iter", "HR iterations.");
    args.AddOption(&param.tmop.n_h_iter, "-nh", "--n_h_iter", "H iterations.");
    args.AddOption(&param.tmop.mesh_node_ordering, "-mno", "--mesh_node_ordering", "Mesh node ordering.");
    args.AddOption(&param.tmop.barrier_type, "-btype", "--barrier-type", "Barrier type.");
    args.AddOption(&param.tmop.worst_case_type, "-wctype", "--worst-case-type", "Worst case type.");

    args.Parse();
    
    if (!args.Good())
    {
        if (myid == 0) args.PrintUsage(std::cout);
        MPI_Finalize();
        exit(0);
    }
    
    // Read configuration file using CLI11
    get_input_parameters(input_parameter_file, param);
    
    if (myid == 0) args.PrintOptions(std::cout);

    param.tmop.mesh_poly_deg = param.mesh.order_v;
    param.tmop.quad_order = 2 * param.mesh.order_v - 1;

    if (param.sim.max_tsteps > -1)
        param.sim.t_final = 1.0e38;
    if (param.sim.year) {
        param.sim.t_final = param.sim.t_final * YEAR2SEC;
        if (myid == 0) std::cout << "Use years in output instead of seconds is true" << std::endl;
    } else {
        if (myid == 0) std::cout << "Use seconds in output instead of years is true" << std::endl;
    }
}

// ============================================================================
// CLI11-based configuration file parser
// ============================================================================
static void get_input_parameters(const char* filename, Param& p)
{
    // Handle help request
    if (std::strncmp(filename, "-h", 3) == 0 || std::strncmp(filename, "--help", 7) == 0) {
        std::cout << "Usage: laghost -i <config.toml>\n"
                  << "Config file format: TOML with sections [sim], [solver], [mesh], [bc], [mat], [tmop]\n";
        std::exit(0);
    }

    // Check file exists
    std::ifstream f(filename);
    if (!f.good()) {
        std::cerr << "Error: Cannot open configuration file '" << filename << "'\n";
        std::exit(1);
    }
    f.close();

    CLI::App app{"Laghost configuration parser"};
    
    // Temporary string holders for material arrays
    std::string rho_str, lambda_str, mu_str, tension_cutoff_str;
    std::string cohesion0_str, cohesion1_str;
    std::string friction_angle0_str, friction_angle1_str;
    std::string dilation_angle0_str, dilation_angle1_str;
    std::string alpha0_str, alpha1_str, plastic_viscosity_str;
    
    // ========== [sim] section ==========
    app.add_option("sim.problem", p.sim.problem)->default_val(1);
    app.add_option("sim.dim", p.sim.dim)->default_val(3);
    app.add_option("sim.t_final", p.sim.t_final)->default_val(1.0);
    app.add_option("sim.max_tsteps", p.sim.max_tsteps)->default_val(-1);
    app.add_option("sim.year", p.sim.year)->default_val(false);
    app.add_option("sim.visualization", p.sim.visualization)->default_val(false);
    app.add_option("sim.vis_steps", p.sim.vis_steps)->default_val(1000);
    app.add_option("sim.visit", p.sim.visit)->default_val(false);
    app.add_option("sim.paraview", p.sim.paraview)->default_val(true);
    app.add_option("sim.gfprint", p.sim.gfprint)->default_val(false);
    app.add_option("sim.basename", p.sim.basename)->default_val("results/Laghost");
    app.add_option("sim.device", p.sim.device)->default_val("cpu");
    app.add_option("sim.dev", p.sim.dev)->default_val(0);
    app.add_option("sim.check", p.sim.check)->default_val(false);
    app.add_option("sim.mem_usage", p.sim.mem_usage)->default_val(false);
    app.add_option("sim.fom", p.sim.fom)->default_val(false);
    app.add_option("sim.gpu_aware_mpi", p.sim.gpu_aware_mpi)->default_val(false);
    
    // ========== [solver] section ==========
    app.add_option("solver.ode_solver_type", p.solver.ode_solver_type)->default_val(7);
    app.add_option("solver.cfl", p.solver.cfl)->default_val(0.5);
    app.add_option("solver.cg_tol", p.solver.cg_tol)->default_val(1.0e-10);
    app.add_option("solver.ftz_tol", p.solver.ftz_tol)->default_val(0.0);
    app.add_option("solver.cg_max_iter", p.solver.cg_max_iter)->default_val(300);
    app.add_option("solver.p_assembly", p.solver.p_assembly)->default_val(false);
    app.add_option("solver.impose_visc", p.solver.impose_visc)->default_val(true);
    
    // ========== [control] section ==========
    app.add_option("control.lithostatic", p.control.lithostatic)->default_val(true);
    app.add_option("control.atmospheric", p.control.atmospheric)->default_val(false);
    app.add_option("control.init_dt", p.control.init_dt)->default_val(1.0);
    app.add_option("control.mscale", p.control.mscale)->default_val(5.0e5);
    app.add_option("control.gravity", p.control.gravity)->default_val(10.0);
    app.add_option("control.thickness", p.control.thickness)->default_val(10.0e3);
    app.add_option("control.mass_bal", p.control.mass_bal)->default_val(false);
    app.add_option("control.dyn_damping", p.control.dyn_damping)->default_val(true);
    app.add_option("control.dyn_factor", p.control.dyn_factor)->default_val(0.8);
    app.add_option("control.max_vbc_val", p.control.max_vbc_val)->default_val(3.1709791983764588e-12);
    
    // ========== [mesh] section ==========
    app.add_option("mesh.mesh_file", p.mesh.mesh_file)->default_val("default");
    app.add_option("mesh.rs_levels", p.mesh.rs_levels)->default_val(2);
    app.add_option("mesh.rp_levels", p.mesh.rp_levels)->default_val(0);
    app.add_option("mesh.partition_type", p.mesh.partition_type)->default_val(0);
    app.add_option("mesh.order_v", p.mesh.order_v)->default_val(2);
    app.add_option("mesh.order_e", p.mesh.order_e)->default_val(1);
    app.add_option("mesh.order_q", p.mesh.order_q)->default_val(-1);
    app.add_option("mesh.local_refinement", p.mesh.local_refinement)->default_val(false);
    app.add_option("mesh.l2_basis", p.mesh.l2_basis)->default_val(1);
    
    // ========== [bc] section ==========
    app.add_option("bc.vbc_unit", p.bc.vbc_unit)->default_val("m/s");
    app.add_option("bc.vbc_factor", p.bc.vbc_factor)->default_val(1.0);
    app.add_option("bc.vbc_x0", p.bc.vbc_x0)->default_val(1);
    app.add_option("bc.vbc_x0_val0", p.bc.vbc_x0_val0)->default_val(0.0);
    app.add_option("bc.vbc_x0_val1", p.bc.vbc_x0_val1)->default_val(0.0);
    app.add_option("bc.vbc_x0_val2", p.bc.vbc_x0_val2)->default_val(0.0);
    app.add_option("bc.vbc_x1", p.bc.vbc_x1)->default_val(1);
    app.add_option("bc.vbc_x1_val0", p.bc.vbc_x1_val0)->default_val(0.0);
    app.add_option("bc.vbc_x1_val1", p.bc.vbc_x1_val1)->default_val(0.0);
    app.add_option("bc.vbc_x1_val2", p.bc.vbc_x1_val2)->default_val(0.0);
    app.add_option("bc.vbc_z0", p.bc.vbc_z0)->default_val(1);
    app.add_option("bc.vbc_z0_val0", p.bc.vbc_z0_val0)->default_val(0.0);
    app.add_option("bc.vbc_z0_val1", p.bc.vbc_z0_val1)->default_val(0.0);
    app.add_option("bc.vbc_z0_val2", p.bc.vbc_z0_val2)->default_val(0.0);
    app.add_option("bc.vbc_z1", p.bc.vbc_z1)->default_val(1);
    app.add_option("bc.vbc_z1_val0", p.bc.vbc_z1_val0)->default_val(0.0);
    app.add_option("bc.vbc_z1_val1", p.bc.vbc_z1_val1)->default_val(0.0);
    app.add_option("bc.vbc_z1_val2", p.bc.vbc_z1_val2)->default_val(0.0);
    app.add_option("bc.vbc_y0", p.bc.vbc_y0)->default_val(1);
    app.add_option("bc.vbc_y0_val0", p.bc.vbc_y0_val0)->default_val(0.0);
    app.add_option("bc.vbc_y0_val1", p.bc.vbc_y0_val1)->default_val(0.0);
    app.add_option("bc.vbc_y0_val2", p.bc.vbc_y0_val2)->default_val(0.0);
    app.add_option("bc.vbc_y1", p.bc.vbc_y1)->default_val(1);
    app.add_option("bc.vbc_y1_val0", p.bc.vbc_y1_val0)->default_val(0.0);
    app.add_option("bc.vbc_y1_val1", p.bc.vbc_y1_val1)->default_val(0.0);
    app.add_option("bc.vbc_y1_val2", p.bc.vbc_y1_val2)->default_val(0.0);
    app.add_option("bc.winkler_foundation", p.bc.winkler_foundation)->default_val(true);
    app.add_option("bc.winkler_flat", p.bc.winkler_flat)->default_val(false);
    app.add_option("bc.winkler_rho", p.bc.winkler_rho)->default_val(3200.0);
    app.add_option("bc.surf_proc", p.bc.surf_proc)->default_val(false);
    app.add_option("bc.surf_diff", p.bc.surf_diff)->default_val(1.0e-6);
    app.add_option("bc.surf_alpha", p.bc.surf_alpha)->default_val(0.0);
    app.add_option("bc.base_proc", p.bc.base_proc)->default_val(false);
    app.add_option("bc.base_diff", p.bc.base_diff)->default_val(1.0e-6);
    app.add_option("bc.base_alpha", p.bc.base_alpha)->default_val(0.0);
    
    // ========== [mat] section ==========
    app.add_option("mat.plastic", p.mat.plastic)->default_val(true);
    app.add_option("mat.viscoplastic", p.mat.viscoplastic)->default_val(true);
    app.add_option("mat.nmat", p.mat.nmat)->default_val(1);
    // Material arrays read as strings and parsed later
    app.add_option("mat.rho", rho_str)->default_val("[2700.0]");
    app.add_option("mat.lambda", lambda_str)->default_val("[3.0e10]");
    app.add_option("mat.mu", mu_str)->default_val("[3.0e10]");
    app.add_option("mat.tension_cutoff", tension_cutoff_str)->default_val("[0.0]");
    app.add_option("mat.cohesion0", cohesion0_str)->default_val("[4.0e7]");
    app.add_option("mat.cohesion1", cohesion1_str)->default_val("[4.0e6]");
    app.add_option("mat.friction_angle0", friction_angle0_str)->default_val("[30.0]");
    app.add_option("mat.friction_angle1", friction_angle1_str)->default_val("[15.0]");
    app.add_option("mat.dilation_angle0", dilation_angle0_str)->default_val("[0.0]");
    app.add_option("mat.dilation_angle1", dilation_angle1_str)->default_val("[0.0]");
    app.add_option("mat.alpha0", alpha0_str)->default_val("[0.0]");
    app.add_option("mat.alpha1", alpha1_str)->default_val("[0.5]");
    app.add_option("mat.plastic_viscosity", plastic_viscosity_str)->default_val("[1.0]");
    app.add_option("mat.weak_rad", p.mat.weak_rad)->default_val(1.0e3);
    app.add_option("mat.weak_x", p.mat.weak_x)->default_val(50.0e3);
    app.add_option("mat.weak_y", p.mat.weak_y)->default_val(2.00e3);
    app.add_option("mat.weak_z", p.mat.weak_z)->default_val(0.00e3);
    app.add_option("mat.ini_alpha", p.mat.ini_alpha)->default_val(0.5);
    
    // ========== [tmop] section ==========
    app.add_option("tmop.tmop", p.tmop.tmop)->default_val(false);
    app.add_option("tmop.amr", p.tmop.amr)->default_val(false);
    app.add_option("tmop.ale", p.tmop.ale)->default_val(1.0);
    app.add_option("tmop.remesh_steps", p.tmop.remesh_steps)->default_val(50000);
    app.add_option("tmop.mesh_poly_deg", p.tmop.mesh_poly_deg)->default_val(2);
    app.add_option("tmop.jitter", p.tmop.jitter)->default_val(0.0);
    app.add_option("tmop.metric_id", p.tmop.metric_id)->default_val(2);
    app.add_option("tmop.target_id", p.tmop.target_id)->default_val(1);
    app.add_option("tmop.lim_const", p.tmop.lim_const)->default_val(0.0);
    app.add_option("tmop.adapt_lim_const", p.tmop.adapt_lim_const)->default_val(0.0);
    app.add_option("tmop.quad_type", p.tmop.quad_type)->default_val(1);
    app.add_option("tmop.quad_order", p.tmop.quad_order)->default_val(8);
    app.add_option("tmop.solver_type", p.tmop.solver_type)->default_val(0);
    app.add_option("tmop.solver_iter", p.tmop.solver_iter)->default_val(20);
    app.add_option("tmop.solver_rtol", p.tmop.solver_rtol)->default_val(1e-10);
    app.add_option("tmop.solver_art_type", p.tmop.solver_art_type)->default_val(0);
    app.add_option("tmop.lin_solver", p.tmop.lin_solver)->default_val(2);
    app.add_option("tmop.max_lin_iter", p.tmop.max_lin_iter)->default_val(100);
    app.add_option("tmop.move_bnd", p.tmop.move_bnd)->default_val(false);
    app.add_option("tmop.combomet", p.tmop.combomet)->default_val(0);
    app.add_option("tmop.bal_expl_combo", p.tmop.bal_expl_combo)->default_val(false);
    app.add_option("tmop.hradaptivity", p.tmop.hradaptivity)->default_val(false);
    app.add_option("tmop.h_metric_id", p.tmop.h_metric_id)->default_val(-1);
    app.add_option("tmop.normalization", p.tmop.normalization)->default_val(false);
    app.add_option("tmop.verbosity_level", p.tmop.verbosity_level)->default_val(0);
    app.add_option("tmop.fdscheme", p.tmop.fdscheme)->default_val(false);
    app.add_option("tmop.adapt_eval", p.tmop.adapt_eval)->default_val(0);
    app.add_option("tmop.exactaction", p.tmop.exactaction)->default_val(false);
    app.add_option("tmop.n_hr_iter", p.tmop.n_hr_iter)->default_val(5);
    app.add_option("tmop.n_h_iter", p.tmop.n_h_iter)->default_val(1);
    app.add_option("tmop.mesh_node_ordering", p.tmop.mesh_node_ordering)->default_val(0);
    app.add_option("tmop.barrier_type", p.tmop.barrier_type)->default_val(0);
    app.add_option("tmop.worst_case_type", p.tmop.worst_case_type)->default_val(0);
    app.add_option("tmop.tmop_cond_num", p.tmop.tmop_cond_num)->default_val(0.5);
    
    // Allow extra keys in config file (for forward compatibility)
    app.allow_config_extras(CLI::config_extras_mode::ignore);
    
    // Parse config file
    try {
        app.set_config("--config", filename, "Configuration file", true);
        app.parse("--config " + std::string(filename));
    } catch (const CLI::ParseError &e) {
        std::cerr << "Error parsing config file '" << filename << "': " << e.what() << std::endl;
        std::exit(1);
    }
    
    // Parse material arrays from strings
    if (p.mat.nmat < 1) {
        std::cerr << "Error: mat.nmat must be greater than 0.\n";
        std::exit(1);
    }
    
    parse_array_string<double>(rho_str, "mat.rho", p.mat.rho, p.mat.nmat);
    parse_array_string<double>(lambda_str, "mat.lambda", p.mat.lambda, p.mat.nmat);
    parse_array_string<double>(mu_str, "mat.mu", p.mat.mu, p.mat.nmat);
    parse_array_string<double>(tension_cutoff_str, "mat.tension_cutoff", p.mat.tension_cutoff, p.mat.nmat);
    parse_array_string<double>(cohesion0_str, "mat.cohesion0", p.mat.cohesion0, p.mat.nmat);
    parse_array_string<double>(cohesion1_str, "mat.cohesion1", p.mat.cohesion1, p.mat.nmat);
    parse_array_string<double>(friction_angle0_str, "mat.friction_angle0", p.mat.friction_angle0, p.mat.nmat);
    parse_array_string<double>(friction_angle1_str, "mat.friction_angle1", p.mat.friction_angle1, p.mat.nmat);
    parse_array_string<double>(dilation_angle0_str, "mat.dilation_angle0", p.mat.dilation_angle0, p.mat.nmat);
    parse_array_string<double>(dilation_angle1_str, "mat.dilation_angle1", p.mat.dilation_angle1, p.mat.nmat);
    parse_array_string<double>(alpha0_str, "mat.alpha0", p.mat.alpha0, p.mat.nmat);
    parse_array_string<double>(alpha1_str, "mat.alpha1", p.mat.alpha1, p.mat.nmat);
    parse_array_string<double>(plastic_viscosity_str, "mat.plastic_viscosity", p.mat.plastic_viscosity, p.mat.nmat);
    
    validate_parameters(p);
}

// ============================================================================
// Validate and post-process parameters
// ============================================================================
static void validate_parameters(Param &p)
{
    std::cout << "Checking consistency of input parameters...\n";

    // Convert velocity BC units
    switch (bc_unit_map[p.bc.vbc_unit]) {
        case 1:  // cm/yr to m/s
            p.bc.vbc_factor = 1.0e-2 / YEAR2SEC;
            break;
        case 2:  // mm/yr to m/s
            p.bc.vbc_factor = 1.0e-3 / YEAR2SEC;
            break;
        case 3:  // cm/s to m/s
            p.bc.vbc_factor = 1.0e-2;
            break;
        default:  // already m/s
            break;
    }
}
