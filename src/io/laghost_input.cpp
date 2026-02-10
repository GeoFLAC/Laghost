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
static void parse_material_arrays_from_file(const char* filename, Param &p);

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

// Helper to extract the value string from a line like "key = value  # comment"
static std::string extract_value(const std::string &line) {
    auto eq_pos = line.find('=');
    if (eq_pos == std::string::npos) return "";
    std::string val = line.substr(eq_pos + 1);
    // Strip inline comment (but not inside brackets)
    auto hash_pos = val.find('#');
    if (hash_pos != std::string::npos) {
        val = val.substr(0, hash_pos);
    }
    // Trim whitespace
    auto start = val.find_first_not_of(" \t");
    auto end = val.find_last_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    return val.substr(start, end - start + 1);
}

// Parse material arrays directly from config file, bypassing CLI11's comma-splitting
static void parse_material_arrays_from_file(const char* filename, Param &p)
{
    std::ifstream f(filename);
    if (!f.good()) {
        std::cerr << "Error: Cannot re-open config file '" << filename << "' for array parsing\n";
        std::exit(1);
    }
    
    // Map of property name -> pointer to the mfem::Vector in Param
    struct ArrayProp { const char* name; Vector* vec; std::string default_val; };
    std::vector<ArrayProp> props = {
        {"rho",               &p.mat.rho,               "[2700.0]"},
        {"lambda",            &p.mat.lambda,             "[3.0e10]"},
        {"mu",                &p.mat.mu,                 "[3.0e10]"},
        {"tension_cutoff",    &p.mat.tension_cutoff,     "[0.0]"},
        {"cohesion0",         &p.mat.cohesion0,          "[4.0e7]"},
        {"cohesion1",         &p.mat.cohesion1,          "[4.0e6]"},
        {"friction_angle0",   &p.mat.friction_angle0,    "[30.0]"},
        {"friction_angle1",   &p.mat.friction_angle1,    "[15.0]"},
        {"dilation_angle0",   &p.mat.dilation_angle0,    "[0.0]"},
        {"dilation_angle1",   &p.mat.dilation_angle1,    "[0.0]"},
        {"alpha0",            &p.mat.alpha0,             "[0.0]"},
        {"alpha1",            &p.mat.alpha1,             "[0.5]"},
        {"plastic_viscosity", &p.mat.plastic_viscosity,  "[1.0]"},
    };
    
    // Track which properties we found
    std::map<std::string, std::string> found_values;
    
    std::string line;
    bool in_mat_section = false;
    while (std::getline(f, line)) {
        // Trim leading whitespace
        auto first_char = line.find_first_not_of(" \t");
        if (first_char == std::string::npos) continue;
        
        // Check for section headers
        if (line[first_char] == '[') {
            auto end_bracket = line.find(']', first_char);
            if (end_bracket != std::string::npos) {
                std::string section = line.substr(first_char + 1, end_bracket - first_char - 1);
                in_mat_section = (section == "mat");
            }
            continue;
        }
        
        if (!in_mat_section) continue;
        if (line[first_char] == '#' || line[first_char] == ';') continue;
        
        // Extract key
        auto eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;
        std::string key = line.substr(first_char, eq_pos - first_char);
        // Trim trailing whitespace from key
        auto key_end = key.find_last_not_of(" \t");
        if (key_end != std::string::npos) key = key.substr(0, key_end + 1);
        
        // Check if this is one of our array properties
        for (auto& prop : props) {
            if (key == prop.name) {
                found_values[prop.name] = extract_value(line);
                break;
            }
        }
    }
    
    // Parse each property (use default if not found in file)
    for (auto& prop : props) {
        auto it = found_values.find(prop.name);
        std::string val_str = (it != found_values.end()) ? it->second : prop.default_val;
        parse_array_string<double>(val_str, prop.name, *prop.vec, p.mat.nmat);
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

    // Parse command-line first to handle -h/--help and validate args
    args.Parse();
    
    if (!args.Good())
    {
        if (myid == 0) args.PrintUsage(std::cout);
        MPI_Finalize();
        exit(0);
    }
    
    // Now load configuration file - this OVERWRITES the MFEM defaults
    get_input_parameters(input_parameter_file, param);
    
    // Print MFEM options (these show command-line args, not config values)
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
    
    // ========== [sim] section ==========
    auto sim = app.add_subcommand("sim");
    sim->allow_config_extras(CLI::config_extras_mode::ignore);
    sim->add_option("problem", p.sim.problem)->default_val(1);
    sim->add_option("dim", p.sim.dim)->default_val(3);
    sim->add_option("t_final", p.sim.t_final)->default_val(1.0);
    sim->add_option("max_tsteps", p.sim.max_tsteps)->default_val(-1);
    sim->add_option("year", p.sim.year)->default_val(false);
    sim->add_option("visualization", p.sim.visualization)->default_val(false);
    sim->add_option("vis_steps", p.sim.vis_steps)->default_val(1000);
    sim->add_option("visit", p.sim.visit)->default_val(false);
    sim->add_option("paraview", p.sim.paraview)->default_val(true);
    sim->add_option("gfprint", p.sim.gfprint)->default_val(false);
    sim->add_option("basename", p.sim.basename)->default_val("results/Laghost");
    sim->add_option("device", p.sim.device)->default_val("cpu");
    sim->add_option("dev", p.sim.dev)->default_val(0);
    sim->add_option("check", p.sim.check)->default_val(false);
    sim->add_option("mem_usage", p.sim.mem_usage)->default_val(false);
    sim->add_option("fom", p.sim.fom)->default_val(false);
    sim->add_option("gpu_aware_mpi", p.sim.gpu_aware_mpi)->default_val(false);
    
    // ========== [solver] section ==========
    auto solver = app.add_subcommand("solver");
    solver->allow_config_extras(CLI::config_extras_mode::ignore);
    solver->add_option("ode_solver_type", p.solver.ode_solver_type)->default_val(7);
    solver->add_option("cfl", p.solver.cfl)->default_val(0.5);
    solver->add_option("cg_tol", p.solver.cg_tol)->default_val(1.0e-10);
    solver->add_option("ftz_tol", p.solver.ftz_tol)->default_val(0.0);
    solver->add_option("cg_max_iter", p.solver.cg_max_iter)->default_val(300);
    solver->add_option("p_assembly", p.solver.p_assembly)->default_val(false);
    solver->add_option("impose_visc", p.solver.impose_visc)->default_val(true);
    
    // ========== [control] section ==========
    auto control = app.add_subcommand("control");
    control->allow_config_extras(CLI::config_extras_mode::ignore);
    control->add_option("lithostatic", p.control.lithostatic)->default_val(true);
    control->add_option("atmospheric", p.control.atmospheric)->default_val(false);
    control->add_option("init_dt", p.control.init_dt)->default_val(1.0);
    control->add_option("mscale", p.control.mscale)->default_val(5.0e5);
    control->add_option("gravity", p.control.gravity)->default_val(10.0);
    control->add_option("thickness", p.control.thickness)->default_val(10.0e3);
    control->add_option("mass_bal", p.control.mass_bal)->default_val(false);
    control->add_option("dyn_damping", p.control.dyn_damping)->default_val(true);
    control->add_option("dyn_factor", p.control.dyn_factor)->default_val(0.8);
    control->add_option("max_vbc_val", p.control.max_vbc_val)->default_val(3.1709791983764588e-12);
    
    // ========== [mesh] section ==========
    auto mesh = app.add_subcommand("mesh");
    mesh->allow_config_extras(CLI::config_extras_mode::ignore);
    mesh->add_option("mesh_file", p.mesh.mesh_file)->default_val("default");
    mesh->add_option("rs_levels", p.mesh.rs_levels)->default_val(2);
    mesh->add_option("rp_levels", p.mesh.rp_levels)->default_val(0);
    mesh->add_option("partition_type", p.mesh.partition_type)->default_val(0);
    mesh->add_option("order_v", p.mesh.order_v)->default_val(2);
    mesh->add_option("order_e", p.mesh.order_e)->default_val(1);
    mesh->add_option("order_q", p.mesh.order_q)->default_val(-1);
    mesh->add_option("local_refinement", p.mesh.local_refinement)->default_val(false);
    mesh->add_option("l2_basis", p.mesh.l2_basis)->default_val(1);
    
    // ========== [bc] section ==========
    auto bc = app.add_subcommand("bc");
    bc->allow_config_extras(CLI::config_extras_mode::ignore);
    bc->add_option("vbc_unit", p.bc.vbc_unit)->default_val("m/s");
    bc->add_option("vbc_factor", p.bc.vbc_factor)->default_val(1.0);
    bc->add_option("vbc_x0", p.bc.vbc_x0)->default_val(1);
    bc->add_option("vbc_x0_val0", p.bc.vbc_x0_val0)->default_val(0.0);
    bc->add_option("vbc_x0_val1", p.bc.vbc_x0_val1)->default_val(0.0);
    bc->add_option("vbc_x0_val2", p.bc.vbc_x0_val2)->default_val(0.0);
    bc->add_option("vbc_x1", p.bc.vbc_x1)->default_val(1);
    bc->add_option("vbc_x1_val0", p.bc.vbc_x1_val0)->default_val(0.0);
    bc->add_option("vbc_x1_val1", p.bc.vbc_x1_val1)->default_val(0.0);
    bc->add_option("vbc_x1_val2", p.bc.vbc_x1_val2)->default_val(0.0);
    bc->add_option("vbc_z0", p.bc.vbc_z0)->default_val(1);
    bc->add_option("vbc_z0_val0", p.bc.vbc_z0_val0)->default_val(0.0);
    bc->add_option("vbc_z0_val1", p.bc.vbc_z0_val1)->default_val(0.0);
    bc->add_option("vbc_z0_val2", p.bc.vbc_z0_val2)->default_val(0.0);
    bc->add_option("vbc_z1", p.bc.vbc_z1)->default_val(1);
    bc->add_option("vbc_z1_val0", p.bc.vbc_z1_val0)->default_val(0.0);
    bc->add_option("vbc_z1_val1", p.bc.vbc_z1_val1)->default_val(0.0);
    bc->add_option("vbc_z1_val2", p.bc.vbc_z1_val2)->default_val(0.0);
    bc->add_option("vbc_y0", p.bc.vbc_y0)->default_val(1);
    bc->add_option("vbc_y0_val0", p.bc.vbc_y0_val0)->default_val(0.0);
    bc->add_option("vbc_y0_val1", p.bc.vbc_y0_val1)->default_val(0.0);
    bc->add_option("vbc_y0_val2", p.bc.vbc_y0_val2)->default_val(0.0);
    bc->add_option("vbc_y1", p.bc.vbc_y1)->default_val(1);
    bc->add_option("vbc_y1_val0", p.bc.vbc_y1_val0)->default_val(0.0);
    bc->add_option("vbc_y1_val1", p.bc.vbc_y1_val1)->default_val(0.0);
    bc->add_option("vbc_y1_val2", p.bc.vbc_y1_val2)->default_val(0.0);
    bc->add_option("winkler_foundation", p.bc.winkler_foundation)->default_val(true);
    bc->add_option("winkler_flat", p.bc.winkler_flat)->default_val(false);
    bc->add_option("winkler_rho", p.bc.winkler_rho)->default_val(3200.0);
    bc->add_option("surf_proc", p.bc.surf_proc)->default_val(false);
    bc->add_option("surf_diff", p.bc.surf_diff)->default_val(1.0e-6);
    bc->add_option("surf_alpha", p.bc.surf_alpha)->default_val(0.0);
    bc->add_option("base_proc", p.bc.base_proc)->default_val(false);
    bc->add_option("base_diff", p.bc.base_diff)->default_val(1.0e-6);
    bc->add_option("base_alpha", p.bc.base_alpha)->default_val(0.0);
    
    // ========== [mat] section ==========
    auto mat = app.add_subcommand("mat");
    mat->allow_config_extras(CLI::config_extras_mode::ignore);
    mat->add_option("plastic", p.mat.plastic)->default_val(true);
    mat->add_option("viscoplastic", p.mat.viscoplastic)->default_val(true);
    mat->add_option("nmat", p.mat.nmat)->default_val(1);
    // Material arrays are NOT registered with CLI11 (comma-splitting breaks them).
    // They are parsed manually from the config file after CLI11 processes scalars.
    mat->add_option("weak_rad", p.mat.weak_rad)->default_val(1.0e3);
    mat->add_option("weak_x", p.mat.weak_x)->default_val(50.0e3);
    mat->add_option("weak_y", p.mat.weak_y)->default_val(2.00e3);
    mat->add_option("weak_z", p.mat.weak_z)->default_val(0.00e3);
    mat->add_option("ini_alpha", p.mat.ini_alpha)->default_val(0.5);
    
    // ========== [tmop] section ==========
    auto tmop = app.add_subcommand("tmop");
    tmop->allow_config_extras(CLI::config_extras_mode::ignore);
    tmop->add_option("tmop", p.tmop.tmop)->default_val(false);
    tmop->add_option("amr", p.tmop.amr)->default_val(false);
    tmop->add_option("ale", p.tmop.ale)->default_val(1.0);
    tmop->add_option("remesh_steps", p.tmop.remesh_steps)->default_val(50000);
    tmop->add_option("mesh_poly_deg", p.tmop.mesh_poly_deg)->default_val(2);
    tmop->add_option("jitter", p.tmop.jitter)->default_val(0.0);
    tmop->add_option("metric_id", p.tmop.metric_id)->default_val(2);
    tmop->add_option("target_id", p.tmop.target_id)->default_val(1);
    tmop->add_option("lim_const", p.tmop.lim_const)->default_val(0.0);
    tmop->add_option("adapt_lim_const", p.tmop.adapt_lim_const)->default_val(0.0);
    tmop->add_option("quad_type", p.tmop.quad_type)->default_val(1);
    tmop->add_option("quad_order", p.tmop.quad_order)->default_val(8);
    tmop->add_option("solver_type", p.tmop.solver_type)->default_val(0);
    tmop->add_option("solver_iter", p.tmop.solver_iter)->default_val(20);
    tmop->add_option("solver_rtol", p.tmop.solver_rtol)->default_val(1e-10);
    tmop->add_option("solver_art_type", p.tmop.solver_art_type)->default_val(0);
    tmop->add_option("lin_solver", p.tmop.lin_solver)->default_val(2);
    tmop->add_option("max_lin_iter", p.tmop.max_lin_iter)->default_val(100);
    tmop->add_option("move_bnd", p.tmop.move_bnd)->default_val(false);
    tmop->add_option("combomet", p.tmop.combomet)->default_val(0);
    tmop->add_option("bal_expl_combo", p.tmop.bal_expl_combo)->default_val(false);
    tmop->add_option("hradaptivity", p.tmop.hradaptivity)->default_val(false);
    tmop->add_option("h_metric_id", p.tmop.h_metric_id)->default_val(-1);
    tmop->add_option("normalization", p.tmop.normalization)->default_val(false);
    tmop->add_option("verbosity_level", p.tmop.verbosity_level)->default_val(0);
    tmop->add_option("fdscheme", p.tmop.fdscheme)->default_val(false);
    tmop->add_option("adapt_eval", p.tmop.adapt_eval)->default_val(0);
    tmop->add_option("exactaction", p.tmop.exactaction)->default_val(false);
    tmop->add_option("n_hr_iter", p.tmop.n_hr_iter)->default_val(5);
    tmop->add_option("n_h_iter", p.tmop.n_h_iter)->default_val(1);
    tmop->add_option("mesh_node_ordering", p.tmop.mesh_node_ordering)->default_val(0);
    tmop->add_option("barrier_type", p.tmop.barrier_type)->default_val(0);
    tmop->add_option("worst_case_type", p.tmop.worst_case_type)->default_val(0);
    tmop->add_option("tmop_cond_num", p.tmop.tmop_cond_num)->default_val(0.5);
    
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
    
    // Parse material arrays directly from config file (CLI11 can't handle bracket arrays)
    if (p.mat.nmat < 1) {
        std::cerr << "Error: mat.nmat must be greater than 0.\n";
        std::exit(1);
    }
    parse_material_arrays_from_file(filename, p);
    
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
