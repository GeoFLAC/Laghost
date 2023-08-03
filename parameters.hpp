#ifndef LAGHOST_PARAMETERS_HPP
#define LAGHOST_PARAMETERS_HPP

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "constants.hpp"
#include "array2d.hpp"

//
// Structures for input parameters
//
struct Sim {
    int         problem;
    int         dim;
    double      t_final;
    int         max_tsteps;
    bool        year;
    bool        visualization;
    int         vis_steps;
    bool        visit;
    bool        paraview;
    bool        gfprint;
    std::string basename; 
    std::string device;
    int         dev;
    bool        check;
    bool        mem_usage;
    bool        fom;
    bool        gpu_aware_mpi;
};

struct Solver {
    int    ode_solver_type;
    double cfl;
    double cg_tol;
    double ftz_tol;
    int    cg_max_iter;
    bool   p_assembly;
    bool   impose_visc;
};

struct Control {
    bool   winkler_foundation;
    bool   winkler_flat;
    bool   lithostatic;
    double init_dt;
    double mscale;
    double gravity; // magnitude 
    double thickness; // meter 
    double winkler_rho; // Density of substratum
};

struct Mesh_param {
    std::string mesh_file;
    int         rs_levels;
    int         rp_levels;
    int         partition_type;
    int         order_v;
    int         order_e;
    int         order_q;
    bool        local_refinement;
};

struct Mat {
    bool   plastic;
    bool   viscoplastic;
    double lambda;
    double mu;
    double weak_rad;
    double weak_x;
    double weak_y;
    double weak_z;
    double ini_pls;
    double tension_cutoff;
    double cohesion0;
    double cohesion1;
    double friction_angle;
    double dilation_angle;
    double pls0;
    double pls1;
    double plastic_viscosity;
};

struct TMOP {
    bool   tmop;
    bool   amr;
    int    remesh_steps;
    int    mesh_poly_deg;
    double jitter;
    int    metric_id;
    int    target_id;
    double lim_const;
    double adapt_lim_const;
    int    quad_type;
    int    quad_order;
    int    solver_type;
    int    solver_iter;
    double solver_rtol;
    int    solver_art_type;
    int    lin_solver;
    int    max_lin_iter;
    bool   move_bnd;
    int    combomet;
    bool   bal_expl_combo;
    bool   hradaptivity;
    int    h_metric_id;
    bool   normalization;
    int    verbosity_level;
    bool   fdscheme;
    int    adapt_eval;
    bool   exactaction;
    int    n_hr_iter;
    int    n_h_iter;
    int    mesh_node_ordering;
    int    barrier_type;
    int    worst_case_type;
};

struct Param {
    Sim sim;
    Solver solver;
    Mesh_param mesh;
    Control control;
    Mat mat;
    TMOP tmop;
};

#endif
