#ifndef MFEM_LAGHOST_PARAMETERS
#define MFEM_LAGHOST_PARAMETERS

#include "mfem.hpp"
#include <map>
#include <string>
#include "../core/laghost_constants.hpp"

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

struct SolverParams {
    int    ode_solver_type;
    double cfl;
    double cg_tol;
    double ftz_tol;
    int    cg_max_iter;
    bool   p_assembly;
    bool   impose_visc;
};

struct Control {
    // double flat_rate;
    bool   lithostatic;
    bool   atmospheric;
    double init_dt;
    double mscale;
    double gravity; // magnitude 
    double thickness; // meter 
    bool   mass_bal;
    bool   dyn_damping;
    double dyn_factor;
    double max_vbc_val;
};

struct Meshing {
    std::string mesh_file;
    int         rs_levels;
    int         rp_levels;
    int         partition_type;
    int         order_v;
    int         order_e;
    int         order_q;
    bool        local_refinement;
    int         l2_basis;
};

struct BC {
    // int_vec bc_ids;
    std::string vbc_unit;
    double vbc_factor;
    // std::vector<int>    bc_id;
    // std::vector<double> bc_vx;
    // std::vector<double> bc_vy;
    // std::vector<double> bc_vz;
    int vbc_x0;
    int vbc_x1;
    int vbc_z0;
    int vbc_z1;
    int vbc_y0;
    int vbc_y1;

    double vbc_x0_val0;
    double vbc_x0_val1;
    double vbc_x0_val2;
    double vbc_x1_val0;
    double vbc_x1_val1;
    double vbc_x1_val2;
    double vbc_z0_val0;
    double vbc_z0_val1;
    double vbc_z0_val2;
    double vbc_z1_val0;
    double vbc_z1_val1;
    double vbc_z1_val2;
    double vbc_y0_val0;
    double vbc_y0_val1;
    double vbc_y0_val2;
    double vbc_y1_val0;
    double vbc_y1_val1;
    double vbc_y1_val2;
    
    bool   winkler_foundation;
    bool   winkler_flat;
    double winkler_rho; // Density of substratum
    bool   surf_proc;
    double surf_diff;
    double surf_alpha;
    bool   base_proc;
    double base_diff;
    double base_alpha;
};

struct Mat {
    bool   plastic;
    bool   viscoplastic;
    int    nmat;
    mfem::Vector rho;
    mfem::Vector lambda;
    mfem::Vector mu;
    // std::string weak_rad;
    // std::string weak_x;
    // std::string weak_y;
    // std::string weak_z;
    // std::string ini_alpha;
    mfem::Vector tension_cutoff;
    mfem::Vector cohesion0;
    mfem::Vector cohesion1;
    mfem::Vector friction_angle0;
    mfem::Vector friction_angle1;
    mfem::Vector dilation_angle0;
    mfem::Vector dilation_angle1;
    mfem::Vector alpha0;  // Internal variable at start of weakening
    mfem::Vector alpha1;  // Internal variable at end of weakening
    mfem::Vector plastic_viscosity;

    // double lambda;
    // double mu;
    double weak_rad;
    double weak_x;
    double weak_y;
    double weak_z;
    double ini_alpha;  // Initial internal variable value
    // double tension_cutoff;
    // double cohesion0;
    // double cohesion1;
    // double friction_angle;
    // double dilation_angle;
    // double alpha0;
    // double alpha1;
    // double plastic_viscosity;
};

struct TMOP {
    bool   tmop;
    bool   amr;
    double ale;
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
    double tmop_cond_num;
};

struct Param {
    Sim sim;
    SolverParams solver;
    BC bc;
    Meshing mesh;
    Control control;
    Mat mat;
    TMOP tmop;
};
#endif // MFEM_LAGHOST_PARAMETERS
