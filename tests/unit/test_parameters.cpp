#include <gtest/gtest.h>
#include "../../src/io/laghost_parameters.hpp"
#include "../../src/core/laghost_constants.hpp"

// Test parameter structures initialization and validation
class ParametersTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize default parameters
        param = Param{};
    }
    
    Param param;
};

TEST_F(ParametersTest, DefaultSimParametersInitialization) {
    // Test default Sim structure initialization
    Sim sim = {};
    
    EXPECT_EQ(sim.problem, 0);
    EXPECT_EQ(sim.dim, 0);
    EXPECT_EQ(sim.t_final, 0.0);
    EXPECT_EQ(sim.max_tsteps, 0);
    EXPECT_FALSE(sim.year);
    EXPECT_FALSE(sim.visualization);
    EXPECT_EQ(sim.vis_steps, 0);
    EXPECT_FALSE(sim.visit);
    EXPECT_FALSE(sim.paraview);
    EXPECT_FALSE(sim.gfprint);
    EXPECT_TRUE(sim.basename.empty());
    EXPECT_TRUE(sim.device.empty());
    EXPECT_EQ(sim.dev, 0);
    EXPECT_FALSE(sim.check);
    EXPECT_FALSE(sim.mem_usage);
    EXPECT_FALSE(sim.fom);
    EXPECT_FALSE(sim.gpu_aware_mpi);
}

TEST_F(ParametersTest, DefaultSolverParametersInitialization) {
    SolverParams solver = {};
    
    EXPECT_EQ(solver.ode_solver_type, 0);
    EXPECT_EQ(solver.cfl, 0.0);
    EXPECT_EQ(solver.cg_tol, 0.0);
    EXPECT_EQ(solver.ftz_tol, 0.0);
    EXPECT_EQ(solver.cg_max_iter, 0);
    EXPECT_FALSE(solver.p_assembly);
    EXPECT_FALSE(solver.impose_visc);
}

TEST_F(ParametersTest, DefaultControlParametersInitialization) {
    Control control = {};
    
    EXPECT_FALSE(control.pseudo_transient);
    EXPECT_EQ(control.transient_num, 0);
    EXPECT_FALSE(control.lithostatic);
    EXPECT_EQ(control.init_dt, 0.0);
    EXPECT_EQ(control.mscale, 0.0);
    EXPECT_EQ(control.gravity, 0.0);
    EXPECT_EQ(control.thickness, 0.0);
    EXPECT_FALSE(control.mass_bal);
    EXPECT_FALSE(control.dyn_damping);
    EXPECT_EQ(control.dyn_factor, 0.0);
    EXPECT_EQ(control.max_vbc_val, 0.0);
}

TEST_F(ParametersTest, DefaultMeshingParametersInitialization) {
    Meshing mesh = {};
    
    EXPECT_TRUE(mesh.mesh_file.empty());
    EXPECT_EQ(mesh.rs_levels, 0);
    EXPECT_EQ(mesh.rp_levels, 0);
    EXPECT_EQ(mesh.partition_type, 0);
    EXPECT_EQ(mesh.order_v, 0);
    EXPECT_EQ(mesh.order_e, 0);
    EXPECT_EQ(mesh.order_q, 0);
    EXPECT_FALSE(mesh.local_refinement);
    EXPECT_EQ(mesh.l2_basis, 0);
}

TEST_F(ParametersTest, DefaultBCParametersInitialization) {
    BC bc = {};
    
    EXPECT_TRUE(bc.bc_unit.empty());
    EXPECT_TRUE(bc.bc_ids.empty());
    EXPECT_TRUE(bc.bc_vxs.empty());
    EXPECT_TRUE(bc.bc_vys.empty());
    EXPECT_TRUE(bc.bc_vzs.empty());
    EXPECT_EQ(bc.vel_unit, 0.0);
    EXPECT_FALSE(bc.winkler_foundation);
    EXPECT_FALSE(bc.winkler_flat);
    EXPECT_EQ(bc.winkler_rho, 0.0);
    EXPECT_FALSE(bc.surf_proc);
    EXPECT_EQ(bc.surf_diff, 0.0);
    EXPECT_EQ(bc.surf_alpha, 0.0);
    EXPECT_FALSE(bc.base_proc);
    EXPECT_EQ(bc.base_diff, 0.0);
    EXPECT_EQ(bc.base_alpha, 0.0);
}

TEST_F(ParametersTest, DefaultMatParametersInitialization) {
    Mat mat = {};
    
    EXPECT_FALSE(mat.plastic);
    EXPECT_FALSE(mat.viscoplastic);
    EXPECT_TRUE(mat.rho.empty());
    EXPECT_TRUE(mat.lambda.empty());
    EXPECT_TRUE(mat.mu.empty());
    EXPECT_TRUE(mat.tension_cutoff.empty());
    EXPECT_TRUE(mat.cohesion0.empty());
    EXPECT_TRUE(mat.cohesion1.empty());
    EXPECT_TRUE(mat.friction_angle0.empty());
    EXPECT_TRUE(mat.friction_angle1.empty());
    EXPECT_TRUE(mat.dilation_angle0.empty());
    EXPECT_TRUE(mat.dilation_angle1.empty());
    EXPECT_TRUE(mat.pls0.empty());
    EXPECT_TRUE(mat.pls1.empty());
    EXPECT_TRUE(mat.plastic_viscosity.empty());
    EXPECT_EQ(mat.weak_rad, 0.0);
    EXPECT_EQ(mat.weak_x, 0.0);
    EXPECT_EQ(mat.weak_y, 0.0);
    EXPECT_EQ(mat.weak_z, 0.0);
    EXPECT_EQ(mat.ini_pls, 0.0);
}

TEST_F(ParametersTest, DefaultTMOPParametersInitialization) {
    TMOP tmop = {};
    
    EXPECT_FALSE(tmop.tmop);
    EXPECT_FALSE(tmop.amr);
    EXPECT_EQ(tmop.ale, 0.0);
    EXPECT_EQ(tmop.remesh_steps, 0);
    EXPECT_EQ(tmop.mesh_poly_deg, 0);
    EXPECT_EQ(tmop.jitter, 0.0);
    EXPECT_EQ(tmop.metric_id, 0);
    EXPECT_EQ(tmop.target_id, 0);
    EXPECT_EQ(tmop.lim_const, 0.0);
    EXPECT_EQ(tmop.adapt_lim_const, 0.0);
    EXPECT_EQ(tmop.quad_type, 0);
    EXPECT_EQ(tmop.quad_order, 0);
    EXPECT_EQ(tmop.solver_type, 0);
    EXPECT_EQ(tmop.solver_iter, 0);
    EXPECT_EQ(tmop.solver_rtol, 0.0);
    EXPECT_EQ(tmop.solver_art_type, 0);
    EXPECT_EQ(tmop.lin_solver, 0);
    EXPECT_EQ(tmop.max_lin_iter, 0);
    EXPECT_FALSE(tmop.move_bnd);
    EXPECT_EQ(tmop.combomet, 0);
    EXPECT_FALSE(tmop.bal_expl_combo);
    EXPECT_FALSE(tmop.hradaptivity);
    EXPECT_EQ(tmop.h_metric_id, 0);
    EXPECT_FALSE(tmop.normalization);
    EXPECT_EQ(tmop.verbosity_level, 0);
    EXPECT_FALSE(tmop.fdscheme);
    EXPECT_EQ(tmop.adapt_eval, 0);
    EXPECT_FALSE(tmop.exactaction);
    EXPECT_EQ(tmop.n_hr_iter, 0);
    EXPECT_EQ(tmop.n_h_iter, 0);
    EXPECT_EQ(tmop.mesh_node_ordering, 0);
    EXPECT_EQ(tmop.barrier_type, 0);
    EXPECT_EQ(tmop.worst_case_type, 0);
    EXPECT_EQ(tmop.tmop_cond_num, 0.0);
}

TEST_F(ParametersTest, ConstantsValidation) {
    // Test mathematical constants
    EXPECT_NEAR(M_PI, 3.1415926535897932384626433832795, 1e-15);
    EXPECT_NEAR(YEAR2SEC, 365.2422 * 86400, 1e-6);
    EXPECT_NEAR(DEG2RAD, M_PI / 180, 1e-15);
    
    // Test dimension-dependent constants
    #ifdef THREED
    EXPECT_EQ(NDIMS, 3);
    EXPECT_EQ(NODES_PER_ELEM, 4);
    EXPECT_EQ(NSTR, 6);
    EXPECT_EQ(FACETS_PER_ELEM, 4);
    EXPECT_EQ(NODES_PER_FACET, 3);
    #else
    EXPECT_EQ(NDIMS, 2);
    EXPECT_EQ(NODES_PER_ELEM, 3);
    EXPECT_EQ(NSTR, 3);
    EXPECT_EQ(FACETS_PER_ELEM, 3);
    EXPECT_EQ(NODES_PER_FACET, 2);
    #endif
}

TEST_F(ParametersTest, BoundaryFlagsValidation) {
    // Test boundary flag constants
    EXPECT_EQ(BOUNDX0, 1);
    EXPECT_EQ(BOUNDX1, 2);
    EXPECT_EQ(BOUNDY0, 4);
    EXPECT_EQ(BOUNDY1, 8);
    EXPECT_EQ(BOUNDZ0, 16);
    EXPECT_EQ(BOUNDZ1, 32);
    EXPECT_EQ(BOUNDN0, 64);
    EXPECT_EQ(BOUNDN1, 128);
    EXPECT_EQ(BOUNDN2, 256);
    EXPECT_EQ(BOUNDN3, 512);
    
    // Test combined boundary flag
    uint expected_bound_any = BOUNDX0 | BOUNDX1 | BOUNDY0 | BOUNDY1 | 
                             BOUNDZ0 | BOUNDZ1 | BOUNDN0 | BOUNDN1 | 
                             BOUNDN2 | BOUNDN3;
    EXPECT_EQ(BOUND_ANY, expected_bound_any);
}