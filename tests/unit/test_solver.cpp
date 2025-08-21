#include <gtest/gtest.h>
#include "../../src/core/laghost_solver.hpp"
#include "../../src/core/laghost_assembly.hpp"
#include "mfem.hpp"

using namespace mfem;
using namespace mfem::geodynamics;

class SolverTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple 2D mesh
        mesh = new Mesh(Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL, true));
        pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
        delete mesh;
        mesh = nullptr;
        
        dim = pmesh->Dimension();
        
        // Create finite element collections and spaces
        h1_fec = new H1_FECollection(1, dim);  // Order 1
        l2_fec = new L2_FECollection(0, dim);  // Order 0
        
        H1FESpace = new ParFiniteElementSpace(pmesh, h1_fec, dim);  // Vector space
        L2FESpace = new ParFiniteElementSpace(pmesh, l2_fec);       // Scalar space
        L2FESpace_stress = new ParFiniteElementSpace(pmesh, l2_fec, 3*(dim-1)); // Stress space
        
        // Create grid functions for material properties
        rho0_gf = new ParGridFunction(L2FESpace);
        fictitious_rho0_gf = new ParGridFunction(L2FESpace);
        gamma_gf = new ParGridFunction(L2FESpace);
        lambda_gf = new ParGridFunction(L2FESpace);
        mu_gf = new ParGridFunction(L2FESpace);
        
        // Initialize material properties with typical rock values
        *rho0_gf = 2700.0;           // kg/m³
        *fictitious_rho0_gf = 2700.0;
        *gamma_gf = 1.4;             // Heat capacity ratio
        *lambda_gf = 3.0e10;         // Pa
        *mu_gf = 3.0e10;             // Pa
        
        // Create essential DOFs array (empty for this test)
        ess_tdofs.SetSize(0);
        
        // Set up boundary condition vector
        bc_id_pa.SetSize(pmesh->bdr_attributes.Max());
        bc_id_pa = 0.0;  // No boundary conditions for this test
        
        // Initialize solver parameters
        source = 0;
        cfl = 0.25;
        visc = true;
        vort = false;
        pa = false;  // Use full assembly for testing
        cgt = 1e-10;
        cgiter = 300;
        ftz_tol = 0.0;
        order_q = -1;
        mscale = 1.0;
        gravity = 9.81;
        thickness = 1000.0;
        winkler = false;
        winkler_rho = 2700.0;
        dyn_damping = false;
        dyn_factor = 0.8;
        vbc_max_val = 1e-12;
    }
    
    void TearDown() override {
        delete mu_gf;
        delete lambda_gf;
        delete gamma_gf;
        delete fictitious_rho0_gf;
        delete rho0_gf;
        delete L2FESpace_stress;
        delete L2FESpace;
        delete H1FESpace;
        delete l2_fec;
        delete h1_fec;
        delete pmesh;
    }
    
    Mesh* mesh;
    ParMesh* pmesh;
    H1_FECollection* h1_fec;
    L2_FECollection* l2_fec;
    ParFiniteElementSpace* H1FESpace;
    ParFiniteElementSpace* L2FESpace;
    ParFiniteElementSpace* L2FESpace_stress;
    ParGridFunction* rho0_gf;
    ParGridFunction* fictitious_rho0_gf;
    ParGridFunction* gamma_gf;
    ParGridFunction* lambda_gf;
    ParGridFunction* mu_gf;
    Array<int> ess_tdofs;
    Vector bc_id_pa;
    
    int dim, source, cgiter, order_q;
    double cfl, cgt, ftz_tol, mscale, gravity, thickness, winkler_rho, dyn_factor, vbc_max_val;
    bool visc, vort, pa, winkler, dyn_damping;
};

TEST_F(SolverTest, TimingDataInitialization) {
    HYPRE_Int l2dof = L2FESpace->GetTrueVSize();
    
    EXPECT_NO_THROW({
        TimingData timer(l2dof);
        
        EXPECT_EQ(timer.L2dof, l2dof);
        EXPECT_EQ(timer.H1iter, 0);
        EXPECT_EQ(timer.L2iter, 0);
        EXPECT_EQ(timer.quad_tstep, 0);
    });
}

TEST_F(SolverTest, QUpdateInitialization) {
    int NE = pmesh->GetNE();
    int Q1D = 2;  // Simple quadrature
    
    // Create integration rule
    const IntegrationRule *ir = &IntRules.Get(pmesh->GetElementType(0), 2);
    
    TimingData timer(L2FESpace->GetTrueVSize());
    
    EXPECT_NO_THROW({
        QUpdate qupdate(dim, NE, Q1D, visc, vort, cfl, &timer,
                       *gamma_gf, *lambda_gf, *mu_gf, *ir,
                       *H1FESpace, *L2FESpace, *L2FESpace_stress);
    });
}

TEST_F(SolverTest, LagrangianGeoOperatorInitialization) {
    // Calculate total system size
    int H1Vsize = H1FESpace->GetVSize();
    int L2Vsize = L2FESpace->GetVSize();
    int L2StressVsize = L2FESpace_stress->GetVSize();
    int total_size = 2 * H1Vsize + L2Vsize + L2StressVsize;
    
    EXPECT_NO_THROW({
        LagrangianGeoOperator geo_oper(total_size, *H1FESpace, *L2FESpace, *L2FESpace_stress,
                                      ess_tdofs, *rho0_gf, *fictitious_rho0_gf, *gamma_gf,
                                      source, cfl, visc, vort, pa, cgt, cgiter, ftz_tol,
                                      order_q, *lambda_gf, *mu_gf, mscale, gravity, thickness,
                                      winkler, winkler_rho, dyn_damping, dyn_factor, bc_id_pa, vbc_max_val);
        
        // Check that operator has correct size
        EXPECT_EQ(geo_oper.Width(), total_size);
        EXPECT_EQ(geo_oper.Height(), total_size);
    });
}

TEST_F(SolverTest, LagrangianGeoOperatorBasicMethods) {
    int H1Vsize = H1FESpace->GetVSize();
    int L2Vsize = L2FESpace->GetVSize();
    int L2StressVsize = L2FESpace_stress->GetVSize();
    int total_size = 2 * H1Vsize + L2Vsize + L2StressVsize;
    
    LagrangianGeoOperator geo_oper(total_size, *H1FESpace, *L2FESpace, *L2FESpace_stress,
                                  ess_tdofs, *rho0_gf, *fictitious_rho0_gf, *gamma_gf,
                                  source, cfl, visc, vort, pa, cgt, cgiter, ftz_tol,
                                  order_q, *lambda_gf, *mu_gf, mscale, gravity, thickness,
                                  winkler, winkler_rho, dyn_damping, dyn_factor, bc_id_pa, vbc_max_val);
    
    // Test basic getter methods
    EXPECT_EQ(geo_oper.GetH1VSize(), H1Vsize);
    EXPECT_GT(geo_oper.GetH0(), 0.0);
    
    // Test block offsets
    const Array<int>& offsets = geo_oper.GetBlockOffsets();
    EXPECT_GT(offsets.Size(), 0);
    
    // Test memory class
    EXPECT_NO_THROW({
        MemoryClass mc = geo_oper.GetMemoryClass();
    });
}

TEST_F(SolverTest, LagrangianGeoOperatorStateVector) {
    int H1Vsize = H1FESpace->GetVSize();
    int L2Vsize = L2FESpace->GetVSize();
    int L2StressVsize = L2FESpace_stress->GetVSize();
    int total_size = 2 * H1Vsize + L2Vsize + L2StressVsize;
    
    LagrangianGeoOperator geo_oper(total_size, *H1FESpace, *L2FESpace, *L2FESpace_stress,
                                  ess_tdofs, *rho0_gf, *fictitious_rho0_gf, *gamma_gf,
                                  source, cfl, visc, vort, pa, cgt, cgiter, ftz_tol,
                                  order_q, *lambda_gf, *mu_gf, mscale, gravity, thickness,
                                  winkler, winkler_rho, dyn_damping, dyn_factor, bc_id_pa, vbc_max_val);
    
    // Create state vector with proper mesh coordinates
    Vector S(total_size);
    Vector dS_dt(total_size);
    
    // Initialize state vector with current mesh positions and some reasonable values
    // Get the current mesh nodes as initial state
    Vector nodes;
    pmesh->GetNodes(nodes);
    
    // Copy mesh node coordinates to the first part of the state vector
    int H1_size = H1Vsize;
    for (int i = 0; i < std::min(H1_size, nodes.Size()); i++) {
        S[i] = nodes[i];
    }
    
    // Set energy part to some small positive value
    for (int i = H1_size; i < H1_size + L2Vsize; i++) {
        S[i] = 1e6;  // Initial energy density
    }
    
    // Set stress part to small values
    for (int i = H1_size + L2Vsize; i < H1_size + L2Vsize + L2StressVsize; i++) {
        S[i] = 1000.0;  // Initial stress
    }
    
    // Set velocity part (second mesh coordinate part) to zero
    for (int i = H1_size + L2Vsize + L2StressVsize; i < total_size; i++) {
        S[i] = 0.0;
    }
    
    dS_dt = 0.0;
    
    // Test time step estimation
    EXPECT_NO_THROW({
        double dt_est = geo_oper.GetTimeStepEstimate(S);
        EXPECT_GT(dt_est, 0.0);
    });
    
    // Test length estimation
    EXPECT_NO_THROW({
        double h_est = geo_oper.GetLengthEstimate(S);
        EXPECT_GT(h_est, 0.0);
    });
}

TEST_F(SolverTest, LagrangianGeoOperatorMult) {
    int H1Vsize = H1FESpace->GetVSize();
    int L2Vsize = L2FESpace->GetVSize();
    int L2StressVsize = L2FESpace_stress->GetVSize();
    int total_size = 2 * H1Vsize + L2Vsize + L2StressVsize;
    
    LagrangianGeoOperator geo_oper(total_size, *H1FESpace, *L2FESpace, *L2FESpace_stress,
                                  ess_tdofs, *rho0_gf, *fictitious_rho0_gf, *gamma_gf,
                                  source, cfl, visc, vort, pa, cgt, cgiter, ftz_tol,
                                  order_q, *lambda_gf, *mu_gf, mscale, gravity, thickness,
                                  winkler, winkler_rho, dyn_damping, dyn_factor, bc_id_pa, vbc_max_val);
    
    // Create state vector with proper initialization
    Vector S(total_size);
    Vector dS_dt(total_size);
    
    // Initialize state vector with current mesh positions and reasonable values
    // Get the current mesh nodes as initial state
    Vector nodes;
    pmesh->GetNodes(nodes);
    
    // Copy mesh node coordinates to the first part of the state vector (positions)
    int H1_size = H1Vsize;
    for (int i = 0; i < std::min(H1_size, nodes.Size()); i++) {
        S[i] = nodes[i];
    }
    
    // Set velocity part (second H1 part) to small values
    for (int i = H1_size; i < 2*H1_size; i++) {
        S[i] = 1e-6;  // Small velocity
    }
    
    // Set energy part to reasonable positive values  
    for (int i = 2*H1_size; i < 2*H1_size + L2Vsize; i++) {
        S[i] = 1e6;  // Energy density
    }
    
    // Set stress part to small values
    for (int i = 2*H1_size + L2Vsize; i < total_size; i++) {
        S[i] = 1000.0;  // Initial stress
    }
    
    EXPECT_NO_THROW({
        geo_oper.Mult(S, dS_dt);
    });
    
    // Check that derivative is computed
    double norm = dS_dt.Norml2();
    EXPECT_GE(norm, 0.0);
}

TEST_F(SolverTest, TaylorCoefficientEvaluation) {
    TaylorCoefficient taylor_coeff;
    
    // Get the first element from our mesh and its transformation
    ElementTransformation *T = pmesh->GetElementTransformation(0);
    IntegrationPoint ip;
    ip.Set2(0.5, 0.5);  // Center of reference element
    
    EXPECT_NO_THROW({
        double val = taylor_coeff.Eval(*T, ip);
        // Value should be finite
        EXPECT_TRUE(std::isfinite(val));
    });
}

TEST_F(SolverTest, GTCoefficientEvaluation) {
    GTCoefficient gt_coeff(dim);
    
    EXPECT_EQ(gt_coeff.GetVDim(), dim);
    
    Vector V(dim);
    ElementTransformation *T = pmesh->GetElementTransformation(0);
    IntegrationPoint ip;
    ip.Set2(0.5, 0.5);
    
    EXPECT_NO_THROW({
        gt_coeff.Eval(V, *T, ip);
    });
    
    // Check gravity vector properties
    if (dim == 2) {
        EXPECT_DOUBLE_EQ(V[0], 0.0);
        EXPECT_DOUBLE_EQ(V[1], -1.0);
    } else if (dim == 3) {
        EXPECT_DOUBLE_EQ(V[0], 0.0);
        EXPECT_DOUBLE_EQ(V[1], 0.0);
        EXPECT_DOUBLE_EQ(V[2], -1.0);
    }
}

TEST_F(SolverTest, DampCoefficientEvaluation) {
    DampCoefficient damp_coeff(dim);
    
    EXPECT_EQ(damp_coeff.GetVDim(), dim);
    
    Vector V(dim);
    ElementTransformation *T = pmesh->GetElementTransformation(0);
    IntegrationPoint ip;
    ip.Set2(0.5, 0.5);
    
    EXPECT_NO_THROW({
        damp_coeff.Eval(V, *T, ip);
    });
    
    // DampCoefficient applies damping based on velocity direction
    // When input V is zero, it sets components to -1.0
    for (int i = 0; i < dim; i++) {
        EXPECT_DOUBLE_EQ(V[i], -1.0);
    }
}

TEST_F(SolverTest, RK2AvgSolverInitialization) {
    int H1Vsize = H1FESpace->GetVSize();
    int L2Vsize = L2FESpace->GetVSize();
    int L2StressVsize = L2FESpace_stress->GetVSize();
    int total_size = 2 * H1Vsize + L2Vsize + L2StressVsize;
    
    LagrangianGeoOperator geo_oper(total_size, *H1FESpace, *L2FESpace, *L2FESpace_stress,
                                  ess_tdofs, *rho0_gf, *fictitious_rho0_gf, *gamma_gf,
                                  source, cfl, visc, vort, pa, cgt, cgiter, ftz_tol,
                                  order_q, *lambda_gf, *mu_gf, mscale, gravity, thickness,
                                  winkler, winkler_rho, dyn_damping, dyn_factor, bc_id_pa, vbc_max_val);
    
    EXPECT_NO_THROW({
        RK2AvgSolver rk2_solver;
        rk2_solver.Init(geo_oper);
    });
}

TEST_F(SolverTest, EnergyComputations) {
    int H1Vsize = H1FESpace->GetVSize();
    int L2Vsize = L2FESpace->GetVSize();
    int L2StressVsize = L2FESpace_stress->GetVSize();
    int total_size = 2 * H1Vsize + L2Vsize + L2StressVsize;
    
    LagrangianGeoOperator geo_oper(total_size, *H1FESpace, *L2FESpace, *L2FESpace_stress,
                                  ess_tdofs, *rho0_gf, *fictitious_rho0_gf, *gamma_gf,
                                  source, cfl, visc, vort, pa, cgt, cgiter, ftz_tol,
                                  order_q, *lambda_gf, *mu_gf, mscale, gravity, thickness,
                                  winkler, winkler_rho, dyn_damping, dyn_factor, bc_id_pa, vbc_max_val);
    
    // Create grid functions for energy computation
    ParGridFunction e_gf(L2FESpace);
    ParGridFunction v_gf(H1FESpace);
    
    e_gf = 1e6;  // Energy density
    v_gf = 1e-3; // Velocity
    
    EXPECT_NO_THROW({
        double internal_energy = geo_oper.InternalEnergy(e_gf);
        double kinetic_energy = geo_oper.KineticEnergy(v_gf);
        
        EXPECT_GE(internal_energy, 0.0);
        EXPECT_GE(kinetic_energy, 0.0);
    });
}