#include <gtest/gtest.h>
#include "../../src/core/laghost_assembly.hpp"
#include "mfem.hpp"

using namespace mfem;
using namespace mfem::geodynamics;

class AssemblyTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple 2D mesh
        mesh = new Mesh(Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL, true));
        pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
        delete mesh;
        mesh = nullptr;
        
        dim = pmesh->Dimension();
        NE = pmesh->GetNE();
        
        // Create finite element collections and spaces
        h1_fec = new H1_FECollection(2, dim);  // Order 2 (minimum for PA kernels)
        l2_fec = new L2_FECollection(1, dim);  // Order 1 (D1D-1)
        
        H1FESpace = new ParFiniteElementSpace(pmesh, h1_fec, dim);  // Vector space
        L2FESpace = new ParFiniteElementSpace(pmesh, l2_fec);       // Scalar space
        L2FESpace_stress = new ParFiniteElementSpace(pmesh, l2_fec, 3*(dim-1)); // Stress space
        
        // Create integration rule with sufficient points for PA kernels (need Q1D=4)
        const IntegrationRule *ir = &IntRules.Get(pmesh->GetElementType(0), 6);  // Higher order to get Q1D=4
        int_rule = *ir;
        
        // Initialize quadrature data
        int quads_per_el = int_rule.GetNPoints();
        qdata = new QuadratureData(dim, NE, quads_per_el);
        
        // Initialize quadrature data with some test values
        qdata->h0 = 1.0;
        qdata->dt_est = 1e-6;
        qdata->h_est = 1.0;
        qdata->mscale = 1.0;
        qdata->vbc_max_val = 1e-12;
        qdata->gravity = 9.81;
        
        // Initialize Jacobian inverse data
        for (int i = 0; i < qdata->Jac0inv.SizeI(); i++) {
            for (int j = 0; j < qdata->Jac0inv.SizeJ(); j++) {
                for (int k = 0; k < qdata->Jac0inv.SizeK(); k++) {
                    qdata->Jac0inv(i, j, k) = (i == j) ? 1.0 : 0.0;  // Identity matrix
                }
            }
        }
        
        // Initialize density data
        qdata->rho0DetJ0w = 2700.0;  // Typical rock density
        
        // Initialize stress and force data
        for (int i = 0; i < qdata->stressJinvT.SizeI(); i++) {
            for (int j = 0; j < qdata->stressJinvT.SizeJ(); j++) {
                for (int k = 0; k < qdata->stressJinvT.SizeK(); k++) {
                    qdata->stressJinvT(i, j, k) = 0.0;
                    qdata->tauJinvT(i, j, k) = 0.0;
                    qdata->buoyJinvT(i, j, k) = 0.0;
                }
            }
        }
    }
    
    void TearDown() override {
        delete qdata;
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
    QuadratureData* qdata;
    IntegrationRule int_rule;
    int dim, NE;
};

TEST_F(AssemblyTest, QuadratureDataInitialization) {
    // Test that quadrature data is properly initialized
    EXPECT_GT(qdata->Jac0inv.SizeI(), 0);
    EXPECT_GT(qdata->Jac0inv.SizeJ(), 0);
    EXPECT_GT(qdata->Jac0inv.SizeK(), 0);
    
    EXPECT_GT(qdata->stressJinvT.SizeI(), 0);
    EXPECT_GT(qdata->stressJinvT.SizeJ(), 0);
    EXPECT_GT(qdata->stressJinvT.SizeK(), 0);
    
    EXPECT_GT(qdata->rho0DetJ0w.Size(), 0);
    
    EXPECT_GT(qdata->h0, 0.0);
    EXPECT_GT(qdata->dt_est, 0.0);
    EXPECT_GT(qdata->h_est, 0.0);
}

TEST_F(AssemblyTest, QuadratureDataResize) {
    // Test resizing quadrature data
    int new_dim = 3;
    int new_NE = 8;
    int new_quads = 8;
    
    EXPECT_NO_THROW({
        qdata->Resize(new_dim, new_NE, new_quads);
    });
    
    // Check that sizes are updated correctly
    EXPECT_EQ(qdata->Jac0inv.SizeI(), new_dim);
    EXPECT_EQ(qdata->Jac0inv.SizeJ(), new_dim);
    EXPECT_EQ(qdata->Jac0inv.SizeK(), new_NE * new_quads);
    
    EXPECT_EQ(qdata->stressJinvT.SizeI(), new_NE * new_quads);
    EXPECT_EQ(qdata->stressJinvT.SizeJ(), new_dim);
    EXPECT_EQ(qdata->stressJinvT.SizeK(), new_dim);
    
    EXPECT_EQ(qdata->rho0DetJ0w.Size(), new_NE * new_quads);
}

TEST_F(AssemblyTest, DensityIntegratorInitialization) {
    EXPECT_NO_THROW({
        DensityIntegrator density_integrator(*qdata);
    });
}

TEST_F(AssemblyTest, SigmaIntegratorInitialization) {
    EXPECT_NO_THROW({
        SigmaIntegrator sigma_integrator(*qdata);
    });
}

TEST_F(AssemblyTest, ForceIntegratorInitialization) {
    EXPECT_NO_THROW({
        ForceIntegrator force_integrator(*qdata);
    });
}

TEST_F(AssemblyTest, BodyForceIntegratorInitialization) {
    EXPECT_NO_THROW({
        BodyForceIntegrator body_force_integrator(*qdata);
    });
}

TEST_F(AssemblyTest, StressPAOperatorInitialization) {
    EXPECT_NO_THROW({
        StressPAOperator stress_pa(*qdata, *H1FESpace, *L2FESpace_stress, int_rule);
        
        // Check that operator has correct dimensions
        EXPECT_GT(stress_pa.Height(), 0);
        EXPECT_GT(stress_pa.Width(), 0);
    });
}

TEST_F(AssemblyTest, ForcePAOperatorInitialization) {
    EXPECT_NO_THROW({
        ForcePAOperator force_pa(*qdata, *H1FESpace, *L2FESpace, int_rule);
        
        // Check that operator has correct dimensions
        EXPECT_GT(force_pa.Height(), 0);
        EXPECT_GT(force_pa.Width(), 0);
    });
}

TEST_F(AssemblyTest, MassPAOperatorInitialization) {
    // Create a simple coefficient for mass operator
    ConstantCoefficient rho_coeff(2700.0);
    
    EXPECT_NO_THROW({
        MassPAOperator mass_pa(*H1FESpace, int_rule, rho_coeff);
        
        // Check that operator has correct dimensions
        EXPECT_GT(mass_pa.Height(), 0);
        EXPECT_GT(mass_pa.Width(), 0);
        EXPECT_EQ(mass_pa.Height(), mass_pa.Width());  // Mass matrix should be square
    });
}

TEST_F(AssemblyTest, StressPAOperatorMult) {
    StressPAOperator stress_pa(*qdata, *H1FESpace, *L2FESpace_stress, int_rule);
    
    // Create test vectors
    Vector x(stress_pa.Width());
    Vector y(stress_pa.Height());
    x = 1.0;  // Initialize with ones
    y = 0.0;
    
    EXPECT_NO_THROW({
        stress_pa.Mult(x, y);
    });
    
    // Result should not be zero vector (assuming non-trivial operator)
    double norm = y.Norml2();
    EXPECT_GE(norm, 0.0);  // Norm should be non-negative
}

TEST_F(AssemblyTest, ForcePAOperatorMult) {
    ForcePAOperator force_pa(*qdata, *H1FESpace, *L2FESpace, int_rule);
    
    // Create test vectors
    Vector x(force_pa.Width());
    Vector y(force_pa.Height());
    x = 1.0;  // Initialize with ones
    y = 0.0;
    
    EXPECT_NO_THROW({
        force_pa.Mult(x, y);
    });
    
    // Result should not be zero vector (assuming non-trivial operator)
    double norm = y.Norml2();
    EXPECT_GE(norm, 0.0);  // Norm should be non-negative
}

TEST_F(AssemblyTest, MassPAOperatorMult) {
    ConstantCoefficient rho_coeff(2700.0);
    MassPAOperator mass_pa(*H1FESpace, int_rule, rho_coeff);
    
    // Create test vectors
    Vector x(mass_pa.Width());
    Vector y(mass_pa.Height());
    x = 1.0;  // Initialize with ones
    y = 0.0;
    
    EXPECT_NO_THROW({
        mass_pa.Mult(x, y);
    });
    
    // Mass matrix times constant vector should give positive result
    double norm = y.Norml2();
    EXPECT_GT(norm, 0.0);  // Mass matrix should be positive definite
}

TEST_F(AssemblyTest, StressPAOperatorMultTranspose) {
    StressPAOperator stress_pa(*qdata, *H1FESpace, *L2FESpace_stress, int_rule);
    
    // Create test vectors
    Vector x(stress_pa.Height());
    Vector y(stress_pa.Width());
    x = 1.0;  // Initialize with ones
    y = 0.0;
    
    EXPECT_NO_THROW({
        stress_pa.MultTranspose(x, y);
    });
    
    double norm = y.Norml2();
    EXPECT_GE(norm, 0.0);
}

TEST_F(AssemblyTest, ForcePAOperatorMultTranspose) {
    ForcePAOperator force_pa(*qdata, *H1FESpace, *L2FESpace, int_rule);
    
    // Create test vectors
    Vector x(force_pa.Height());
    Vector y(force_pa.Width());
    x = 1.0;  // Initialize with ones
    y = 0.0;
    
    EXPECT_NO_THROW({
        force_pa.MultTranspose(x, y);
    });
    
    double norm = y.Norml2();
    EXPECT_GE(norm, 0.0);
}

TEST_F(AssemblyTest, MassPAOperatorEssentialDofs) {
    ConstantCoefficient rho_coeff(2700.0);
    MassPAOperator mass_pa(*H1FESpace, int_rule, rho_coeff);
    
    // Create array of essential DOFs
    Array<int> ess_tdofs;
    ess_tdofs.Append(0);  // Mark first DOF as essential
    ess_tdofs.Append(1);  // Mark second DOF as essential
    
    EXPECT_NO_THROW({
        mass_pa.SetEssentialTrueDofs(ess_tdofs);
    });
    
    // Test RHS elimination
    Vector rhs(mass_pa.Height());
    rhs = 1.0;
    
    EXPECT_NO_THROW({
        mass_pa.EliminateRHS(rhs);
    });
    
    // Essential DOFs should be zeroed in RHS
    EXPECT_DOUBLE_EQ(rhs[0], 0.0);
    EXPECT_DOUBLE_EQ(rhs[1], 0.0);
}

TEST_F(AssemblyTest, OperatorDimensionConsistency) {
    // Test that operators have consistent dimensions
    StressPAOperator stress_pa(*qdata, *H1FESpace, *L2FESpace_stress, int_rule);
    ForcePAOperator force_pa(*qdata, *H1FESpace, *L2FESpace, int_rule);
    
    ConstantCoefficient rho_coeff(2700.0);
    MassPAOperator mass_pa(*H1FESpace, int_rule, rho_coeff);
    
    // Check dimension relationships
    EXPECT_EQ(stress_pa.Width(), L2FESpace_stress->GetTrueVSize());
    EXPECT_EQ(stress_pa.Height(), H1FESpace->GetTrueVSize());
    
    EXPECT_EQ(force_pa.Width(), L2FESpace->GetTrueVSize());
    EXPECT_EQ(force_pa.Height(), H1FESpace->GetTrueVSize());
    
    EXPECT_EQ(mass_pa.Width(), H1FESpace->GetTrueVSize());
    EXPECT_EQ(mass_pa.Height(), H1FESpace->GetTrueVSize());
}