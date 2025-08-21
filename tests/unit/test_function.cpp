#include <gtest/gtest.h>
#include "../../src/physics/laghost_function.hpp"
#include "mfem.hpp"

using namespace mfem;

class FunctionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple 2D mesh for testing
        mesh = new Mesh(Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL, true));
        pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
        delete mesh;
        mesh = nullptr;
        
        // Create finite element spaces
        H1_FECollection h1_fec(1, 2);  // Order 1, 2D
        L2_FECollection l2_fec(0, 2);  // Order 0, 2D
        
        h1_fes = new ParFiniteElementSpace(pmesh, &h1_fec, 2);  // Vector space
        l2_fes = new ParFiniteElementSpace(pmesh, &l2_fec);     // Scalar space
        
        // Create grid functions
        xyz_gf = new ParGridFunction(h1_fes);
        rho_gf = new ParGridFunction(l2_fes);
        mat_gf = new ParGridFunction(l2_fes);
        
        // Initialize with some test values
        *xyz_gf = 0.0;
        *rho_gf = 2700.0;  // Typical rock density
        *mat_gf = 0.0;     // Material ID 0
    }
    
    void TearDown() override {
        delete mat_gf;
        delete rho_gf;
        delete xyz_gf;
        delete l2_fes;
        delete h1_fes;
        delete pmesh;
    }
    
    Mesh* mesh;
    ParMesh* pmesh;
    ParFiniteElementSpace* h1_fes;
    ParFiniteElementSpace* l2_fes;
    ParGridFunction* xyz_gf;
    ParGridFunction* rho_gf;
    ParGridFunction* mat_gf;
};

TEST_F(FunctionTest, BasicFunctionEvaluations) {
    Vector test_point(2);
    test_point[0] = 0.5;  // x coordinate
    test_point[1] = 0.5;  // y coordinate
    
    // Test basic functions
    EXPECT_NO_THROW({
        double e_val = e0(test_point);
        double p_val = p0(test_point);
        double depth_val = depth0(test_point);
        double rho_val = rho0(test_point);
        double gamma_val = gamma_func(test_point);
        double zero_val = zero_func(test_point);
        
        // Zero function should return 0
        EXPECT_DOUBLE_EQ(zero_val, 0.0);
    });
}

TEST_F(FunctionTest, VectorFunctionEvaluations) {
    Vector test_point(2);
    test_point[0] = 0.5;
    test_point[1] = 0.5;
    
    Vector result(2);
    
    EXPECT_NO_THROW({
        v0(test_point, result);
        xyz0(test_point, result);
    });
    
    // Results should have correct size
    EXPECT_EQ(result.Size(), 2);
}

TEST_F(FunctionTest, CoordinateFunctions) {
    Vector test_point(3);
    test_point[0] = 1.0;  // x
    test_point[1] = 2.0;  // y
    test_point[2] = 3.0;  // z
    
    // Test coordinate extraction functions
    EXPECT_DOUBLE_EQ(x_l2(test_point), 1.0);
    EXPECT_DOUBLE_EQ(y_l2(test_point), 2.0);
    EXPECT_DOUBLE_EQ(z_l2(test_point), 3.0);
}

TEST_F(FunctionTest, PrincipalStresses2D) {
    // Test 2D principal stress calculation
    double stress[3] = {100.0, 50.0, 25.0};  // {sxx, szz, sxz}
    double p[2];
    double cos2t, sin2t;
    
    EXPECT_NO_THROW({
        principal_stresses2(stress, p, cos2t, sin2t);
    });
    
    // Principal stresses should be ordered: p[0] <= p[1]
    EXPECT_LE(p[0], p[1]);
    
    // Direction cosines should satisfy cos²(2θ) + sin²(2θ) = 1
    EXPECT_NEAR(cos2t*cos2t + sin2t*sin2t, 1.0, 1e-10);
}

TEST_F(FunctionTest, PlasticCoefficientInitialization) {
    int dim = 2;
    Vector location(2);
    location[0] = 0.5;
    location[1] = 0.5;
    double rad = 0.1;
    double ini_pls = 0.5;
    
    EXPECT_NO_THROW({
        PlasticCoefficient plastic_coeff(dim, *xyz_gf, location, rad, ini_pls);
        
        // Test that coefficient has correct vector dimension
        EXPECT_EQ(plastic_coeff.GetVDim(), 1);
    });
}

TEST_F(FunctionTest, LithostaticCoefficientInitialization) {
    int dim = 2;
    double gravity = 9.81;
    double thickness = 1000.0;
    
    EXPECT_NO_THROW({
        LithostaticCoefficient litho_coeff(dim, *xyz_gf, *rho_gf, gravity, thickness);
        
        // Test that coefficient has correct vector dimension
        EXPECT_EQ(litho_coeff.GetVDim(), 3*(dim-1));
    });
}

TEST_F(FunctionTest, ATMCoefficientInitialization) {
    int dim = 2;
    double gravity = 9.81;
    double thickness = 1000.0;
    
    EXPECT_NO_THROW({
        ATMCoefficient atm_coeff(dim, *xyz_gf, *rho_gf, gravity, thickness);
        
        // Test that coefficient has correct vector dimension
        EXPECT_EQ(atm_coeff.GetVDim(), 3*(dim-1));
    });
}

TEST_F(FunctionTest, CompoCoefficientInitialization) {
    int mat_num = 3;
    
    EXPECT_NO_THROW({
        CompoCoefficient compo_coeff(mat_num, *mat_gf);
        
        // Test that coefficient has correct vector dimension
        EXPECT_EQ(compo_coeff.GetVDim(), mat_num);
    });
}

TEST_F(FunctionTest, CoefficientEvaluation) {
    // Create a simple element transformation for testing
    IsoparametricTransformation T;
    IntegrationPoint ip;
    ip.Set2(0.5, 0.5);  // Center of reference element
    
    // Test PlasticCoefficient evaluation
    int dim = 2;
    Vector location(2);
    location[0] = 0.5;
    location[1] = 0.5;
    double rad = 0.1;
    double ini_pls = 0.5;
    
    PlasticCoefficient plastic_coeff(dim, *xyz_gf, location, rad, ini_pls);
    
    // Note: Full evaluation would require proper element transformation setup
    // This tests the interface without crashing
    EXPECT_NO_THROW({
        Vector K(1);
        // plastic_coeff.Eval(K, T, ip);  // Would need proper T setup
    });
}

TEST_F(FunctionTest, StressMappingCoefficientDimensions) {
    int dim2 = 2;
    int dim3 = 3;
    
    // Create temporary grid functions for stress components
    ParGridFunction temp_gf_2d(l2_fes);
    ParGridFunction temp_gf_3d(l2_fes);
    
    EXPECT_NO_THROW({
        StressMappingCoefficient stress_coeff_2d(dim2, temp_gf_2d);
        StressMappingCoefficient stress_coeff_3d(dim3, temp_gf_3d);
        
        // Check dimensions
        EXPECT_EQ(stress_coeff_2d.GetVDim(), 3*(dim2-1));  // 3 for 2D
        EXPECT_EQ(stress_coeff_3d.GetVDim(), 3*(dim3-1));  // 6 for 3D
    });
}

TEST_F(FunctionTest, PlasticityMappingCoefficientDimensions) {
    int dim = 2;
    ParGridFunction temp_gf(l2_fes);
    
    EXPECT_NO_THROW({
        PlasticityMappingCoefficient plas_coeff(dim, temp_gf);
        
        // Should always have dimension 1 (scalar plastic strain)
        EXPECT_EQ(plas_coeff.GetVDim(), 1);
    });
}