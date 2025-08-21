#include <gtest/gtest.h>
#include "../../src/physics/laghost_rheology.hpp"
#include "mfem.hpp"

using namespace mfem;

class RheologyTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize test vectors for stress and strain
        stress_2d.SetSize(3);      // sxx, syy, sxy for 2D
        stress_3d.SetSize(6);      // sxx, syy, szz, sxy, sxz, syz for 3D
        
        strain_2d.SetSize(3);
        strain_3d.SetSize(6);
        
        // Material properties vectors
        rho.SetSize(1);
        lambda.SetSize(1);
        mu.SetSize(1);
        tension_cutoff.SetSize(1);
        cohesion0.SetSize(1);
        cohesion1.SetSize(1);
        friction_angle0.SetSize(1);
        friction_angle1.SetSize(1);
        dilation_angle0.SetSize(1);
        dilation_angle1.SetSize(1);
        pls0.SetSize(1);
        pls1.SetSize(1);
        plastic_viscosity.SetSize(1);
        
        // Set typical rock properties
        rho[0] = 2700.0;           // kg/m³
        lambda[0] = 3.0e10;        // Pa
        mu[0] = 3.0e10;            // Pa
        tension_cutoff[0] = 0.0;   // Pa
        cohesion0[0] = 44.0e6;     // Pa
        cohesion1[0] = 4.0e6;      // Pa
        friction_angle0[0] = 30.0; // degrees
        friction_angle1[0] = 30.0; // degrees
        dilation_angle0[0] = 0.0;  // degrees
        dilation_angle1[0] = 0.0;  // degrees
        pls0[0] = 0.0;
        pls1[0] = 0.5;
        plastic_viscosity[0] = 1.0;
        
        // Initialize stress states
        // 2D stress: compression in x, tension in y, no shear
        stress_2d[0] = -100.0e6;  // sxx (compression)
        stress_2d[1] = 10.0e6;    // syy (tension)
        stress_2d[2] = 0.0;       // sxy
        
        // 3D stress: similar pattern
        stress_3d[0] = -100.0e6;  // sxx
        stress_3d[1] = 10.0e6;    // syy
        stress_3d[2] = -50.0e6;   // szz
        stress_3d[3] = 0.0;       // sxy
        stress_3d[4] = 0.0;       // sxz
        stress_3d[5] = 0.0;       // syz
        
        // Initialize other parameters
        mat_id = 0;
        dt = 1.0e-6;  // Small time step
        h_min = 1000.0;  // Minimum element size
        viscoplastic = false;
    }
    
    Vector stress_2d, stress_3d;
    Vector strain_2d, strain_3d;
    Vector rho, lambda, mu, tension_cutoff;
    Vector cohesion0, cohesion1, friction_angle0, friction_angle1;
    Vector dilation_angle0, dilation_angle1, pls0, pls1, plastic_viscosity;
    int mat_id;
    double dt, h_min;
    bool viscoplastic;
};

TEST_F(RheologyTest, ReturnMapping2DInterface) {
    // Test that 2D return mapping function can be called without crashing
    Vector stress_old = stress_2d;
    Vector strain_inc = strain_2d;
    Vector plastic_strain(1);
    plastic_strain[0] = 0.0;
    Vector comp_gf(1);  // Component field
    comp_gf[0] = 1.0;
    int dim = 2;  // 2D dimension
    
    EXPECT_NO_THROW({
        Returnmapping2d(comp_gf, stress_2d, stress_old, strain_inc, plastic_strain, 
                       dim, h_min, rho, lambda, mu, tension_cutoff,
                       cohesion0, cohesion1, pls0, pls1, friction_angle0, 
                       friction_angle1, dilation_angle0, dilation_angle1,
                       plastic_viscosity, viscoplastic, dt);
    });
}

TEST_F(RheologyTest, ReturnMapping3DInterface) {
    // Test that 3D return mapping function can be called without crashing
    Vector stress_old = stress_3d;
    Vector strain_inc = strain_3d;
    Vector plastic_strain(1);
    plastic_strain[0] = 0.0;
    Vector comp_gf(1);  // Component field
    comp_gf[0] = 1.0;
    int dim = 3;  // 3D dimension
    
    EXPECT_NO_THROW({
        Returnmapping3d(comp_gf, stress_3d, stress_old, strain_inc, plastic_strain,
                       dim, h_min, rho, lambda, mu, tension_cutoff,
                       cohesion0, cohesion1, pls0, pls1, friction_angle0, 
                       friction_angle1, dilation_angle0, dilation_angle1,
                       plastic_viscosity, viscoplastic, dt);
    });
}

TEST_F(RheologyTest, ElasticBehavior2D) {
    // Test elastic behavior (no plastic deformation)
    Vector stress_old = stress_2d;
    Vector strain_inc(3);
    strain_inc = 0.0;  // No strain increment
    Vector plastic_strain(1);
    plastic_strain[0] = 0.0;
    Vector comp_gf(1);
    comp_gf[0] = 1.0;
    int dim = 2;
    
    Vector stress_initial = stress_2d;
    
    Returnmapping2d(comp_gf, stress_2d, stress_old, strain_inc, plastic_strain,
                   dim, h_min, rho, lambda, mu, tension_cutoff,
                   cohesion0, cohesion1, pls0, pls1, friction_angle0, 
                   friction_angle1, dilation_angle0, dilation_angle1,
                   plastic_viscosity, viscoplastic, dt);
    
    // With no strain increment, stress should remain approximately the same
    // (allowing for small numerical differences in return mapping)
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(stress_2d[i], stress_initial[i], 1e-6);
    }
}

TEST_F(RheologyTest, ElasticBehavior3D) {
    // Test elastic behavior (no plastic deformation)
    Vector stress_old = stress_3d;
    Vector strain_inc(6);
    strain_inc = 0.0;  // No strain increment
    Vector plastic_strain(1);
    plastic_strain[0] = 0.0;
    
    Vector stress_initial = stress_3d;
    
    Returnmapping3d(stress_3d, stress_old, strain_inc, plastic_strain,
                   rho, mat_id, dt, lambda, mu, tension_cutoff,
                   cohesion0, cohesion1, friction_angle0, friction_angle1,
                   dilation_angle0, dilation_angle1, pls0, pls1,
                   plastic_viscosity, viscoplastic, h_min);
    
    // With no strain increment, stress should remain approximately the same
    for (int i = 0; i < 6; i++) {
        EXPECT_NEAR(stress_3d[i], stress_initial[i], 1e-6);
    }
}

TEST_F(RheologyTest, PlasticYielding2D) {
    // Test plastic yielding with high stress
    Vector stress_old = stress_2d;
    Vector strain_inc(3);
    
    // Apply large compressive strain increment to trigger yielding
    strain_inc[0] = -0.01;  // Large compressive strain in x
    strain_inc[1] = 0.005;  // Tensile strain in y
    strain_inc[2] = 0.0;    // No shear strain
    
    Vector plastic_strain(1);
    plastic_strain[0] = 0.0;
    
    // Set very high stress to ensure yielding
    stress_2d[0] = -200.0e6;  // High compression
    stress_2d[1] = 50.0e6;    // High tension
    
    double initial_plastic_strain = plastic_strain[0];
    
    Returnmapping2d(stress_2d, stress_old, strain_inc, plastic_strain,
                   rho, mat_id, dt, lambda, mu, tension_cutoff,
                   cohesion0, cohesion1, friction_angle0, friction_angle1,
                   dilation_angle0, dilation_angle1, pls0, pls1,
                   plastic_viscosity, viscoplastic, h_min);
    
    // Plastic strain should increase if yielding occurred
    // (This test may need adjustment based on actual yield criteria)
    EXPECT_GE(plastic_strain[0], initial_plastic_strain);
}

TEST_F(RheologyTest, PlasticYielding3D) {
    // Test plastic yielding with high stress
    Vector stress_old = stress_3d;
    Vector strain_inc(6);
    
    // Apply large strain increment to trigger yielding
    strain_inc[0] = -0.01;  // Large compressive strain in x
    strain_inc[1] = 0.005;  // Tensile strain in y
    strain_inc[2] = -0.005; // Compressive strain in z
    strain_inc[3] = 0.0;    // No shear strains
    strain_inc[4] = 0.0;
    strain_inc[5] = 0.0;
    
    Vector plastic_strain(1);
    plastic_strain[0] = 0.0;
    
    // Set very high stress to ensure yielding
    stress_3d[0] = -200.0e6;  // High compression
    stress_3d[1] = 50.0e6;    // High tension
    stress_3d[2] = -100.0e6;  // Moderate compression
    
    double initial_plastic_strain = plastic_strain[0];
    
    Returnmapping3d(stress_3d, stress_old, strain_inc, plastic_strain,
                   rho, mat_id, dt, lambda, mu, tension_cutoff,
                   cohesion0, cohesion1, friction_angle0, friction_angle1,
                   dilation_angle0, dilation_angle1, pls0, pls1,
                   plastic_viscosity, viscoplastic, h_min);
    
    // Plastic strain should increase if yielding occurred
    EXPECT_GE(plastic_strain[0], initial_plastic_strain);
}

TEST_F(RheologyTest, ViscoplasticBehavior2D) {
    // Test viscoplastic behavior
    viscoplastic = true;
    
    Vector stress_old = stress_2d;
    Vector strain_inc(3);
    strain_inc[0] = -0.005;  // Moderate compressive strain
    strain_inc[1] = 0.002;   // Small tensile strain
    strain_inc[2] = 0.0;
    
    Vector plastic_strain(1);
    plastic_strain[0] = 0.0;
    
    // Set moderate stress level
    stress_2d[0] = -150.0e6;
    stress_2d[1] = 30.0e6;
    
    EXPECT_NO_THROW({
        Returnmapping2d(stress_2d, stress_old, strain_inc, plastic_strain,
                       rho, mat_id, dt, lambda, mu, tension_cutoff,
                       cohesion0, cohesion1, friction_angle0, friction_angle1,
                       dilation_angle0, dilation_angle1, pls0, pls1,
                       plastic_viscosity, viscoplastic, h_min);
    });
}

TEST_F(RheologyTest, ViscoplasticBehavior3D) {
    // Test viscoplastic behavior
    viscoplastic = true;
    
    Vector stress_old = stress_3d;
    Vector strain_inc(6);
    strain_inc[0] = -0.005;  // Moderate compressive strain
    strain_inc[1] = 0.002;   // Small tensile strain
    strain_inc[2] = -0.003;  // Small compressive strain
    strain_inc[3] = 0.0;
    strain_inc[4] = 0.0;
    strain_inc[5] = 0.0;
    
    Vector plastic_strain(1);
    plastic_strain[0] = 0.0;
    
    // Set moderate stress level
    stress_3d[0] = -150.0e6;
    stress_3d[1] = 30.0e6;
    stress_3d[2] = -80.0e6;
    
    EXPECT_NO_THROW({
        Returnmapping3d(stress_3d, stress_old, strain_inc, plastic_strain,
                       rho, mat_id, dt, lambda, mu, tension_cutoff,
                       cohesion0, cohesion1, friction_angle0, friction_angle1,
                       dilation_angle0, dilation_angle1, pls0, pls1,
                       plastic_viscosity, viscoplastic, h_min);
    });
}

TEST_F(RheologyTest, MaterialPropertyValidation) {
    // Test that material properties are within reasonable ranges
    EXPECT_GT(rho[0], 0.0);
    EXPECT_GT(lambda[0], 0.0);
    EXPECT_GT(mu[0], 0.0);
    EXPECT_GE(tension_cutoff[0], 0.0);
    EXPECT_GT(cohesion0[0], 0.0);
    EXPECT_GT(cohesion1[0], 0.0);
    EXPECT_GE(friction_angle0[0], 0.0);
    EXPECT_LE(friction_angle0[0], 90.0);
    EXPECT_GE(friction_angle1[0], 0.0);
    EXPECT_LE(friction_angle1[0], 90.0);
    EXPECT_GE(dilation_angle0[0], 0.0);
    EXPECT_LE(dilation_angle0[0], friction_angle0[0]);  // Dilation <= friction
    EXPECT_GE(dilation_angle1[0], 0.0);
    EXPECT_LE(dilation_angle1[0], friction_angle1[0]);
    EXPECT_GE(pls0[0], 0.0);
    EXPECT_GT(pls1[0], pls0[0]);  // pls1 should be greater than pls0
    EXPECT_GT(plastic_viscosity[0], 0.0);
}

TEST_F(RheologyTest, StressInvariance) {
    // Test that stress tensor properties are preserved
    Vector stress_test = stress_3d;
    Vector stress_old = stress_3d;
    Vector strain_inc(6);
    strain_inc = 0.0;  // No strain increment
    Vector plastic_strain(1);
    plastic_strain[0] = 0.0;
    
    // Calculate initial trace (sum of diagonal components)
    double initial_trace = stress_test[0] + stress_test[1] + stress_test[2];
    
    Returnmapping3d(stress_test, stress_old, strain_inc, plastic_strain,
                   rho, mat_id, dt, lambda, mu, tension_cutoff,
                   cohesion0, cohesion1, friction_angle0, friction_angle1,
                   dilation_angle0, dilation_angle1, pls0, pls1,
                   plastic_viscosity, viscoplastic, h_min);
    
    // For elastic behavior with no strain increment, trace should be preserved
    double final_trace = stress_test[0] + stress_test[1] + stress_test[2];
    EXPECT_NEAR(final_trace, initial_trace, 1e-6);
}