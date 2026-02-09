// Copyright (c) 2024, Laghost Authors
// SPDX-License-Identifier: BSD-3-Clause
//
// Constitutive Model Interface
// Abstract base class for all material constitutive models

#ifndef LAGHOST_CONSTITUTIVE_MODEL_HPP
#define LAGHOST_CONSTITUTIVE_MODEL_HPP

#include "mfem.hpp"
#include <string>

namespace mfem
{
namespace geodynamics
{

/// Material properties for a single material
struct MaterialProperties
{
   double rho;               // Density
   double lambda;            // First Lamé parameter
   double mu;                // Shear modulus (second Lamé parameter)
   double tension_cutoff;    // Tensile strength
   double cohesion0;         // Initial cohesion
   double cohesion1;         // Final cohesion (after weakening)
   double friction_angle0;   // Initial friction angle (degrees)
   double friction_angle1;   // Final friction angle (degrees)
   double dilation_angle0;   // Initial dilation angle (degrees)
   double dilation_angle1;   // Final dilation angle (degrees)
   double alpha0;            // Strain at start of weakening
   double alpha1;            // Strain at end of weakening
   double viscosity;         // Plastic viscosity (for viscoplastic models)
};

/// State variables at a quadrature point
struct MaterialState
{
   DenseMatrix stress;       // Cauchy stress tensor
   DenseMatrix stress_old;   // Stress at previous time step
   double alpha;             // Internal hardening variable (accumulated plastic strain)
   int material_id;          // Material identifier
   
   MaterialState(int dim = 3) : stress(dim), stress_old(dim), alpha(0.0), material_id(0) {}
};

/// Abstract base class for constitutive models
class ConstitutiveModel
{
public:
   virtual ~ConstitutiveModel() = default;
   
   /// Update stress and internal variable α at a single quadrature point
   /// @param strain_rate Strain rate tensor (symmetric)
   /// @param spin Spin tensor (antisymmetric part of velocity gradient)
   /// @param state Material state (stress, α) - modified in place
   /// @param props Material properties
   /// @param dt Time step
   virtual void Update(
      const DenseMatrix &strain_rate,
      const DenseMatrix &spin,
      MaterialState &state,
      const MaterialProperties &props,
      double dt
   ) = 0;
   
   /// Batch update for all quadrature points (CPU version)
   /// Default implementation loops over points; override for efficiency
   virtual void UpdateBatch(
      const Vector &strain_rate_flat,   // [6 * nqp] for 3D
      const Vector &spin_flat,          // [3 * nqp] for 3D
      Vector &stress_flat,              // [6 * nqp] for 3D (in/out)
      Vector &alpha,                    // [nqp] (in/out)
      const Vector &props_flat,         // Material properties per qp
      const Array<int> &mat_ids,        // Material ID per qp
      double dt,
      int nqp,
      int dim
   );
   
   /// Device kernel interface for GPU/SYCL (flattened data)
   /// Override this for GPU-compatible implementations
   virtual void UpdateDevice(
      const double *d_strain_rate,
      const double *d_spin,
      double *d_stress,
      double *d_alpha,
      const double *d_props,
      const int *d_mat_ids,
      double dt,
      int nqp,
      int dim
   );
   
   /// Name of the constitutive model
   virtual std::string Name() const = 0;
   
   /// Whether this model requires temperature field
   virtual bool RequiresTemperature() const { return false; }
   
   /// Number of internal state variables (beyond stress)
   virtual int NumStateVariables() const { return 1; }  // Default: just α
   
   /// Get state variable name for output (e.g., "accumulated_plastic_strain")
   virtual std::string StateVariableName(int i) const 
   { 
      if (i == 0) return "accumulated_plastic_strain";
      return "unknown";
   }
};

} // namespace geodynamics
} // namespace mfem

#endif // LAGHOST_CONSTITUTIVE_MODEL_HPP
