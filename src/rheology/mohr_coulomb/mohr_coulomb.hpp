// Copyright (c) 2024, Laghost Authors
// SPDX-License-Identifier: BSD-3-Clause
//
// Mohr-Coulomb Constitutive Model
// Strain-weakening elasto-plastic model with Mohr-Coulomb yield criterion

#ifndef LAGHOST_MOHR_COULOMB_HPP
#define LAGHOST_MOHR_COULOMB_HPP

#include "../constitutive_model.hpp"
#include "../rheology_factory.hpp"

namespace mfem
{
namespace geodynamics
{

/// Mohr-Coulomb elasto-plastic model with strain weakening
class MohrCoulombModel : public ConstitutiveModel
{
public:
   MohrCoulombModel() = default;
   
   /// Update stress and internal variable α
   void Update(
      const DenseMatrix &strain_rate,
      const DenseMatrix &spin,
      MaterialState &state,
      const MaterialProperties &props,
      double dt
   ) override;
   
   std::string Name() const override { return "mohr_coulomb"; }
   
   bool RequiresTemperature() const override { return false; }
   
private:
   /// Compute strain-weakened material parameters
   void ComputeWeakenedParams(
      double alpha,
      const MaterialProperties &props,
      double &cohesion,
      double &friction_angle,
      double &dilation_angle
   ) const;
   
   /// Apply Jaumann stress rate correction
   void ApplyJaumannCorrection(
      DenseMatrix &stress,
      const DenseMatrix &spin,
      double dt
   ) const;
   
   /// Return mapping for shear yield
   double ReturnMappingShear(
      double sig1, double sig3,
      double cohesion, double friction_angle, double dilation_angle,
      double lambda, double mu,
      double &dsig1, double &dsig3
   ) const;
   
   /// Return mapping for tensile yield
   double ReturnMappingTensile(
      double sig3, double tension_cutoff,
      double lambda, double mu,
      double &dsig3
   ) const;
};

// Auto-register with factory
REGISTER_RHEOLOGY(MohrCoulombModel, "mohr_coulomb")

} // namespace geodynamics
} // namespace mfem

#endif // LAGHOST_MOHR_COULOMB_HPP
