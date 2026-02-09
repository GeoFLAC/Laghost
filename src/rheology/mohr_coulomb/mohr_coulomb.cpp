// Copyright (c) 2024, Laghost Authors
// SPDX-License-Identifier: BSD-3-Clause
//
// Mohr-Coulomb Constitutive Model - Implementation
// Placeholder implementation to be filled with extracted code from laghost_rheology.cpp

#include "mohr_coulomb.hpp"
#include <cmath>

namespace mfem
{
namespace geodynamics
{

void MohrCoulombModel::Update(
   const DenseMatrix &strain_rate,
   const DenseMatrix &spin,
   MaterialState &state,
   const MaterialProperties &props,
   double dt)
{
   const int dim = strain_rate.Height();
   const double DEG2RAD = M_PI / 180.0;
   
   // Get weakened parameters based on α
   double cohesion, friction_angle, dilation_angle;
   ComputeWeakenedParams(state.alpha, props, cohesion, friction_angle, dilation_angle);
   
   // Elastic parameters
   const double lambda = props.lambda;
   const double mu = props.mu;
   
   // Compute elastic trial stress increment
   double trace_eps = 0.0;
   for (int i = 0; i < dim; i++)
   {
      trace_eps += strain_rate(i, i);
   }
   
   // σ_trial = σ_old + λ tr(ε̇) I dt + 2μ ε̇ dt
   DenseMatrix stress_trial(state.stress);
   for (int i = 0; i < dim; i++)
   {
      for (int j = 0; j < dim; j++)
      {
         stress_trial(i, j) += 2.0 * mu * strain_rate(i, j) * dt;
         if (i == j)
         {
            stress_trial(i, j) += lambda * trace_eps * dt;
         }
      }
   }
   
   // Principal stress decomposition
   // TODO: Extract eigenvalue computation from laghost_rheology.cpp
   // For now, use MFEM's CalcEigenvalues for symmetric matrices
   
   Vector eig_vals(dim);
   DenseMatrix eig_vecs(dim);
   
   if (dim == 2)
   {
      // 2D eigenvalue problem
      stress_trial.CalcEigenvalues(eig_vals.GetData(), eig_vecs.GetData());
   }
   else
   {
      // 3D eigenvalue problem
      stress_trial.CalcEigenvalues(eig_vals.GetData(), eig_vecs.GetData());
   }
   
   // Sort eigenvalues: sig1 >= sig2 >= sig3 (compressive positive convention)
   // Note: MFEM returns in ascending order, so we may need to reorder
   double sig1 = eig_vals[dim - 1];  // Most compressive
   double sig3 = eig_vals[0];        // Least compressive (most tensile)
   
   // Yield functions
   double N_phi = (1.0 + sin(friction_angle * DEG2RAD)) / 
                  (1.0 - sin(friction_angle * DEG2RAD));
   double fs = sig1 - N_phi * sig3 + 2.0 * cohesion * sqrt(N_phi);  // Shear
   double ft = sig3 - props.tension_cutoff;                           // Tensile
   
   double dalpha = 0.0;
   
   if (fs > 0.0 || ft > 0.0)
   {
      // Plastic correction needed
      // Determine yield mode using the "hybrid" criterion
      double fh_threshold = sig3 - props.tension_cutoff + 
                            (sqrt(N_phi) + N_phi) * 
                            (sig1 - N_phi * props.tension_cutoff + 
                             2.0 * cohesion * sqrt(N_phi));
      
      if (ft > 0.0 && fh_threshold < 0.0)
      {
         // Tensile failure
         dalpha = ReturnMappingTensile(sig3, props.tension_cutoff, lambda, mu, sig3);
         // Update sig3 in eigenvalue array
         eig_vals[0] = sig3;
      }
      else if (fs > 0.0)
      {
         // Shear failure
         double dsig1, dsig3;
         dalpha = ReturnMappingShear(sig1, sig3, cohesion, friction_angle, 
                                     dilation_angle, lambda, mu, dsig1, dsig3);
         eig_vals[dim - 1] = sig1 - dsig1;
         eig_vals[0] = sig3 - dsig3;
      }
   }
   
   // Reconstruct stress from corrected principal values
   // σ = V * diag(λ) * V^T
   DenseMatrix diag_sig(dim);
   diag_sig = 0.0;
   for (int i = 0; i < dim; i++)
   {
      diag_sig(i, i) = eig_vals[i];
   }
   
   DenseMatrix temp(dim);
   Mult(eig_vecs, diag_sig, temp);
   MultAtB(eig_vecs, temp, state.stress);
   
   // Apply Jaumann correction for objectivity
   ApplyJaumannCorrection(state.stress, spin, dt);
   
   // Update internal variable
   state.alpha += dalpha;
}

void MohrCoulombModel::ComputeWeakenedParams(
   double alpha,
   const MaterialProperties &props,
   double &cohesion,
   double &friction_angle,
   double &dilation_angle) const
{
   // Linear strain weakening between alpha0 and alpha1
   double f = 0.0;
   if (alpha >= props.alpha1)
   {
      f = 1.0;
   }
   else if (alpha > props.alpha0)
   {
      f = (alpha - props.alpha0) / (props.alpha1 - props.alpha0);
   }
   
   cohesion = props.cohesion0 + f * (props.cohesion1 - props.cohesion0);
   friction_angle = props.friction_angle0 + f * (props.friction_angle1 - props.friction_angle0);
   dilation_angle = props.dilation_angle0 + f * (props.dilation_angle1 - props.dilation_angle0);
}

void MohrCoulombModel::ApplyJaumannCorrection(
   DenseMatrix &stress,
   const DenseMatrix &spin,
   double dt) const
{
   // Jaumann rate: σ̌ = σ + (σ·ω - ω·σ) dt
   const int dim = stress.Height();
   DenseMatrix sw(dim), ws(dim);
   
   Mult(stress, spin, sw);
   Mult(spin, stress, ws);
   
   for (int i = 0; i < dim; i++)
   {
      for (int j = 0; j < dim; j++)
      {
         stress(i, j) += dt * (sw(i, j) - ws(i, j));
      }
   }
}

double MohrCoulombModel::ReturnMappingShear(
   double sig1, double sig3,
   double cohesion, double friction_angle, double dilation_angle,
   double lambda, double mu,
   double &dsig1, double &dsig3) const
{
   const double DEG2RAD = M_PI / 180.0;
   
   double N_phi = (1.0 + sin(friction_angle * DEG2RAD)) / 
                  (1.0 - sin(friction_angle * DEG2RAD));
   double N_psi = (1.0 + sin(dilation_angle * DEG2RAD)) / 
                  (1.0 - sin(dilation_angle * DEG2RAD));
   
   // Yield function value
   double fs = sig1 - N_phi * sig3 + 2.0 * cohesion * sqrt(N_phi);
   
   // Plastic multiplier
   double denom = (lambda + 2.0 * mu) * (1.0 + N_psi * N_phi) + 
                  2.0 * (lambda * (1.0 - N_psi) + mu * (1.0 + N_psi));
   double beta = fs / denom;
   
   // Stress corrections
   dsig1 = beta * ((lambda + 2.0 * mu) + lambda * (1.0 - N_psi));
   dsig3 = beta * (lambda * (1.0 + N_psi) - (lambda + 2.0 * mu) * N_psi);
   
   // Plastic strain increment (second invariant)
   double dalpha = beta * sqrt((1.0 + N_psi * N_psi) / 2.0);
   
   return dalpha;
}

double MohrCoulombModel::ReturnMappingTensile(
   double sig3, double tension_cutoff,
   double lambda, double mu,
   double &dsig3_out) const
{
   // Tensile yield: project sig3 back to tension_cutoff
   double ft = sig3 - tension_cutoff;
   
   // For associated flow in tensile mode
   double beta = ft / (lambda + 2.0 * mu);
   dsig3_out = tension_cutoff;
   
   return beta;
}

} // namespace geodynamics
} // namespace mfem
