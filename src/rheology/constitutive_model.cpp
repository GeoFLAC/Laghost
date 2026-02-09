// Copyright (c) 2024, Laghost Authors
// SPDX-License-Identifier: BSD-3-Clause
//
// Constitutive Model - Default CPU implementations

#include "constitutive_model.hpp"

namespace mfem
{
namespace geodynamics
{

void ConstitutiveModel::UpdateBatch(
   const Vector &strain_rate_flat,
   const Vector &spin_flat,
   Vector &stress_flat,
   Vector &alpha,
   const Vector &props_flat,
   const Array<int> &mat_ids,
   double dt,
   int nqp,
   int dim)
{
   // Default CPU implementation: loop over quadrature points
   const int stress_size = (dim == 2) ? 3 : 6;  // Voigt notation
   const int spin_size = (dim == 2) ? 1 : 3;
   
   // Number of material properties per qp
   const int props_per_qp = sizeof(MaterialProperties) / sizeof(double);
   
   DenseMatrix strain_rate(dim), spin_mat(dim), stress_mat(dim);
   
   for (int i = 0; i < nqp; i++)
   {
      // Unpack strain rate from Voigt to tensor
      if (dim == 2)
      {
         strain_rate(0,0) = strain_rate_flat[i*stress_size + 0];
         strain_rate(1,1) = strain_rate_flat[i*stress_size + 1];
         strain_rate(0,1) = strain_rate_flat[i*stress_size + 2] * 0.5;
         strain_rate(1,0) = strain_rate(0,1);
      }
      else // dim == 3
      {
         strain_rate(0,0) = strain_rate_flat[i*stress_size + 0];
         strain_rate(1,1) = strain_rate_flat[i*stress_size + 1];
         strain_rate(2,2) = strain_rate_flat[i*stress_size + 2];
         strain_rate(0,1) = strain_rate_flat[i*stress_size + 3] * 0.5;
         strain_rate(1,0) = strain_rate(0,1);
         strain_rate(0,2) = strain_rate_flat[i*stress_size + 4] * 0.5;
         strain_rate(2,0) = strain_rate(0,2);
         strain_rate(1,2) = strain_rate_flat[i*stress_size + 5] * 0.5;
         strain_rate(2,1) = strain_rate(1,2);
      }
      
      // Unpack spin
      spin_mat = 0.0;
      if (dim == 2)
      {
         spin_mat(0,1) = spin_flat[i];
         spin_mat(1,0) = -spin_flat[i];
      }
      else
      {
         spin_mat(0,1) = spin_flat[i*spin_size + 0];
         spin_mat(1,0) = -spin_mat(0,1);
         spin_mat(0,2) = spin_flat[i*spin_size + 1];
         spin_mat(2,0) = -spin_mat(0,2);
         spin_mat(1,2) = spin_flat[i*spin_size + 2];
         spin_mat(2,1) = -spin_mat(1,2);
      }
      
      // Unpack stress
      if (dim == 2)
      {
         stress_mat(0,0) = stress_flat[i*stress_size + 0];
         stress_mat(1,1) = stress_flat[i*stress_size + 1];
         stress_mat(0,1) = stress_flat[i*stress_size + 2];
         stress_mat(1,0) = stress_mat(0,1);
      }
      else
      {
         stress_mat(0,0) = stress_flat[i*stress_size + 0];
         stress_mat(1,1) = stress_flat[i*stress_size + 1];
         stress_mat(2,2) = stress_flat[i*stress_size + 2];
         stress_mat(0,1) = stress_flat[i*stress_size + 3];
         stress_mat(1,0) = stress_mat(0,1);
         stress_mat(0,2) = stress_flat[i*stress_size + 4];
         stress_mat(2,0) = stress_mat(0,2);
         stress_mat(1,2) = stress_flat[i*stress_size + 5];
         stress_mat(2,1) = stress_mat(1,2);
      }
      
      // Build material properties struct
      MaterialProperties props;
      const double *p = props_flat.GetData() + i * props_per_qp;
      props.rho = p[0];
      props.lambda = p[1];
      props.mu = p[2];
      props.tension_cutoff = p[3];
      props.cohesion0 = p[4];
      props.cohesion1 = p[5];
      props.friction_angle0 = p[6];
      props.friction_angle1 = p[7];
      props.dilation_angle0 = p[8];
      props.dilation_angle1 = p[9];
      props.alpha0 = p[10];
      props.alpha1 = p[11];
      props.viscosity = p[12];
      
      // Build state
      MaterialState state(dim);
      state.stress = stress_mat;
      state.alpha = alpha[i];
      state.material_id = mat_ids[i];
      
      // Call single-point update
      Update(strain_rate, spin_mat, state, props, dt);
      
      // Pack stress back
      if (dim == 2)
      {
         stress_flat[i*stress_size + 0] = state.stress(0,0);
         stress_flat[i*stress_size + 1] = state.stress(1,1);
         stress_flat[i*stress_size + 2] = state.stress(0,1);
      }
      else
      {
         stress_flat[i*stress_size + 0] = state.stress(0,0);
         stress_flat[i*stress_size + 1] = state.stress(1,1);
         stress_flat[i*stress_size + 2] = state.stress(2,2);
         stress_flat[i*stress_size + 3] = state.stress(0,1);
         stress_flat[i*stress_size + 4] = state.stress(0,2);
         stress_flat[i*stress_size + 5] = state.stress(1,2);
      }
      
      alpha[i] = state.alpha;
   }
}

void ConstitutiveModel::UpdateDevice(
   const double *d_strain_rate,
   const double *d_spin,
   double *d_stress,
   double *d_alpha,
   const double *d_props,
   const int *d_mat_ids,
   double dt,
   int nqp,
   int dim)
{
   // Default: error - subclass must implement for GPU
   MFEM_ABORT("GPU kernel not implemented for this constitutive model. "
              "Override UpdateDevice() for GPU support.");
}

} // namespace geodynamics
} // namespace mfem
