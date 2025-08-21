// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

#ifndef MFEM_LAGHOST_SOLVER
#define MFEM_LAGHOST_SOLVER

#include "mfem.hpp"
#include "laghost_assembly.hpp"
#include "laghost_parameters.hpp"

#ifdef MFEM_USE_MPI

namespace mfem
{

namespace geodynamics
{

/// Visualize the given parallel grid function, using a GLVis server on the
/// specified host and port. Set the visualization window title, and optionally,
/// its geometry.
void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    ParGridFunction &gf, const char *title,
                    int x = 0, int y = 0, int w = 400, int h = 400,
                    bool vec = false);

struct TimingData
{
   // Total times for all major computations:
   // CG solves (H1 and L2) / force RHS assemblies / quadrature computations.
   StopWatch sw_cgH1, sw_cgL2, sw_force, sw_qdata;

   // Store the number of dofs of the corresponding local CG
   const HYPRE_Int L2dof;

   // These accumulate the total processed dofs or quad points:
   // #(CG iterations) for the L2 CG solve.
   // #quads * #(RK sub steps) for the quadrature data computations.
   HYPRE_Int H1iter, L2iter;
   HYPRE_Int quad_tstep;

   TimingData(const HYPRE_Int l2d) :
      L2dof(l2d), H1iter(0), L2iter(0), quad_tstep(0) { }
};

class QUpdate
{
private:
   const int dim, vdim, NQ, NE, Q1D;
   const bool use_viscosity, use_vorticity;
   const double cfl;
   double max_vel_q;
   double mscale_q; 
   double gravity_q;
   TimingData *timer;
   const IntegrationRule &ir;
   ParFiniteElementSpace &H1, &L2, &L2_stress;
   const Operator *H1R;
   Vector q_dt_est, q_e, e_vec, q_dx, q_dv, q_sig, q_h_est;
   const QuadratureInterpolator *q1,*q2,*q3;
   const ParGridFunction &gamma_gf;
   const ParGridFunction &lambda_gf;
   const ParGridFunction &mu_gf;
public:
   QUpdate(const int d, const int ne, const int q1d,
           const bool visc, const bool vort,
           const double cfl, TimingData *t,
           const ParGridFunction &gamma_gf,
           const ParGridFunction &lambda_gf,
           const ParGridFunction &mu_gf,
           const IntegrationRule &ir,
           ParFiniteElementSpace &h1, ParFiniteElementSpace &l2, ParFiniteElementSpace &l2_stress):
      dim(d), vdim(h1.GetVDim()),
      NQ(ir.GetNPoints()), NE(ne), Q1D(q1d),
      use_viscosity(visc), use_vorticity(vort), cfl(cfl),
      timer(t), ir(ir), H1(h1), L2(l2), L2_stress(l2_stress),
      H1R(H1.GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC)),
      q_dt_est(NE*NQ),
      q_h_est(NE*NQ),
      q_e(NE*NQ),
      q_sig(NE*NQ*3*(vdim-1)),
      e_vec(NQ*NE*vdim),
      q_dx(NQ*NE*vdim*vdim),
      q_dv(NQ*NE*vdim*vdim),
      q1(H1.GetQuadratureInterpolator(ir)),
      q2(L2.GetQuadratureInterpolator(ir)),
      q3(L2_stress.GetQuadratureInterpolator(ir)),
      gamma_gf(gamma_gf),
      lambda_gf(lambda_gf),
      mu_gf(mu_gf) { }
   
   // QUpdate(const int d, const int ne, const int q1d,
   //         const bool visc, const bool vort,
   //         const double dt,
   //         const double cfl, TimingData *t,
   //         const ParGridFunction &gamma_gf,
   //         const ParGridFunction &lambda_gf,
   //         const ParGridFunction &mu_gf,
   //         const IntegrationRule &ir,
   //         ParFiniteElementSpace &h1, ParFiniteElementSpace &l2, ParFiniteElementSpace &l2_stress):
   //    dim(d), vdim(h1.GetVDim()),
   //    NQ(ir.GetNPoints()), NE(ne), Q1D(q1d),
   //    use_viscosity(visc), use_vorticity(vort), cfl(cfl),
   //    timer(t), ir(ir), H1(h1), L2(l2), L2_stress(l2_stress),
   //    H1R(H1.GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC)),
   //    q_dt_est(NE*NQ),
   //    q_h_est(NE*NQ),
   //    q_e(NE*NQ),
   //    q_sig(NE*NQ*3*(vdim-1)),
   //    e_vec(NQ*NE*vdim),
   //    q_dx(NQ*NE*vdim*vdim),
   //    q_dv(NQ*NE*vdim*vdim),
   //    q1(H1.GetQuadratureInterpolator(ir)),
   //    q2(L2.GetQuadratureInterpolator(ir)),
   //    q3(L2_stress.GetQuadratureInterpolator(ir)),
   //    gamma_gf(gamma_gf),
   //    lambda_gf(lambda_gf),
   //    mu_gf(mu_gf) { }

   void UpdateQuadratureData(const Vector &S, QuadratureData &qdata);
   // void UpdateQuadratureData(const Vector &S, QuadratureData &qdata, const double dt);
};

// Given a solutions state (x, v, e), this class performs all necessary
// computations to evaluate the new slopes (dx_dt, dv_dt, de_dt).
class LagrangianGeoOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace &H1, &L2, &L2_stress;
   mutable ParFiniteElementSpace H1c;
   ParMesh *pmesh;
   // FE spaces local and global sizes
   // const int H1Vsize;
   // const int H1TVSize;
   // const HYPRE_Int H1GTVSize;
   // const int L2Vsize;
   // const int L2TVSize;
   // const HYPRE_Int L2GTVSize;
   mutable int H1Vsize;
   mutable int H1TVSize;
   mutable HYPRE_Int H1GTVSize;
   mutable int L2Vsize;
   mutable int L2TVSize;
   mutable HYPRE_Int L2GTVSize;
   mutable Array<int> block_offsets;
   // Reference to the current mesh configuration.
   mutable ParGridFunction x_gf;
   // ParGridFunction x0_gf; // copy of initial mesh position

   ParGridFunction &rho0_gf;
   ParGridFunction &fictitious_rho0_gf;
   // ParGridFunction x0_gf; // copy of initial mesh position

   const Array<int> &ess_tdofs;
   // const int dim, NE, l2dofs_cnt, l2_stress_dofs_cnt, h1dofs_cnt, source_type;
   const int dim, l2dofs_cnt, l2_stress_dofs_cnt, h1dofs_cnt, source_type;
   mutable int NE;
   const double cfl;
   const bool use_viscosity, use_vorticity, p_assembly, winkler_foundation, dyn_damping;
   const double cg_rel_tol;
   mutable double mass_scale, grav_mag, thickness, winkler_rho, dyn_factor, vbc_max_val;
   const int cg_max_iter;
   const double ftz_tol;
   ParGridFunction &gamma_gf;
   ParGridFunction &lambda_gf;
   ParGridFunction &mu_gf;
   Vector bc_id_pa;
   // mutable Vector tension_cutoff, cohesion, friction_angle, dilation_angle;
   // Velocity mass matrix and local inverses of the energy mass matrices. These
   // are constant in time, due to the pointwise mass conservation property.
   mutable ParBilinearForm Mv;
   mutable SparseMatrix Mv_spmat_copy;
   mutable ParBilinearForm fic_Mv;
   mutable SparseMatrix fic_Mv_spmat_copy;

   mutable DenseTensor Me, Me_inv;

   // 
   GridFunctionCoefficient rho0_coeff; // TODO: remove when Mv update improved
   GridFunctionCoefficient scale_rho0_coeff; // TODO: remove when fic_Mv update improved
   
   // DenseTensor Me, Me_inv;
   // Integration rule for all assemblies.
   const IntegrationRule &ir;
   // Data associated with each quadrature point in the mesh.
   // These values are recomputed at each time step.
   const int Q1D;
   mutable QuadratureData qdata;
   mutable bool qdata_is_current, forcemat_is_assembled;
   mutable bool gmat_is_assembled;
   // Force matrix that combines the kinematic and thermodynamic spaces. It is
   // assembled in each time step and then it is used to compute the final
   // right-hand sides for momentum and specific internal energy.
   mutable MixedBilinearForm Force;
   // G matrix is in thermodynamic spaces.
   
   mutable LinearForm Body_Force;
   // mutable LinearForm *Body_Force;
   // Same as above, but done through partial assembly.
   ForcePAOperator *ForcePA;
   StressPAOperator *StressPA; // partial assembly for stress rate, slee
   // Mass matrices done through partial assembly:
   // velocity (coupled H1 assembly) and energy (local L2 assemblies).
   MassPAOperator *VMassPA, *EMassPA;
   OperatorJacobiSmoother *VMassPA_Jprec;
   // Linear solver for energy.
   CGSolver CG_VMass, CG_EMass;

   mutable Vector zone_max_visc, zone_vgrad;

   mutable TimingData timer;
   mutable QUpdate *qupdate;
   mutable Vector X, B, one, rhs, e_rhs, s_rhs;
   // mutable Vector e_rhs;
   mutable ParGridFunction rhs_c_gf, dvc_gf;
   mutable Array<int> c_tdofs[3];

   virtual void ComputeMaterialProperties(int nvalues, const double gamma[],
                                          const double rho[], const double e[],
                                          double p[], double cs[], double pmod[], double mscale) const
   {
      for (int v = 0; v < nvalues; v++)
      {
         p[v]  = rho[v] * e[v];
         // cs[v] = sqrt(pmod[v]/(rho[v]*mscale));
         cs[v] = sqrt(pmod[v]/(rho[v]));
         // cs[v] = sqrt(pmod[v]/rho[v]);
         // p[v]  = (gamma[v] - 1.0) * rho[v] * e[v];
         // cs[v] = sqrt(gamma[v] * (gamma[v]-1.0) * e[v]);
      }
   }

   void UpdateQuadratureData(const Vector &S) const;
   // void UpdateQuadratureData(const Vector &S, const double dt) const;
   void AssembleForceMatrix() const;
   void AssembleSigmaMatrix() const;

public:
   LagrangianGeoOperator(  const int size,
                           ParFiniteElementSpace &h1_fes,
                           ParFiniteElementSpace &l2_fes,
                           ParFiniteElementSpace &l2_stress_fes,
                           const Array<int> &ess_tdofs,
                           // Coefficient &rho0_coeff,
                           // Coefficient &scale_rho0_coeff,
                           ParGridFunction &rho0_gf,
                           ParGridFunction &fictitious_rho0_gf,
                           ParGridFunction &gamma_gf,
                           const int source,
                           const bool visc, const bool vort, 
                           ParGridFunction &lambda_gf, ParGridFunction &mu_gf,
                           const Param &param, 
                           const double _vbc_max_val);
   ~LagrangianGeoOperator();

   // Solve for dx_dt, dv_dt and de_dt.
   virtual void Mult(const Vector &S, Vector &dS_dt) const;
   // virtual void Mult(const Vector &S, Vector &dS_dt, const double dt) const;

   virtual MemoryClass GetMemoryClass() const
   { return Device::GetMemoryClass(); }

   void SolveVelocity(const Vector &S, Vector &dS_dt) const;
   void SolveEnergy(const Vector &S, const Vector &v, Vector &dS_dt) const;
   void SolveStress(const Vector &S, Vector &dS_dt) const;

   // void SolveVelocity(const Vector &S, Vector &dS_dt, const double dt) const;
   // void SolveEnergy(const Vector &S, const Vector &v, Vector &dS_dt, const double dt) const;
   // void SolveStress(const Vector &S, Vector &dS_dt, const double dt) const;

   // void RadialReturn(const Vector &S, Vector &dS_dt, const double dt) const;
   void UpdateMesh(const Vector &S) const;
   // void test_function(const Vector &S, Vector &_test) const;
   void Getdamping(const Vector &S, Vector &_v_damping) const;
   void Getdamping_comp(const Vector &S, const int &comp, Vector &_v_damping) const;
   void Winkler(const Vector &S, Vector &_winkler, double &_thickness) const;
   
   // Calls UpdateQuadratureData to compute the new qdata.dt_estimate.
   // double GetTimeStepEstimate(const Vector &S) const;
   double GetTimeStepEstimate(const Vector &S) const;
   double GetLengthEstimate(const Vector &S) const;
   // double GetTimeStepEstimate(const Vector &S, const double dt) const;
   // double GetLengthEstimate(const Vector &S, const double dt) const;
   // double GetTimeStepEstimate(const Vector &S, const double dt, bool IamRoot) const;
   void ResetTimeStepEstimate() const;
   void ResetQuadratureData() const { qdata_is_current = false; }

   // The density values, which are stored only at some quadrature points,
   // are projected as a ParGridFunction.
   void ComputeDensity(ParGridFunction &rho) const;
   double InternalEnergy(const ParGridFunction &e) const;
   double KineticEnergy(const ParGridFunction &v) const;

   int GetH1VSize() const { return H1.GetVSize(); }
   const Array<int> &GetBlockOffsets() const 
   { 
      return block_offsets; 
   }

   // Update all internal data on mesh change.
   void TMOPUpdate(const Vector &S, bool quick);
   void GetErrorEstimates(ParGridFunction &e_gf, Vector & errors);

   // void TMOPUpdate(const Vector &S, Coefficient &rho0_coeff, bool quick);

   void SetH0(double h0) { qdata.h0 = h0; }
   double GetH0() const { return qdata.h0; }

   Vector& GetZoneMaxVisc() { return zone_max_visc; }
   Vector& GetZoneVGrad() { return zone_vgrad; }

   void PrintTimingData(bool IamRoot, int steps, const bool fom) const;
};

// TaylorCoefficient used in the 2D Taylor-Green problem.
class TaylorCoefficient : public Coefficient
{
public:
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
   {
      Vector x(2);
      T.Transform(ip, x);
      return 3.0 / 8.0 * M_PI * ( cos(3.0*M_PI*x(0)) * cos(M_PI*x(1)) -
                                  cos(M_PI*x(0))     * cos(3.0*M_PI*x(1)) );
   }
};

// Acceleration source coefficient used in the 2D Rayleigh-Taylor problem.
class GTCoefficient : public VectorCoefficient
{
public:
   GTCoefficient(int dim) : VectorCoefficient(dim) { }
   using VectorCoefficient::Eval;
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip)
   {
      if(V.Size() == 2)
      {
         V = 0.0; V(1) = -1.0;
      }
      else if(V.Size() == 3)
      {
         V = 0.0; V(2) = -1.0; 
      }
   }
};

// gravity coefficient
class DampCoefficient : public VectorCoefficient
{
public:
   DampCoefficient(int dim) : VectorCoefficient(dim) { }
   using VectorCoefficient::Eval;
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip)
   {
      if(V.Size() == 2)
      {
         V = 0.0; 
         V(0) = -1.0*copysignl(1.0, V(0));
         V(1) = -1.0*copysignl(1.0, V(1));
      }
      else if(V.Size() == 3)
      {
         V = 0.0; 
         V(0) = -1.0*copysignl(1.0, V(0));
         V(1) = -1.0*copysignl(1.0, V(1)); 
         V(2) = -1.0*copysignl(1.0, V(2));
      }
   }
};

} // namespace geodynamics

class GeoODESolver : public ODESolver
{
protected:
   geodynamics::LagrangianGeoOperator *geo_oper;
public:
   GeoODESolver() : geo_oper(NULL) { }
   virtual void Init(TimeDependentOperator&);
   virtual void Step(Vector&, double&, double&)
   { MFEM_ABORT("Time stepping is undefined."); }
};

class RK2AvgSolver : public GeoODESolver
{
protected:
   Vector V;
   BlockVector dS_dt, S0;
public:
   RK2AvgSolver() { }
   virtual void Init(TimeDependentOperator &_f);
   virtual void Step(Vector &S, double &t, double &dt);
};

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_LAGHOS
