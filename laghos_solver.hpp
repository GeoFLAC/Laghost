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

#ifndef MFEM_LAGHOS_SOLVER
#define MFEM_LAGHOS_SOLVER

#include "mfem.hpp"
#include "laghos_assembly.hpp"

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
   TimingData *timer;
   const IntegrationRule &ir;
   ParFiniteElementSpace &H1, &L2;
   const Operator *H1R;
   Vector q_dt_est, q_e, e_vec, q_dx, q_dv;
   const QuadratureInterpolator *q1,*q2;
   const ParGridFunction &gamma_gf;
public:
   QUpdate(const int d, const int ne, const int q1d,
           const bool visc, const bool vort,
           const double cfl, TimingData *t,
           const ParGridFunction &gamma_gf,
           const IntegrationRule &ir,
           ParFiniteElementSpace &h1, ParFiniteElementSpace &l2):
      dim(d), vdim(h1.GetVDim()),
      NQ(ir.GetNPoints()), NE(ne), Q1D(q1d),
      use_viscosity(visc), use_vorticity(vort), cfl(cfl),
      timer(t), ir(ir), H1(h1), L2(l2),
      H1R(H1.GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC)),
      q_dt_est(NE*NQ),
      q_e(NE*NQ),
      e_vec(NQ*NE*vdim),
      q_dx(NQ*NE*vdim*vdim),
      q_dv(NQ*NE*vdim*vdim),
      q1(H1.GetQuadratureInterpolator(ir)),
      q2(L2.GetQuadratureInterpolator(ir)),
      gamma_gf(gamma_gf) { }
   
   QUpdate(const int d, const int ne, const int q1d,
           const bool visc, const bool vort,
           const double dt,
           const double cfl, TimingData *t,
           const ParGridFunction &gamma_gf,
           const IntegrationRule &ir,
           ParFiniteElementSpace &h1, ParFiniteElementSpace &l2):
      dim(d), vdim(h1.GetVDim()),
      NQ(ir.GetNPoints()), NE(ne), Q1D(q1d),
      use_viscosity(visc), use_vorticity(vort), cfl(cfl),
      timer(t), ir(ir), H1(h1), L2(l2),
      H1R(H1.GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC)),
      q_dt_est(NE*NQ),
      q_e(NE*NQ),
      e_vec(NQ*NE*vdim),
      q_dx(NQ*NE*vdim*vdim),
      q_dv(NQ*NE*vdim*vdim),
      q1(H1.GetQuadratureInterpolator(ir)),
      q2(L2.GetQuadratureInterpolator(ir)),
      gamma_gf(gamma_gf) { }

   void UpdateQuadratureData(const Vector &S, QuadratureData &qdata);
   void UpdateQuadratureData(const Vector &S, QuadratureData &qdata, const double dt);
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
   const int H1Vsize;
   const int H1TVSize;
   const HYPRE_Int H1GTVSize;
   const int L2Vsize;
   const int L2TVSize;
   const HYPRE_Int L2GTVSize;
   Array<int> block_offsets;
   // Reference to the current mesh configuration.
   mutable ParGridFunction x_gf;
   const Array<int> &ess_tdofs;
   const int dim, NE, l2dofs_cnt, l2_stress_dofs_cnt, h1dofs_cnt, source_type;
   const double cfl;
   const bool use_viscosity, use_vorticity, p_assembly, winkler_foundation;
   const double cg_rel_tol;
   mutable double mass_scale, grav_mag, thickness;
   const int cg_max_iter;
   const double ftz_tol;
   const ParGridFunction &gamma_gf;
   const ParGridFunction &lambda_gf;
   const ParGridFunction &mu_gf;
   mutable Vector tension_cutoff, cohesion, friction_angle, dilation_angle;
   // Velocity mass matrix and local inverses of the energy mass matrices. These
   // are constant in time, due to the pointwise mass conservation property.
   mutable ParBilinearForm Mv;
   SparseMatrix Mv_spmat_copy;
   DenseTensor Me, Me_inv;
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
   // Same as above, but done through partial assembly.
   ForcePAOperator *ForcePA;
   ForcePAOperator *SigmaPA;
   // Mass matrices done through partial assembly:
   // velocity (coupled H1 assembly) and energy (local L2 assemblies).
   MassPAOperator *VMassPA, *EMassPA;
   OperatorJacobiSmoother *VMassPA_Jprec;
   // Linear solver for energy.
   CGSolver CG_VMass, CG_EMass;
   mutable TimingData timer;
   mutable QUpdate *qupdate;
   mutable Vector X, B, one, rhs, e_rhs;
   mutable ParGridFunction rhs_c_gf, dvc_gf;
   mutable Array<int> c_tdofs[3];

   virtual void ComputeMaterialProperties(int nvalues, const double gamma[],
                                          const double rho[], const double e[],
                                          double p[], double cs[], double pmod[], double mscale) const
   {
      for (int v = 0; v < nvalues; v++)
      {
         p[v]  = rho[v] * e[v];
         cs[v] = sqrt(pmod[v]/(rho[v]*mscale));
         // cs[v] = sqrt(pmod[v]/rho[v]);
         // p[v]  = (gamma[v] - 1.0) * rho[v] * e[v];
         // cs[v] = sqrt(gamma[v] * (gamma[v]-1.0) * e[v]);
      }
   }

   void UpdateQuadratureData(const Vector &S) const;
   void UpdateQuadratureData(const Vector &S, const double dt) const;
   void AssembleForceMatrix() const;
   void AssembleSigmaMatrix() const;

public:
   LagrangianGeoOperator(const int size,
                           ParFiniteElementSpace &h1_fes,
                           ParFiniteElementSpace &l2_fes,
                           ParFiniteElementSpace &l2_stress_fes,
                           const Array<int> &ess_tdofs,
                           Coefficient &rho0_coeff,
                           Coefficient &scale_rho0_coeff,
                           ParGridFunction &rho0_gf,
                           ParGridFunction &gamma_gf,
                           const int source,
                           const double cfl,
                           const bool visc, const bool vort, const bool pa,
                           const double cgt, const int cgiter, double ftz_tol,
                           const int order_q,
                           ParGridFunction &lambda_gf, ParGridFunction &mu_gf, double mscale, const double gravity, const double _thickness,
                           Vector _lambda, Vector _mu, Vector _tension_cutoff, Vector _cohesion, Vector _friction_angle, Vector _dilation_angle, const bool winkler);
   ~LagrangianGeoOperator();


   // Solve for dx_dt, dv_dt and de_dt.
   // virtual void Mult(const Vector &S, Vector &dS_dt) const;
   virtual void Mult(const Vector &S, Vector &dS_dt, const double dt) const;

   virtual MemoryClass GetMemoryClass() const
   { return Device::GetMemoryClass(); }

   // void SolveVelocity(const Vector &S, Vector &dS_dt) const;
   // void SolveEnergy(const Vector &S, const Vector &v, Vector &dS_dt) const;
   
   void SolveVelocity(const Vector &S, Vector &dS_dt, const double dt) const;
   // void SolveVelocity(const Vector &S, Vector &dS_dt, const double dt, const double t) const;
   void SolveEnergy(const Vector &S, const Vector &v, Vector &dS_dt, const double dt) const;
   void SolveStress(const Vector &S, Vector &dS_dt, const double dt) const;
   // void RadialReturn(const Vector &S, Vector &dS_dt, const double dt) const;
   void UpdateMesh(const Vector &S) const;
   // void test_function(const Vector &S, Vector &_test) const;
   void Getdamping(const Vector &S, Vector &_v_damping) const;
   void Winkler(const Vector &S, Vector &_winkler, double &_thickness) const;
   
   // Calls UpdateQuadratureData to compute the new qdata.dt_estimate.
   // double GetTimeStepEstimate(const Vector &S) const;
   double GetTimeStepEstimate(const Vector &S, const double dt) const;
   double GetTimeStepEstimate(const Vector &S, const double dt, bool IamRoot) const;
   void ResetTimeStepEstimate() const;
   void ResetQuadratureData() const { qdata_is_current = false; }

   // The density values, which are stored only at some quadrature points,
   // are projected as a ParGridFunction.
   void ComputeDensity(ParGridFunction &rho) const;
   double InternalEnergy(const ParGridFunction &e) const;
   double KineticEnergy(const ParGridFunction &v) const;

   int GetH1VSize() const { return H1.GetVSize(); }
   const Array<int> &GetBlockOffsets() const { return block_offsets; }

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
class RTCoefficient : public VectorCoefficient
{
public:
   RTCoefficient(int dim) : VectorCoefficient(dim) { }
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
