
#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "remhos_ho.hpp"
#include "remhos_lo.hpp"
#include "remhos_fct.hpp"
#include "remhos_mono.hpp"
#include "remhos_tools.hpp"
#include "remhos_sync.hpp"

// #include <unistd.h>

using namespace std;
using namespace mfem;
using std::cout;
using std::endl;

enum class HOSolverType {None, Neumann, CG, LocalInverse};
enum class FCTSolverType {None, FluxBased, ClipScale,
                          NonlinearPenalty, FCTProject};
enum class LOSolverType {None,    DiscrUpwind,    DiscrUpwindPrec,
                         ResDist, ResDistSubcell, MassBased};
enum class MonolithicSolverType {None, ResDistMono, ResDistMonoSubcell};

enum class TimeStepControl {FixedTimeStep, LOBoundsError};

// Choice for the problem setup. The fluid velocity, initial condition and
// inflow boundary condition are chosen based on this parameter.
int problem_num;

// 0 is standard transport.
// 1 is standard remap (mesh moves, solution is fixed).
int exec_mode;

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v);

// Initial condition
double u0_function(const Vector &x);
double s0_function(const Vector &x);

// Inflow boundary condition
double inflow_function(const Vector &x);

// Mesh bounding box
Vector bb_min, bb_max;


// void attach_debugger(bool condition) {
//     if (!condition) return;
//     bool attached = false;

//     // also write PID to a file
//     ofstream os("/tmp/mpi_debug.pid");
//     os << getpid() << endl;
//     os.close();

//     cout << "Waiting for debugger to be attached, PID: "
//         << getpid() << endl;
//     while (!attached) sleep(1);
// }

class AdvectionOperator : public TimeDependentOperator
{
private:
   BilinearForm &Mbf, &ml;
   ParBilinearForm &Kbf;
   ParBilinearForm &M_HO, &K_HO;
   Vector &lumpedM;

   Vector start_mesh_pos, start_submesh_pos;
   GridFunction &mesh_pos, *submesh_pos, &mesh_vel, &submesh_vel;

   mutable ParGridFunction x_gf;

   double dt;
   TimeStepControl dt_control;
   mutable double dt_est;
   Assembly &asmbl;

   LowOrderMethod &lom;
   DofInfo &dofs;

   HOSolver *ho_solver;
   LOSolver *lo_solver;
   FCTSolver *fct_solver;
   MonolithicSolver *mono_solver;

   void UpdateTimeStepEstimate(const Vector &x, const Vector &dx,
                               const Vector &x_min, const Vector &x_max) const;

public:
   AdvectionOperator(int size, BilinearForm &Mbf_, BilinearForm &_ml,
                     Vector &_lumpedM,
                     ParBilinearForm &Kbf_,
                     ParBilinearForm &M_HO_, ParBilinearForm &K_HO_,
                     GridFunction &pos, GridFunction *sub_pos,
                     GridFunction &vel, GridFunction &sub_vel,
                     Assembly &_asmbl, LowOrderMethod &_lom, DofInfo &_dofs,
                     HOSolver *hos, LOSolver *los, FCTSolver *fct,
                     MonolithicSolver *mos);

   virtual void Mult(const Vector &x, Vector &y) const;

   void SetTimeStepControl(TimeStepControl tsc)
   {
      if (tsc == TimeStepControl::LOBoundsError)
      {
         MFEM_VERIFY(lo_solver,
                     "The selected time step control requires a LO solver.");
      }
      dt_control = tsc;
   }
   void SetDt(double _dt) { dt = _dt; dt_est = dt; }
   double GetTimeStepEstimate() { return dt_est; }

   void SetRemapStartPos(const Vector &m_pos, const Vector &sm_pos)
   {
      start_mesh_pos    = m_pos;
      start_submesh_pos = sm_pos;
   }

   virtual ~AdvectionOperator() { }
};


#ifdef MFEM_USE_MPI

void Remapping(ParMesh *pmesh, ParGridFunction &x, ParGridFunction &x_mod, ParGridFunction &u_gf, int &mesh_order, int &order, bool &pa, bool &ncmesh) 
{
   MPI_Comm comm = pmesh->GetComm();
   int num_procs, myid;
   MPI_Comm_size(comm, &num_procs);
   MPI_Comm_rank(comm, &myid);

   int ode_solver_type = 3;

   HOSolverType ho_type           = HOSolverType::LocalInverse;
   LOSolverType lo_type           = LOSolverType::DiscrUpwind;
   FCTSolverType fct_type         = FCTSolverType::FluxBased;
   MonolithicSolverType mono_type = MonolithicSolverType::None;

   // HOSolverType ho_type           = HOSolverType::LocalInverse;
   // LOSolverType lo_type           = LOSolverType::MassBased;
   // FCTSolverType fct_type         = FCTSolverType::FluxBased;
   // MonolithicSolverType mono_type = MonolithicSolverType::None;

   // enum class HOSolverType {None, Neumann, CG, LocalInverse};
   // enum class FCTSolverType {None, FluxBased, ClipScale,
   //                        NonlinearPenalty, FCTProject};
   // enum class LOSolverType {None,    DiscrUpwind,    DiscrUpwindPrec,
   //                       ResDist, ResDistSubcell, MassBased};
   // enum class MonolithicSolverType {None, ResDistMono, ResDistMonoSubcell};

   // enum class TimeStepControl {FixedTimeStep, LOBoundsError};

   if(ncmesh)
   {
      if(myid==0){std::cout << "ncmesh" << std::endl;}
      ho_type           = HOSolverType::LocalInverse;
      lo_type           = LOSolverType::None;
      fct_type          = FCTSolverType::None;
      // lo_type           = LOSolverType::MassBased;
      // fct_type          = FCTSolverType::FluxBased;
   }
   

   int bounds_type = 0;
   // bool pa = false;
   pa = false;
   bool next_gen_full = false;
   int smth_ind_type = 0;
   double t_final = 1.0;
   TimeStepControl dt_control = TimeStepControl::FixedTimeStep;
   // double dt = 0.01;
   double dt = 0.005;
   bool visualization = false;
   bool visit = false;
   bool verify_bounds = false;
   bool product_sync = false;
   int vis_steps = 100;
   // const char *device_config = "cpu";
   int precision = 16;
   int dim = pmesh->Dimension();
   // cout.precision(precision);

   exec_mode = 1; // Always standard remap
   problem_num = 10;
   // if (myid == 0) {std::cout << "standard remap (mesh moves, solution is fixed)" <<std::endl;}
  
   // Only standard assembly in 1D (some mfem functions just abort in 1D).
   if ((pa || next_gen_full) && dim == 1)
   {
      MFEM_WARNING("Disabling PA / FA for 1D.");
      pa = false;
      next_gen_full = false;
   }

   const int NE  = pmesh->GetNE();

   // Define the ODE solver used for time integration. Several explicit
   // Runge-Kutta methods are available.
   ODESolver *ode_solver = NULL;
   switch (ode_solver_type)
   {
      case 1: ode_solver = new ForwardEulerSolver; break;
      case 2: ode_solver = new RK2Solver(1.0); break;
      case 3: ode_solver = new RK3SSPSolver; break;
      case 4:
         if (myid == 0) { MFEM_WARNING("RK4 may violate the bounds."); }
         ode_solver = new RK4Solver; break;
      case 6:
         if (myid == 0) { MFEM_WARNING("RK6 may violate the bounds."); }
         ode_solver = new RK6Solver; break;
      default:
         cout << "Unknown ODE solver type for remapping: " << ode_solver_type << '\n';
         delete pmesh; MPI_Finalize();
   }

   // Check if the input mesh is periodic.
   const bool periodic = pmesh->GetNodes() != NULL &&
                         dynamic_cast<const L2_FECollection *>
                         (pmesh->GetNodes()->FESpace()->FEColl()) != NULL;
   if(ncmesh){pmesh->EnsureNCMesh(true);}
   pmesh->SetCurvature(mesh_order, periodic);

   FiniteElementCollection *mesh_fec;
   if (periodic)
   {
      if (myid == 0) {std::cout << "input mesh is periodic" << std::endl;}
      mesh_fec = new L2_FECollection(mesh_order, dim, BasisType::GaussLobatto);
   }
   else
   {
      if (myid == 0) {std::cout << "input mesh is not periodic" << std::endl;}
      mesh_fec = new H1_FECollection(mesh_order, dim, BasisType::GaussLobatto);
   }

   // Current mesh positions.
   ParFiniteElementSpace mesh_pfes(pmesh, mesh_fec, dim);
   // pmesh->GetBoundingBox(bb_min, bb_max, max(order, 1));

   // std::cout << bb_min.Size() << "," << bb_min[0] << "," << bb_min[1] << std::endl;
   // std::cout << bb_max.Size() << "," << bb_max[0] << "," << bb_max[1] << std::endl;
   bb_min.SetSize(dim);
   bb_max.SetSize(dim);
   bb_min[0] =0.0; bb_min[1] =0.0;
   bb_max[0] =1e3; bb_max[1] =1e3;

   
   // Store initial mesh positions.
   Vector x0(x.Size());
   x0 = x;

   // Velocity for the problem. Depending on the execution mode, this is the
   // advective velocity (transport) or mesh velocity (remap).
   VectorFunctionCoefficient velocity(dim, velocity_function);
   
   // // Mesh velocity.
   // // If remap is on, obtain the mesh velocity by moving the mesh to the final
   // // mesh positions, and taking the displacement vector.
   // // The mesh motion resembles a time-dependent deformation, e.g., similar to
   // // a deformation that is obtained by a Lagrangian simulation.
   // GridFunction v_gf(x.FESpace());
   // v_gf = 0.0;
   // VectorGridFunctionCoefficient v_mesh_coeff(&v_gf);
   // if (exec_mode == 1)
   // {
   //    ParGridFunction v(&mesh_pfes);
   //    VectorFunctionCoefficient vcoeff(dim, velocity_function);
   //    v.ProjectCoefficient(vcoeff);

   //    double t = 0.0;
   //    while (t < t_final)
   //    {
   //       t += dt;
   //       // Move the mesh nodes.
   //       x.Add(std::min(dt, t_final-t), v);
   //       // Update the node velocities.
   //       v.ProjectCoefficient(vcoeff);
   //    }

   //    // Pseudotime velocity.
   //    add(x, -1.0, x0, v_gf);

   //    // Return the mesh to the initial configuration.
   //    x = x0;
   // }

   // Mesh velocity.
   // If remap is on, obtain the mesh velocity by moving the mesh to the final
   // mesh positions, and taking the displacement vector.
   // The mesh motion resembles a time-dependent deformation, e.g., similar to
   // a deformation that is obtained by a Lagrangian simulation.
   GridFunction v_gf(x.FESpace());
   VectorGridFunctionCoefficient v_mesh_coeff(&v_gf);
   v_gf  = x; v_gf -= x_mod;
   v_gf.Neg();
   
   // Define the discontinuous DG finite element space of the given
   // polynomial order on the refined mesh.
   const int btype = BasisType::Positive;
   DG_FECollection fec(order, dim, btype);
   ParFiniteElementSpace pfes(pmesh, &fec);
   
   // Check for meaningful combinations of parameters.
   const bool forced_bounds = lo_type   != LOSolverType::None ||
                              mono_type != MonolithicSolverType::None;
   if (forced_bounds)
   {
      MFEM_VERIFY(btype == 2,
                  "Monotonicity treatment requires Bernstein basis.");

      if (order == 0)
      {
         // Disable monotonicity treatment for piecewise constants.
         if (myid == 0)
         { mfem_warning("For -o 0, monotonicity treatment is disabled."); }
         lo_type = LOSolverType::None;
         fct_type = FCTSolverType::None;
         mono_type = MonolithicSolverType::None;
      }
   }
   const bool use_subcell_RD =
      ( lo_type   == LOSolverType::ResDistSubcell ||
        mono_type == MonolithicSolverType::ResDistMonoSubcell );
   if (use_subcell_RD)
   {
      MFEM_VERIFY(order > 1, "Subcell schemes require FE order > 2.");
   }
   if (dt_control == TimeStepControl::LOBoundsError)
   {
      MFEM_VERIFY(bounds_type == 1, "Error: -dtc 1 requires -bt 1.");
   }

   const int prob_size = pfes.GlobalTrueVSize();
   if (myid == 0) { cout << "Number of unknowns: " << prob_size << endl; }

   // Fields related to inflow BC.
   FunctionCoefficient inflow(inflow_function);
   ParGridFunction inflow_gf(&pfes);
   inflow_gf =0.0;
   // inflow_gf.ProjectCoefficient(inflow);

   // Set up the bilinear and linear forms corresponding to the DG
   // discretization.
   ParBilinearForm m(&pfes);
   m.AddDomainIntegrator(new MassIntegrator);

   ParBilinearForm M_HO(&pfes);
   M_HO.AddDomainIntegrator(new MassIntegrator);

   ParBilinearForm k(&pfes);
   ParBilinearForm K_HO(&pfes);

   if (exec_mode == 0)
   {
      k.AddDomainIntegrator(new ConvectionIntegrator(velocity, -1.0));
      K_HO.AddDomainIntegrator(new ConvectionIntegrator(velocity, -1.0));
   }
   else if (exec_mode == 1)
   {
      k.AddDomainIntegrator(new ConvectionIntegrator(v_mesh_coeff)); // mesh velocity
      K_HO.AddDomainIntegrator(new ConvectionIntegrator(v_mesh_coeff)); // mesh velocity
   }

   if (ho_type == HOSolverType::CG ||
       ho_type == HOSolverType::LocalInverse ||
       fct_type == FCTSolverType::FluxBased)
   {
      if (exec_mode == 0)
      {
         DGTraceIntegrator *dgt_i = new DGTraceIntegrator(velocity, 1.0, -0.5);
         DGTraceIntegrator *dgt_b = new DGTraceIntegrator(velocity, 1.0, -0.5);
         K_HO.AddInteriorFaceIntegrator(new TransposeIntegrator(dgt_i));
         K_HO.AddBdrFaceIntegrator(new TransposeIntegrator(dgt_b));
      }
      else if (exec_mode == 1)
      {
         auto dgt_i = new DGTraceIntegrator(v_mesh_coeff, -1.0, -0.5);
         auto dgt_b = new DGTraceIntegrator(v_mesh_coeff, -1.0, -0.5);
         K_HO.AddInteriorFaceIntegrator(new TransposeIntegrator(dgt_i));
         K_HO.AddBdrFaceIntegrator(new TransposeIntegrator(dgt_b));
      }

      K_HO.KeepNbrBlock(true);
   }

   if (pa)
   {
      M_HO.SetAssemblyLevel(AssemblyLevel::PARTIAL);
      K_HO.SetAssemblyLevel(AssemblyLevel::PARTIAL);
   }

   if (next_gen_full)
   {
      K_HO.SetAssemblyLevel(AssemblyLevel::FULL);
   }

   
   M_HO.Assemble();
   // std::cout <<"K_HO.Assemble(0);" << std::endl;
   K_HO.Assemble(0);
   // std::cout <<"K_HO.Assemble(1);" << std::endl;

   if (pa == false)
   {
      M_HO.Finalize();
      K_HO.Finalize(0);
   }

   // Compute the lumped mass matrix.
   Vector lumpedM;
   ParBilinearForm ml(&pfes);
   ml.AddDomainIntegrator(new LumpedIntegrator(new MassIntegrator));
   ml.Assemble();
   ml.Finalize();
   ml.SpMat().GetDiag(lumpedM);

   m.Assemble();
   m.Finalize();
   int skip_zeros = 0;
   k.Assemble(skip_zeros);
   k.Finalize(skip_zeros);

   // Store topological dof data.
   DofInfo dofs(pfes, bounds_type);

   // Precompute data required for high and low order schemes. This could be put
   // into a separate routine. I am using a struct now because the various
   // schemes require quite different information.
   LowOrderMethod lom;
   lom.subcell_scheme = use_subcell_RD;

   lom.pk = NULL;
   if (lo_type == LOSolverType::DiscrUpwind)
   {
      lom.smap = SparseMatrix_Build_smap(k.SpMat());
      lom.D = k.SpMat();

      if (exec_mode == 0)
      {
         ComputeDiscreteUpwindingMatrix(k.SpMat(), lom.smap, lom.D);
      }
   }
   else if (lo_type == LOSolverType::DiscrUpwindPrec)
   {
      lom.pk = new ParBilinearForm(&pfes);
      if (exec_mode == 0)
      {
         lom.pk->AddDomainIntegrator(
            new PrecondConvectionIntegrator(velocity, -1.0) );
      }
      else if (exec_mode == 1)
      {
         lom.pk->AddDomainIntegrator(
            new PrecondConvectionIntegrator(v_mesh_coeff) );
      }
      lom.pk->Assemble(skip_zeros);
      lom.pk->Finalize(skip_zeros);

      lom.smap = SparseMatrix_Build_smap(lom.pk->SpMat());
      lom.D = lom.pk->SpMat();

      if (exec_mode == 0)
      {
         ComputeDiscreteUpwindingMatrix(lom.pk->SpMat(), lom.smap, lom.D);
      }
   }
   if (exec_mode == 1) { lom.coef = &v_mesh_coeff; }
   else                { lom.coef = &velocity; }

   // Face integration rule.
   const FaceElementTransformations *ft = pmesh->GetFaceElementTransformations(0);
   const int el_order = pfes.GetFE(0)->GetOrder();
   int ft_order = ft->Elem1->OrderW() + 2 * el_order;
   if (pfes.GetFE(0)->Space() == FunctionSpace::Pk) { ft_order++; }
   lom.irF = &IntRules.Get(ft->FaceGeom, ft_order);

   DG_FECollection fec0(0, dim, btype);
   DG_FECollection fec1(1, dim, btype);

   ParMesh *subcell_mesh = NULL;
   lom.SubFes0 = NULL;
   lom.SubFes1 = NULL;
   FiniteElementCollection *fec_sub = NULL;
   ParFiniteElementSpace *pfes_sub = NULL;;
   ParGridFunction *xsub = NULL;
   ParGridFunction v_sub_gf;
   VectorGridFunctionCoefficient v_sub_coef;
   Vector x0_sub;

   if (use_subcell_RD)
   {
      // Get a uniformly refined mesh.
      const int btype = BasisType::ClosedUniform;
      subcell_mesh = new ParMesh(ParMesh::MakeRefined(*pmesh, order, btype));

      // Check if the mesh is periodic.
      const L2_FECollection *L2_coll = dynamic_cast<const L2_FECollection *>
                                       (pmesh->GetNodes()->FESpace()->FEColl());
      if (L2_coll == NULL)
      {
         // Standard non-periodic mesh.
         // Note that the fine mesh is always linear.
         fec_sub = new H1_FECollection(1, dim, BasisType::ClosedUniform);
         pfes_sub = new ParFiniteElementSpace(subcell_mesh, fec_sub, dim);
         xsub = new ParGridFunction(pfes_sub);
         subcell_mesh->SetCurvature(1);
         subcell_mesh->SetNodalGridFunction(xsub);
      }
      else
      {
         // Periodic mesh - the node positions must be corrected after the call
         // to the above Mesh constructor. Note that the fine mesh is always
         // linear.
         const bool disc_nodes = true;
         subcell_mesh->SetCurvature(1, disc_nodes);

         fec_sub = new L2_FECollection(1, dim, BasisType::ClosedUniform);
         pfes_sub = new ParFiniteElementSpace(subcell_mesh, fec_sub, dim);
         xsub = new ParGridFunction(pfes_sub);
         subcell_mesh->SetNodalGridFunction(xsub);

         GridFunction *coarse = pmesh->GetNodes();
         InterpolationGridTransfer transf(*coarse->FESpace(), *pfes_sub);
         transf.ForwardOperator().Mult(*coarse, *xsub);
      }

      lom.SubFes0 = new FiniteElementSpace(subcell_mesh, &fec0);
      lom.SubFes1 = new FiniteElementSpace(subcell_mesh, &fec1);

      // Submesh velocity.
      v_sub_gf.SetSpace(pfes_sub);
      v_sub_gf.ProjectCoefficient(velocity);
      // v_sub_gf.ProjectCoefficient(v_mesh_coeff);

      // Zero it out on boundaries (not moving boundaries).
      Array<int> ess_bdr, ess_vdofs;
      if (subcell_mesh->bdr_attributes.Size() > 0)
      {
         ess_bdr.SetSize(subcell_mesh->bdr_attributes.Max());
      }
      ess_bdr = 1;
      xsub->ParFESpace()->GetEssentialVDofs(ess_bdr, ess_vdofs);
      for (int i = 0; i < ess_vdofs.Size(); i++)
      {
         if (ess_vdofs[i] == -1) { v_sub_gf(i) = 0.0; }
      }
      v_sub_coef.SetGridFunction(&v_sub_gf);

      // Store initial submesh positions.
      x0_sub = *xsub;

      // Integrator on the submesh.
      if (exec_mode == 0)
      {
         lom.subcellCoeff = &velocity;
         lom.VolumeTerms = new MixedConvectionIntegrator(velocity, -1.0);
      }
      else if (exec_mode == 1)
      {
         lom.subcellCoeff = &v_sub_coef;
         lom.VolumeTerms = new MixedConvectionIntegrator(v_sub_coef);
      }
   }
   else { subcell_mesh = pmesh; }

   Assembly asmbl(dofs, lom, inflow_gf, pfes, subcell_mesh, exec_mode);

   // Setup the initial conditions.
   const int vsize = pfes.GetVSize();
   Array<int> offset((product_sync) ? 3 : 2);
   for (int i = 0; i < offset.Size(); i++) { offset[i] = i*vsize; }
   BlockVector S(offset, Device::GetMemoryType());
   // Primary scalar field is u.
   ParGridFunction u(&pfes);
   u.MakeRef(&pfes, S, offset[0]);
   // FunctionCoefficient u0(u0_function);
   // u.ProjectCoefficient(u0);
   u = u_gf;
   u.SyncAliasMemory(S);
   // For the case of product remap, we also solve for s and u_s.
   ParGridFunction s, us;
   Array<bool> u_bool_el, u_bool_dofs;
   if (product_sync)
   {
      s.SetSpace(&pfes);
      ComputeBoolIndicators(pmesh->GetNE(), u, u_bool_el, u_bool_dofs);
      BoolFunctionCoefficient sc(s0_function, u_bool_el);
      s.ProjectCoefficient(sc);

      us.MakeRef(&pfes, S, offset[1]);
      double *h_us = us.HostWrite();
      const double *h_u = u.HostRead();
      const double *h_s = s.HostRead();
      // Simple - we don't target conservation at initialization.
      for (int i = 0; i < s.Size(); i++) { h_us[i] = h_u[i] * h_s[i]; }
      us.SyncAliasMemory(S);
   }

   // Smoothness indicator.
   SmoothnessIndicator *smth_indicator = NULL;
   if (smth_ind_type)
   {
      smth_indicator = new SmoothnessIndicator(smth_ind_type, *subcell_mesh,
                                               pfes, u, dofs);
   }

   // Setup of the high-order solver (if any).
   HOSolver *ho_solver = NULL;
   if (ho_type == HOSolverType::Neumann)
   {
      ho_solver = new NeumannHOSolver(pfes, m, k, lumpedM, asmbl);
   }
   else if (ho_type == HOSolverType::CG)
   {
      ho_solver = new CGHOSolver(pfes, M_HO, K_HO);
   }
   else if (ho_type == HOSolverType::LocalInverse)
   {
      ho_solver = new LocalInverseHOSolver(pfes, M_HO, K_HO);
   }

   // Setup the low order solver (if any).
   LOSolver *lo_solver = NULL;
   Array<int> lo_smap;
   const bool time_dep = (exec_mode == 0) ? false : true;
   if (lo_type == LOSolverType::DiscrUpwind)
   {
      lo_smap = SparseMatrix_Build_smap(k.SpMat());
      lo_solver = new DiscreteUpwind(pfes, k.SpMat(), lo_smap,
                                     lumpedM, asmbl, time_dep);
   }
   else if (lo_type == LOSolverType::DiscrUpwindPrec)
   {
      lo_smap = SparseMatrix_Build_smap(lom.pk->SpMat());
      lo_solver = new DiscreteUpwind(pfes, lom.pk->SpMat(), lo_smap,
                                     lumpedM, asmbl, time_dep);
   }
   else if (lo_type == LOSolverType::ResDist)
   {
      const bool subcell_scheme = false;
      if (pa)
      {
         lo_solver = new PAResidualDistribution(pfes, k, asmbl, lumpedM,
                                                subcell_scheme, time_dep);
         if (exec_mode == 0)
         {
            const PAResidualDistribution *RD_ptr =
               dynamic_cast<const PAResidualDistribution*>(lo_solver);
            RD_ptr->SampleVelocity(FaceType::Interior);
            RD_ptr->SampleVelocity(FaceType::Boundary);
            RD_ptr->SetupPA(FaceType::Interior);
            RD_ptr->SetupPA(FaceType::Boundary);
         }
      }
      else
      {
         lo_solver = new ResidualDistribution(pfes, k, asmbl, lumpedM,
                                              subcell_scheme, time_dep);
      }
   }
   else if (lo_type == LOSolverType::ResDistSubcell)
   {
      const bool subcell_scheme = true;
      if (pa)
      {
         lo_solver = new PAResidualDistributionSubcell(pfes, k, asmbl, lumpedM,
                                                       subcell_scheme, time_dep);
         if (exec_mode == 0)
         {
            const PAResidualDistributionSubcell *RD_ptr =
               dynamic_cast<const PAResidualDistributionSubcell*>(lo_solver);
            RD_ptr->SampleVelocity(FaceType::Interior);
            RD_ptr->SampleVelocity(FaceType::Boundary);
            RD_ptr->SetupPA(FaceType::Interior);
            RD_ptr->SetupPA(FaceType::Boundary);
         }
      }
      else
      {
         lo_solver = new ResidualDistribution(pfes, k, asmbl, lumpedM,
                                              subcell_scheme, time_dep);
      }
   }
   else if (lo_type == LOSolverType::MassBased)
   {
      MFEM_VERIFY(ho_solver != nullptr,
                  "Mass-Based LO solver requires a choice of a HO solver.");
      lo_solver = new MassBasedAvg(pfes, *ho_solver,
                                   (exec_mode == 1) ? &v_gf : nullptr);
   }

   // Setup of the monolithic solver (if any).
   MonolithicSolver *mono_solver = NULL;
   bool mass_lim = (problem_num != 6 && problem_num != 7) ? true : false;
   if (mono_type == MonolithicSolverType::ResDistMono)
   {
      const bool subcell_scheme = false;
      mono_solver = new MonoRDSolver(pfes, k.SpMat(), m.SpMat(), lumpedM,
                                     asmbl, smth_indicator, velocity,
                                     subcell_scheme, time_dep, mass_lim);
   }
   else if (mono_type == MonolithicSolverType::ResDistMonoSubcell)
   {
      const bool subcell_scheme = true;
      mono_solver = new MonoRDSolver(pfes, k.SpMat(), m.SpMat(), lumpedM,
                                     asmbl, smth_indicator, velocity,
                                     subcell_scheme, time_dep, mass_lim);
   }

   // Print the starting meshes and initial condition.
   if(visualization)
   {
      ofstream meshHO("meshHO_init.mesh");
      meshHO.precision(precision);
      pmesh->PrintAsOne(meshHO);
      if (subcell_mesh)
      {
         ofstream meshLO("meshLO_init.mesh");
         meshLO.precision(precision);
         subcell_mesh->PrintAsOne(meshLO);
      }
      ofstream sltn("sltn_init.gf");
      sltn.precision(precision);
      u.SaveAsOne(sltn);
   }

   // Create data collection for solution output: either VisItDataCollection for
   // ASCII data files, or SidreDataCollection for binary data files.
   DataCollection *dc = NULL;
   if (visit)
   {
      dc = new VisItDataCollection("Remhos", pmesh);
      dc->SetPrecision(precision);
      dc->RegisterField("solution", &u);
      dc->SetCycle(0);
      dc->SetTime(0.0);
      dc->Save();
   }
   
   socketstream sout, vis_s, vis_us;
   char vishost[] = "localhost";
   int  visport   = 19916;
   if (visualization)
   {
      // Make sure all MPI ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh->GetComm());

      sout.precision(8);
      vis_s.precision(8);
      vis_us.precision(8);

      int Wx = 0, Wy = 0; // window position
      const int Ww = 400, Wh = 400; // window size
      u.HostRead();
      s.HostRead();
      VisualizeField(sout, vishost, visport, u, "Solution u", Wx, Wy, Ww, Wh);
      if (product_sync)
      {
         VisualizeField(vis_s, vishost, visport, s, "Solution s",
                        Wx + Ww, Wy, Ww, Wh);
         VisualizeField(vis_us, vishost, visport, us, "Solution u_s",
                        Wx + 2*Ww, Wy, Ww, Wh);
      }
   }

   // Record the initial mass.
   Vector masses(lumpedM);
   const double mass0_u_loc = lumpedM * u;
   double mass0_u, mass0_us;
   MPI_Allreduce(&mass0_u_loc, &mass0_u, 1, MPI_DOUBLE, MPI_SUM, comm);
   if (product_sync)
   {
      const double mass0_us_loc = lumpedM * us;
      MPI_Allreduce(&mass0_us_loc, &mass0_us, 1, MPI_DOUBLE, MPI_SUM, comm);
   }

   // Setup of the FCT solver (if any).
   Array<int> K_HO_smap;
   FCTSolver *fct_solver = NULL;
   if (fct_type == FCTSolverType::FluxBased)
   {
      MFEM_VERIFY(pa == false, "Flux-based FCT and PA are incompatible.");
      K_HO.SpMat().HostReadI();
      K_HO.SpMat().HostReadJ();
      K_HO.SpMat().HostReadData();
      K_HO_smap = SparseMatrix_Build_smap(K_HO.SpMat());
      const int fct_iterations = 1;
      fct_solver = new FluxBasedFCT(pfes, smth_indicator, dt, K_HO.SpMat(),
                                    K_HO_smap, M_HO.SpMat(), fct_iterations);
   }
   else if (fct_type == FCTSolverType::ClipScale)
   {
      fct_solver = new ClipScaleSolver(pfes, smth_indicator, dt);
   }
   else if (fct_type == FCTSolverType::NonlinearPenalty)
   {
      fct_solver = new NonlinearPenaltySolver(pfes, smth_indicator, dt);
   }
   else if (fct_type == FCTSolverType::FCTProject)
   {
      fct_solver = new ElementFCTProjection(pfes, dt);
   }

   AdvectionOperator adv(S.Size(), m, ml, lumpedM, k, M_HO, K_HO,
                         x, xsub, v_gf, v_sub_gf, asmbl, lom, dofs,
                         ho_solver, lo_solver, fct_solver, mono_solver);

   double t = 0.0;
   adv.SetTime(t);
   adv.SetTimeStepControl(dt_control);
   ode_solver->Init(adv);

   double umin, umax;
   GetMinMax(u, umin, umax);

   if (exec_mode == 1)
   {
      adv.SetRemapStartPos(x0, x0_sub);

      // For remap, the pseudo-time always evolves from 0 to 1.
      t_final = 1.0;
   }

   ParGridFunction res = u;
   double residual = 0.0;
   double s_min_glob = numeric_limits<double>::infinity(),
          s_max_glob = -numeric_limits<double>::infinity();

   // Time-integration (loop over the time iterations, ti, with a time-step dt).
   bool done = false;
   BlockVector Sold(S);
   int ti_total = 0, ti = 0;

   while (done == false)
   {
      double dt_real = min(dt, t_final - t);

      // This also resets the time step estimate when automatic dt is on.
      adv.SetDt(dt_real);
      if (lo_solver)  { lo_solver->UpdateTimeStep(dt_real); }
      if (fct_solver) { fct_solver->UpdateTimeStep(dt_real); }

      if (product_sync)
      {
         ComputeMinMaxS(pmesh->GetNE(), us, u, s_min_glob, s_max_glob);
#ifdef REMHOS_FCT_PRODUCT_DEBUG
         if (myid == 0)
         {
            std::cout << "   --- Full time step" << std::endl;
            std::cout << "   in:  ";
            std::cout << std::scientific << std::setprecision(5);
            std::cout << "min_s: " << s_min_glob
                      << "; max_s: " << s_max_glob << std::endl;
         }
#endif
      }

      Sold = S;
      ode_solver->Step(S, t, dt_real);
      ti++;
      ti_total++;

      if (dt_control != TimeStepControl::FixedTimeStep)
      {
         double dt_est = adv.GetTimeStepEstimate();
         if (dt_est < dt_real)
         {
            // Repeat with the proper time step.
            if (myid == 0)
            {
               cout << "Repeat / decrease dt: "
                    << dt_real << " --> " << 0.85 * dt << endl;
            }
            ti--;
            t -= dt_real;
            S  = Sold;
            dt = 0.85 * dt;
            if (dt < 1e-12) { MFEM_ABORT("The time step crashed!"); }
            continue;
         }
         else if (dt_est > 1.25 * dt_real) { dt *= 1.02; }
      }

      // S has been modified, update the alias
      u.SyncMemory(S);
      if (product_sync)
      {
         us.SyncMemory(S);

         // It is known that RK time integrators with more than 1 stage may
         // cause violation of the lower bounds for us.
         // The lower bound is corrected, causing small conservation error.
         // Correction can also be done with localized bounds for s, but for
         // now we have implemented only the minimum global bound.
         u.HostRead();
         us.HostReadWrite();
         const int s = u.Size();
         Array<bool> active_elem, active_dofs;
         ComputeBoolIndicators(NE, u, active_elem, active_dofs);
         for (int i = 0; i < s; i++)
         {
            if (active_dofs[i] == false) { continue; }

            double us_min = u(i) * s_min_glob;
            if (us(i) < us_min) { us(i) = us_min; }
         }

#ifdef REMHOS_FCT_PRODUCT_DEBUG
         ComputeMinMaxS(NE, us, u, s_min_glob, s_max_glob);
         if (myid == 0)
         {
            std::cout << "   out: ";
            std::cout << std::scientific << std::setprecision(5);
            std::cout << "min_s: " << s_min_glob
                      << "; max_s: " << s_max_glob << std::endl;
         }
#endif
      }

      // Monotonicity check for debug purposes mainly.
      if (verify_bounds && forced_bounds && smth_indicator == NULL)
      {
         double umin_new, umax_new;
         GetMinMax(u, umin_new, umax_new);
         if (problem_num % 10 != 6 && problem_num % 10 != 7)
         {
            if (myid == 0)
            {
               MFEM_VERIFY(umin_new > umin - 1e-12,
                           "Undershoot of " << umin - umin_new);
               MFEM_VERIFY(umax_new < umax + 1e-12,
                           "Overshoot of " << umax_new - umax);
            }
            umin = umin_new;
            umax = umax_new;
         }
         else
         {
            if (myid == 0)
            {
               MFEM_VERIFY(umin_new > 0.0 - 1e-12,
                           "Undershoot of " << 0.0 - umin_new);
               MFEM_VERIFY(umax_new < 1.0 + 1e-12,
                           "Overshoot of " << umax_new - 1.0);
            }
         }
      }

      if (exec_mode == 1)
      {
         x0.HostReadWrite(); v_sub_gf.HostReadWrite();
         x.HostReadWrite();
         add(x0, t, v_gf, x);
         if (use_subcell_RD)
         {
            x0_sub.HostReadWrite(); v_sub_gf.HostReadWrite();
            MFEM_VERIFY(xsub != NULL, "Subcell mesh not defined!");
            xsub->HostReadWrite();
            add(x0_sub, t, v_sub_gf, *xsub);
         }
      }

      if (problem_num != 6 && problem_num != 7 && problem_num != 8)
      {
         done = (t >= t_final - 1.e-8*dt);
      }
      else
      {
         // Steady state problems - stop at convergence.
         double res_loc = 0.;
         lumpedM.HostReadWrite(); u.HostReadWrite(); res.HostReadWrite();
         for (int i = 0; i < res.Size(); i++)
         {
            res_loc += pow( (lumpedM(i) * u(i) / dt) - (lumpedM(i) * res(i) / dt), 2. );
         }
         MPI_Allreduce(&res_loc, &residual, 1, MPI_DOUBLE, MPI_SUM, comm);

         residual = sqrt(residual);
         if (residual < 1.e-12 && t >= 1.) { done = true; u = res; }
         else { res = u; }
      }


      // if (done || ti % vis_steps == 0)
      if (vis_steps == 0)
      {
         if (myid == 0)
         {
            cout << "time step: " << ti << ", time: " << t
                 << ", dt: " << dt << ", residual: " << residual << endl;
         }

         if (visualization)
         {
            int Wx = 0, Wy = 0; // window position
            int Ww = 400, Wh = 400; // window size
            VisualizeField(sout, vishost, visport, u, "Solution",
                           Wx, Wy, Ww, Wh);
            if (product_sync)
            {
               // Recompute s = u_s / u.
               ComputeRatio(pmesh->GetNE(), us, u, s, u_bool_el, u_bool_dofs);
               VisualizeField(vis_s, vishost, visport, s, "Solution s",
                              Wx + Ww, Wy, Ww, Wh);
               VisualizeField(vis_us, vishost, visport, us, "Solution u_s",
                              Wx + 2*Ww, Wy, Ww, Wh);
            }
         }

         if (visit)
         {
            dc->SetCycle(ti);
            dc->SetTime(t);
            dc->Save();
         }
      }
   }

   if (dt_control != TimeStepControl::FixedTimeStep && myid == 0)
   {
      cout << "Total time steps: " << ti_total
           << " (" << ti_total-ti << " repeated)." << endl;
   }

   // Print the final meshes and solution.
   if (visualization)
   {
      // Print the final meshes and solution.
      {
         ofstream meshHO("meshHO_final.mesh");
         meshHO.precision(precision);
         pmesh->PrintAsOne(meshHO);
         if (subcell_mesh)
         {
            ofstream meshLO("meshLO_final.mesh");
            meshLO.precision(precision);
            subcell_mesh->PrintAsOne(meshLO);
         }
         ofstream sltn("sltn_final.gf");
         sltn.precision(precision);
         u.SaveAsOne(sltn);
      }
   }

   // Check for mass conservation.
   double mass_u_loc = 0.0, mass_us_loc = 0.0;
   if (exec_mode == 1)
   {
      ml.BilinearForm::operator=(0.0);
      ml.Assemble();
      lumpedM.HostRead();
      ml.SpMat().GetDiag(lumpedM);
      mass_u_loc = lumpedM * u;
      if (product_sync) { mass_us_loc = lumpedM * us; }
   }
   else
   {
      mass_u_loc = masses * u;
      if (product_sync) { mass_us_loc = masses * us; }
   }
   double mass_u, mass_us = 0.0, s_max = 0.0;
   MPI_Allreduce(&mass_u_loc, &mass_u, 1, MPI_DOUBLE, MPI_SUM, comm);
   const double umax_loc = u.Max();
   MPI_Allreduce(&umax_loc, &umax, 1, MPI_DOUBLE, MPI_MAX, comm);
   if (product_sync)
   {
      ComputeRatio(pmesh->GetNE(), us, u, s, u_bool_el, u_bool_dofs);
      const double s_max_loc = s.Max();
      MPI_Allreduce(&mass_us_loc, &mass_us, 1, MPI_DOUBLE, MPI_SUM, comm);
      MPI_Allreduce(&s_max_loc, &s_max, 1, MPI_DOUBLE, MPI_MAX, comm);
   }
   if (myid == 0)
   {
      cout << setprecision(10)
           << "Final mass u:  " << mass_u << endl
           << "Max value u:   " << umax << endl << setprecision(6)
           << "Mass loss u:   " << abs(mass0_u - mass_u) << endl;
      if (product_sync)
      {
         cout << setprecision(10)
              << "Final mass us: " << mass_us << endl
              << "Max value s:   " << s_max << endl << setprecision(6)
              << "Mass loss us:  " << abs(mass0_us - mass_us) << endl;
      }
   }
   u_gf = u; 

   // // Compute errors, if the initial condition is equal to the final solution
   // if (problem_num == 4) // solid body rotation
   // {
   //    double err = u.ComputeLpError(1., u0);
   //    if (myid == 0) { cout << "L1-error: " << err << "." << endl; }
   // }
   // else if (problem_num == 7)
   // {
   //    FunctionCoefficient u_ex(inflow_function);
   //    double e1 = u.ComputeLpError(1., u_ex);
   //    double e2 = u.ComputeLpError(2., u_ex);
   //    double eInf = u.ComputeLpError(numeric_limits<double>::infinity(), u_ex);
   //    if (myid == 0)
   //    {
   //       cout << "L1-error: " << e1 << "." << endl;

   //       // write output
   //       ofstream file("errors.txt", ios_base::app);

   //       if (!file)
   //       {
   //          MFEM_ABORT("Error opening file.");
   //       }
   //       else
   //       {
   //          ostringstream strs;
   //          strs << e1 << " " << e2 << " " << eInf << "\n";
   //          string str = strs.str();
   //          file << str;
   //          file.close();
   //       }
   //    }
   // }

   // if (smth_indicator)
   // {
   //    // Print the values of the smoothness indicator.
   //    ParGridFunction si_val;
   //    smth_indicator->ComputeSmoothnessIndicator(u, si_val);
   //    {
   //       ofstream smth("si_final.gf");
   //       smth.precision(precision);
   //       si_val.SaveAsOne(smth);
   //    }
   // }

   // Free the used memory.
   delete mono_solver;
   delete fct_solver;
   delete smth_indicator;
   delete lo_solver;
   delete ho_solver;

   delete ode_solver;
   delete mesh_fec;
   delete lom.pk;
   delete dc;

   if (use_subcell_RD)
   {
      delete subcell_mesh;
      delete fec_sub;
      delete pfes_sub;
      delete xsub;
      delete lom.SubFes0;
      delete lom.SubFes1;
      delete lom.VolumeTerms;
   }

   // return 0;
}  

AdvectionOperator::AdvectionOperator(int size, BilinearForm &Mbf_,
                                     BilinearForm &_ml, Vector &_lumpedM,
                                     ParBilinearForm &Kbf_,
                                     ParBilinearForm &M_HO_, ParBilinearForm &K_HO_,
                                     GridFunction &pos, GridFunction *sub_pos,
                                     GridFunction &vel, GridFunction &sub_vel,
                                     Assembly &_asmbl,
                                     LowOrderMethod &_lom, DofInfo &_dofs,
                                     HOSolver *hos, LOSolver *los, FCTSolver *fct,
                                     MonolithicSolver *mos) :
   TimeDependentOperator(size), Mbf(Mbf_), ml(_ml), Kbf(Kbf_),
   M_HO(M_HO_), K_HO(K_HO_),
   lumpedM(_lumpedM),
   start_mesh_pos(pos.Size()), start_submesh_pos(sub_vel.Size()),
   mesh_pos(pos), submesh_pos(sub_pos),
   mesh_vel(vel), submesh_vel(sub_vel),
   x_gf(Kbf.ParFESpace()),
   asmbl(_asmbl), lom(_lom), dofs(_dofs),
   ho_solver(hos), lo_solver(los), fct_solver(fct), mono_solver(mos) { }

void AdvectionOperator::Mult(const Vector &X, Vector &Y) const
{
   if (exec_mode == 1)
   {
      // Move the mesh positions.
      const double t = GetTime();
      add(start_mesh_pos, t, mesh_vel, mesh_pos);
      if (submesh_pos)
      {
         add(start_submesh_pos, t, submesh_vel, *submesh_pos);
      }
      // Reset precomputed geometric data.
      Mbf.FESpace()->GetMesh()->DeleteGeometricFactors();

      // Reassemble on the new mesh. Element contributions.
      // Currently needed to have the sparse matrices used by the LO methods.
      Mbf.BilinearForm::operator=(0.0);
      Mbf.Assemble();
      Kbf.BilinearForm::operator=(0.0);
      Kbf.Assemble(0);
      ml.BilinearForm::operator=(0.0);
      ml.Assemble();
      lumpedM.HostReadWrite();
      ml.SpMat().GetDiag(lumpedM);

      M_HO.BilinearForm::operator=(0.0);
      M_HO.Assemble();
      K_HO.BilinearForm::operator=(0.0);
      K_HO.Assemble(0);

      if (lom.pk)
      {
         lom.pk->BilinearForm::operator=(0.0);
         lom.pk->Assemble();
      }

      // Face contributions.
      asmbl.bdrInt = 0.;
      Mesh *mesh = M_HO.FESpace()->GetMesh();
      const int dim = mesh->Dimension(), ne = mesh->GetNE();
      Array<int> bdrs, orientation;
      FaceElementTransformations *Trans;

      if (auto RD_ptr = dynamic_cast<const PAResidualDistribution*>(lo_solver))
      {
         RD_ptr->SampleVelocity(FaceType::Interior);
         RD_ptr->SampleVelocity(FaceType::Boundary);
         RD_ptr->SetupPA(FaceType::Interior);
         RD_ptr->SetupPA(FaceType::Boundary);
      }
      else
      {
         for (int k = 0; k < ne; k++)
         {
            if (dim == 1)      { mesh->GetElementVertices(k, bdrs); }
            else if (dim == 2) { mesh->GetElementEdges(k, bdrs, orientation); }
            else if (dim == 3) { mesh->GetElementFaces(k, bdrs, orientation); }

            for (int i = 0; i < dofs.numBdrs; i++)
            {
               Trans = mesh->GetFaceElementTransformations(bdrs[i]);
               asmbl.ComputeFluxTerms(k, i, Trans, lom);
            }
         }
      }
   }

   const int size = Kbf.ParFESpace()->GetVSize();
   const int NE   = Kbf.ParFESpace()->GetNE();

   // Needed because X and Y are allocated on the host by the ODESolver.
   X.Read(); Y.Read();

   Vector u, d_u;
   Vector* xptr = const_cast<Vector*>(&X);
   u.MakeRef(*xptr, 0, size);
   d_u.MakeRef(Y, 0, size);
   Vector du_HO(u.Size()), du_LO(u.Size());

   x_gf = u;
   x_gf.ExchangeFaceNbrData();

   if (mono_solver) { mono_solver->CalcSolution(u, d_u); }
   else if (fct_solver)
   {
      MFEM_VERIFY(ho_solver && lo_solver, "FCT requires HO and LO solvers.");

      lo_solver->CalcLOSolution(u, du_LO);
      ho_solver->CalcHOSolution(u, du_HO);

      dofs.ComputeElementsMinMax(u, dofs.xe_min, dofs.xe_max, NULL, NULL);
      dofs.ComputeBounds(dofs.xe_min, dofs.xe_max, dofs.xi_min, dofs.xi_max);
      fct_solver->CalcFCTSolution(x_gf, lumpedM, du_HO, du_LO,
                                  dofs.xi_min, dofs.xi_max, d_u);

      if (dt_control == TimeStepControl::LOBoundsError)
      {
         UpdateTimeStepEstimate(u, du_LO, dofs.xi_min, dofs.xi_max);
      }
   }
   else if (lo_solver)
   {
      lo_solver->CalcLOSolution(u, d_u);

      if (dt_control == TimeStepControl::LOBoundsError)
      {
         dofs.ComputeElementsMinMax(u, dofs.xe_min, dofs.xe_max, NULL, NULL);
         dofs.ComputeBounds(dofs.xe_min, dofs.xe_max, dofs.xi_min, dofs.xi_max);
         UpdateTimeStepEstimate(u, d_u, dofs.xi_min, dofs.xi_max);
      }
   }
   // The HO option must be last, since some LO solvers use the HO. Then if the
   // user only wants to run LO, this order will give him the LO solution.
   else if (ho_solver) { ho_solver->CalcHOSolution(u, d_u); }
   else { MFEM_ABORT("No solver was chosen."); }

   d_u.SyncAliasMemory(Y);

   // Remap the product field, if there is a product field.
   if (X.Size() > size)
   {
      MFEM_VERIFY(exec_mode == 1, "Products are processed only in remap mode.");
      MFEM_VERIFY(dt_control == TimeStepControl::FixedTimeStep,
                  "Automatic time step is not implemented for product remap.");

      Vector us, d_us;
      us.MakeRef(*xptr, size, size);
      d_us.MakeRef(Y, size, size);

      x_gf = us;
      x_gf.ExchangeFaceNbrData();

      if (mono_solver) { mono_solver->CalcSolution(us, d_us); }
      else if (fct_solver)
      {
         MFEM_VERIFY(ho_solver && lo_solver, "FCT requires HO and LO solvers.");

         Vector d_us_HO(us.Size()), d_us_LO;
         if (fct_solver->NeedsLOProductInput())
         {
            d_us_LO.SetSize(us.Size());
            lo_solver->CalcLOSolution(us, d_us_LO);
         }
         ho_solver->CalcHOSolution(us, d_us_HO);

         // Compute the ratio s = us_old / u_old, and old active dofs.
         Vector s(size);
         Array<bool> s_bool_el, s_bool_dofs;
         ComputeRatio(NE, us, u, s, s_bool_el, s_bool_dofs);
#ifdef REMHOS_FCT_PRODUCT_DEBUG
         const int myid = x_gf.ParFESpace()->GetMyRank();
         if (myid == 0) { std::cout << "      --- RK stage" << std::endl; }
         std::cout << "      in:  ";
         ComputeMinMaxS(s, s_bool_dofs, myid);
#endif

         // Bounds for s, based on the old values (and old active dofs).
         // This doesn't consider s values from the old inactive dofs, because
         // there were no bounds restriction on them at the previous time step.
         dofs.ComputeElementsMinMax(s, dofs.xe_min, dofs.xe_max,
                                    &s_bool_el, &s_bool_dofs);
         dofs.ComputeBounds(dofs.xe_min, dofs.xe_max,
                            dofs.xi_min, dofs.xi_max, &s_bool_el);

         // Evolve u and get the new active dofs.
         Vector u_new(size);
         add(1.0, u, dt, d_u, u_new);
         Array<bool> s_bool_el_new, s_bool_dofs_new;
         ComputeBoolIndicators(NE, u_new, s_bool_el_new, s_bool_dofs_new);

         fct_solver->CalcFCTProduct(x_gf, lumpedM, d_us_HO, d_us_LO,
                                    dofs.xi_min, dofs.xi_max,
                                    u_new,
                                    s_bool_el_new, s_bool_dofs_new, d_us);

#ifdef REMHOS_FCT_PRODUCT_DEBUG
         Vector us_new(size);
         add(1.0, us, dt, d_us, us_new);
         std::cout << "      out: ";
         ComputeMinMaxS(NE, us_new, u_new, myid);
#endif
      }
      else if (lo_solver) { lo_solver->CalcLOSolution(us, d_us); }
      else if (ho_solver) { ho_solver->CalcHOSolution(us, d_us); }
      else { MFEM_ABORT("No solver was chosen."); }

      d_us.SyncAliasMemory(Y);
   }
}

void AdvectionOperator::UpdateTimeStepEstimate(const Vector &x,
                                               const Vector &dx,
                                               const Vector &x_min,
                                               const Vector &x_max) const
{
   if (dt_control == TimeStepControl::FixedTimeStep) { return; }

   // x_min <= x + dt * dx <= x_max.
   int n = x.Size();
   const double eps = 1e-12;
   double dt = numeric_limits<double>::infinity();

   for (int i = 0; i < n; i++)
   {
      if (dx(i) > eps)
      {
         dt = fmin(dt, (x_max(i) - x(i)) / dx(i) );
      }
      else if (dx(i) < -eps)
      {
         dt = fmin(dt, (x_min(i) - x(i)) / dx(i) );
      }
   }

   MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN,
                 Kbf.ParFESpace()->GetComm());

   dt_est = fmin(dt_est, dt);
}

// // Velocity coefficient
// // Velocity coefficient
// void velocity_function(const Vector &x, Vector &v)
// {
//    int dim = x.Size();

//    // map to the reference [-1,1] domain
//    Vector X(dim);
//    for (int i = 0; i < dim; i++)
//    {
//       double center = (bb_min[i] + bb_max[i]) * 0.5;
//       X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);

//       // std::cout << i << ", bb_min: " << bb_min[i] << ", bb_max: " << bb_max[i] << std::endl;
//       // std::cout << "X(i) ," << i << ", " << X(i) << std::endl;
//    }

//    int ProbExec = problem_num % 20;

//    switch (ProbExec)
//    {
//       case 0:
//       {
//          // Translations in 1D, 2D, and 3D
//          switch (dim)
//          {
//             case 1: v(0) = 1.0; break;
//             case 2: v(0) = sqrt(2./3.); v(1) = sqrt(1./3.); break;
//             case 3: v(0) = sqrt(3./6.); v(1) = sqrt(2./6.); v(2) = sqrt(1./6.);
//                break;
//          }
//          break;
//       }
//       case 1:
//       case 2:
//       case 4:
//       {
//          // Clockwise rotation in 2D around the origin
//          const double w = M_PI/2;
//          switch (dim)
//          {
//             case 1: v(0) = 1.0; break;
//             case 2: v(0) = -w*X(1); v(1) = w*X(0); break;
//             case 3: v(0) = -w*X(1); v(1) = w*X(0); v(2) = 0.0; break;
//          }
//          break;
//       }
//       case 3:
//       {
//          // Clockwise twisting rotation in 2D around the origin
//          const double w = M_PI/2;
//          double d = max((X(0)+1.)*(1.-X(0)),0.) * max((X(1)+1.)*(1.-X(1)),0.);
//          d = d*d;
//          switch (dim)
//          {
//             case 1: v(0) = 1.0; break;
//             case 2: v(0) = d*w*X(1); v(1) = -d*w*X(0); break;
//             case 3: v(0) = d*w*X(1); v(1) = -d*w*X(0); v(2) = 0.0; break;
//          }
//          break;
//       }
//       case 5:
//       {
//          switch (dim)
//          {
//             case 1: v(0) = 1.0; break;
//             case 2: v(0) = 1.0; v(1) = 1.0; break;
//             case 3: v(0) = 1.0; v(1) = 1.0; v(2) = 1.0; break;
//          }
//          break;
//       }
//       case 6:
//       case 7:
//       {
//          switch (dim)
//          {
//             case 1: v(0) = 1.0; break;
//             case 2: v(0) = x(1); v(1) = -x(0); break;
//             case 3: v(0) = x(1); v(1) = -x(0); v(2) = 0.0; break;
//          }
//          break;
//       }
//       case 11:
//       {
//          // Gresho deformation used for mesh motion in remap tests.
//          const double r = sqrt(x(0)*x(0) + x(1)*x(1));
//          if (r < 0.2)
//          {
//             v(0) =  5.0 * x(1);
//             v(1) = -5.0 * x(0);
//          }
//          else if (r < 0.4)
//          {
//             v(0) =  2.0 * x(1) / r - 5.0 * x(1);
//             v(1) = -2.0 * x(0) / r + 5.0 * x(0);
//          }
//          else { v = 0.0; }
//          break;
//       }
//       case 10:
//       case 12:
//       case 13:
//       case 14:
//       case 15:
//       case 16:
//       case 17:
//       {
//          // Taylor-Green deformation used for mesh motion in remap tests.

//          // Map [-1,1] to [0,1].
//          for (int d = 0; d < dim; d++) { X(d) = X(d) * 0.5 + 0.5; }

//          if (dim == 1) { MFEM_ABORT("Not implemented."); }
//          // Clockwise rotation in 2D around the origin
//          const double w = M_PI/2;
//          // v(0) = (x(0)*cos(45*M_PI/180) - x(1)*sin(45*M_PI/180))/1000; 
//          // v(1) = (x(0)*sin(45*M_PI/180) + x(1)*cos(45*M_PI/180))/1000; 

//          v(0) = -w*X(1)*100; 
//          v(1) =  w*X(0)*100;
      
//          if (dim == 3)
//          {
//             v(0) *= 0.1*cos(M_PI*X(2));
//             v(1) *= 0.1*cos(M_PI*X(2));
//             v(2) = 0.0;
//          }
//          break;
//       }
//    }
// }

void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();

   // map to the reference [-1,1] domain
   Vector X(dim);
   for (int i = 0; i < dim; i++)
   {
      double center = (bb_min[i] + bb_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
   }

   // v(0) = v(0);
   // v(1) = v(1);
   // if(dim ==3) {v(2) = v(2);}
   // v(2) = 0.0;
   v *= 1.0;
}

double inflow_function(const Vector &x)
{
   return 0.0;
}

double s0_function(const Vector &x)
{
   return 1.0;
}


#endif // MFEM_USE_MPI
