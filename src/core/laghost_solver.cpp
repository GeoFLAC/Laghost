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

#include "general/forall.hpp"
#include "laghost_solver.hpp"
#include "linalg/kernels.hpp"
#include <unordered_map>
#include <cmath>

#ifdef MFEM_USE_MPI

namespace mfem
{

namespace geodynamics
{

void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    ParGridFunction &gf, const char *title,
                    int x, int y, int w, int h, bool vec)
{
   gf.HostRead();
   ParMesh &pmesh = *gf.ParFESpace()->GetParMesh();
   MPI_Comm comm = pmesh.GetComm();

   int num_procs, myid;
   MPI_Comm_size(comm, &num_procs);
   MPI_Comm_rank(comm, &myid);

   bool newly_opened = false;
   int connection_failed;

   do
   {
      if (myid == 0)
      {
         if (!sock.is_open() || !sock)
         {
            sock.open(vishost, visport);
            sock.precision(8);
            newly_opened = true;
         }
         sock << "solution\n";
      }

      pmesh.PrintAsOne(sock);
      gf.SaveAsOne(sock);

      if (myid == 0 && newly_opened)
      {
         const char* keys = (gf.FESpace()->GetMesh()->Dimension() == 2)
                            ? "mAcRjl" : "mmaaAcl";

         sock << "window_title '" << title << "'\n"
              << "window_geometry "
              << x << " " << y << " " << w << " " << h << "\n"
              << "keys " << keys;
         if ( vec ) { sock << "vvv"; }
         sock << std::endl;
      }

      if (myid == 0)
      {
         connection_failed = !sock && !newly_opened;
      }
      MPI_Bcast(&connection_failed, 1, MPI_INT, 0, comm);
   }
   while (connection_failed);
}

static void Rho0DetJ0Vol(const int dim, const int NE,
                         const IntegrationRule &ir,
                         ParMesh *pmesh,
                         ParFiniteElementSpace &L2,
                         const ParGridFunction &rho0,
                         QuadratureData &qdata,
                         double &volume);

LagrangianGeoOperator::LagrangianGeoOperator(const int size,
                                                 ParFiniteElementSpace &h1,
                                                 ParFiniteElementSpace &l2,
                                                 ParFiniteElementSpace &l2_stress,
                                                 const Array<int> &ess_tdofs,
                                                //  Coefficient &_rho0_coeff,
                                                //  Coefficient &_scale_rho0_coeff,
                                                 ParGridFunction &rho0_gf,
                                                 ParGridFunction &fictitious_rho0_gf,
                                                 ParGridFunction &gamma_gf,
                                                 const int source,
                                                 const double cfl,
                                                 const bool visc,
                                                 const bool vort,
                                                 const bool p_assembly,
                                                 const double cgt,
                                                 const int cgiter,
                                                 double ftz,
                                                 const int oq,
                                                 ParGridFunction &lambda_gf, ParGridFunction &mu_gf, double mscale, const double gravity, const double _thickness, 
                                                 const bool winkler_foundation, const double _winkler_rho, const bool dyn_damping, const double _dyn_factor, Vector _bc_id_pa, const double _vbc_max_val) : // -0-
   TimeDependentOperator(size),
   H1(h1), L2(l2), L2_stress(l2_stress), H1c(H1.GetParMesh(), H1.FEColl(), 1),
   pmesh(H1.GetParMesh()),
   H1Vsize(H1.GetVSize()),
   H1TVSize(H1.TrueVSize()),
   H1GTVSize(H1.GlobalTrueVSize()),
   L2Vsize(L2.GetVSize()),
   L2TVSize(L2.TrueVSize()),
   L2GTVSize(L2.GlobalTrueVSize()),
   block_offsets(5),
   x_gf(&H1),
   ess_tdofs(ess_tdofs),
   dim(pmesh->Dimension()),
   NE(pmesh->GetNE()),
   l2dofs_cnt(L2.GetFE(0)->GetDof()),
   l2_stress_dofs_cnt(L2_stress.GetFE(0)->GetDof()),
   h1dofs_cnt(H1.GetFE(0)->GetDof()),
   source_type(source), cfl(cfl),
   mass_scale(mscale), grav_mag(gravity), thickness(_thickness), winkler_rho(_winkler_rho), dyn_factor(_dyn_factor), vbc_max_val(_vbc_max_val),
   use_viscosity(visc),
   use_vorticity(vort),
   p_assembly(p_assembly),
   winkler_foundation(winkler_foundation),
   dyn_damping(dyn_damping),
   cg_rel_tol(cgt), cg_max_iter(cgiter),ftz_tol(ftz),
   rho0_gf(rho0_gf),
   fictitious_rho0_gf(fictitious_rho0_gf),
   gamma_gf(gamma_gf),
   lambda_gf(lambda_gf),
   mu_gf(mu_gf),
   bc_id_pa(_bc_id_pa),
   // tension_cutoff(_tension_cutoff.Size()),
   // cohesion(_cohesion.Size()),
   // friction_angle(_friction_angle.Size()),
   // dilation_angle(_dilation_angle.Size()),
   Mv(&H1), Mv_spmat_copy(),
   fic_Mv(&H1), fic_Mv_spmat_copy(),
   Me(l2dofs_cnt, l2dofs_cnt, NE),
   Me_inv(l2dofs_cnt, l2dofs_cnt, NE),
   rho0_coeff(&rho0_gf),
   scale_rho0_coeff(&fictitious_rho0_gf),
   ir(IntRules.Get(pmesh->GetElementBaseGeometry(0),
                   (oq > 0) ? oq : 3 * H1.GetOrder(0) + L2.GetOrder(0) - 1)),
   Q1D(int(floor(0.7 + pow(ir.GetNPoints(), 1.0 / dim)))),
   qdata(dim, NE, ir.GetNPoints()),
   qdata_is_current(false),
   forcemat_is_assembled(false),
   gmat_is_assembled(false),
   Force(&l2, &h1),
   // Force(&L2, &H1),
   // Body_Force(nullptr),
   Body_Force(&H1),
   ForcePA(nullptr), VMassPA(nullptr), EMassPA(nullptr), StressPA(nullptr),
   VMassPA_Jprec(nullptr),
   CG_VMass(H1.GetParMesh()->GetComm()),
   CG_EMass(L2.GetParMesh()->GetComm()),
   timer(p_assembly ? L2TVSize : 1),
   qupdate(nullptr),
   X(H1c.GetTrueVSize()),
   B(H1c.GetTrueVSize()),
   one(L2Vsize),
   rhs(H1Vsize),
   e_rhs(L2Vsize),
   s_rhs(L2Vsize), // rhs for stress vector
   rhs_c_gf(&H1c),
   dvc_gf(&H1c)
{
   // If you add block vector, you should add offset 
   block_offsets[0] = 0;
   block_offsets[1] = block_offsets[0] + H1.GetVSize();
   block_offsets[2] = block_offsets[1] + H1.GetVSize();
   block_offsets[3] = block_offsets[2] + L2.GetVSize();
   block_offsets[4] = block_offsets[3] + L2.GetVSize()*3*(dim-1);
   // block_offsets[5] = block_offsets[4] + H1.GetVSize();
   one.UseDevice(true);
   one = 1.0;
   qdata.mscale  = mscale;
   qdata.vbc_max_val  = vbc_max_val;
   qdata.gravity = gravity;
   qdata.h_est   = 1e+38;

   // for (int i = 0; i < pmesh->attributes.Max(); i++)
   // {
   //    if(_tension_cutoff[i] == 0)
   //    {
   //       tension_cutoff[i] = _cohesion[i]/tan(M_PI*_friction_angle[i]/180.0); // Tension limit
   //    }
   //    else
   //    {
   //       tension_cutoff[i] = _tension_cutoff[i];
   //    }
      
   //    cohesion[i] = _cohesion[i];
   //    friction_angle[i] = M_PI*_friction_angle[i]/180.0;
   //    dilation_angle[i] = M_PI*_dilation_angle[i]/180.0;
   // }
   
   if (p_assembly)
   {
      qupdate = new QUpdate(dim, NE, Q1D, visc, vort, cfl, &timer, gamma_gf, lambda_gf, mu_gf, ir, H1, L2, L2_stress); // -1-
      // qupdate = new QUpdate(dim, NE, Q1D, visc, vort, cfl,
      //                       &timer, gamma_gf, ir, H1, L2);                            
      ForcePA = new ForcePAOperator(qdata, H1, L2, ir);
      StressPA = new StressPAOperator(qdata, H1, L2, ir); // stress rate operator, slee
      VMassPA = new MassPAOperator(H1c, ir, rho0_coeff);
      // VMassPA = new MassPAOperator(H1c, ir, scale_rho0_coeff);
      EMassPA = new MassPAOperator(L2, ir, rho0_coeff);
      // Inside the above constructors for mass, there is reordering of the mesh
      // nodes which is performed on the host. Since the mesh nodes are a
      // subvector, so we need to sync with the rest of the base vector (which
      // is assumed to be in the memory space used by the mfem::Device).
      H1.GetParMesh()->GetNodes()->ReadWrite();
      // Attributes 1/2/3 correspond to fixed-x/y/z boundaries, i.e.,
      // we must enforce v_x/y/z = 0 for the velocity components.
      const int bdr_attr_max = H1.GetMesh()->bdr_attributes.Max();
      Array<int> ess_bdr(bdr_attr_max);
   
      for (int c = 0; c < dim; c++)
      {
         ess_bdr = 0;

         if(dim == 2)
         {
            if(c == 0)
            {
               for (int i = 0; i < bc_id_pa.Size(); ++i) {if(bc_id_pa[i]==1 || bc_id_pa[i]==3){ess_bdr[i] = 1;}}
               H1c.GetEssentialTrueDofs(ess_bdr, c_tdofs[c]);
               c_tdofs[c].Read();
            }
            else if(c == 1)
            {
               for (int i = 0; i < bc_id_pa.Size(); ++i) {if(bc_id_pa[i]==2 || bc_id_pa[i]==3){ess_bdr[i] = 1;}}
               H1c.GetEssentialTrueDofs(ess_bdr, c_tdofs[c]);
               c_tdofs[c].Read();
            }
         }
         else
         {
            if(c == 0)
            {
               for (int i = 0; i < bc_id_pa.Size(); ++i) {if(bc_id_pa[i]==1 || bc_id_pa[i]==4 || bc_id_pa[i]==5 || bc_id_pa[i]==6){ess_bdr[i] = 1;}}
               H1c.GetEssentialTrueDofs(ess_bdr, c_tdofs[c]);
               c_tdofs[c].Read();
            }
            else if(c == 1)
            {
               for (int i = 0; i < bc_id_pa.Size(); ++i) {if(bc_id_pa[i]==2 || bc_id_pa[i]==4 || bc_id_pa[i]==5 || bc_id_pa[i]==7){ess_bdr[i] = 1;}}
               H1c.GetEssentialTrueDofs(ess_bdr, c_tdofs[c]);
               c_tdofs[c].Read();
            }
            else
            {
               for (int i = 0; i < bc_id_pa.Size(); ++i) {if(bc_id_pa[i]==3 || bc_id_pa[i]==4 || bc_id_pa[i]==6 || bc_id_pa[i]==7){ess_bdr[i] = 1;}}
               H1c.GetEssentialTrueDofs(ess_bdr, c_tdofs[c]);
               c_tdofs[c].Read();
            }
            
         }
      }
      
      X.UseDevice(true);
      B.UseDevice(true);
      rhs.UseDevice(true);
      e_rhs.UseDevice(true);
      s_rhs.UseDevice(true); // for gpu calculation, slee

      // Standard local (full) assembly and inversion for energy mass matrices.
      // 'Me' is used in the computation of stress
      // std::cout << "Standard local assembly and inversion for energy mass matrices" << std::endl;
      MassIntegrator mi(rho0_coeff, &ir);
      for (int e = 0; e < NE; e++)
      {
         DenseMatrixInverse inv(&Me(e));
         const FiniteElement &fe = *L2.GetFE(e);
         ElementTransformation &Tr = *L2.GetElementTransformation(e);
         mi.AssembleElementMatrix(fe, Tr, Me(e));
         inv.Factor();
         inv.GetInverseMatrix(Me_inv(e));
      }

   }
   else
   {
      // Standard local assembly and inversion for energy mass matrices.
      // 'Me' is used in the computation of the internal energy
      // which is used twice: once at the start and once at the end of the run.
      MassIntegrator mi(rho0_coeff, &ir);
      for (int e = 0; e < NE; e++)
      {
         DenseMatrixInverse inv(&Me(e));
         const FiniteElement &fe = *L2.GetFE(e);
         ElementTransformation &Tr = *L2.GetElementTransformation(e);
         mi.AssembleElementMatrix(fe, Tr, Me(e));
         inv.Factor();
         inv.GetInverseMatrix(Me_inv(e));
      }
      
      // Standard assembly for the velocity mass matrix.

      VectorMassIntegrator *vmi = new VectorMassIntegrator(rho0_coeff, &ir);
      Mv.AddDomainIntegrator(vmi);
      Mv.Assemble();
      Mv_spmat_copy = Mv.SpMat();

      VectorMassIntegrator *vmi_scale = new VectorMassIntegrator(scale_rho0_coeff, &ir);
      fic_Mv.AddDomainIntegrator(vmi_scale);
      fic_Mv.Assemble();
      fic_Mv_spmat_copy = Mv.SpMat();
   }

   // Values of rho0DetJ0 and Jac0inv at all quadrature points.
   // Initial local mesh size (assumes all mesh elements are the same).
   int Ne, ne = NE;
   double Volume, vol = 0.0;
   if (dim > 1 && p_assembly)
   {
      Rho0DetJ0Vol(dim, NE, ir, pmesh, L2, rho0_gf, qdata, vol);
   }
   else
   {
      const int NQ = ir.GetNPoints();
      Vector rho_vals(NQ);
      for (int e = 0; e < NE; e++)
      {
         rho0_gf.GetValues(e, ir, rho_vals);
         ElementTransformation &Tr = *H1.GetElementTransformation(e);
         for (int q = 0; q < NQ; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            Tr.SetIntPoint(&ip);
            DenseMatrixInverse Jinv(Tr.Jacobian());
            Jinv.GetInverseMatrix(qdata.Jac0inv(e*NQ + q));
            const double rho0DetJ0 = Tr.Weight() * rho_vals(q);
            qdata.rho0DetJ0w(e*NQ + q) = rho0DetJ0 * ir.IntPoint(q).weight;
         }
      }
      for (int e = 0; e < NE; e++) { vol += pmesh->GetElementVolume(e); }
   }
   MPI_Allreduce(&vol, &Volume, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
   MPI_Allreduce(&ne, &Ne, 1, MPI_INT, MPI_SUM, pmesh->GetComm());
   switch (pmesh->GetElementBaseGeometry(0))
   {
      case Geometry::SEGMENT: qdata.h0 = Volume / Ne; break;
      case Geometry::SQUARE: qdata.h0 = sqrt(Volume / Ne); break;
      case Geometry::TRIANGLE: qdata.h0 = sqrt(2.0 * Volume / Ne); break;
      case Geometry::CUBE: qdata.h0 = pow(Volume / Ne, 1./3.); break;
      case Geometry::TETRAHEDRON: qdata.h0 = pow(6.0 * Volume / Ne, 1./3.); break;
      default: MFEM_ABORT("Unknown zone type!");
   }
   qdata.h0 /= (double) H1.GetOrder(0);

   if (p_assembly)
   {
      // Setup the preconditioner of the velocity mass operator.
      // BC are handled by the VMassPA, so ess_tdofs here can be empty.
      Array<int> empty_tdofs;
      VMassPA_Jprec = new OperatorJacobiSmoother(VMassPA->GetBF(), empty_tdofs);
      CG_VMass.SetPreconditioner(*VMassPA_Jprec);

      CG_VMass.SetOperator(*VMassPA);
      CG_VMass.SetRelTol(cg_rel_tol);
      CG_VMass.SetAbsTol(0.0);
      CG_VMass.SetMaxIter(cg_max_iter);
      CG_VMass.SetPrintLevel(-1);

      CG_EMass.SetOperator(*EMassPA);
      CG_EMass.iterative_mode = false;
      CG_EMass.SetRelTol(cg_rel_tol);
      CG_EMass.SetAbsTol(0.0);
      CG_EMass.SetMaxIter(cg_max_iter);
      CG_EMass.SetPrintLevel(-1);
   }
   else
   {
      ForceIntegrator *fi = new ForceIntegrator(qdata);
      fi->SetIntRule(&ir);
      Force.AddDomainIntegrator(fi);
      // Make a dummy assembly to figure out the sparsity.
      Force.Assemble(0);
      Force.Finalize(0);
      // Standarad fully aseembly for body force integarator
      BodyForceIntegrator *bi = new BodyForceIntegrator(qdata);
      bi->SetIntRule(&ir);
      Body_Force.AddDomainIntegrator(bi);
      Body_Force.Assemble();
   }
}

LagrangianGeoOperator::~LagrangianGeoOperator()
{
   delete qupdate;
   if (p_assembly)
   {
      delete EMassPA;
      delete VMassPA;
      delete VMassPA_Jprec;
      delete ForcePA;
      delete StressPA;
      // delete BodyFPA;
   }
}

void LagrangianGeoOperator::Mult(const Vector &S, Vector &dS_dt) const
{
   // Make sure that the mesh positions correspond to the ones in S. This is
   // needed only because some mfem time integrators don't update the solution
   // vector at every intermediate stage (hence they don't change the mesh).
   UpdateMesh(S); // comment out for testing pseudo transient loop

   // The monolithic BlockVector stores the unknown fields as follows:
   // (Position, Velocity, Specific Internal Energy).
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction v;
   const int VsizeH1 = H1.GetVSize();
   v.MakeRef(&H1, *sptr, VsizeH1);
   // Set dx_dt = v (explicit).
   ParGridFunction dx; // comment out for testing pseudo transient loop
   dx.MakeRef(&H1, dS_dt, 0); // comment out for testing pseudo transient loop
   dx = v; // comment out for testing pseudo transient loop
   
   SolveVelocity(S, dS_dt);
   SolveEnergy(S, v, dS_dt);
   SolveStress(S, dS_dt);
   qdata_is_current = false;
}


void LagrangianGeoOperator::SolveVelocity(const Vector &S, Vector &dS_dt) const
{
   UpdateQuadratureData(S);
   AssembleForceMatrix();

   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction v;
   v.MakeRef(&H1, *sptr, H1.GetVSize());
   
   Vector vel_mag(H1.GetVSize()/dim);
   for (int i = 0; i < H1.GetVSize()/dim; i++)
   {
      if(dim == 2){vel_mag[i] = sqrt(pow(v[i], 2) + pow(v[i+H1.GetVSize()/dim], 2));}
      else{vel_mag[i] = sqrt(pow(v[i], 2) + pow(v[i+H1.GetVSize()/dim], 2) + pow(v[i+2*H1.GetVSize()/dim], 2));}
   }

   double bulkm = lambda_gf.Max() + 2 * mu_gf.Max(); 
   double denm  = rho0_gf.Max();

   // double pseudo_speed = global_max_vel * mass_scale;
   double pseudo_speed = vbc_max_val * mass_scale;
   double mfactor = (pseudo_speed * pseudo_speed) / bulkm; 
   if(mfactor > 1.0){mfactor = 1.0;} // if mass scaling is greater than unity, mass scaling is off.
   if(mass_scale == 1.0){mfactor = 1.0;} // no mass scaling
   mfactor = mfactor * denm; // p_assembly requires small cfl < 0.5 (emperically)
   // The monolithic BlockVector stores the unknown fields as follows:
   // (Position, Velocity, Specific Internal Energy).
   ParGridFunction dv;
   dv.MakeRef(&H1, dS_dt, H1.GetVSize());
   dv = 0.0;

   // AMR requires resize of vectors
   // Solve for velocity.
   one.SetSize(L2.GetVSize());
   rhs.SetSize(H1.GetVSize());
   // Vector B, X; In AMR B and X should be redefined.
   one = 1.0; rhs = 0.0;
   B.HostRead();

   // Body Force vector (F = 1 * g)
   ParGridFunction accel_src_gf;
   accel_src_gf.SetSpace(&H1);
   GTCoefficient accel_coeff(dim);
   accel_src_gf.ProjectCoefficient(accel_coeff);
   accel_src_gf.Read();
   accel_src_gf *= grav_mag;

   if (p_assembly)
   {
      // Damping vector (F = -1 * sign(v))
      ParGridFunction sign_gf;
      sign_gf.SetSpace(&H1);
      sign_gf = 0.0; sign_gf += v;
      
      timer.sw_force.Start();
      ForcePA->Mult(one, rhs); // F*1
      timer.sw_force.Stop();
      rhs.Neg(); // -F

      if(winkler_foundation)
      {
         // ParGridFunction x;
         // x.MakeRef(&H1, *sptr, H1.GetVSize()*0);
         // Vector WC(H1.GetVSize()); WC.UseDevice(true); 
         // Vector WD(H1.GetVSize()); WD.UseDevice(true); 
         // Vector WK(H1.GetVSize()); WK.UseDevice(true); 
         // WC = x; WD = thickness; WC -= WD; WC.Neg(); 

         Array<int> nbc_bdr(pmesh->bdr_attributes.Max());   
         nbc_bdr = 0; nbc_bdr[2] = 1; // bottome boundary
         LinearForm winkler(&H1);

         VectorArrayCoefficient winkler_load(dim);
         for (int i = 0; i < dim-1; i++)
         {
            winkler_load.Set(i, new ConstantCoefficient(0.0));
         }

         Vector pull_force(pmesh->bdr_attributes.Max());
         pull_force = 0.0;
         pull_force(2) = 1.0;
         winkler_load.Set(dim-1, new PWConstCoefficient(pull_force));

         winkler.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(winkler_load), nbc_bdr);
         winkler.Assemble();

         // Applying winkler      
         Winkler(S, winkler, thickness);
         rhs.Add(winkler_rho*grav_mag,  winkler);

         // Multvv(WC, winkler, WK);
         // rhs.Add(winkler_rho*grav_mag,  WK);
      }

      // Partial assembly solve for each velocity component
      const int size = H1c.GetVSize();
      const Operator *Pconf = H1c.GetProlongationMatrix();
      for (int c = 0; c < dim; c++)
      {
         dvc_gf.MakeRef(&H1c, dS_dt, H1Vsize + c*size);
         rhs_c_gf.MakeRef(&H1c, rhs, c*size);

         if (Pconf) { Pconf->MultTranspose(rhs_c_gf, B); }
         else { B = rhs_c_gf; }

         if(c == dim -1)
         {
            ParGridFunction accel_comp;
            accel_comp.MakeRef(&H1c, accel_src_gf, c*size);
            Vector AC;
            accel_comp.GetTrueDofs(AC);
            Vector BA(AC.Size());
            VMassPA->MultFull(AC, BA);
            B += BA; // -(F + rho*g)
         }
         
         // // Applying damping for all forces such internal, external, and body
         // if(dyn_damping)
         // {
         //    ParGridFunction sign_comp;
         //    sign_comp.MakeRef(&H1c, sign_gf, c*size);
         //    Vector AC; sign_comp.GetTrueDofs(AC);
         //    AC.Signcopy(); AC.Neg();
         //    Vector damping(B.Size()); Vector sdamping(B.Size());  
         //    damping.UseDevice(true); sdamping.UseDevice(true);
         //    damping = 0.0; damping += B; sdamping = 0.0;
         //    damping.Abs(); damping *=dyn_factor;
         //    Multvv(AC, damping, sdamping);
         //    damping.HostReadWrite(); sdamping.HostReadWrite();
         //    B += sdamping;
            
         //    // Getdamping_comp(S, c, damping);
         //    // B += damping;
         // }

         // Applying damping for all forces such internal, external, and body
         if(dyn_damping)
         {
            Vector damping(B.Size());
            damping=0.0;
            damping.Add(1.0, B);
            Getdamping_comp(S, c, damping);
            B.Add(dyn_factor, damping);
         }
         
         H1c.GetRestrictionMatrix()->Mult(dvc_gf, X);
         VMassPA->SetEssentialTrueDofs(c_tdofs[c]);
         VMassPA->EliminateRHS(B);

         timer.sw_cgH1.Start();
         CG_VMass.Mult(B, X);
         timer.sw_cgH1.Stop();
         timer.H1iter += CG_VMass.GetNumIterations();
         X*=mfactor;
         if (Pconf) { Pconf->Mult(X, dvc_gf); }
         else { dvc_gf = X; }
         // We need to sync the subvector 'dvc_gf' with its base vector
         // because it may have been moved to a different memory space.
         dvc_gf.GetMemory().SyncAlias(dS_dt.GetMemory(), dvc_gf.Size());
      }
   }

   else // full assembly
   {
      timer.sw_force.Start();
      Force.Mult(one, rhs);
      timer.sw_force.Stop();
      rhs.Neg();
      // rhs += *v_source; 

      // // Applying body force
      // Body_Force.Update();
      // Body_Force.Assemble();
      // rhs.Add(1.0,  Body_Force);

      Vector rhs_accel(rhs.Size());
      Mv_spmat_copy.Mult(accel_src_gf, rhs_accel);
      rhs += rhs_accel; 

      if(winkler_foundation)
      {
         Array<int> nbc_bdr(pmesh->bdr_attributes.Max());   
         nbc_bdr = 0; nbc_bdr[2] = 1; // bottome boundary
         LinearForm winkler(&H1);

         VectorArrayCoefficient winkler_load(dim);
         for (int i = 0; i < dim-1; i++)
         {
            winkler_load.Set(i, new ConstantCoefficient(0.0));
         }

         Vector pull_force(pmesh->bdr_attributes.Max());
         pull_force = 0.0;
         pull_force(2) = 1.0;
         winkler_load.Set(dim-1, new PWConstCoefficient(pull_force));

         winkler.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(winkler_load), nbc_bdr);
         winkler.Assemble();

         // Applying winkler      
         Winkler(S, winkler, thickness);
         rhs.Add(winkler_rho*grav_mag,  winkler);
      }

      HypreParMatrix A;
      // Mv.FormLinearSystem(ess_tdofs, dv, rhs, A, X, B);
      fic_Mv.FormLinearSystem(ess_tdofs, dv, rhs, A, X, B);

      // Applying damping for all forces such internal, external, and body
      if(dyn_damping)
      {
         Vector damping(B.Size());
         damping=0.0;
         damping.Add(1.0, B);
         Getdamping(S, damping);
         B.Add(dyn_factor, damping);
         //
      }

      CGSolver cg(H1.GetParMesh()->GetComm());
      HypreSmoother prec;
      // prec.SetType(HypreSmoother::Jacobi, 1); // l1-scaled Jacobi 
      prec.SetType(HypreSmoother::Chebyshev);  // Chebyshev
      cg.SetPreconditioner(prec);
      cg.SetOperator(A);
      cg.SetRelTol(cg_rel_tol);
      cg.SetAbsTol(0.0);
      cg.SetMaxIter(cg_max_iter);
      cg.SetPrintLevel(-1);
      timer.sw_cgH1.Start();
      cg.Mult(B, X);
      timer.sw_cgH1.Stop();
      // X*=mfactor;
      timer.H1iter += cg.GetNumIterations();
      // Mv.RecoverFEMSolution(X, rhs, dv);
      fic_Mv.RecoverFEMSolution(X, rhs, dv);
   }
}

void LagrangianGeoOperator::SolveEnergy(const Vector &S, const Vector &v, Vector &dS_dt) const
{
   // UpdateQuadratureData(S, dt);
   UpdateQuadratureData(S);
   AssembleForceMatrix();

   // The monolithic BlockVector stores the unknown fields as follows:
   // (Position, Velocity, Specific Internal Energy).
   ParGridFunction de;
   de.MakeRef(&L2, dS_dt, H1.GetVSize()*2);
   de = 0.0;

   // resize for e_rhs
   /*
   AMR requires resize of vectors
   e_rhs.SetSize(L2.GetVSize());
   e_rhs = 0.0;
   */
   
   // Solve for energy, assemble the energy source if such exists.
   LinearForm *e_source = nullptr;
   if (source_type == 1) // 2D Taylor-Green.
   {
      // Needed since the Assemble() defaults to PA.
      L2.GetMesh()->DeleteGeometricFactors();
      e_source = new LinearForm(&L2);
      TaylorCoefficient coeff;
      DomainLFIntegrator *d = new DomainLFIntegrator(coeff, &ir);
      e_source->AddDomainIntegrator(d);
      e_source->Assemble();
   }

   // TODO
   // MassIntegrator mi(&ir);
   // Me.SetSize(l2dofs_cnt, l2dofs_cnt, NE);
   // Me_inv.SetSize(l2dofs_cnt, l2dofs_cnt, NE);
   // MassIntegrator mi(rho0_coeff, &ir);
      // for (int e = 0; e < NE; e++)
      // {
      //    DenseMatrixInverse inv(&Me(e));
      //    const FiniteElement &fe = *L2.GetFE(e);
      //    ElementTransformation &Tr = *L2.GetElementTransformation(e);
      //    mi.AssembleElementMatrix(fe, Tr, Me(e));
      //    inv.Factor();
      //    inv.GetInverseMatrix(Me_inv(e));
      // }
      
   Array<int> l2dofs;
   if (p_assembly)
   {
      timer.sw_force.Start();
      ForcePA->MultTranspose(v, e_rhs);
      timer.sw_force.Stop();
      if (e_source) { e_rhs += *e_source; }
      timer.sw_cgL2.Start();
      CG_EMass.Mult(e_rhs, de);
      timer.sw_cgL2.Stop();
      const HYPRE_Int cg_num_iter = CG_EMass.GetNumIterations();
      timer.L2iter += (cg_num_iter==0) ? 1 : cg_num_iter;
      // Move the memory location of the subvector 'de' to the memory
      // location of the base vector 'dS_dt'.
      de.GetMemory().SyncAlias(dS_dt.GetMemory(), de.Size());
   }
   else // not p_assembly
   {
      timer.sw_force.Start();
      Force.MultTranspose(v, e_rhs);
      timer.sw_force.Stop();

      if (e_source) { e_rhs += *e_source; }
      Vector loc_rhs(l2dofs_cnt), loc_de(l2dofs_cnt);
      for (int e = 0; e < NE; e++)
      {
         L2.GetElementDofs(e, l2dofs);
         e_rhs.GetSubVector(l2dofs, loc_rhs);
         timer.sw_cgL2.Start();
         Me_inv(e).Mult(loc_rhs, loc_de);
         timer.sw_cgL2.Stop();
         timer.L2iter += 1;
         de.SetSubVector(l2dofs, loc_de);
      }
   }
   delete e_source;
}

void LagrangianGeoOperator::SolveStress(const Vector &S, Vector &dS_dt) const
{
   UpdateQuadratureData(S);
   
   if (p_assembly)
   {
      ParGridFunction v(&H1);
      v = 1.0; // dummy velocity for partial assembly
      // Stress rate devided by mass matrix
      for (int i = 0; i < 3*(dim - 1); i++)
      {
         // Partial assembly solve for each stress component
         timer.sw_force.Start();
         const int comp = i;
         StressPA->MultTranspose(v, s_rhs, comp);
         timer.sw_force.Stop();
         
         ParGridFunction ds;
         ds.MakeRef(&L2, dS_dt, H1Vsize*2 + L2Vsize + L2Vsize*i);
         ds = 0.0;
         timer.sw_cgL2.Start();
         CG_EMass.Mult(s_rhs, ds);
         timer.sw_cgL2.Stop();
         const HYPRE_Int cg_num_iter = CG_EMass.GetNumIterations();
         timer.L2iter += (cg_num_iter==0) ? 1 : cg_num_iter;
         // Move the memory location of the subvector 'ds' to the memory
         // location of the base vector 'dS_dt'.
         ds.GetMemory().SyncAlias(dS_dt.GetMemory(), ds.Size());
      }
   }
   else
   {
      ParGridFunction dsig;
      dsig.MakeRef(&L2_stress, dS_dt, H1.GetVSize()*2 + L2.GetVSize());
      int NED = NE*l2_stress_dofs_cnt;
      int dim2 = dim*dim;
      if(dim == 2)
      {
         Vector sub_rhs1(l2_stress_dofs_cnt), sub_rhs2(l2_stress_dofs_cnt), sub_rhs3(l2_stress_dofs_cnt);
         Vector loc_dsig1(l2_stress_dofs_cnt), loc_dsig2(l2_stress_dofs_cnt), loc_dsig3(l2_stress_dofs_cnt);
         loc_dsig1=0.0; loc_dsig2=0.0; loc_dsig3=0.0;
         Array<int> offset(4); 
         offset[0] = 0;
         offset[1] = offset[0] + l2_stress_dofs_cnt;
         offset[2] = offset[1] + l2_stress_dofs_cnt;
         offset[3] = offset[2] + l2_stress_dofs_cnt;

         BlockVector loc_rhs(offset, Device::GetMemoryType());

         sub_rhs1.MakeRef(loc_rhs, 0*l2_stress_dofs_cnt);
         sub_rhs2.MakeRef(loc_rhs, 1*l2_stress_dofs_cnt);
         sub_rhs3.MakeRef(loc_rhs, 2*l2_stress_dofs_cnt);
         loc_rhs = 0.0;
   
         SigmaIntegrator gi(qdata);
         gi.SetIntRule(&ir);
         Array<int> dof_loc1(l2_stress_dofs_cnt);
         Array<int> dof_loc2(l2_stress_dofs_cnt);
         Array<int> dof_loc3(l2_stress_dofs_cnt);
      
         for (int e = 0; e < NE; e++)
         {

            L2_stress.GetElementDofs(e, dof_loc1);
            L2_stress.GetElementDofs(e, dof_loc2);
            L2_stress.GetElementDofs(e, dof_loc3);
         
            const FiniteElement &fe = *L2_stress.GetFE(e);
            ElementTransformation &eltr = *L2_stress.GetElementTransformation(e);
            gi.AssembleRHSElementVect(fe, eltr, loc_rhs);

            Me_inv(e).Mult(sub_rhs1, loc_dsig1);
            Me_inv(e).Mult(sub_rhs2, loc_dsig2);
            Me_inv(e).Mult(sub_rhs3, loc_dsig3);
   
            for (int i = 0; i < dof_loc1.Size(); i++)
            {
               dof_loc1[i] = i + (e+0*NE)*dof_loc1.Size(); 
               dof_loc2[i] = i + (e+1*NE)*dof_loc1.Size();
               dof_loc3[i] = i + (e+2*NE)*dof_loc1.Size(); 
            }

            dsig.SetSubVector(dof_loc1, loc_dsig1);
            dsig.SetSubVector(dof_loc2, loc_dsig2);
            dsig.SetSubVector(dof_loc3, loc_dsig3);
         }
      }
      else if(dim == 3)
      {
         Vector sub_rhs1(l2_stress_dofs_cnt), sub_rhs2(l2_stress_dofs_cnt), sub_rhs3(l2_stress_dofs_cnt);
         Vector sub_rhs4(l2_stress_dofs_cnt), sub_rhs5(l2_stress_dofs_cnt), sub_rhs6(l2_stress_dofs_cnt);

         Vector loc_dsig1(l2_stress_dofs_cnt), loc_dsig2(l2_stress_dofs_cnt), loc_dsig3(l2_stress_dofs_cnt);
         Vector loc_dsig4(l2_stress_dofs_cnt), loc_dsig5(l2_stress_dofs_cnt), loc_dsig6(l2_stress_dofs_cnt);

         loc_dsig1=0.0; loc_dsig2=0.0; loc_dsig3=0.0; loc_dsig4=0.0; loc_dsig5=0.0; loc_dsig6=0.0;
         Array<int> offset(7); 
         offset[0] = 0;
         offset[1] = offset[0] + l2_stress_dofs_cnt;
         offset[2] = offset[1] + l2_stress_dofs_cnt;
         offset[3] = offset[2] + l2_stress_dofs_cnt;
         offset[4] = offset[3] + l2_stress_dofs_cnt;
         offset[5] = offset[4] + l2_stress_dofs_cnt;
         offset[6] = offset[5] + l2_stress_dofs_cnt;
         
         BlockVector loc_rhs(offset, Device::GetMemoryType());

         sub_rhs1.MakeRef(loc_rhs, 0*l2_stress_dofs_cnt);
         sub_rhs2.MakeRef(loc_rhs, 1*l2_stress_dofs_cnt);
         sub_rhs3.MakeRef(loc_rhs, 2*l2_stress_dofs_cnt);
         sub_rhs4.MakeRef(loc_rhs, 3*l2_stress_dofs_cnt);
         sub_rhs5.MakeRef(loc_rhs, 4*l2_stress_dofs_cnt);
         sub_rhs6.MakeRef(loc_rhs, 5*l2_stress_dofs_cnt);
         
         loc_rhs = 0.0;
   
         SigmaIntegrator gi(qdata);
         gi.SetIntRule(&ir);
         Array<int> dof_loc1(l2_stress_dofs_cnt);
         Array<int> dof_loc2(l2_stress_dofs_cnt);
         Array<int> dof_loc3(l2_stress_dofs_cnt);
         Array<int> dof_loc4(l2_stress_dofs_cnt);
         Array<int> dof_loc5(l2_stress_dofs_cnt);
         Array<int> dof_loc6(l2_stress_dofs_cnt);
         
         for (int e = 0; e < NE; e++)
         {

            L2_stress.GetElementDofs(e, dof_loc1);
            L2_stress.GetElementDofs(e, dof_loc2);
            L2_stress.GetElementDofs(e, dof_loc3);
            L2_stress.GetElementDofs(e, dof_loc4);
            L2_stress.GetElementDofs(e, dof_loc5);
            L2_stress.GetElementDofs(e, dof_loc6);
            
            const FiniteElement &fe = *L2_stress.GetFE(e);
            ElementTransformation &eltr = *L2_stress.GetElementTransformation(e);
            gi.AssembleRHSElementVect(fe, eltr, loc_rhs);

            Me_inv(e).Mult(sub_rhs1, loc_dsig1);
            Me_inv(e).Mult(sub_rhs2, loc_dsig2);
            Me_inv(e).Mult(sub_rhs3, loc_dsig3);
            Me_inv(e).Mult(sub_rhs4, loc_dsig4);
            Me_inv(e).Mult(sub_rhs5, loc_dsig5);
            Me_inv(e).Mult(sub_rhs6, loc_dsig6);

            for (int i = 0; i < dof_loc1.Size(); i++)
            {
               dof_loc1[i] = i + (e+0*NE)*dof_loc1.Size(); 
               dof_loc2[i] = i + (e+1*NE)*dof_loc1.Size();
               dof_loc3[i] = i + (e+2*NE)*dof_loc1.Size(); 
               dof_loc4[i] = i + (e+3*NE)*dof_loc1.Size(); 
               dof_loc5[i] = i + (e+4*NE)*dof_loc1.Size(); 
               dof_loc6[i] = i + (e+5*NE)*dof_loc1.Size(); 
            }

            dsig.SetSubVector(dof_loc1, loc_dsig1);
            dsig.SetSubVector(dof_loc2, loc_dsig2);
            dsig.SetSubVector(dof_loc3, loc_dsig3);
            dsig.SetSubVector(dof_loc4, loc_dsig4);
            dsig.SetSubVector(dof_loc5, loc_dsig5);
            dsig.SetSubVector(dof_loc6, loc_dsig6);
         }
      }
   }
}

void LagrangianGeoOperator::UpdateMesh(const Vector &S) const
{
   Vector* sptr = const_cast<Vector*>(&S);
   x_gf.MakeRef(&H1, *sptr, 0);
   H1.GetParMesh()->NewNodes(x_gf, false);
}

void LagrangianGeoOperator::Getdamping(const Vector &S, Vector &_v_damping) const
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction v;
   v.MakeRef(&H1, *sptr, H1.GetVSize());
   v.SetTrueVector();
   Vector vt = v.GetTrueVector();

   for( int i = 0; i < _v_damping.Size(); i++ )
   {  
      _v_damping[i]=-1*copysignl(1.0, vt[i])*fabs(_v_damping[i]); //
   }
}

void LagrangianGeoOperator::Getdamping_comp(const Vector &S, const int &comp, Vector &_v_damping) const
{
   // damping by each component in p_assembly
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction v;
   v.MakeRef(&H1, *sptr, H1.GetVSize());
   v.SetTrueVector();
   Vector vt = v.GetTrueVector();
   for( int i = 0; i < _v_damping.Size(); i++ )
   {
      _v_damping[i]=-1*copysignl(1.0, vt[i+_v_damping.Size()*comp])*fabs(_v_damping[i]); //
   }
}


void LagrangianGeoOperator::Winkler(const Vector &S, Vector &_winkler, double &_thickness) const
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x;
   // ParGridFunction x, x_ini;
   x.MakeRef(&H1, *sptr, H1.GetVSize()*0);
   double bottom = 0.0;
   // x_ini.MakeRef(&H1, *sptr, H1.GetVSize()*2 + L2.GetVSize() + L2.GetVSize()*3*(dim-1));
   for( int i = 0; i < _winkler.Size(); i++ )
   {
      
      // _winkler[i]=-1.0*(x[i]-(bottom + _thickness))*_winkler[i]; //
      _winkler[i]= (_thickness - x[i])*_winkler[i];
      
      // p = var.compensation_pressure -
      //               (var.mat->rho(e) + param.bc.winkler_delta_rho) *
      //               param.control.gravity * (zcenter + param.mesh.zlength);
   }
}

double LagrangianGeoOperator::GetTimeStepEstimate(const Vector &S) const
{
   UpdateMesh(S);
   UpdateQuadratureData(S);
   double glob_dt_est;
   const MPI_Comm comm = H1.GetParMesh()->GetComm();
   MPI_Allreduce(&qdata.dt_est, &glob_dt_est, 1, MPI_DOUBLE, MPI_MIN, comm);
   return glob_dt_est;
}

double LagrangianGeoOperator::GetLengthEstimate(const Vector &S) const
{
   UpdateMesh(S);
   UpdateQuadratureData(S);
   double glob_h_est;
   const MPI_Comm comm = H1.GetParMesh()->GetComm();
   MPI_Allreduce(&qdata.h_est, &glob_h_est, 1, MPI_DOUBLE, MPI_MIN, comm);
   return glob_h_est;
}

void LagrangianGeoOperator::ResetTimeStepEstimate() const
{
   qdata.dt_est = std::numeric_limits<double>::infinity();
}

void LagrangianGeoOperator::ComputeDensity(ParGridFunction &rho) const
{
   rho.SetSpace(&L2);
   DenseMatrix Mrho(l2dofs_cnt);
   Vector rhs(l2dofs_cnt), rho_z(l2dofs_cnt);
   Array<int> dofs(l2dofs_cnt);
   DenseMatrixInverse inv(&Mrho);
   MassIntegrator mi(&ir);
   DensityIntegrator di(qdata);
   di.SetIntRule(&ir);

   for (int e = 0; e < NE; e++)
   {
      const FiniteElement &fe = *L2.GetFE(e);
      ElementTransformation &eltr = *L2.GetElementTransformation(e);
      di.AssembleRHSElementVect(fe, eltr, rhs);
      mi.AssembleElementMatrix(fe, eltr, Mrho);
      inv.Factor();
      inv.Mult(rhs, rho_z);
      L2.GetElementDofs(e, dofs);
      rho.SetSubVector(dofs, rho_z);
   }
}

double ComputeVolumeIntegral(const ParFiniteElementSpace &pfes,
                             const int DIM, const int NE, const int NQ,
                             const int Q1D, const int VDIM, const double norm,
                             const Vector& mass, const Vector& f)
{
   MFEM_VERIFY(pfes.GetNE() > 0, "Empty local mesh should have been handled!");
   MFEM_VERIFY(DIM==1 || DIM==2 || DIM==3, "Unsuported dimension!");
   const bool use_tensors = UsesTensorBasis(pfes);
   const int QX = use_tensors ? Q1D : NQ,
             QY = use_tensors ? Q1D :  1,
             QZ = use_tensors ? Q1D :  1;

   auto f_vals = mfem::Reshape(f.Read(), VDIM, NQ, NE);
   Vector integrand(NE*NQ);
   auto I = Reshape(integrand.Write(), NQ, NE);

   if (DIM == 1)
   {
      MFEM_FORALL(e, NE,
      {
         for (int q = 0; q < NQ; ++q)
         {
            double vmag = 0;
            for (int k = 0; k < VDIM; k++)
            {
               vmag += pow(f_vals(k,q,e), norm);
            }
            I(q,e) = vmag;
         }
      });
   }
   else if (DIM == 2)
   {
      MFEM_FORALL_2D(e, NE, QX, QY, 1,
      {
         MFEM_FOREACH_THREAD(qy,y,QY)
         {
            MFEM_FOREACH_THREAD(qx,x,QX)
            {
               const int q = qx + qy * QX;
               double vmag = 0.0;
               for (int k = 0; k < VDIM; k++)
               {
                  vmag += pow(f_vals(k, q, e), norm);
               }
               I(q, e) = vmag;
            }
         }
      });
   }
   else if (DIM == 3)
   {
      MFEM_FORALL_3D(e, NE, QX, QY, QZ,
      {
         MFEM_FOREACH_THREAD(qz,z,QZ)
         {
            MFEM_FOREACH_THREAD(qy,y,QY)
            {
               MFEM_FOREACH_THREAD(qx,x,QX)
               {
                  const int q = qx + (qy + qz * QY) * QX;
                  double vmag = 0;
                  for (int k = 0; k < VDIM; k++)
                  {
                     vmag += pow(f_vals(k, q, e), norm);
                  }
                  I(q, e) = vmag;
               }
            }
         }
      });
   }
   return integrand * mass;

}
double LagrangianGeoOperator::InternalEnergy(const ParGridFunction &gf) const
{
   double glob_ie = 0.0, internal_energy = 0.0;

   if (L2.GetNE() > 0) // UsesTensorBasis does not handle empty local mesh
   {
      auto L2ordering =
         UsesTensorBasis(L2) ?
         ElementDofOrdering::LEXICOGRAPHIC : ElementDofOrdering::NATIVE;
      // get the restriction and interpolator objects
      auto L2qi = L2.GetQuadratureInterpolator(ir);
      L2qi->SetOutputLayout(QVectorLayout::byVDIM);
      auto L2r = L2.GetElementRestriction(L2ordering);
      const int NQ = ir.GetNPoints();
      const int ND = L2.GetFE(0)->GetDof();
      Vector e_vec(NE*ND), q_val(NE*NQ);
      // Get internal energy at the quadrature points
      L2r->Mult(gf, e_vec);
      L2qi->Values(e_vec, q_val);
      internal_energy =
         ComputeVolumeIntegral(L2, dim, NE, NQ, Q1D,  1, 1.0, qdata.rho0DetJ0w, q_val);
   }

   MPI_Allreduce(&internal_energy, &glob_ie, 1, MPI_DOUBLE, MPI_SUM,
                 L2.GetParMesh()->GetComm());

   return glob_ie;
}

double LagrangianGeoOperator::KineticEnergy(const ParGridFunction &v) const
{
   double glob_ke = 0.0, kinetic_energy = 0.0;

   if (H1.GetNE() > 0) // UsesTensorBasis does not handle empty local mesh
   {
      auto H1ordering =
         UsesTensorBasis(H1) ?
         ElementDofOrdering::LEXICOGRAPHIC : ElementDofOrdering::NATIVE;
      // get the restriction and interpolator objects
      auto h1_interpolator = H1.GetQuadratureInterpolator(ir);
      h1_interpolator->SetOutputLayout(QVectorLayout::byVDIM);
      auto H1r = H1.GetElementRestriction(H1ordering);
      const int NQ = ir.GetNPoints();
      const int ND = H1.GetFE(0)->GetDof();
      Vector e_vec(dim*NE*ND), q_val(dim*NE*NQ);
      // Get internal energy at the quadrature points
      H1r->Mult(v, e_vec);
      h1_interpolator->Values(e_vec, q_val);
      // Get the IE, initial weighted mass
      kinetic_energy =
         ComputeVolumeIntegral(H1, dim, NE, NQ, Q1D, dim, 2.0, qdata.rho0DetJ0w, q_val);
   }

   MPI_Allreduce(&kinetic_energy, &glob_ke, 1, MPI_DOUBLE, MPI_SUM,
                 H1.GetParMesh()->GetComm());

   return 0.5*glob_ke;
}

void LagrangianGeoOperator::PrintTimingData(bool IamRoot, int steps,
                                              const bool fom) const
{
   const MPI_Comm com = H1.GetComm();
   double my_rt[5], T[5];
   my_rt[0] = timer.sw_cgH1.RealTime();
   my_rt[1] = timer.sw_cgL2.RealTime();
   my_rt[2] = timer.sw_force.RealTime();
   my_rt[3] = timer.sw_qdata.RealTime();
   my_rt[4] = my_rt[0] + my_rt[2] + my_rt[3];
   MPI_Reduce(my_rt, T, 5, MPI_DOUBLE, MPI_MAX, 0, com);

   HYPRE_Int mydata[3], alldata[3];
   mydata[0] = timer.L2dof * timer.L2iter;
   mydata[1] = timer.quad_tstep;
   mydata[2] = NE;
   MPI_Reduce(mydata, alldata, 3, HYPRE_MPI_INT, MPI_SUM, 0, com);

   if (IamRoot)
   {
      using namespace std;
      // FOM = (FOM1 * T1 + FOM2 * T2 + FOM3 * T3) / (T1 + T2 + T3)
      const HYPRE_Int H1iter = p_assembly ? (timer.H1iter/dim) : timer.H1iter;
      const double FOM1 = 1e-6 * H1GTVSize * H1iter / T[0];
      const double FOM2 = 1e-6 * steps * (H1GTVSize + L2GTVSize) / T[2];
      const double FOM3 = 1e-6 * alldata[1] * ir.GetNPoints() / T[3];
      const double FOM = (FOM1 * T[0] + FOM2 * T[2] + FOM3 * T[3]) / T[4];
      const double FOM0 = 1e-6 * steps * (H1GTVSize + L2GTVSize) / T[4];
      cout << endl;
      cout << "CG (H1) total time: " << T[0] << endl;
      cout << "CG (H1) rate (megadofs x cg_iterations / second): "
           << FOM1 << endl;
      cout << endl;
      cout << "CG (L2) total time: " << T[1] << endl;
      cout << "CG (L2) rate (megadofs x cg_iterations / second): "
           << 1e-6 * alldata[0] / T[1] << endl;
      cout << endl;
      cout << "Forces total time: " << T[2] << endl;
      cout << "Forces rate (megadofs x timesteps / second): "
           << FOM2 << endl;
      cout << endl;
      cout << "UpdateQuadData total time: " << T[3] << endl;
      cout << "UpdateQuadData rate (megaquads x timesteps / second): "
           << FOM3 << endl;
      cout << endl;
      cout << "Major kernels total time (seconds): " << T[4] << endl;
      cout << "Major kernels total rate (megadofs x time steps / second): "
           << FOM << endl;
      if (!fom) { return; }
      const int QPT = ir.GetNPoints();
      const HYPRE_Int GNZones = alldata[2];
      const long ndofs = 2*H1GTVSize + L2GTVSize + QPT*GNZones;
      cout << endl;
      cout << "| Ranks " << "| Zones   "
           << "| H1 dofs " << "| L2 dofs "
           << "| QP "      << "| N dofs   "
           << "| FOM0   "
           << "| FOM1   " << "| T1   "
           << "| FOM2   " << "| T2   "
           << "| FOM3   " << "| T3   "
           << "| FOM    " << "| TT   "
           << "|" << endl;
      cout << setprecision(3);
      cout << "| " << setw(6) << H1.GetNRanks()
           << "| " << setw(8) << GNZones
           << "| " << setw(8) << H1GTVSize
           << "| " << setw(8) << L2GTVSize
           << "| " << setw(3) << QPT
           << "| " << setw(9) << ndofs
           << "| " << setw(7) << FOM0
           << "| " << setw(7) << FOM1
           << "| " << setw(5) << T[0]
           << "| " << setw(7) << FOM2
           << "| " << setw(5) << T[2]
           << "| " << setw(7) << FOM3
           << "| " << setw(5) << T[3]
           << "| " << setw(7) << FOM
           << "| " << setw(5) << T[4]
           << "| " << endl;
   }
}

// Smooth transition between 0 and 1 for x in [-eps, eps].
MFEM_HOST_DEVICE inline double smooth_step_01(double x, double eps)
{
   const double y = (x + eps) / (2.0 * eps);
   if (y < 0.0) { return 0.0; }
   if (y > 1.0) { return 1.0; }
   return (3.0 - 2.0 * y) * y * y;
}

// void LagrangianGeoOperator::UpdateQuadratureData(const Vector &S, const double dt) const
void LagrangianGeoOperator::UpdateQuadratureData(const Vector &S) const
{
   if (qdata_is_current) { return; }

   qdata_is_current = true;
   forcemat_is_assembled = false;
   gmat_is_assembled = false;

   // if (dim > 1 && p_assembly) { return qupdate->UpdateQuadratureData(S, qdata, dt); }
   if (dim > 1 && p_assembly) { return qupdate->UpdateQuadratureData(S, qdata); }
   // This code is only for the 1D/FA mode
   timer.sw_qdata.Start();
   const int nqp = ir.GetNPoints();
   ParGridFunction x, v, e, sig;
   Vector* sptr = const_cast<Vector*>(&S);
   x.MakeRef(&H1, *sptr, 0);
   v.MakeRef(&H1, *sptr, H1.GetVSize());
   e.MakeRef(&L2, *sptr, 2*H1.GetVSize());
   sig.MakeRef(&L2_stress, *sptr, 2*H1.GetVSize()+L2.GetVSize());

   Vector e_vals;
   Vector sxx, syy, szz;
   Vector sxy, sxz, syz;

   DenseMatrix Jpi(dim), sgrad_v(dim), Jinv(dim), stress(dim), stressJiT(dim);
   DenseMatrix spin(dim), srate(dim);
   DenseMatrix tau0(dim), tau1(dim);
   DenseMatrix old_sig(dim);
   DenseMatrix crot1(dim), crot2(dim);
   DenseMatrix buoy(dim), buoyJiT(dim);

   bool elasticity = true;
   double lame1{1.0};
   double lame2{1.0};
   // double max_vel{0.0};
   double mscale{1.0};
   double grav{10.0};
   int    index = 0;
   mscale = qdata.mscale;
   grav   = -1.0*qdata.gravity;
   
   Vector vel_mag(H1.GetVSize()/dim);
   for (int i = 0; i < H1.GetVSize()/dim; i++)
   {
      if(dim == 2){vel_mag[i] = sqrt(pow(v[i], 2) + pow(v[i+H1.GetVSize()/dim], 2));}
      else{vel_mag[i] = sqrt(pow(v[i], 2) + pow(v[i+H1.GetVSize()/dim], 2) + pow(v[i+2*H1.GetVSize()/dim], 2));}
   }

   double cut_off_vel   = 1.0/86400.0/365.25/1000.0/100.0; // 0.01 mm/yr
   double local_max_vel = std::max(cut_off_vel, vel_mag.Max());
   double global_max_vel;
   MPI_Reduce(&local_max_vel, &global_max_vel, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   // Broadcast the global_max from process 0 to all other processes
   MPI_Bcast(&global_max_vel, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


   // Batched computations are needed, because geodynamic codes usually
   // involve expensive computations of material properties. Although this
   // miniapp uses simple EOS equations, we still want to represent the batched
   // cycle structure.
   int nzones_batch = 3;
   const int nbatches =  NE / nzones_batch + 1; // +1 for the remainder.
   int nqp_batch = nqp * nzones_batch;
   double *gamma_b = new double[nqp_batch],
   *rho_b = new double[nqp_batch],
   *e_b   = new double[nqp_batch],
   *p_b   = new double[nqp_batch],
   *cs_b  = new double[nqp_batch];
   double *lambda_b = new double[nqp_batch], 
   *mu_b = new double[nqp_batch],
   *pmod_b = new double[nqp_batch];
   // Jacobians of reference->physical transformations for all quadrature points
   // in the batch.
   DenseTensor *Jpr_b = new DenseTensor[nzones_batch];

   for (int b = 0; b < nbatches; b++)
   {
      int z_id = b * nzones_batch; // Global index over zones.
      // The last batch might not be full.
      if (z_id == NE) { break; }
      else if (z_id + nzones_batch > NE)
      {
         nzones_batch = NE - z_id;
         nqp_batch    = nqp * nzones_batch;
      }

      double min_detJ = std::numeric_limits<double>::infinity();
      for (int z = 0; z < nzones_batch; z++)
      {
         ElementTransformation *T = H1.GetElementTransformation(z_id);
         Jpr_b[z].SetSize(dim, dim, nqp);
         e.GetValues(z_id, ir, e_vals);
         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            T->SetIntPoint(&ip);
            Jpr_b[z](q) = T->Jacobian();
            const double detJ = Jpr_b[z](q).Det();
            min_detJ = fmin(min_detJ, detJ);
            const int idx = z * nqp + q;
            // Assuming piecewise constant gamma that moves with the mesh.
            gamma_b[idx] = gamma_gf(z_id);
            rho_b[idx] = qdata.rho0DetJ0w(z_id*nqp + q) / detJ / ip.weight;
            e_b[idx] = fmax(0.0, e_vals(q));
            lambda_b[idx] = lambda_gf(z_id);
            mu_b[idx] = mu_gf(z_id);
            pmod_b[idx] = lambda_b[idx] + 2*mu_b[idx];
         }
         ++z_id;
      }
      // Batched computation of material properties.
      ComputeMaterialProperties(nqp_batch, gamma_b, rho_b, e_b, p_b, cs_b, pmod_b, mscale);

      z_id -= nzones_batch;
      for (int z = 0; z < nzones_batch; z++)
      {
         // qvol =0.0;
         ElementTransformation *T = H1.GetElementTransformation(z_id);
         // ElementTransformation &eltr = *L2_stress.GetElementTransformation(z_id);

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            T->SetIntPoint(&ip);
            // Note that the Jacobian was already computed above. We've chosen
            // not to store the Jacobians for all batched quadrature points.
            const DenseMatrix &Jpr = Jpr_b[z](q);
            CalcInverse(Jpr, Jinv);
            const double detJ = Jpr.Det(), rho = 1.0*rho_b[z*nqp + q],
                         p = p_b[z*nqp + q], sound_speed = cs_b[z*nqp + q];
            
            // const double detJ = Jpr.Det(), rho = rho_b[z*nqp + q],
            //              p = p_b[z*nqp + q], sound_speed = cs_b[z*nqp + q];
            // index = gamma_b[z*nqp + q];
            lame1 = lambda_b[z*nqp + q];
            lame2 = mu_b[z*nqp + q];
            stress = 0.0; tau0 = 0.0; tau1 = 0.0; old_sig = 0.0;
            buoy   = 0.0; buoyJiT = 0.0; 
            buoy(dim-1, dim-1) = 1.0*grav;
            // plastic_str = 0.0; plastic_sig = 0.0;

            for (int d = 0; d < dim; d++) { stress(d, d) = 0; }
            // for (int d = 0; d < dim; d++) { stress(d, d) = -p; }
            for (int d = 0; d < dim; d++) { tau1(d, d) = 1.0;} // Identity matrix
 
            double visc_coeff = 0.0;

            
            if (elasticity)
            {
               // Compression-based length scale at the point. The first
               // eigenvector of the symmetric velocity gradient gives the
               // direction of maximal compression. This is used to define the
               // relative change of the initial length scale.
               v.GetVectorGradient(*T, sgrad_v);

               // Get deviatoric stress components at quadrature points
               if(dim == 2)
               {
                  sig.GetValues(z_id, ir, sxx, 1);
                  sig.GetValues(z_id, ir, syy, 2);
                  sig.GetValues(z_id, ir, sxy, 3);
                  old_sig(0,0) = sxx(q) ; old_sig(0,1) = sxy(q); old_sig(1,0) = sxy(q); old_sig(1,1) = syy(q);
               }
               else if(dim == 3)
               {

                  sig.GetValues(z_id, ir, sxx, 1);
                  sig.GetValues(z_id, ir, syy, 2);
                  sig.GetValues(z_id, ir, szz, 3);
                  sig.GetValues(z_id, ir, sxy, 4);
                  sig.GetValues(z_id, ir, sxz, 5);
                  sig.GetValues(z_id, ir, syz, 6);
                                    
                  old_sig(0,0) = sxx(q) ; old_sig(0,1) = sxy(q); old_sig(0,2) = sxz(q); 
                  old_sig(1,0) = sxy(q) ; old_sig(1,1) = syy(q); old_sig(1,2) = syz(q);
                  old_sig(2,0) = sxz(q) ; old_sig(2,1) = syz(q); old_sig(2,2) = szz(q);
               }

               double vorticity_coeff = 1.0;
               if (use_vorticity)
               {
                  const double grad_norm = sgrad_v.FNorm();
                  const double div_v = fabs(sgrad_v.Trace());
                  vorticity_coeff = (grad_norm > 0.0) ? div_v / grad_norm : 1.0;
               }
               // copy sgrad_v to srate and then symmetrize srate.
               srate = sgrad_v;
               srate.Symmetrize();
               // copy sgrad_v to spin and subtract srate from it.
               spin = sgrad_v;
               spin.Add(-1.0, srate);

               double eig_val_data[3], eig_vec_data[9];
               double eig_sig_var[3], eig_sig_vec[9];
                              
               if (dim==1)
               {
                  eig_val_data[0] = sgrad_v(0, 0);
                  eig_vec_data[0] = 1.;
               }
               else { sgrad_v.CalcEigenvalues(eig_val_data, eig_vec_data); }
               Vector compr_dir(eig_vec_data, dim);
               // Computes the initial->physical transformation Jacobian.
               mfem::Mult(Jpr, qdata.Jac0inv(z_id*nqp + q), Jpi);
               Vector ph_dir(dim); Jpi.Mult(compr_dir, ph_dir);
               // Change of the initial mesh size in the compression direction.
               const double h = qdata.h0 * ph_dir.Norml2() /
                                compr_dir.Norml2();
               // Measure of maximal compression.
               const double mu = eig_val_data[0];
               visc_coeff = 2.0 * rho * h * h * fabs(mu);
               // The following represents a "smooth" version of the statement
               // "if (mu < 0) visc_coeff += 0.5 rho h sound_speed".  Note that
               // eps must be scaled appropriately if a different unit system is
               // being used.
               const double eps = 1e-12;
               visc_coeff += 0.5 * rho * h * sound_speed * vorticity_coeff *
                             (1.0 - smooth_step_01(mu - 2.0 * eps, eps));

               if(use_viscosity){stress.Add(1.0*visc_coeff, sgrad_v);} // applying visc_coeff
               stress.Add(1.0, old_sig); // previous total stress

               // #if
               tau0.Set(2*lame2,srate); // 2 * G * eij
               tau1.Set(lame1*srate.Trace(), tau1); // lambda*tr(strain rate)*Identity
               tau0.Add(1.0, tau1); // stress rate for s_ij grid function.
               // tau0.Add(-1.0, tau1); // deviatoric stress rate for s_ij grid function.
               
               mfem::Mult(old_sig, spin, crot1);
               mfem::Mult(spin, old_sig, crot2);

               // Juamman stress increment
               tau0.Add(1.0,  crot1); 
               tau0.Add(-1.0, crot2); 
               // tau0.Add(-1.0, plastic_sig); // Applying plastic stress correction.
               // #ifend
            }
            // Time step estimate at the point. Here the more relevant length
            // scale is related to the actual mesh deformation; we use the min
            // singular value of the ref->physical Jacobian. In addition, the
            // time step estimate should be aware of the presence of shocks.
            
            const double h_min =
               Jpr.CalcSingularvalue(dim-1) / (double) H1.GetOrder(0);

            // // using initial size
            // const double h_min =
            //     1/qdata.Jac0inv(z_id*nqp + q).CalcSingularvalue(dim-1) / (double) H1.GetOrder(0);

            // const double inv_dt = sound_speed / h_min;
            // const double inv_dt = sound_speed / h_min +
                                 //  2.5 * visc_coeff / rho / h_min / h_min;
            // const double inv_dt = sound_speed / h_min +
            //                       1.0 * visc_coeff / rho / h_min / h_min;


            const double vel_max = fmax(qdata.vbc_max_val, global_max_vel);
            const double inv_dt  = (vel_max * qdata.mscale) / h_min;
            
            if (min_detJ < 0.0)
            {
               // This will force repetition of the step with smaller dt.
               // qdata.dt_est = 0.0;
               MFEM_ABORT("Negative min_detJ!");
            }
            else
            {
               if (inv_dt>0.0)
               {
                  qdata.dt_est = fmin(qdata.dt_est, cfl*(1.0/inv_dt));
                  qdata.h_est  = fmin(qdata.h_est, h_min);
               }
            }
            // Quadrature data for partial assembly of the force operator.
            MultABt(stress, Jinv, stressJiT);
            stressJiT *= ir.IntPoint(q).weight * detJ;
            tau0 *= ir.IntPoint(q).weight * detJ * (rho / 1.0);
            buoy *= ir.IntPoint(q).weight * detJ * (rho / 1.0);
            
            for (int vd = 0 ; vd < dim; vd++)
            {
               for (int gd = 0; gd < dim; gd++)
               {
                  const int offset = z_id*nqp + q + nqp*NE*(gd + vd*dim);
                  
                  qdata.stressJinvT(vd)(z_id*nqp + q, gd) = stressJiT(vd, gd);
                  qdata.tauJinvT(vd)(z_id*nqp + q, gd)    = tau0(vd, gd);
                  qdata.buoyJinvT(vd)(z_id*nqp + q, gd)   = buoy(vd, gd);
                  // qdata.buoyJinvT(vd)(z_id*nqp + q, gd)   = buoyJiT(vd, gd);
               }
            }
         }
         ++z_id;
      }
   }
   delete [] gamma_b;
   delete [] rho_b;
   delete [] e_b;
   delete [] p_b;
   delete [] cs_b;
   delete [] Jpr_b;
   delete [] lambda_b;
   delete [] mu_b;
   delete [] pmod_b;
   timer.sw_qdata.Stop();
   timer.quad_tstep += NE;
}

/// Trace of a square matrix
template<int H, int W, typename T>
MFEM_HOST_DEVICE inline
double Trace(const T * __restrict__ data)
{
   double t = 0.0;
   for (int i = 0; i < W; i++) { t += data[i+i*H]; }
   return t;
}

template<int H, int W, typename T>
MFEM_HOST_DEVICE static inline
void SFNorm(double &scale_factor, double &scaled_fnorm2,
            const T * __restrict__ data)
{
   int i;
   constexpr int hw = H * W;
   T max_norm = 0.0, entry, fnorm2;

   for (i = 0; i < hw; i++)
   {
      entry = fabs(data[i]);
      if (entry > max_norm)
      {
         max_norm = entry;
      }
   }

   if (max_norm == 0.0)
   {
      scale_factor = scaled_fnorm2 = 0.0;
      return;
   }

   fnorm2 = 0.0;
   for (i = 0; i < hw; i++)
   {
      entry = data[i] / max_norm;
      fnorm2 += entry * entry;
   }

   scale_factor = max_norm;
   scaled_fnorm2 = fnorm2;
}

/// Compute the Frobenius norm of the matrix
template<int H, int W, typename T>
MFEM_HOST_DEVICE inline
double FNorm(const T * __restrict__ data)
{
   double s, n2;
   SFNorm<H,W>(s, n2, data);
   return s*sqrt(n2);
}

template<int DIM> MFEM_HOST_DEVICE static inline
void QUpdateBody(const int NE, const int e,
                 const int NQ, const int q,
                 const bool use_viscosity,
                 const bool use_vorticity,
                 const double h0,
                 const double h1order,
                 const double cfl,
                 const double infinity,
                 double max_vel_q,
                 double mscale_q, 
                 double gravity_q,
                 double* __restrict__ Jinv,
                 double* __restrict__ stress,
                 double* __restrict__ tau0,
                 double* __restrict__ tau1,
                 double* __restrict__ tau2,
                 double* __restrict__ tau0_Jit,
                 double* __restrict__ tau1_Jit,
                 double* __restrict__ tau2_Jit,
                 double* __restrict__ sgrad_v,
                 double* __restrict__ spin,
                 double* __restrict__ eig_val_data,
                 double* __restrict__ eig_vec_data,
                 double* __restrict__ sig_val_data,
                 double* __restrict__ sig_vec_data,
                 double* __restrict__ compr_dir,
                 double* __restrict__ Jpi,
                 double* __restrict__ ph_dir,
                 double* __restrict__ stressJiT,
                 const double* __restrict__ d_gamma,
                 const double* __restrict__ d_lambda,
                 const double* __restrict__ d_mu,
                 const double* __restrict__ d_weights,
                 const double* __restrict__ d_Jacobians,
                 const double* __restrict__ d_rho0DetJ0w,
                 const double* __restrict__ d_e_quads,
                 const double* __restrict__ d_sig_quads,
                 const double* __restrict__ d_grad_v_ext,
                 const double* __restrict__ d_Jac0inv,
                 double *d_dt_est,
                 double *d_h_est,
                 double *d_stressJinvT,
                 double *d_tauJinvT,
                 double *d_buoyJinvT) // -6-
{
   constexpr int DIM2 = DIM*DIM;
   double min_detJ = infinity;

   const int eq = e * NQ + q;
   const double gamma = d_gamma[e];
   double lame1 = d_lambda[e];
   double lame2 = d_mu[e];
   const double weight =  d_weights[q];
   const double inv_weight = 1. / weight;
   const double *J = d_Jacobians + DIM2*(NQ*e + q); // q_dx -> Jacobians -> d_Jacobians
   const double detJ = kernels::Det<DIM>(J);
   min_detJ = fmin(min_detJ, detJ);
   kernels::CalcInverse<DIM>(J, Jinv);
   const double R = inv_weight * d_rho0DetJ0w[eq] / detJ;
   const double E = fmax(0.0, d_e_quads[eq]);
   const double P = (gamma - 1.0) * R * E;
   // const double S = sqrt(gamma * (gamma - 1.0) * E);
   double old_sig[DIM2], crot1[DIM2], crot2[DIM2];
   double buoy[DIM2];
   double sxx{0.0}, syy{0.0}, szz{0.0};
   double sxy{0.0}, sxz{0.0}, syz{0.0};
   double grav{-10.0}, mscale{1e16};
   const double S = sqrt((lame1+2*lame2)/(R * mscale));
   bool elasticity = true;
   bool plasticity = true;

   for (int k = 0; k < DIM2; k++) 
      { 
         stress[k] = 0.0; tau0[k] = 0.0; tau1[k] = 0.0; tau2[k] = 0.0; tau1_Jit[k] = 0.0; tau2_Jit[k] = 0.0; 
         crot1[k] = 0.0; crot2[k] = 0.0; buoy[k] = 0.0;
      }
   for (int d = 0; d < DIM; d++) { tau1[d*DIM+d] = 1.0;} // I matrix
   
   if (DIM == 2) {buoy[3] = -1*gravity_q*1.0;} 
   else {buoy[8] = -1*gravity_q*1.0;}
   // for (int d = 0; d < DIM; d++) { stress[d*DIM+d] = -P; }
   double visc_coeff = 0.0;
   // double d_lambda{1e-1};
   // double d_mu{1e-1};
   // DenseMatrix old_sig(DIM);
   // DenseMatrix crot1(DIM), crot2(DIM);
   
   if(DIM == 2)
   {  
      sxx = d_sig_quads[NQ*e + q + 0*NQ*NE];
      syy = d_sig_quads[NQ*e + q + 1*NQ*NE];
      sxy = d_sig_quads[NQ*e + q + 2*NQ*NE];
      old_sig[0] = sxx; old_sig[1] = sxy; old_sig[2] = sxy; old_sig[3] = syy;
   }
   else if(DIM == 3)
   {
      sxx = d_sig_quads[NQ*e + q + 0*NQ*NE];
      syy = d_sig_quads[NQ*e + q + 1*NQ*NE];
      szz = d_sig_quads[NQ*e + q + 2*NQ*NE];
      sxy = d_sig_quads[NQ*e + q + 3*NQ*NE];
      sxz = d_sig_quads[NQ*e + q + 4*NQ*NE];
      syz = d_sig_quads[NQ*e + q + 5*NQ*NE];

      old_sig[0] = sxx; old_sig[1] = sxy; old_sig[2] = sxz; 
      old_sig[3] = sxy; old_sig[4] = syy; old_sig[5] = syz;
      old_sig[6] = sxz; old_sig[7] = syz; old_sig[8] = szz;
   }


   if (plasticity)
   {
      kernels::CalcEigenvalues<DIM>(sgrad_v, sig_val_data, sig_vec_data);
   }

   if (elasticity)
   {
      // Compression-based length scale at the point. The first
      // eigenvector of the symmetric velocity gradient gives the
      // direction of maximal compression. This is used to define the
      // relative change of the initial length scale.
      const double *dV = d_grad_v_ext + DIM2*(NQ*e + q);
      kernels::Mult(DIM, DIM, DIM, dV, Jinv, sgrad_v); // sgrad_v = dV*Jinv (J^-1)
      for (int k = 0; k < DIM2; k++) { spin[k] = sgrad_v[k]; }

      double vorticity_coeff = 1.0;
      if (use_vorticity)
      {
         const double grad_norm = FNorm<DIM,DIM>(sgrad_v);
         const double div_v = fabs(Trace<DIM,DIM>(sgrad_v));
         vorticity_coeff = (grad_norm > 0.0) ? div_v / grad_norm : 1.0;
      }

      // sgrad_v is now (symmetric) strain rate
      kernels::Symmetrize(DIM, sgrad_v);
      // 2nd spin = 1st spin (= the original sgrad_v) - 1.0 * sgrad_v (symmetrized)
      kernels::Add(DIM, DIM, -1.0, spin, sgrad_v, spin);

      if (DIM == 1)
      {
         eig_val_data[0] = sgrad_v[0];
         eig_vec_data[0] = 1.;
      }
      else
      {
         kernels::CalcEigenvalues<DIM>(sgrad_v, eig_val_data, eig_vec_data);
      }
      for (int k=0; k<DIM; k++) { compr_dir[k] = eig_vec_data[k]; }
      // Computes the initial->physical transformation Jacobian.
      kernels::Mult(DIM, DIM, DIM, J, d_Jac0inv + eq*DIM*DIM, Jpi);
      kernels::Mult(DIM, DIM, Jpi, compr_dir, ph_dir);
      // Change of the initial mesh size in the compression direction.
      const double ph_dir_nl2 = kernels::Norml2(DIM, ph_dir);
      const double compr_dir_nl2 = kernels::Norml2(DIM, compr_dir);
      const double H = h0 * ph_dir_nl2 / compr_dir_nl2;
      // Measure of maximal compression.
      const double mu = eig_val_data[0];
      visc_coeff = 2.0 * R * H * H * fabs(mu); // R = inv_weight * d_rho0DetJ0w[eq] / detJ;
      // The following represents a "smooth" version of the statement
      // "if (mu < 0) visc_coeff += 0.5 rho h sound_speed".  Note that
      // eps must be scaled appropriately if a different unit system is
      // being used.
      const double eps = 1e-12;
      visc_coeff += 0.5 * R * H  * S * vorticity_coeff *
                    (1.0 - smooth_step_01(mu-2.0*eps, eps));
      // stress = DIM x DIM matrix
      // sgrad_v = DIM x DIM matrix
      // stress(2nd) = stress(first) + visc_coeff * sgrad_v
      if(use_viscosity){kernels::Add(DIM, DIM, 1.0*visc_coeff, stress, sgrad_v, stress);}

      // for (int k = 0; k < DIM2; k++) {stress[k] = old_sig[k];}
      kernels::Add(DIM, DIM, 1.0, stress, old_sig, stress); // previous total stress

      // stress deviator tensor
      const double div_vel = Trace<DIM,DIM>(sgrad_v); 
      kernels::Set(DIM, DIM, 2.0*lame2, sgrad_v, tau0); // 2 * G * eij
      
      kernels::Set(DIM, DIM, lame1*div_vel, tau1, tau1); // lambda*tr(strain rate)*Identity
      kernels::Add(DIM, DIM, 1.0, tau0, tau1, tau0); // stress rate tensor
      
      kernels::Mult(DIM, DIM, DIM, old_sig, spin, crot1); // corotational term; sigma*w
      kernels::Mult(DIM, DIM, DIM, spin, old_sig, crot2); // corotational term; -w*sigma

      kernels::Add(DIM, DIM,  1.0, tau0, crot1, tau0);
      kernels::Add(DIM, DIM, -1.0, tau0, crot2, tau0);
   }

   // Time step estimate at the point. Here the more relevant length
   // scale is related to the actual mesh deformation; we use the min
   // singular value of the ref->physical Jacobian. In addition, the
   // time step estimate should be aware of the presence of shocks.
   const double sv = kernels::CalcSingularvalue<DIM>(J, DIM - 1);
   const double h_min = sv / h1order;
   const double ih_min = 1. / h_min;
   const double irho_ih_min_sq = ih_min * ih_min / R ;
   // const double idt = S * ih_min + 2.5 * visc_coeff / R / h_min / h_min;
   // const double idt = S * ih_min + 1.0 * visc_coeff * irho_ih_min_sq;
   const double idt = (max_vel_q * mscale_q) / h_min;
   if (min_detJ < 0.0)
   {
      // This will force repetition of the step with smaller dt.
      d_dt_est[eq] = 0.0;
   }
   else
   {
      if (idt > 0.0)
      {
         const double cfl_inv_dt = cfl / idt;
         d_dt_est[eq] = fmin(d_dt_est[eq], cfl_inv_dt);
         d_h_est[eq] = fmin(d_h_est[eq], h_min);
      }
   }
   // Quadrature data for partial assembly of the force operator.
   kernels::MultABt(DIM, DIM, DIM, stress, Jinv, stressJiT);
   for (int k = 0; k < DIM2; k++) { stressJiT[k] *= weight * detJ; }
   for (int k = 0; k < DIM2; k++) { tau0[k] *= R * weight * detJ; }
   for (int k = 0; k < DIM2; k++) { buoy[k] *= R * weight * detJ; }
   
   for (int vd = 0 ; vd < DIM; vd++)
   {
      for (int gd = 0; gd < DIM; gd++)
      {
         const int offset = eq + NQ*NE*(gd + vd*DIM);
         d_stressJinvT[offset] = stressJiT[vd + gd*DIM];
         d_tauJinvT[offset] = tau0[vd + gd*DIM];
         d_buoyJinvT[offset] = buoy[vd + gd*DIM];
      }
   }
}

static void Rho0DetJ0Vol(const int dim, const int NE,
                         const IntegrationRule &ir,
                         ParMesh *pmesh,
                         ParFiniteElementSpace &L2,
                         const ParGridFunction &rho0,
                         QuadratureData &qdata,
                         double &volume)
{
   const int NQ = ir.GetNPoints();
   const int Q1D = IntRules.Get(Geometry::SEGMENT,ir.GetOrder()).GetNPoints();
   const int flags = GeometricFactors::JACOBIANS|GeometricFactors::DETERMINANTS;
   const GeometricFactors *geom = pmesh->GetGeometricFactors(ir, flags);
   Vector rho0Q(NQ*NE);
   rho0Q.UseDevice(true);
   Vector j, detj;
   const QuadratureInterpolator *qi = L2.GetQuadratureInterpolator(ir);
   qi->Mult(rho0, QuadratureInterpolator::VALUES, rho0Q, j, detj);
   const auto W = ir.GetWeights().Read();
   const auto R = Reshape(rho0Q.Read(), NQ, NE);
   const auto J = Reshape(geom->J.Read(), NQ, dim, dim, NE);
   const auto detJ = Reshape(geom->detJ.Read(), NQ, NE);
   auto V = Reshape(qdata.rho0DetJ0w.Write(), NQ, NE);
   Memory<double> &Jinv_m = qdata.Jac0inv.GetMemory();
   const MemoryClass mc = Device::GetMemoryClass();
   const int Ji_total_size = qdata.Jac0inv.TotalSize();
   auto invJ = Reshape(Jinv_m.Write(mc, Ji_total_size), dim, dim, NQ, NE);
   Vector vol(NE*NQ), one(NE*NQ);
   auto A = Reshape(vol.Write(), NQ, NE);
   auto O = Reshape(one.Write(), NQ, NE);
   MFEM_ASSERT(dim==2 || dim==3, "");
   if (dim==2)
   {
      MFEM_FORALL_2D(e, NE, Q1D, Q1D, 1,
      {
         MFEM_FOREACH_THREAD(qy,y,Q1D)
         {
            MFEM_FOREACH_THREAD(qx,x,Q1D)
            {
               const int q = qx + qy * Q1D;
               const double J11 = J(q,0,0,e);
               const double J12 = J(q,1,0,e);
               const double J21 = J(q,0,1,e);
               const double J22 = J(q,1,1,e);
               const double det = detJ(q,e);
               V(q,e) =  W[q] * R(q,e) * det;
               const double r_idetJ = 1.0 / det;
               invJ(0,0,q,e) =  J22 * r_idetJ;
               invJ(1,0,q,e) = -J12 * r_idetJ;
               invJ(0,1,q,e) = -J21 * r_idetJ;
               invJ(1,1,q,e) =  J11 * r_idetJ;
               A(q,e) = W[q] * det;
               O(q,e) = 1.0;
            }
         }
      });
   }
   else
   {
      MFEM_FORALL_3D(e, NE, Q1D, Q1D, Q1D,
      {
         MFEM_FOREACH_THREAD(qz,z,Q1D)
         {
            MFEM_FOREACH_THREAD(qy,y,Q1D)
            {
               MFEM_FOREACH_THREAD(qx,x,Q1D)
               {
                  const int q = qx + (qy + qz * Q1D) * Q1D;
                  const double J11 = J(q,0,0,e), J12 = J(q,0,1,e), J13 = J(q,0,2,e);
                  const double J21 = J(q,1,0,e), J22 = J(q,1,1,e), J23 = J(q,1,2,e);
                  const double J31 = J(q,2,0,e), J32 = J(q,2,1,e), J33 = J(q,2,2,e);
                  const double det = detJ(q,e);
                  V(q,e) = W[q] * R(q,e) * det;
                  const double r_idetJ = 1.0 / det;
                  invJ(0,0,q,e) = r_idetJ * ((J22 * J33)-(J23 * J32));
                  invJ(1,0,q,e) = r_idetJ * ((J32 * J13)-(J33 * J12));
                  invJ(2,0,q,e) = r_idetJ * ((J12 * J23)-(J13 * J22));
                  invJ(0,1,q,e) = r_idetJ * ((J23 * J31)-(J21 * J33));
                  invJ(1,1,q,e) = r_idetJ * ((J33 * J11)-(J31 * J13));
                  invJ(2,1,q,e) = r_idetJ * ((J13 * J21)-(J11 * J23));
                  invJ(0,2,q,e) = r_idetJ * ((J21 * J32)-(J22 * J31));
                  invJ(1,2,q,e) = r_idetJ * ((J31 * J12)-(J32 * J11));
                  invJ(2,2,q,e) = r_idetJ * ((J11 * J22)-(J12 * J21));
                  A(q,e) = W[q] * det;
                  O(q,e) = 1.0;
               }
            }
         }
      });
   }
   qdata.rho0DetJ0w.HostRead();
   volume = vol * one;
}

// dt
template<int DIM, int Q1D> static inline
void QKernel(const int NE, const int NQ,
             const bool use_viscosity,
             const bool use_vorticity,
             const double h0,
             const double h1order,
             const double cfl,
             const double infinity,
             double max_vel_q,
             double mscale_q, 
             double gravity_q,
             const ParGridFunction &gamma_gf,
             const ParGridFunction &lambda_gf,
             const ParGridFunction &mu_gf,
             const Array<double> &weights,
             const Vector &Jacobians,
             const Vector &rho0DetJ0w,
             const Vector &e_quads,
             const Vector &sig_quads,
             const Vector &grad_v_ext,
             const DenseTensor &Jac0inv,
             Vector &dt_est, Vector &h_est,
             DenseTensor &stressJinvT, DenseTensor &tauJinvT, DenseTensor &buoyJinvT) //-4-
{
   constexpr int DIM2 = DIM*DIM;
   const auto d_gamma = gamma_gf.Read();
   const auto d_lambda = lambda_gf.Read();
   const auto d_mu = mu_gf.Read();
   // auto d_sigma = sigma_gf.ReadWrite();
   const auto d_weights = weights.Read();
   const auto d_Jacobians = Jacobians.Read();
   const auto d_rho0DetJ0w = rho0DetJ0w.Read();
   const auto d_e_quads = e_quads.Read();
   const auto d_sig_quads = sig_quads.Read();
   const auto d_grad_v_ext = grad_v_ext.Read();
   const auto d_Jac0inv = Read(Jac0inv.GetMemory(), Jac0inv.TotalSize());
   auto d_dt_est = dt_est.ReadWrite();
   auto d_h_est = h_est.ReadWrite();
   auto d_stressJinvT = Write(stressJinvT.GetMemory(), stressJinvT.TotalSize());
   auto d_tauJinvT = Write(tauJinvT.GetMemory(), tauJinvT.TotalSize());
   auto d_buoyJinvT = Write(buoyJinvT.GetMemory(), buoyJinvT.TotalSize());
   // auto d_tauJinvT = tauJinvT.ReadWrite();
   if (DIM == 2)
   {
      MFEM_FORALL_2D(e, NE, Q1D, Q1D, 1,
      {
         double Jinv[DIM2];
         double stress[DIM2];
         double tau0[DIM2];
         double tau1[DIM2];
         double tau2[DIM2];
         double tau0_Jit[DIM2];
         double tau1_Jit[DIM2];
         double tau2_Jit[DIM2];
         double sgrad_v[DIM2];
         double spin[DIM2];
         double eig_val_data[3];
         double eig_vec_data[9];
         double sig_val_data[2];
         double sig_vec_data[4];
         double compr_dir[DIM];
         double Jpi[DIM2];
         double ph_dir[DIM];
         double stressJiT[DIM2];
         MFEM_FOREACH_THREAD(qx,x,Q1D)
         {
            MFEM_FOREACH_THREAD(qy,y,Q1D)
            {
               QUpdateBody<DIM>(NE, e, NQ, qx + qy * Q1D,
                                use_viscosity, use_vorticity, h0, h1order, cfl, infinity, max_vel_q, mscale_q, gravity_q,
                                Jinv, stress, tau0, tau1, tau2, tau0_Jit, tau1_Jit, tau2_Jit,
                                sgrad_v, spin, eig_val_data, eig_vec_data, sig_val_data, sig_vec_data,
                                compr_dir, Jpi, ph_dir, stressJiT,
                                d_gamma, d_lambda, d_mu, d_weights, d_Jacobians, d_rho0DetJ0w,
                                d_e_quads, d_sig_quads, d_grad_v_ext, d_Jac0inv,
                                d_dt_est,  d_h_est, d_stressJinvT, d_tauJinvT, d_buoyJinvT); //-5a-
            }
         }
         MFEM_SYNC_THREAD;
      });
   }
   if (DIM == 3)
   {
      MFEM_FORALL_3D(e, NE, Q1D, Q1D, Q1D,
      {
         double Jinv[DIM2];
         double stress[DIM2];
         double tau0[DIM2];
         double tau1[DIM2];
         double tau2[DIM2];
         double tau0_Jit[DIM2];
         double tau1_Jit[DIM2];
         double tau2_Jit[DIM2];
         double sgrad_v[DIM2];
         double spin[DIM2];
         double eig_val_data[3];
         double eig_vec_data[9];
         double sig_val_data[3];
         double sig_vec_data[9];
         double compr_dir[DIM];
         double Jpi[DIM2];
         double ph_dir[DIM];
         double stressJiT[DIM2];
         MFEM_FOREACH_THREAD(qx,x,Q1D)
         {
            MFEM_FOREACH_THREAD(qy,y,Q1D)
            {
               MFEM_FOREACH_THREAD(qz,z,Q1D)
               {
                  QUpdateBody<DIM>(NE, e, NQ, qx + Q1D * (qy + qz * Q1D),
                                   use_viscosity, use_vorticity, h0, h1order, cfl, infinity, max_vel_q, mscale_q, gravity_q,
                                   Jinv, stress, tau0, tau1, tau2, tau0_Jit, tau1_Jit, tau2_Jit,
                                   sgrad_v, spin, eig_val_data, eig_vec_data, sig_val_data, sig_vec_data,
                                   compr_dir, Jpi, ph_dir, stressJiT,
                                   d_gamma, d_lambda, d_mu, d_weights, d_Jacobians, d_rho0DetJ0w,
                                   d_e_quads, d_sig_quads, d_grad_v_ext, d_Jac0inv,
                                   d_dt_est, d_h_est, d_stressJinvT, d_tauJinvT, d_buoyJinvT); //-5b-
               }
            }
         }
         MFEM_SYNC_THREAD;
      });
   }
}

void QUpdate::UpdateQuadratureData(const Vector &S, QuadratureData &qdata)
{
   timer->sw_qdata.Start();
   Vector* S_p = const_cast<Vector*>(&S);
   const int H1_size = H1.GetVSize();
   const int L2_size = L2.GetVSize();
   const double h1order = (double) H1.GetOrder(0);
   const double infinity = std::numeric_limits<double>::infinity();
   ParGridFunction x, v, e, sig;
   x.MakeRef(&H1,*S_p, 0);
   H1R->Mult(x, e_vec);
   q1->SetOutputLayout(QVectorLayout::byVDIM);
   q1->Derivatives(e_vec, q_dx); // q_dx -> Jacobians -> d_Jacobians
   v.MakeRef(&H1,*S_p, H1_size);
   H1R->Mult(v, e_vec); // v * unit vector
   q1->Derivatives(e_vec, q_dv); // q_dv -> grad_v_ext -> d_grad_v_ext
   e.MakeRef(&L2, *S_p, 2*H1_size);
   q2->SetOutputLayout(QVectorLayout::byVDIM);
   q2->Values(e, q_e); // q_e -> e_quads -> d_e_quads
   sig.MakeRef(&L2_stress, *S_p, 2*H1_size+L2_size);
   // q3->SetOutputLayout(QVectorLayout::byVDIM);
   q3->SetOutputLayout(QVectorLayout::byNODES);
   q3->Values(sig, q_sig); // q_sig -> e_quads -> d_e_quads
   q_dt_est = qdata.dt_est;
   q_h_est = qdata.h_est;

   Vector vel_mag(H1.GetVSize()/dim);
   for (int i = 0; i < H1.GetVSize()/dim; i++)
   {
      if(dim == 2){vel_mag[i] = sqrt(pow(v[i], 2) + pow(v[i+H1.GetVSize()/dim], 2));}
      else{vel_mag[i] = sqrt(pow(v[i], 2) + pow(v[i+H1.GetVSize()/dim], 2) + pow(v[i+2*H1.GetVSize()/dim], 2));}
   }

   double cut_off_vel   = 1.0/86400.0/365.25/1000.0/100.0; // 0.01 mm/yr
   double local_max_vel = std::max(cut_off_vel, vel_mag.Max());
   double global_max_vel;
   MPI_Reduce(&local_max_vel, &global_max_vel, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   // Broadcast the global_max from process 0 to all other processes
   MPI_Bcast(&global_max_vel, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
   double max_vel_q = fmax(qdata.vbc_max_val, global_max_vel);

   const int id = (dim << 4) | Q1D;
   typedef void (*fQKernel)(const int NE, const int NQ,
                            const bool use_viscosity,
                            const bool use_vorticity,
                            const double h0, const double h1order,
                            const double cfl, const double infinity, double max_vel_q, double mscale_q, double gravity_q,
                            const ParGridFunction &gamma_gf,
                            const ParGridFunction &lambda_gf,
                            const ParGridFunction &mu_gf,
                            const Array<double> &weights,
                            const Vector &Jacobians, const Vector &rho0DetJ0w,
                            const Vector &e_quads, const Vector &sig_quads, const Vector &grad_v_ext,
                            const DenseTensor &Jac0inv,
                            Vector &dt_est, Vector &h_est, 
                            DenseTensor &stressJinvT, DenseTensor &tauJinvT, DenseTensor &buoyJinvT); // -2-
   static std::unordered_map<int, fQKernel> qupdate =
   {
      // 2D.
      {0x24,&QKernel<2,4>}, {0x26,&QKernel<2,6>},
      {0x28,&QKernel<2,8>}, {0x2A,&QKernel<2,10>}, 
      {0x2C,&QKernel<2,12>},{0x2E,&QKernel<2,14>},
      {0x30,&QKernel<2,16>},
      // 3D.
      {0x34,&QKernel<3,4>}, {0x36,&QKernel<3,6>}, {0x38,&QKernel<3,8>}
   };
   if (!qupdate[id])
   {
      mfem::out << "Unknown kernel 0x" << std::hex << id << std::endl;
      MFEM_ABORT("Unknown kernel");
   }

   qupdate[id](NE, NQ, use_viscosity, use_vorticity, qdata.h0, h1order, 
               cfl, infinity, qdata.vbc_max_val, qdata.mscale, qdata.gravity, gamma_gf, lambda_gf, mu_gf, ir.GetWeights(), q_dx,
               // dt, cfl, infinity, qdata.vbc_max_val, qdata.mscale, qdata.gravity, gamma_gf, lambda_gf, mu_gf, ir.GetWeights(), q_dx,
               qdata.rho0DetJ0w, q_e, q_sig, q_dv, qdata.Jac0inv, q_dt_est, q_h_est, qdata.stressJinvT, qdata.tauJinvT, qdata.buoyJinvT); // -3-
   qdata.dt_est = q_dt_est.Min();
   qdata.h_est = q_h_est.Min();
   timer->sw_qdata.Stop();
   timer->quad_tstep += NE;
}

void LagrangianGeoOperator::AssembleForceMatrix() const
{
   if (forcemat_is_assembled || p_assembly) { return; }
   Force = 0.0;
   timer.sw_force.Start();
   Force.Assemble();
   timer.sw_force.Stop();
   forcemat_is_assembled = true;
}

void LagrangianGeoOperator::AssembleSigmaMatrix() const
{
   // not used anymore
}

void LagrangianGeoOperator::GetErrorEstimates(ParGridFunction &e_gf, Vector & errors)
{
   int order = L2.GetOrder(0);
   // Space for the discontinuous (original) flux
   ConstantCoefficient one(1.0);
   DiffusionIntegrator flux_integrator(one);
   L2_FECollection flux_fec(order, pmesh->Dimension());
   // ND_FECollection flux_fec(order_, pmesh_->Dimension());
   ParFiniteElementSpace flux_fes(pmesh, &flux_fec, pmesh->SpaceDimension());

   // Space for the smoothed (conforming) flux
   int norm_p = 1;
   RT_FECollection smooth_flux_fec(order-1, pmesh->Dimension());
   ParFiniteElementSpace smooth_flux_fes(pmesh, &smooth_flux_fec);

   L2ZZErrorEstimator(flux_integrator, e_gf,
                      smooth_flux_fes, flux_fes, errors, norm_p);
}


// void LagrangianGeoOperator::TMOPUpdate(const Vector &S, bool quick)
// {
   
//    ParMesh *pmesh = H1.GetParMesh();
//    width = height = S.Size();

//    Vector* sptr = const_cast<Vector*>(&S);
//    ParGridFunction x_ini;
//    x_ini.MakeRef(&H1, *sptr, H1.GetVSize()*2 + L2.GetVSize() + L2.GetVSize()*3*(dim-1));
   
//    // Element number update
//    NE = pmesh->GetNE();
//    const int NQ = ir.GetNPoints();

//    if (quick) { return; }
   
//    // go back to initial mesh configuration temporarily
//    int own_nodes = 0;
//    GridFunction *x_gf = &x_ini;
//    pmesh->SwapNodes(x_gf, own_nodes);

//    // update mass matrix
//    Mv.Update();
//    Mv.Assemble();
//    // Mv_spmat_copy = Mv.SpMat();

//    // update force matrix
//    Force.Update();
//    Force.Assemble();
//    Force = 0.0;
   
//    // // Standard local assembly and inversion for energy mass matrices.
//    // // 'Me' is used in the computation of the internal energy
//    // // which is used twice: once at the start and once at the end of the run.
//    {
//       Me.SetSize(l2dofs_cnt, l2dofs_cnt, NE);
//       Me_inv.SetSize(l2dofs_cnt, l2dofs_cnt, NE);
//       MassIntegrator mi(rho0_coeff, &ir);
//       for (int e = 0; e < NE; e++)
//       {
//          DenseMatrixInverse inv(&Me(e));
//          const FiniteElement &fe = *L2.GetFE(e);
//          ElementTransformation &Tr = *L2.GetElementTransformation(e);
//          mi.AssembleElementMatrix(fe, Tr, Me(e));
//          inv.Factor();
//          inv.GetInverseMatrix(Me_inv(e));
//       }
//    }
   
//    // resize quadrature data and make sure 'stressJinvT' will be recomputed
//    qdata.Resize(dim, NE, NQ);
//    qdata_is_current = false;
//    // swap back to deformed mesh configuration
//    pmesh->SwapNodes(x_gf, own_nodes);
//    // Force.Assemble();

//    block_offsets[0] = 0;
//    block_offsets[1] = block_offsets[0] + H1.GetVSize();
//    block_offsets[2] = block_offsets[1] + H1.GetVSize();
//    block_offsets[3] = block_offsets[2] + L2.GetVSize();
//    block_offsets[4] = block_offsets[3] + L2.GetVSize()*3*(dim-1);
//    // block_offsets[5] = block_offsets[4] + H1.GetVSize();
// }

void LagrangianGeoOperator::TMOPUpdate(const Vector &S, bool quick)
{
   
   // Element number update
   NE = pmesh->GetNE();
   const int NQ = ir.GetNPoints();

   if (quick) { return; }
   
   // update mass matrix
   Mv.Update();
   Mv.Assemble();
   Mv_spmat_copy = Mv.SpMat();

   // update mass matrix
   fic_Mv.Update();
   fic_Mv.Assemble();
   fic_Mv_spmat_copy = Mv.SpMat();

   {
      // Me.SetSize(l2dofs_cnt, l2dofs_cnt, NE);
      // Me_inv.SetSize(l2dofs_cnt, l2dofs_cnt, NE);
      MassIntegrator mi(rho0_coeff, &ir);
      // MassIntegrator mi(scale_rho0_coeff, &ir);
      for (int e = 0; e < NE; e++)
      {
         DenseMatrixInverse inv(&Me(e));
         const FiniteElement &fe = *L2.GetFE(e);
         ElementTransformation &Tr = *L2.GetElementTransformation(e);
         mi.AssembleElementMatrix(fe, Tr, Me(e));
         inv.Factor();
         inv.GetInverseMatrix(Me_inv(e));
      }
   }

   // update 'rho0DetJ0' and 'Jac0inv' at all quadrature points
   // TODO: remove code duplication
   Vector rho_vals(NQ);
   for (int i = 0; i < NE; i++)
   {
      rho0_gf.GetValues(i, ir, rho_vals);
      ElementTransformation *T = H1.GetElementTransformation(i);
      for (int q = 0; q < NQ; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         T->SetIntPoint(&ip);

         DenseMatrixInverse Jinv(T->Jacobian());
         Jinv.GetInverseMatrix(qdata.Jac0inv(i*NQ + q));

         const double rho0DetJ0 = T->Weight() * rho_vals(q);
         qdata.rho0DetJ0w(i*NQ + q) = rho0DetJ0 *
                                           ir.IntPoint(q).weight;
      }
   }
}


} // namespace geodynamics


void GeoODESolver::Init(TimeDependentOperator &tdop)
{
   ODESolver::Init(tdop);
   geo_oper = dynamic_cast<geodynamics::LagrangianGeoOperator *>(f);
   MFEM_VERIFY(geo_oper, "GeoSolvers expect LagrangianGeoOperator.");
}

void RK2AvgSolver::Init(TimeDependentOperator &tdop)
{
   GeoODESolver::Init(tdop);

   // If you add block vector, you should add offset 
   const Array<int> &block_offsets = geo_oper->GetBlockOffsets();

   V.SetSize(block_offsets[1], mem_type);
   V.UseDevice(true);
   dS_dt.Update(block_offsets, mem_type);
   dS_dt = 0.0;
   S0.Update(block_offsets, mem_type);
}

void RK2AvgSolver::Step(Vector &S, double &t, double &dt)
{
   // The monolithic BlockVector stores the unknown fields as follows:
   // (Position, Velocity, Specific Internal Energy).
   S0.Vector::operator=(S);
   Vector &v0 = S0.GetBlock(1);
   Vector &dx_dt = dS_dt.GetBlock(0);
   Vector &dv_dt = dS_dt.GetBlock(1);

   // In each sub-step:
   // - Update the global state Vector S.
   // - Compute dv_dt using S.
   // - Update V using dv_dt.
   // - Compute de_dt and dx_dt using S and V.

   // -- 1.
   // S is S0.
   geo_oper->UpdateMesh(S);
   geo_oper->SolveVelocity(S, dS_dt);
   // V = v0 + 0.5 * dt * dv_dt;
   add(v0, 0.5 * dt, dv_dt, V);
   geo_oper->SolveEnergy(S, V, dS_dt);
   geo_oper->SolveStress(S, dS_dt);
   dx_dt = V;

   // -- 2.
   // S = S0 + 0.5 * dt * dS_dt;
   add(S0, 0.5 * dt, dS_dt, S);
   geo_oper->ResetQuadratureData();
   geo_oper->UpdateMesh(S);
   geo_oper->SolveVelocity(S, dS_dt);
   // V = v0 + 0.5 * dt * dv_dt;
   add(v0, 0.5 * dt, dv_dt, V);
   geo_oper->SolveEnergy(S, V, dS_dt);
   geo_oper->SolveStress(S, dS_dt);
   dx_dt = V;

   // S = S0 + dt * dS_dt.
   add(S0, dt, dS_dt, S);
   geo_oper->ResetQuadratureData();
   t += dt;

   // geo_oper->RadialReturn(S, dS_dt, 1.0);
   // geo_oper->RadialReturn(S, dS_dt, 1.0);
}

} // namespace mfem

#endif // MFEM_USE_MPI
