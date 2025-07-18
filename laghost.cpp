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
//
//                     __                __
//                    / /   ____  ____  / /_  ____  _____
//                   / /   / __ `/ __ `/ __ \/ __ \/ ___/
//                  / /___/ /_/ / /_/ / / / / /_/ (__  )
//                 /_____/\__,_/\__, /_/ /_/\____/____/
//                             /____/
//
//             High-order Lagrangian Geodynamics Solver
//
// Laghos(LAGrangian High-Order Solver) is a miniapp that solves the
// time-dependent Euler equation of compressible gas dynamics in a moving
// Lagrangian frame using unstructured high-order finite element spatial
// discretization and explicit high-order time-stepping. Laghos is based on the
// numerical algorithm described in the following article:
//
//    V. Dobrev, Tz. Kolev and R. Rieben, "High-order curvilinear finite element
//    methods for Lagrangian geodynamics", SIAM Journal on Scientific
//    Computing, (34) 2012, pp. B606–B641, https://doi.org/10.1137/120864672.
//
//                     __                __               __
//                    / /   ____ _____ _/ /_  ____  _____/ /_
//                   / /   / __ `/ __ `/ __ \/ __ \/ ___/ __/
//                  / /___/ /_/ / /_/ / / / / /_/ (__  ) /_
//                 /_____/\__,_/\__, /_/ /_/\____/____/\__/
//                             /____/
//             Lagrangian High-order Solver for Tectonics
//
// Laghost inherits the main structure of LAGHOS. However, it solves
// the dynamic form of the general momenntum balance equation for continnum,
// aquiring quasi-static solution with dynamic relaxation; reasonably large
// time steps with mass scaling. The target applications include long-term
// brittle and ductile deformations of rocks coupled with time-evolving
// thermal state.
//
// -- How to run LAGHOST
// mpirun -np 8 laghost -i ./defaults.cfg

#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>
#include <cmath>
#include <memory>
// #include "laghost_parameters.hpp"
#include "laghost_solver.hpp"
#include "laghost_rheology.hpp"
#include "laghost_function.hpp"
#include "laghost_input.hpp"
#include "laghost_tmop.hpp"
#include "laghost_remhos.hpp"

using std::cout;
using std::endl;
using namespace mfem;

// Forward declarations for our new functions
struct AppState;
void initialize(AppState& appState, int argc, char *argv[]);
void run(AppState& appState);
void finalize(AppState& appState);
void print_progress(const int ti, AppState &appState);


// Choice for the problem setup.
static int problem, dim;
static long GetMaxRssMB();
static void display_banner(std::ostream&);
static void Checks(const int ti, const double norm, int &checks);

class ConductionOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

   ParBilinearForm *M;
   ParBilinearForm *K;

   HypreParMatrix Mmat;
   HypreParMatrix Kmat;
   HypreParMatrix *T; // T = M + dt K
   double current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver T_solver;    // Implicit solver for T = M + dt K
   HypreSmoother T_prec; // Preconditioner for the implicit solver

   double alpha, kappa;

   mutable Vector z; // auxiliary vector

public:
   ConductionOperator(ParFiniteElementSpace &f, double alpha, double kappa,
                      const Vector &u);

   virtual void Mult(const Vector &u, Vector &du_dt) const;
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

   /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
   void SetParameters(const Vector &u);

   virtual ~ConductionOperator();
};

void TMOPUpdate(BlockVector &S, BlockVector &S_old,
               Array<int> &offset,
               ParGridFunction &x_gf,
               ParGridFunction &v_gf,
               ParGridFunction &e_gf,
               ParGridFunction &s_gf,
               ParGridFunction &x_ini_gf,
               ParGridFunction &p_gf,
               ParGridFunction &n_p_gf,
               ParGridFunction &ini_p_gf,
               ParGridFunction &u_gf,
               ParGridFunction &rho0_gf,
               ParGridFunction &lambda0_gf,
               ParGridFunction &mu0_gf,
               ParGridFunction &mat_gf,
               // ParLinearForm &flattening,
               int dim, bool amr);

static void Generate_and_refine_initial_mesh(Mesh *&mesh, Param &param)
{
   if (param.mesh.mesh_file.compare("default") != 0)
      mesh = new Mesh(param.mesh.mesh_file.c_str(), true, true);
   else {
      switch (param.sim.dim) {
         case 1:
               mesh = new Mesh(Mesh::MakeCartesian1D(2));
               if( mesh ) {
                  mesh->GetBdrElement(0)->SetAttribute(1);
                  mesh->GetBdrElement(1)->SetAttribute(1);
               }
               else
                  MFEM_ABORT("Failed to create 1D mesh.");
               break;
         case 2:
               mesh = new Mesh(Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL, true));
               if( mesh ) {
                  const int NBE = mesh->GetNBE();
                  for (int b = 0; b < NBE; b++) {
                     Element *bel = mesh->GetBdrElement(b);
                     const int attr = (b < NBE / 2) ? 2 : 1;
                     std::cout << NBE << "," << b << "," << attr << std::endl;
                     bel->SetAttribute(attr);
                  }
               }
               else
                  MFEM_ABORT("Failed to create 2D mesh.");
               break;
         case 3:
               mesh = new Mesh(Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON, true));
               if( mesh ) {
                  const int NBE = mesh->GetNBE();
                  for (int b = 0; b < NBE; b++) {
                     Element *bel = mesh->GetBdrElement(b);
                     const int attr = (b < NBE / 3) ? 3 : (b < 2 * NBE / 3) ? 1 : 2;
                     bel->SetAttribute(attr);
                  }
               }
               else
                  MFEM_ABORT("Failed to create 3D mesh.");
               break;
         default:
               break;
      }
   }
   dim = mesh->Dimension();
   MFEM_VERIFY( param.sim.dim == dim, "Dimension mismatch.");

   // 1D vs partial assembly sanity check.
   if (param.solver.p_assembly && dim == 1) {
      param.solver.p_assembly = false;
      if (Mpi::Root())
         cout << "Laghos does not support PA in 1D. Switching to FA." << endl;
   }

   // Refine the mesh in serial to increase the resolution.
   for (int lev = 0; lev < param.mesh.rs_levels; lev++) {
      mesh->UniformRefinement();
   }
   if (Mpi::Root())
      cout << "Number of zones in the serial mesh: " << mesh->GetNE() << endl;

   // and then refine locally. Non-conforming elements might be generated at this stage.
   if (param.mesh.local_refinement) {
      mesh->EnsureNCMesh(true);
      Array<int> refs;
      for (int i = 0; i < mesh->GetNE(); i++) {
         if (mesh->GetAttribute(i) >= 2)
               refs.Append(i);
      }
      mesh->GeneralRefinement(refs, 1);
      refs.DeleteAll();

      for (int i = 0; i < mesh->GetNE(); i++) {
         if (mesh->GetAttribute(i) >= 3)
               refs.Append(i);
      }
      mesh->GeneralRefinement(refs, 1);
      refs.DeleteAll();

      mesh->Finalize(true);
   }
}

static void Partition_initial_mesh(std::unique_ptr<ParMesh>& pmesh, Mesh *&mesh, Param &param)
{
   const int num_tasks = Mpi::WorldSize();
   int unit = 1;
   const int dim = mesh->Dimension();
   int *nxyz = new int[dim];
   switch (param.mesh.partition_type)
   {
      case 0:
         for (int d = 0; d < dim; d++)
         {
            nxyz[d] = unit;
         }
         break;
      case 11:
      case 111:
         unit = static_cast<int>(floor(pow(num_tasks, 1.0 / dim) + 1e-2));
         for (int d = 0; d < dim; d++)
         {
            nxyz[d] = unit;
         }
         break;
      case 21: // 2D
         unit = static_cast<int>(floor(pow(num_tasks / 2, 1.0 / 2) + 1e-2));
         nxyz[0] = 2 * unit;
         nxyz[1] = unit;
         break;
      case 31: // 2D
         unit = static_cast<int>(floor(pow(num_tasks / 3, 1.0 / 2) + 1e-2));
         nxyz[0] = 3 * unit;
         nxyz[1] = unit;
         break;
      case 32: // 2D
         unit = static_cast<int>(floor(pow(2 * num_tasks / 3, 1.0 / 2) + 1e-2));
         nxyz[0] = 3 * unit / 2;
         nxyz[1] = unit;
         break;
      case 49: // 2D
         unit = static_cast<int>(floor(pow(9 * num_tasks / 4, 1.0 / 2) + 1e-2));
         nxyz[0] = 4 * unit / 9;
         nxyz[1] = unit;
         break;
      case 51: // 2D
         unit = static_cast<int>(floor(pow(num_tasks / 5, 1.0 / 2) + 1e-2));
         nxyz[0] = 5 * unit;
         nxyz[1] = unit;
         break;
      case 211: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 2, 1.0 / 3) + 1e-2));
         nxyz[0] = 2 * unit;
         nxyz[1] = unit;
         nxyz[2] = unit;
         break;
      case 221: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 4, 1.0 / 3) + 1e-2));
         nxyz[0] = 2 * unit;
         nxyz[1] = 2 * unit;
         nxyz[2] = unit;
         break;
      case 311: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 3, 1.0 / 3) + 1e-2));
         nxyz[0] = 3 * unit;
         nxyz[1] = unit;
         nxyz[2] = unit;
         break;
      case 321: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 6, 1.0 / 3) + 1e-2));
         nxyz[0] = 3 * unit;
         nxyz[1] = 2 * unit;
         nxyz[2] = unit;
         break;
      case 322: // 3D.
         unit = static_cast<int>(floor(pow(2 * num_tasks / 3, 1.0 / 3) + 1e-2));
         nxyz[0] = 3 * unit / 2;
         nxyz[1] = unit;
         nxyz[2] = unit;
         break;
      case 432: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 3, 1.0 / 3) + 1e-2));
         nxyz[0] = 2 * unit;
         nxyz[1] = 3 * unit / 2;
         nxyz[2] = unit;
         break;
      case 511: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 5, 1.0 / 3) + 1e-2));
         nxyz[0] = 5 * unit;
         nxyz[1] = unit;
         nxyz[2] = unit;
         break;
      case 521: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 10, 1.0 / 3) + 1e-2));
         nxyz[0] = 5 * unit;
         nxyz[1] = 2 * unit;
         nxyz[2] = unit;
         break;
      case 522: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 20, 1.0 / 3) + 1e-2));
         nxyz[0] = 5 * unit;
         nxyz[1] = 2 * unit;
         nxyz[2] = 2 * unit;
         break;
      case 911: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 9, 1.0 / 3) + 1e-2));
         nxyz[0] = 9 * unit;
         nxyz[1] = unit;
         nxyz[2] = unit;
         break;
      case 921: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 18, 1.0 / 3) + 1e-2));
         nxyz[0] = 9 * unit;
         nxyz[1] = 2 * unit;
         nxyz[2] = unit;
         break;
      case 922: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 36, 1.0 / 3) + 1e-2));
         nxyz[0] = 9 * unit;
         nxyz[1] = 2 * unit;
         nxyz[2] = 2 * unit;
         break;
      default:
         if ( Mpi::Root() )
            cout << "Unknown partition type: " << param.mesh.partition_type << '\n';
         delete mesh;
         MPI_Finalize();
         MFEM_ABORT("Unknown partition type.");
         break;
   }

   Array<int> cxyz; // Leave undefined. It won't be used.
   int product = 1;
   for (int d = 0; d < dim; d++)
   {
      product *= nxyz[d];
   }
   const bool cartesian_partitioning = (cxyz.Size() > 0) ? true : false;
   if (product == num_tasks || cartesian_partitioning)
   {
      if (cartesian_partitioning)
      {
         int cproduct = 1;
         for (int d = 0; d < dim; d++)
         {
            cproduct *= cxyz[d];
         }
         MFEM_VERIFY(!cartesian_partitioning || cxyz.Size() == dim,
                  "Expected " << mesh->SpaceDimension() << " integers with the "
                  "option --cartesian-partitioning.");
         MFEM_VERIFY(!cartesian_partitioning || num_tasks == cproduct,
                  "Expected cartesian partitioning product to match number of ranks.");
      }
      int *partitioning = cartesian_partitioning ?
                     mesh->CartesianPartitioning(cxyz) :
                     mesh->CartesianPartitioning(nxyz);
      pmesh.reset(new ParMesh(MPI_COMM_WORLD, *mesh, partitioning));
      delete[] partitioning;
   }
   else
   {
      if ( Mpi::Root() )
      {
         cout << "Non-Cartesian partitioning through METIS will be used.\n";
#ifndef MFEM_USE_METIS
         MFEM_ABORT("MFEM was built without METIS.\nAdjust the number of tasks to use a Cartesian split.");
#endif
      }
      pmesh.reset(new ParMesh(MPI_COMM_WORLD, *mesh));
   }
   delete[] nxyz;
   delete mesh;

   // Refine the mesh further in parallel to increase the resolution.
   for(int lev = 0; lev < param.mesh.rp_levels; lev++)
      pmesh->UniformRefinement();
   // pmesh->Rebalance();

   int NE = pmesh->GetNE(), ne_min, ne_max;
   MPI_Reduce(&NE, &ne_min, 1, MPI_INT, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&NE, &ne_max, 1, MPI_INT, MPI_MAX, 0, pmesh->GetComm());
   if( Mpi::Root() )
      cout << "Zones min/max: " << ne_min << " " << ne_max << endl;

}

static void Collect_boundingbox_info( ParMesh *pmesh, Param &param, Vector &bb_center, Vector &bb_length)
{
   const int dim = pmesh->Dimension();
   MFEM_VERIFY( bb_center.Size() == dim && bb_length.Size() == dim, "The size of bb_center and bb_length should be equal to dim.");
   // Mesh bounding box
   Vector bb_min(dim), bb_max(dim);
   pmesh->GetBoundingBox(bb_min, bb_max, max(param.mesh.order_v, 1));
   for( int i = 0; i < dim; i++)
   {
      bb_center[i] = (bb_min[i] + bb_max[i]) * 0.5;
      bb_length[i] = (bb_max[i] - bb_min[i]);
   }
}

static void get_ess_dofs_for_component(ParFiniteElementSpace &H1FESpace, Array<int> &ess_bdr, int component, Array<int> &ess_tdofs, Array<int> &ess_vdofs, Array<int> &vdofs_list)
{
   // Temporary storage for essential dofs.
   // vdofs_list is given because it's used outside this function.
   Array<int> dofs_marker, tdofs_list;

   H1FESpace.GetEssentialTrueDofs(ess_bdr, tdofs_list, component);
   ess_tdofs.Append(tdofs_list);
   H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker, component);
   FiniteElementSpace::MarkerToList(dofs_marker, vdofs_list);
   ess_vdofs.Append(vdofs_list);
}

static void set_vbc_val(ParGridFunction &v_gf, Array<int> &vdofs_list, Param &param, double &max_vbc_val)
{
   double val = param.bc.vbc_x0_val0 * param.bc.vbc_factor;
   for (int i = 0; i < vdofs_list.Size(); i++)
      v_gf[vdofs_list[i]] = val;

   max_vbc_val = std::max(max_vbc_val, val);
}

static void unknown_bc_abort(ParMesh *pmesh, const int &myid, const char *bc_name, const char *bc_symb, const int &bc_attr)
{
   if (myid == 0)
      cout << "Unknown BC type for the "<<bc_name<<" boundary "<<bc_symb<<" (attribute "<<bc_attr<<")" << endl;
   delete pmesh;
   MPI_Finalize();
}

static void set_essential_bc( const int myid, ParMesh *pmesh, ParFiniteElementSpace &H1FESpace, Param &param, Array<int> &ess_bdr, Array<int> &ess_tdofs, Array<int> &ess_vdofs, ParGridFunction &v_gf, double &max_vbc_val)
{
   //
   // MFEM assumes bdr_attribute is sequentially numbered from 1.
   // And ess_bdr's (bdr_attribute-1)-th element should be non-zero to set the essential bc on that boundary.
   // Example: Boundary attributes for 2D, 1 to 4, corresponding to internal designation, x0, x1, z0 and z1.
   //   ^ z (x1)
   //   |     4 (z1, top)     Notation: bdy id for meshing (internal designation, geographic meaning)
   //   +-----------------+
   //   |                 |
   // 1 | (x0, west)      | 2 (x1, east)
   //   |                 |
   //   +-----------------+--> x (x0)
   //         3 (z0, bottom)
   //   - To set velocity BC on boundaries 1 and 2, ess_bdr[0] = ess_bdr[1] = 1: i.e., ess_bdr = {1, 1, 0, 0}
   //   - To set velocity BC on boundaries 1, 2, and 3, ess_bdr = {1, 1, 1, 0}
   // However, MFEM's FiniteElementSpace::GetEssentialVDofs() can take only one component at a time while
   // looping over all the boundary elements to check the assigned attribute, setting the essential dofs
   // has to be done for each component separately. Also, the input parameter system assumes the essential
   // boundaries to be defined individually: e.g., vbc_x1 = 1 (0 by default), vbc_z0 = 1, etc.
   // So, we work on
   //

   const int dim = pmesh->Dimension();
   int component;
   ess_bdr = 0;
   switch (param.bc.vbc_x0)
   {
      case 0: // All velocity components free. Don't do anything.
         break;
      case 1: // Normal component fixed, shear components free"
      {
         ess_bdr[0] = 1; // x0's attribute is 1, so the 0-th element is set to 1.
         component = 0;  // component 0 is fixed since x-normal.
         break;
      }
      default:
         unknown_bc_abort( pmesh, myid, "western", "x0", 1 );
   }
   switch (param.bc.vbc_x1)
   {
      case 0: // All velocity components free. Don't do anything.
         break;
      case 1:
      {
         ess_bdr[1] = 1; // x1's attribute is 2, so the 1st element is set to 1.
         component = 0;  // component 0 is fixed since x-normal.
         break;
      }
      default:
         unknown_bc_abort( pmesh, myid, "eastern", "x1", 2 );
   }
   // Here we call GetEssentialVDofs() for 0-th component of velocity (i.e., vx).
   {
      Array<int> vdofs_list;
      get_ess_dofs_for_component(H1FESpace, ess_bdr, component, ess_tdofs, ess_vdofs, vdofs_list);
      set_vbc_val(v_gf, vdofs_list, param, max_vbc_val);
   }
   ess_bdr = 0;
   switch (param.bc.vbc_z0)
   {
      case 0: // All velocity components free. Don't do anything.
         break;
      case 1:
      {
         ess_bdr[2] = 1;  // y0's attribute is 3, so the 2nd element is set to 1.
         component = 1;   // since y-normal boundary, component 1 is fixed.
         break;
      }
      default:
         unknown_bc_abort( pmesh, myid, "bottom", "z0", 3 );
   }
   switch (param.bc.vbc_z1)
   {
      case 0: // All velocity components free. Don't do anything.
         break;
      case 1:
      {
         ess_bdr[3] = 1;  // y1's attribute is 4, so the 3rd element is set to 1.
         component = 1;   // since y-normal boundary, component 1 is fixed.
         break;
      }
      default:
         unknown_bc_abort( pmesh, myid, "top", "z1", 4 );
   }
   // Here we call GetEssentialVDofs() for 1st component of velocity (i.e., vz).
   {
      Array<int> vdofs_list;
      get_ess_dofs_for_component(H1FESpace, ess_bdr, component, ess_tdofs, ess_vdofs, vdofs_list);
      set_vbc_val(v_gf, vdofs_list, param, max_vbc_val);
   }
   if( dim == 3 ) {
      //          ^ z (x1)
      //          |      Notation: bdy id for meshing (internal designation, geographic meaning)
      //           +-------------------------+
      //         / |                        /|
      //        /  |                       / |
      //       /   |    4 (z1, top)       /  |
      //      /    |                     /   |
      //     /  1  |     5 (y0, north)  /    |
      //    +----- |------------------+      |
      //    | (x0, |__________________|__2___|___> x (x0)
      //    | west)/                  | (x1, /
      //    |     /                   |east)/
      //    |    /     3 (z0, bottom) |    /
      //    |   /                     |   /
      //    |  /   6 (y1, south)      |  /
      //    | /                       | /
      //    |/                        |/
      //    +--------- ---------------+
      //   / y (x2)
      //
      ess_bdr = 0;
      switch (param.bc.vbc_y0)
      {
         case 0: // All velocity components free. Don't do anything.
            break;
         case 1:
         {
            ess_bdr[4] = 1;  // y0's attribute is 5, so the 4th element is set to 1.
            component = 2;   // since y-normal boundary, the third component (index 2) is fixed.
            break;
         }
         default:
            unknown_bc_abort( pmesh, myid, "northern", "y0", 5 );
      }
      switch (param.bc.vbc_y1)
      {
         case 0: // All velocity components free. Don't do anything.
            break;
         case 1:
         {
            ess_bdr[5] = 1;  // z1's attribute is 6, so the 5th element is set to 1.
            component = 2;   // since y-normal boundary, the third component (index 2) is fixed.
            break;
         }
         default:
            unknown_bc_abort( pmesh, myid, "southern", "y1", 6 );
      }
         // Here we call GetEssentialVDofs() for 0-th component of velocity (i.e., vx).
      {
         Array<int> vdofs_list;
         get_ess_dofs_for_component(H1FESpace, ess_bdr, component, ess_tdofs, ess_vdofs, vdofs_list);
         set_vbc_val(v_gf, vdofs_list, param, max_vbc_val);
      }
   }
}

void initial_velocity(const Vector &x, Vector &v)
{
   // Initial velocity can be set here. For instance,
   //    v(0) =  sin(M_PI*x(0)) * cos(M_PI*x(1));
   //    v(1) = -cos(M_PI*x(0)) * sin(M_PI*x(1));
   //    if (x.Size() == 3)
   //    {
   //       v(0) *= cos(M_PI*x(2));
   //       v(1) *= cos(M_PI*x(2));
   //       v(2) = 0.0;
   //    }
   // For now, initial velocity is uniformly zero.
   v = 0.0;
}

// =============================================================================
// AppState struct to hold all shared variables
// =============================================================================
struct AppState
{
   // Core Simulation Objects
   Param param;
   std::unique_ptr<ParMesh> pmesh;
   std::unique_ptr<geodynamics::LagrangianGeoOperator> geo;
   std::unique_ptr<ODESolver> ode_solver;
   std::unique_ptr<ODESolver> ode_solver_sub;
   std::unique_ptr<ODESolver> ode_solver_sub2;

   // FE Collections and Spaces
   std::unique_ptr<H1_FECollection> H1FEC;
   std::unique_ptr<L2_FECollection> L2FEC;
   std::unique_ptr<L2_FECollection> L2FEC_positive;
   std::unique_ptr<ParFiniteElementSpace> H1FESpace;
   std::unique_ptr<ParFiniteElementSpace> L2FESpace;
   std::unique_ptr<ParFiniteElementSpace> L2FESpace_stress;
   std::unique_ptr<ParFiniteElementSpace> L2FESpace_mat;

   // State Vectors and GridFunctions
   BlockVector S;
   Array<int> offset;
   Array<int> ess_tdofs;

   ParGridFunction x_gf, v_gf, e_gf, s_gf;
   ParGridFunction u_gf, p_gf, n_p_gf, ini_p_gf, s_old_gf, p_gf_old, ini_p_old_gf, x_old_gf;
   ParGridFunction rho0_gf, fictitious_rho0_gf, lambda0_gf, mu0_gf, mat_gf;
   ParGridFunction x_ini_gf, vol_ini_gf, skew_ini_gf;
   ParGridFunction comp_gf, comp_ref_gf;

   // Submesh and Boundary Objects
   std::unique_ptr<ParSubMesh> submesh;
   std::unique_ptr<ParSubMesh> submesh_bottom;
   std::unique_ptr<ParFiniteElementSpace> sub_fespace0;
   std::unique_ptr<ParFiniteElementSpace> sub_fespace1;
   std::unique_ptr<ParFiniteElementSpace> sub_fespace2;
   std::unique_ptr<ParFiniteElementSpace> sub_fespace3;

   ParGridFunction x_top, topo, x_bottom, bottom;
   Vector topo_t, topo_t_old, bottom_t, bottom_t_old;
   std::unique_ptr<ConductionOperator> oper_sub;
   std::unique_ptr<ConductionOperator> oper_sub2;

   // Time and Loop Control
   double t = 0.0, dt = 1.0;
   bool last_step = false;
   int steps = 0;

   // Visualization and I/O
   std::unique_ptr<VisItDataCollection> visit_dc;
   ParaViewDataCollection *pd = nullptr;
   socketstream vis_rho, vis_v, vis_e;
   std::unique_ptr<VectorFunctionCoefficient> v_coeff;

   // Finalization Variables
   double energy_init = 0.0;
   int checks = 0;
};


// =============================================================================
// initialize() function
// =============================================================================
void initialize(AppState& appState, int argc, char *argv[])
{
   Mpi::Init();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   if (Mpi::Root()) { display_banner(cout); }

   OptionsParser args(argc, argv);
   read_and_assign_input_parameters(args, appState.param, myid);

   Device backend;
   backend.Configure(appState.param.sim.device, appState.param.sim.dev);
   if (Mpi::Root()) { backend.Print(); }
   backend.SetGPUAwareMPI(appState.param.sim.gpu_aware_mpi);

   Mesh *mesh = nullptr;
   Generate_and_refine_initial_mesh(mesh, appState.param);
   Partition_initial_mesh(appState.pmesh, mesh, appState.param);

   Vector bb_center(appState.pmesh->Dimension());
   Vector bb_length(appState.pmesh->Dimension());
   Collect_boundingbox_info(appState.pmesh.get(), appState.param, bb_center, bb_length);

   appState.L2FEC.reset(new L2_FECollection(appState.param.mesh.order_e, dim, BasisType::GaussLobatto));
   appState.H1FEC.reset(new H1_FECollection(appState.param.mesh.order_v, dim));
   appState.L2FEC_positive.reset(new L2_FECollection(appState.param.mesh.order_e, dim, BasisType::Positive));

   appState.L2FESpace.reset(new ParFiniteElementSpace(appState.pmesh.get(), appState.L2FEC.get()));
   appState.L2FESpace_stress.reset(new ParFiniteElementSpace(appState.pmesh.get(), appState.L2FEC.get(), 3*(dim-1)));
   appState.H1FESpace.reset(new ParFiniteElementSpace(appState.pmesh.get(), appState.H1FEC.get(), appState.pmesh->Dimension()));
   ParFiniteElementSpace L2FESpace_positive(appState.pmesh.get(), appState.L2FEC_positive.get());

   ODESolver* temp_ode_solver = nullptr;
   switch (appState.param.solver.ode_solver_type)
   {
      case 1: temp_ode_solver = new ForwardEulerSolver; break;
      case 2: temp_ode_solver = new RK2Solver(0.5); break;
      case 3: temp_ode_solver = new RK3SSPSolver; break;
      case 4: temp_ode_solver = new RK4Solver; break;
      case 6: temp_ode_solver = new RK6Solver; break;
      case 7: temp_ode_solver = new RK2AvgSolver; break;
      default:
         if (myid == 0) { cout << "Unknown ODE solver type: " << appState.param.solver.ode_solver_type << '\n'; }
         MPI_Finalize();
         exit(3);
   }
   appState.ode_solver.reset(temp_ode_solver);

   appState.ode_solver_sub.reset(new RK2Solver(0.5));
   appState.ode_solver_sub2.reset(new RK2Solver(0.5));

   appState.L2FESpace->Update();
   appState.L2FESpace_stress->Update();
   appState.H1FESpace->Update();

   const int Vsize_l2 = appState.L2FESpace->GetVSize();
   const int Vsize_h1 = appState.H1FESpace->GetVSize();
   appState.offset.SetSize(5);
   appState.offset[0] = 0;
   appState.offset[1] = appState.offset[0] + Vsize_h1;
   appState.offset[2] = appState.offset[1] + Vsize_h1;
   appState.offset[3] = appState.offset[2] + Vsize_l2;
   appState.offset[4] = appState.offset[3] + Vsize_l2*3*(dim-1);
   appState.S.Update(appState.offset, Device::GetMemoryType());

   appState.x_gf.MakeRef(appState.H1FESpace.get(), appState.S, appState.offset[0]);
   appState.v_gf.MakeRef(appState.H1FESpace.get(), appState.S, appState.offset[1]);
   appState.e_gf.MakeRef(appState.L2FESpace.get(), appState.S, appState.offset[2]);
   appState.s_gf.MakeRef(appState.L2FESpace_stress.get(), appState.S, appState.offset[3]);
   appState.pmesh->SetNodalGridFunction(&appState.x_gf);
   appState.x_gf.SyncAliasMemory(appState.S);

   Array<int> bdr_attrs(1);
   bdr_attrs[0] = 4;
   appState.submesh.reset(new ParSubMesh(ParSubMesh::CreateFromBoundary(*appState.pmesh, bdr_attrs)));
   appState.sub_fespace0.reset(new ParFiniteElementSpace(appState.submesh.get(), appState.H1FEC.get(), appState.pmesh->Dimension()));
   appState.sub_fespace1.reset(new ParFiniteElementSpace(appState.submesh.get(), appState.H1FEC.get()));
   appState.x_top.SetSpace(appState.sub_fespace0.get());
   appState.topo.SetSpace(appState.sub_fespace1.get());
   appState.submesh->SetNodalGridFunction(&appState.x_top);
   for (int i = 0; i < appState.topo.Size(); i++){appState.topo[i] = appState.x_top[i+appState.topo.Size()];}
   appState.topo.GetTrueDofs(appState.topo_t); appState.topo_t_old = appState.topo_t;

   Array<int> bdr_attrs_b(1);
   bdr_attrs_b[0] = 3;
   appState.submesh_bottom.reset(new ParSubMesh(ParSubMesh::CreateFromBoundary(*appState.pmesh, bdr_attrs_b)));
   appState.sub_fespace2.reset(new ParFiniteElementSpace(appState.submesh_bottom.get(), appState.H1FEC.get(), appState.pmesh->Dimension()));
   appState.sub_fespace3.reset(new ParFiniteElementSpace(appState.submesh_bottom.get(), appState.H1FEC.get()));
   appState.x_bottom.SetSpace(appState.sub_fespace2.get());
   appState.bottom.SetSpace(appState.sub_fespace3.get());
   appState.submesh_bottom->SetNodalGridFunction(&appState.x_bottom);
   for (int i = 0; i < appState.bottom.Size(); i++){appState.bottom[i] = appState.x_bottom[i+appState.bottom.Size()];}
   appState.bottom.GetTrueDofs(appState.bottom_t); appState.bottom_t_old = appState.bottom_t;

   appState.oper_sub.reset(new ConductionOperator(*appState.sub_fespace1, appState.param.bc.surf_alpha, appState.param.bc.surf_diff, appState.topo_t));
   appState.oper_sub2.reset(new ConductionOperator(*appState.sub_fespace3, appState.param.bc.base_alpha, appState.param.bc.base_diff, appState.bottom_t));

   ParFiniteElementSpace L2FESpace_xyz(appState.pmesh.get(), appState.L2FEC_positive.get(), dim);
   ParGridFunction xyz_gf_l2(&L2FESpace_xyz);
   VectorFunctionCoefficient xyz_coeff(appState.pmesh->Dimension(), xyz0);
   xyz_gf_l2.ProjectCoefficient(xyz_coeff);

   appState.v_gf = 0.0;
   appState.v_coeff.reset(new VectorFunctionCoefficient(appState.pmesh->Dimension(), initial_velocity));
   appState.v_gf.ProjectCoefficient(*appState.v_coeff);

   double max_vbc_val = appState.param.control.max_vbc_val;
   Array<int> ess_vdofs;
   Array<int> ess_bdr(appState.pmesh->bdr_attributes.Max());
   set_essential_bc( myid, appState.pmesh.get(), *appState.H1FESpace, appState.param, ess_bdr, appState.ess_tdofs, ess_vdofs, appState.v_gf, max_vbc_val );
   appState.ess_tdofs.Read();
   appState.v_gf.SyncAliasMemory(appState.S);

   int num_materials = appState.pmesh->attributes.Max();
   if(num_materials != appState.param.mat.nmat)
   {
      if (myid == 0) {
         cout << __FILE__<<":"<<__LINE__<< endl;
         cout <<"\tThe number of mesh attributes, "<<num_materials<<", are not consistent with the number of materials, "<<appState.param.mat.nmat<<", in the input file."<< endl;
      }
      MPI_Finalize();
      exit(3);
   }
   Vector rho0(appState.pmesh->attributes.Max());
   Vector fictitious_rho0(appState.pmesh->attributes.Max());
   double pseudo_speed =  max_vbc_val * appState.param.control.mscale;
   double pseudo_speed_sqrd =  pseudo_speed * pseudo_speed;
   for (int i = 0; i < appState.pmesh->attributes.Max(); i++) {
      rho0[i] = appState.param.mat.rho[i];
      fictitious_rho0[i] = (appState.param.mat.lambda[i] + 2*appState.param.mat.mu[i]) / pseudo_speed_sqrd;
   }
   PWConstCoefficient rho0_coeff(rho0);
   appState.rho0_gf.SetSpace(appState.L2FESpace.get());
   appState.rho0_gf.ProjectCoefficient(rho0_coeff);

   PWConstCoefficient fictitious_rho0_coeff(fictitious_rho0);
   appState.fictitious_rho0_gf.SetSpace(appState.L2FESpace.get());
   appState.fictitious_rho0_gf.ProjectCoefficient(fictitious_rho0_coeff);

   ParGridFunction l2_e(&L2FESpace_positive);
   if (appState.param.sim.problem == 1)
   {
      DeltaCoefficient e_coeff(0.0, 0.0, 0.0, 0.0);
      l2_e.ProjectCoefficient(e_coeff);
   }
   else
   {
      FunctionCoefficient e_coeff(e0);
      l2_e.ProjectCoefficient(e_coeff);
   }
   appState.e_gf.ProjectGridFunction(l2_e);
   appState.e_gf.SyncAliasMemory(appState.S);

   if(appState.param.mat.lambda.Size() != appState.pmesh->attributes.Max() ||
      appState.param.mat.mu.Size() != appState.pmesh->attributes.Max())
   {
      if (myid == 0){cout << "Material property arrays are not consistent with material IDs." << endl; }
      MPI_Finalize();
      exit(3);
   }
   PWConstCoefficient lambda_func(appState.param.mat.lambda);
   appState.lambda0_gf.SetSpace(appState.L2FESpace.get());
   appState.lambda0_gf.ProjectCoefficient(lambda_func);

   PWConstCoefficient mu_func(appState.param.mat.mu);
   appState.mu0_gf.SetSpace(appState.L2FESpace.get());
   appState.mu0_gf.ProjectCoefficient(mu_func);

   Vector mat(appState.pmesh->attributes.Max());
   for (int i = 0; i < mat.Size(); i++)
      mat[i] = i;
   PWConstCoefficient mat_func(mat);
   appState.mat_gf.SetSpace(appState.L2FESpace.get());
   appState.mat_gf.ProjectCoefficient(mat_func);

   appState.L2FESpace_mat.reset(new ParFiniteElementSpace(appState.pmesh.get(), appState.L2FEC.get(), num_materials));
   appState.comp_gf.SetSpace(appState.L2FESpace_mat.get());
   appState.comp_ref_gf.SetSpace(appState.L2FESpace_mat.get());
   CompoCoefficient comp_coeff(num_materials, appState.mat_gf);
   appState.comp_gf.ProjectCoefficient(comp_coeff);
   appState.comp_ref_gf = appState.comp_gf;

   appState.s_gf=0.0;
   if( appState.param.control.gravity > 0.0 )
   {
      if(appState.param.control.lithostatic)
      {
         LithostaticCoefficient Lithostatic_coeff(dim, xyz_gf_l2, appState.rho0_gf, appState.param.control.gravity, appState.param.control.thickness);
         appState.s_gf.ProjectCoefficient(Lithostatic_coeff);
      }
      else if(appState.param.control.atmospheric){
         ATMCoefficient ATM_coeff(dim, xyz_gf_l2, appState.rho0_gf, appState.param.control.gravity, appState.param.control.thickness);
         appState.s_gf.ProjectCoefficient(ATM_coeff);
      }
   }
   appState.s_gf.SyncAliasMemory(appState.S);
   appState.s_old_gf.SetSpace(appState.L2FESpace_stress.get());
   appState.s_old_gf = appState.s_gf;

   appState.x_ini_gf.SetSpace(appState.H1FESpace.get());
   appState.x_old_gf.SetSpace(appState.H1FESpace.get());
   appState.x_ini_gf = appState.x_gf;
   appState.x_old_gf = 0.0;

   appState.p_gf.SetSpace(appState.L2FESpace.get());
   appState.p_gf_old.SetSpace(appState.L2FESpace.get());
   appState.p_gf = 0.0; appState.p_gf_old = 0.0;
   Vector weak_location(dim);
   if(dim == 2){weak_location[0] = appState.param.mat.weak_x; weak_location[1] = appState.param.mat.weak_y;}
   else if(dim ==3){weak_location[0] = appState.param.mat.weak_x; weak_location[1] = appState.param.mat.weak_y; weak_location[2] = appState.param.mat.weak_z;}
   PlasticCoefficient p_coeff(dim, xyz_gf_l2, weak_location, appState.param.mat.weak_rad, appState.param.mat.ini_pls);
   ParGridFunction l2_p_gf(&L2FESpace_positive);
   l2_p_gf.ProjectCoefficient(p_coeff);
   appState.p_gf.ProjectGridFunction(l2_p_gf);
   appState.p_gf_old = appState.p_gf;
   appState.ini_p_gf.SetSpace(appState.L2FESpace.get());
   appState.ini_p_old_gf.SetSpace(appState.L2FESpace.get());
   appState.n_p_gf.SetSpace(appState.L2FESpace.get());
   appState.ini_p_gf = appState.p_gf; appState.ini_p_old_gf = appState.p_gf;
   appState.n_p_gf = 0.0;

   appState.u_gf.SetSpace(appState.H1FESpace.get());
   appState.u_gf = 0.0;

   int source = 0; bool visc = false, vorticity = false;
   if (appState.param.solver.impose_visc) { visc = true; }

   appState.geo.reset(new geodynamics::LagrangianGeoOperator(appState.S.Size(),
                                       *appState.H1FESpace, *appState.L2FESpace, *appState.L2FESpace_stress, appState.ess_tdofs,
                                       appState.rho0_gf, appState.fictitious_rho0_gf,
                                       appState.mat_gf, source,
                                       visc, vorticity,
                                       appState.lambda0_gf, appState.mu0_gf,
                                       appState.param, max_vbc_val));

   char vishost[] = "localhost";
   int visport = 19916;

   appState.energy_init = appState.geo->InternalEnergy(appState.e_gf) +
                        appState.geo->KineticEnergy(appState.v_gf);

   if (appState.param.sim.visualization)
   {
      MPI_Barrier(appState.pmesh->GetComm());
      appState.vis_rho.precision(8);
      appState.vis_v.precision(8);
      appState.vis_e.precision(8);
      int Wx = 0, Wy = 0;
      const int Ww = 350, Wh = 350;
      int offx = Ww+10;
      if (appState.param.sim.problem != 0 && appState.param.sim.problem != 4)
      {
         geodynamics::VisualizeField(appState.vis_rho, vishost, visport, appState.rho0_gf,
                                       "Density", Wx, Wy, Ww, Wh);
      }
      Wx += offx;
      geodynamics::VisualizeField(appState.vis_v, vishost, visport, appState.v_gf,
                                    "Velocity", Wx, Wy, Ww, Wh);
      Wx += offx;
      geodynamics::VisualizeField(appState.vis_e, vishost, visport, appState.e_gf,
                                    "Specific Internal Energy", Wx, Wy, Ww, Wh);
   }

   if (appState.param.sim.visit)
   {
      appState.visit_dc.reset(new VisItDataCollection(appState.param.sim.basename, appState.pmesh.get()));
      appState.visit_dc->RegisterField("Density",  &appState.rho0_gf);
      appState.visit_dc->RegisterField("Displacement", &appState.u_gf);
      appState.visit_dc->RegisterField("Velocity", &appState.v_gf);
      appState.visit_dc->RegisterField("Specific Internal Energy", &appState.e_gf);
      appState.visit_dc->RegisterField("Stress", &appState.s_gf);
      appState.visit_dc->RegisterField("Plastic Strain", &appState.p_gf);
      appState.visit_dc->RegisterField("Non-inital Plastic Strain", &appState.n_p_gf);
      appState.visit_dc->RegisterField("Composition", &appState.comp_gf);
      appState.visit_dc->RegisterField("Lambda", &appState.lambda0_gf);
      appState.visit_dc->RegisterField("Mu", &appState.mu0_gf);
      appState.visit_dc->SetCycle(0);
      appState.visit_dc->SetTime(0.0);
      appState.visit_dc->Save();
   }

   if (appState.param.sim.paraview)
   {
      appState.pd = new ParaViewDataCollection(appState.param.sim.basename, appState.pmesh.get());
      appState.pd->RegisterField("Density",  &appState.rho0_gf);
      appState.pd->RegisterField("Displacement", &appState.u_gf);
      appState.pd->RegisterField("Velocity", &appState.v_gf);
      appState.pd->RegisterField("Specific Internal Energy", &appState.e_gf);
      appState.pd->RegisterField("Stress", &appState.s_gf);
      appState.pd->RegisterField("Plastic Strain", &appState.p_gf);
      appState.pd->RegisterField("inital Plastic Strain", &appState.ini_p_gf);
      appState.pd->RegisterField("Non-inital Plastic Strain", &appState.n_p_gf);
      // Geometric parameters would need to be added to AppState if used
      // appState.pd->RegisterField("Geometric Parameters", &quality);
      appState.pd->RegisterField("Composition", &appState.comp_gf);
      appState.pd->RegisterField("Lambda", &appState.lambda0_gf);
      appState.pd->RegisterField("Mu", &appState.mu0_gf);
      appState.pd->SetLevelsOfDetail(appState.param.mesh.order_v);
      appState.pd->SetDataFormat(VTKFormat::BINARY);
      appState.pd->SetHighOrderOutput(true);
      appState.pd->SetCycle(0);
      appState.pd->SetTime(0.0);
      appState.pd->Save();
   }

   // Initialize the time integrator.
   appState.ode_solver->Init(*appState.geo);
   appState.geo->ResetTimeStepEstimate();
   appState.dt = appState.geo->GetTimeStepEstimate(appState.S);

   // Initialize submesh ODE solvers
   appState.ode_solver_sub->Init(*appState.oper_sub);
   appState.ode_solver_sub2->Init(*appState.oper_sub2);

   if (Mpi::Root())
   {
      std::cout<<""<<std::endl;
      std::cout<<"simulation starts"<<std::endl;
   }
}


// =============================================================================
// run() function
// =============================================================================
void run(AppState& appState)
{
   int myid = Mpi::WorldRank();
   bool mesh_changed = false;
   BlockVector S_old(appState.offset, Device::GetMemoryType());

   double h_min = appState.geo->GetLengthEstimate(appState.S);
   double t_old, dt_old = 0.0;

   for (int ti = 1; !appState.last_step; ti++)
   {
      if (appState.t + appState.dt >= appState.param.sim.t_final)
      {
         appState.dt = appState.param.sim.t_final - appState.t;
         appState.last_step = true;
      }
      if (appState.steps == appState.param.sim.max_tsteps) { appState.last_step = true; }

      S_old = appState.S;
      t_old = appState.t;
      double year = appState.t/86400/365.25;
      appState.p_gf_old = appState.p_gf; appState.ini_p_old_gf = appState.ini_p_gf; appState.x_old_gf = appState.x_gf;
      appState.geo->ResetTimeStepEstimate();

      if(appState.param.control.pseudo_transient)
      {
         for (int i = 0; i < appState.param.control.transient_num; i++)
         {
            appState.x_gf = appState.x_old_gf;
            appState.s_gf = appState.s_old_gf;
            appState.ode_solver->Step(appState.S, appState.t, appState.dt);
         }
         appState.t = appState.t - appState.dt*(appState.param.control.transient_num-1.0);
      }
      else
      {
         appState.ode_solver->Step(appState.S, appState.t, appState.dt);
      }

      if(appState.param.bc.surf_proc)
      {
         ParSubMesh::Transfer(appState.x_gf, appState.x_top);
         for (int i = 0; i < appState.topo.Size(); i++){appState.topo[i] = appState.x_top[i+appState.topo.Size()];}
         appState.topo.GetTrueDofs(appState.topo_t);
         appState.topo_t_old = appState.topo_t;
         appState.ode_solver_sub->Step(appState.topo_t, appState.t, appState.dt); appState.t = appState.t-appState.dt;
         appState.topo.SetFromTrueDofs(appState.topo_t);
         for (int i = 0; i < appState.topo.Size(); i++){appState.x_top[i+appState.topo.Size()] = appState.topo[i];}
         appState.submesh->NewNodes(appState.x_top, false);
         ParSubMesh::Transfer(appState.x_top, appState.x_gf);
      }

      if(appState.param.bc.winkler_foundation & appState.param.bc.base_proc)
      {
         ParSubMesh::Transfer(appState.x_gf, appState.x_bottom);
         for (int i = 0; i < appState.bottom.Size(); i++){appState.bottom[i] = appState.x_bottom[i+appState.bottom.Size()];}
         appState.bottom.GetTrueDofs(appState.bottom_t);
         appState.bottom_t_old = appState.bottom_t;
         appState.ode_solver_sub2->Step(appState.bottom_t, appState.t, appState.dt); appState.t=appState.t-appState.dt;
         appState.bottom.SetFromTrueDofs(appState.bottom_t);
         if(appState.param.bc.winkler_flat)
         {
            for (int i = 0; i < appState.bottom.Size(); i++){appState.x_bottom[i+appState.bottom.Size()] = 0.0;}
         }
         else
         {
            for (int i = 0; i < appState.bottom.Size(); i++){appState.x_bottom[i+appState.bottom.Size()] = appState.bottom[i];}
         }
         appState.submesh_bottom->NewNodes(appState.x_bottom, false);
         ParSubMesh::Transfer(appState.x_bottom, appState.x_gf);
      }

      if(appState.param.mat.plastic)
      {
         Returnmapping(dim, appState.comp_gf, appState.s_gf, appState.s_old_gf, appState.p_gf, appState.mat_gf, appState.param, h_min, dt_old);
         appState.n_p_gf  = appState.ini_p_gf;
         appState.n_p_gf -= appState.p_gf;
         appState.n_p_gf.Neg();
      }

      appState.steps++;
      dt_old = appState.dt;

      double dt_est = appState.geo->GetTimeStepEstimate(appState.S);
      h_min = appState.geo->GetLengthEstimate(appState.S);

      // NOTE: Remeshing logic is very complex.
      // A full port of the original remeshing logic is required here.
      if (appState.param.tmop.tmop && (ti % appState.param.tmop.remesh_steps == 0))
      {
          if (Mpi::Root()) cout << "*** Remeshing logic placeholder. Full logic from original code is needed. ***" << endl;
          mesh_changed = true;
      }

      if(mesh_changed)
      {
        mesh_changed = false;
      }
      else
      {
         if (dt_est < appState.dt)
         {
            appState.dt  = dt_est;
            if (appState.dt < 1.0E-38)
            {
               if (appState.visit_dc) { appState.visit_dc->SetCycle(ti); appState.visit_dc->SetTime(appState.t); appState.visit_dc->Save(); }
               if (appState.pd) { appState.pd->SetCycle(ti); appState.pd->SetTime(appState.t); appState.pd->Save(); }
               MFEM_ABORT("Aborting. dt became smaller than 1.0e-38 s!");
            }
            appState.t = t_old;
            appState.S = S_old;
            appState.p_gf = appState.p_gf_old; appState.ini_p_gf = appState.ini_p_old_gf;
            appState.geo->ResetQuadratureData();
            if (appState.steps < appState.param.sim.max_tsteps) { appState.last_step = false; }
            ti--; continue;
         }
         else if (dt_est > 1.25 * appState.dt) { appState.dt *= 1.02; }
      }

      appState.x_gf.SyncAliasMemory(appState.S);
      appState.v_gf.SyncAliasMemory(appState.S);
      appState.e_gf.SyncAliasMemory(appState.S);
      appState.s_gf.SyncAliasMemory(appState.S);
      appState.s_old_gf = appState.s_gf;
      appState.u_gf.Add(appState.dt, appState.v_gf);

      if (appState.param.tmop.tmop)
      {
         if(appState.param.control.mass_bal && ti > 1) { appState.geo->TMOPUpdate(appState.S, true); }
         else { appState.geo->TMOPUpdate(appState.S, false); }
      }
      else
      {
         appState.geo->TMOPUpdate(appState.S, false);
      }

      // Print progress.
      const bool print_progress_now = (ti % appState.param.sim.vis_steps == 0) || appState.last_step;
      if (print_progress_now)
      {
         print_progress(ti, appState);
      }

      // Save data at every vis_steps or at the last time step.
      // save data when printing progress for now. save_steps and print_steps will be separated in the future.
      const bool save_data_now = print_progress_now;
      if (save_data_now)
      {
         if (Mpi::Root()) { cout << "Saving data ..." << endl; }
         if (appState.param.sim.paraview)
         {
            appState.pd->SetCycle(ti);
            appState.pd->SetTime(appState.t);
            appState.pd->Save();
         }
         if (appState.param.sim.visit)
         {
            appState.visit_dc->SetCycle(ti);
            appState.visit_dc->SetTime(appState.t);
            appState.visit_dc->Save();
         }
      }
   }
   if (Mpi::Root())
   {
      cout << "Simulation finished." << endl;
      cout << "Total time steps: " << appState.steps << endl;
      cout << "Final time: " << appState.t << endl;
      cout << "Final dt: " << appState.dt << endl;
   }
}

// =============================================================================
// finalize() function
// =============================================================================
void finalize(AppState& appState)
{
   long mem=0, mmax=0, msum=0;

   switch (appState.param.solver.ode_solver_type)
   {
      case 2: appState.steps *= 2; break;
      case 3: appState.steps *= 3; break;
      case 4: appState.steps *= 4; break;
      case 6: appState.steps *= 6; break;
      case 7: appState.steps *= 2;
   }

   appState.geo->PrintTimingData(Mpi::Root(), appState.steps, appState.param.sim.fom);

   if (appState.param.sim.mem_usage)
   {
      mem = GetMaxRssMB();
      MPI_Reduce(&mem, &mmax, 1, MPI_LONG, MPI_MAX, 0, appState.pmesh->GetComm());
      MPI_Reduce(&mem, &msum, 1, MPI_LONG, MPI_SUM, 0, appState.pmesh->GetComm());
   }

   const double energy_final = appState.geo->InternalEnergy(appState.e_gf) +
                               appState.geo->KineticEnergy(appState.v_gf);
   if (Mpi::Root())
   {
      cout << endl;
      cout << "Energy  diff: " << std::scientific << std::setprecision(2)
           << fabs(appState.energy_init - energy_final) << endl;
      if (appState.param.sim.mem_usage)
      {
         cout << "Maximum memory resident set size: "
              << mmax << "/" << msum << " MB" << endl;
      }
   }

   if (appState.param.sim.problem == 0 || appState.param.sim.problem == 4)
   {
      const double error_max = appState.v_gf.ComputeMaxError(*appState.v_coeff),
                   error_l1  = appState.v_gf.ComputeL1Error(*appState.v_coeff),
                   error_l2  = appState.v_gf.ComputeL2Error(*appState.v_coeff);
      if (Mpi::Root())
      {
         cout << "L_inf  error: " << error_max << endl
              << "L_1    error: " << error_l1 << endl
              << "L_2    error: " << error_l2 << endl;
      }
   }

   if (appState.param.sim.visualization)
   {
      appState.vis_v.close();
      appState.vis_e.close();
   }

   delete appState.pd;
}


int main(int argc, char *argv[])
{
   AppState appState;
   initialize(appState, argc, argv);
   run(appState);
   finalize(appState);
   return 0;
}

// ... The rest of the functions (TMOPUpdate, ConductionOperator, etc.) remain unchanged ...
void TMOPUpdate(BlockVector &S, BlockVector &S_old,
               Array<int> &offset,
               ParGridFunction &x_gf,
               ParGridFunction &v_gf,
               ParGridFunction &e_gf,
               ParGridFunction &s_gf,
               ParGridFunction &x_ini_gf,
               ParGridFunction &p_gf,
               ParGridFunction &n_p_gf,
               ParGridFunction &ini_p_gf,
               ParGridFunction &u_gf,
               ParGridFunction &rho0_gf,
               ParGridFunction &lambda0_gf,
               ParGridFunction &mu0_gf,
               ParGridFunction &mat_gf,
               // ParLinearForm &flattening,
               int dim, bool amr)
{
   ParFiniteElementSpace* H1FESpace = x_gf.ParFESpace();
   ParFiniteElementSpace* L2FESpace = e_gf.ParFESpace();
   ParFiniteElementSpace* L2FESpace_stress = s_gf.ParFESpace();

   H1FESpace->Update();
   L2FESpace->Update();
   L2FESpace_stress->Update();

   int Vsize_h1 = H1FESpace->GetVSize();
   int Vsize_l2 = L2FESpace->GetVSize();

   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   offset[2] = offset[1] + Vsize_h1;
   offset[3] = offset[2] + Vsize_l2;
   offset[4] = offset[3] + Vsize_l2*3*(dim-1);
   // offset[5] = offset[4] + Vsize_h1;

   S_old = S;
   S.Update(offset);

   x_gf.Update();
   v_gf.Update();
   e_gf.Update();
   s_gf.Update();
   // x_ini_gf.Update();

   if(amr)
   {
      const Operator* H1Update = H1FESpace->GetUpdateOperator();
      const Operator* L2Update = L2FESpace->GetUpdateOperator();
      const Operator* L2Update_stress = L2FESpace_stress->GetUpdateOperator();

      H1Update->Mult(S_old.GetBlock(0), S.GetBlock(0));
      H1Update->Mult(S_old.GetBlock(1), S.GetBlock(1));
      L2Update->Mult(S_old.GetBlock(2), S.GetBlock(2));
      L2Update_stress->Mult(S_old.GetBlock(3), S.GetBlock(3));
      H1Update->Mult(S_old.GetBlock(4), S.GetBlock(4));
   }

   x_gf.MakeRef(H1FESpace, S, offset[0]);
   v_gf.MakeRef(H1FESpace, S, offset[1]);
   e_gf.MakeRef(L2FESpace, S, offset[2]);
   s_gf.MakeRef(L2FESpace_stress, S, offset[3]);
   // x_ini_gf.MakeRef(H1FESpace, S, offset[4]);
   S_old.Update(offset);

   // Gridfunction update (Non-blcok vector )
   p_gf.Update();
   n_p_gf.Update();
   ini_p_gf.Update();
   u_gf.Update();
   rho0_gf.Update();
   lambda0_gf.Update();
   mu0_gf.Update();
   mat_gf.Update();

   //
   // flattening.Update();
   // flattening.Assemble();

   H1FESpace->UpdatesFinished();
   L2FESpace->UpdatesFinished();
   L2FESpace_stress->UpdatesFinished();
}

ConductionOperator::ConductionOperator(ParFiniteElementSpace &f, double al,
                                       double kap, const Vector &u)
   : TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), M(NULL), K(NULL),
     T(NULL), current_dt(0.0),
     M_solver(f.GetComm()), T_solver(f.GetComm()), z(height)
{
   const double rel_tol = 1e-8;

   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator());
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
   M_prec.SetType(HypreSmoother::Jacobi);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(Mmat);

   alpha = al;
   kappa = kap;

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_solver.SetPreconditioner(T_prec);

   SetParameters(u);
}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
   // Compute:
   //    du_dt = M^{-1}*-Ku
   // for du_dt, where K is linearized by using u from the previous timestep
   Kmat.Mult(u, z);
   z.Neg(); // z = -z
   M_solver.Mult(z, du_dt);
}

void ConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
   // for du_dt, where K is linearized by using u from the previous timestep
   if (!T)
   {
      T = Add(1.0, Mmat, dt, Kmat);
      current_dt = dt;
      T_solver.SetOperator(*T);
   }
   MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
   Kmat.Mult(u, z);
   z.Neg();
   T_solver.Mult(z, du_dt);
}

void ConductionOperator::SetParameters(const Vector &u)
{
   ParGridFunction u_alpha_gf(&fespace);
   u_alpha_gf.SetFromTrueDofs(u);
   for (int i = 0; i < u_alpha_gf.Size(); i++)
   {
      u_alpha_gf(i) = kappa + alpha*u_alpha_gf(i);
   }

   delete K;
   K = new ParBilinearForm(&fespace);

   GridFunctionCoefficient u_coeff(&u_alpha_gf);

   K->AddDomainIntegrator(new DiffusionIntegrator(u_coeff));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);
   delete T;
   T = NULL; // re-compute T on the next ImplicitSolve
}

ConductionOperator::~ConductionOperator()
{
   delete T;
   delete M;
   delete K;
}

static void display_banner(std::ostream &os)
{
   os << endl
      << "       __                __               __    " << endl
      << "      / /   ____ _____ _/ /_  ____  _____/ /_   " << endl
      << "     / /   / __ `/ __ `/ __ \\/ __ \\/ ___/ __/ " << endl
      << "    / /___/ /_/ / /_/ / / / / /_/ (__  ) /_     " << endl
      << "   /_____/\\__,_/\\__, /_/ /_/\\____/____/\\__/ " << endl
      << "               /____/                           " << endl << endl;
}

static long GetMaxRssMB()
{
   struct rusage usage;
   if (getrusage(RUSAGE_SELF, &usage)) { return -1; }
#ifndef __APPLE__
   const long unit = 1024; // kilo
#else
   const long unit = 1024*1024; // mega
#endif
   return usage.ru_maxrss/unit; // mega bytes
}

static void Checks(const int ti, const double nrm, int &chk)
{
   const double eps = 1.e-13;
   //printf("\033[33m%.15e\033[m\n",nrm);

   auto check = [&](int p, int i, const double res)
   {
      auto rerr = [](const double a, const double v, const double eps)
      {
         MFEM_VERIFY(fabs(a) > eps && fabs(v) > eps, "One value is near zero!");
         const double err_a = fabs((a-v)/a);
         const double err_v = fabs((a-v)/v);
         return fmax(err_a, err_v) < eps;
      };
      if (problem == p && ti == i)
      { chk++; MFEM_VERIFY(rerr(nrm, res, eps), "P"<<problem<<", #"<<i); }
   };

   const double it_norms[2][8][2][2] = // dim, problem, {it,norm}
   {
      {
         {{5,  6.546538624534384e+00}, { 27, 7.588576357792927e+00}},
         {{5, 3.508254945225794e+00}, { 15, 2.756444596823211e+00}},
         {{5, 1.020745795651244e+01}, { 59, 1.721590205901898e+01}},
         {{5, 8.000000000000000e+00}, { 16, 8.000000000000000e+00}},
         {{5, 3.446324942352448e+01}, { 18, 3.446844033767240e+01}},
         {{5, 1.030899557252528e+01}, { 36, 1.057362418574309e+01}},
         {{5, 8.039707010835693e+00}, { 36, 8.316970976817373e+00}},
         {{5, 1.514929259650760e+01}, { 25, 1.514931278155159e+01}},
      },
      {
         {{5, 1.198510951452527e+03}, {188, 1.199384410059154e+03}},
         {{5, 1.339163718592566e+01}, { 28, 7.521073677397994e+00}},
         {{5, 2.041491591302486e+01}, {  59, 3.443180411803796e+01}},
         {{5, 1.600000000000000e+01}, { 16, 1.600000000000000e+01}},
         {{5, 6.892649884704898e+01}, { 18, 6.893688067534482e+01}},
         {{5, 2.061984481890964e+01}, { 36, 2.114519664792607e+01}},
         {{5, 1.607988713996459e+01}, { 36, 1.662736010353023e+01}},
         {{5, 3.029858112572883e+01}, { 24, 3.029858832743707e+01}}
      }
   };

   for (int p=0; p<8; p++)
   {
      for (int i=0; i<2; i++)
      {
         const int it = it_norms[dim-2][p][i][0];
         const double norm = it_norms[dim-2][p][i][1];
         check(p, it, norm);
      }
   }
}

void print_progress(const int ti, AppState &appState)
{
   // Get global maximum velocity
   Vector vel_mag(appState.v_gf.Size()/dim);
   int n = appState.v_gf.Size() / dim;
   for (int i = 0; i < n; i++)
   {
      double vx = appState.v_gf(i);
      double vy = appState.v_gf(i + n);
      double vz = (dim == 3) ? appState.v_gf(i + 2 * n) : 0.0;
      vel_mag[i] = std::sqrt(vx * vx + vy * vy + vz * vz);
   }
   double local_max_vel = vel_mag.Max();
   double global_max_vel;
   MPI_Reduce(&local_max_vel, &global_max_vel, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

   if (Mpi::Root())
   {
      double t = appState.t;
      double dt = appState.dt;
      std::string dt_unit = "sec";


      if (appState.param.sim.year)
      {
         const double secs_per_year = 86400.0 * 365.25;
         t /= secs_per_year;
         dt /= secs_per_year;
         dt_unit = "yr";
      }

      cout << std::fixed << std::setprecision(6) << std::scientific;
      cout << "step " << std::setw(5) << ti
           << ", t (" << dt_unit << ") = " << std::setw(5) << t
           << ", dt (" << dt_unit << ") = " << std::setw(5) << dt
           << ", max vel = " << std::setw(5) << global_max_vel;
      cout << endl;
   }
}