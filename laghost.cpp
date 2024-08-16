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
#include "laghost_solver.hpp"
#include "laghost_rheology.hpp"
#include "laghost_function.hpp"
#include "laghost_parameters.hpp"
#include "laghost_input.hpp"
#include "laghost_tmop.hpp"
#include "laghost_remhos.hpp"

using std::cout;
using std::endl;
using namespace mfem;

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

int main(int argc, char *argv[])
{
   // Initialize MPI.
   mfem::MPI_Session mpi(argc, argv);
   const int myid = mpi.WorldRank();

   // Print the banner.
   if (mpi.Root()) { display_banner(cout); }

   // Take care of input parameters given in a file or command line.
   OptionsParser args(argc, argv);
   Param param;
   read_and_assign_input_parameters( args, param, myid );

   Array<int> cxyz; // Leave undefined. It won't be used.
   double init_dt = 1.0;
   double blast_energy = 0.0;
   // double v_unit = 1.0/86400/365.25;
   double v_unit = 1.0;
   // double blast_energy = 1.0e-6;
   double blast_position[] = {0.0, 0.0, 0.0};
   // double blast_energy2 = 0.1;
   // double blast_position2[] = {8.0, 0.5};
   // Mesh bounding box
   Vector bb_min, bb_max;
   double bb_center[] = {0.0, 0.0, 0.0};
   double bb_length[] = {0.0, 0.0, 0.0};

   Vector bb_min2, bb_max2;
   double bb_center2[] = {0.0, 0.0, 0.0};
   double bb_length2[] = {0.0, 0.0, 0.0};
   double stretching_factor[] = {0.0, 0.0, 0.0};
   int num_materials = 1;

   // Remeshing performs using the Target-Matrix Optimization Paradigm (TMOP)
   bool mesh_changed = false;
   bool mesh_control_side = false;
   // int  multi_comp   = 0;

   param.tmop.mesh_poly_deg = param.mesh.order_v;
   param.tmop.quad_order = 2*param.mesh.order_v - 1; // integration order = 2p  - 1
   
   if(param.sim.max_tsteps > -1)
   {
      param.sim.t_final = 1.0e38;
   }   
   if(param.sim.year)
   {
      param.sim.t_final = param.sim.t_final * 86400 * 365.25;
      v_unit = 1.0/86400/365.25;
      
      if (mpi.Root())
      {
         std::cout << "Use years in output instead of seconds is true" << std::endl;
      }
   }
   else
   {
      // init_dt = 1e-1;
      if (mpi.Root())
      {
         std::cout << "Use seconds in output instead of years is true" << std::endl;
      }
      
   }

   // Configure the device from the command line options
   Device backend;
   backend.Configure(param.sim.device, param.sim.dev);
   if (mpi.Root()) { backend.Print(); }
   backend.SetGPUAwareMPI(param.sim.gpu_aware_mpi);

   // On all processors, use the default builtin 1D/2D/3D mesh or read the
   // serial one given on the command line.
   Mesh *mesh;

   if (param.mesh.mesh_file.compare("default") != 0)
   {
      mesh = new Mesh(param.mesh.mesh_file.c_str(), true, true);
   }
   else
   {
      if (param.sim.dim == 1)
      {
         mesh = new Mesh(Mesh::MakeCartesian1D(2));
         mesh->GetBdrElement(0)->SetAttribute(1);
         mesh->GetBdrElement(1)->SetAttribute(1);
      }
      if (param.sim.dim == 2)
      {
         mesh = new Mesh(Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL,
                                               true));
         const int NBE = mesh->GetNBE();
         for (int b = 0; b < NBE; b++)
         {
            Element *bel = mesh->GetBdrElement(b);
            const int attr = (b < NBE/2) ? 2 : 1;
            std::cout<< NBE <<"," << b << "," << attr <<std::endl;
            bel->SetAttribute(attr);
         }
      }
      if (param.sim.dim == 3)
      {
         mesh = new Mesh(Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON,
                                               true));
         const int NBE = mesh->GetNBE();
         for (int b = 0; b < NBE; b++)
         {
            Element *bel = mesh->GetBdrElement(b);
            const int attr = (b < NBE/3) ? 3 : (b < 2*NBE/3) ? 1 : 2;
            bel->SetAttribute(attr);
         }
      }
   }
   dim = mesh->Dimension();

   // 1D vs partial assembly sanity check.
   if (param.solver.p_assembly && dim == 1)
   {
      param.solver.p_assembly = false;
      if (mpi.Root())
      {
         cout << "Laghos does not support PA in 1D. Switching to FA." << endl;
      }
   }

   // Refine the mesh in serial to increase the resolution.
   for (int lev = 0; lev < param.mesh.rs_levels; lev++) { mesh->UniformRefinement(); }
   // mesh->EnsureNCMesh(true);

   if(param.mesh.local_refinement)
   {
      // Local refiement
      mesh->EnsureNCMesh(true);

      Array<int> refs;
      // for (int i = 0; i < mesh->GetNE(); i++)
      // {
      //    if(mesh->GetAttribute(i) >= 2)
      //    {
      //       refs.Append(i);
      //    }
      // }

      // mesh->GeneralRefinement(refs, 1);
      // refs.DeleteAll();

      for (int i = 0; i < mesh->GetNE(); i++)
      {
         if(mesh->GetAttribute(i) >= 2)
         {
            refs.Append(i);
         }
      }

      mesh->GeneralRefinement(refs, 1);
      refs.DeleteAll();

      for (int i = 0; i < mesh->GetNE(); i++)
      {
         if(mesh->GetAttribute(i) >= 3)
         {
            refs.Append(i);
         }
      }

      mesh->GeneralRefinement(refs, 1);
      refs.DeleteAll();

      // for (int i = 0; i < mesh->GetNE(); i++)
      // {
      //    if(mesh->GetAttribute(i) == 3)
      //    {
      //       refs.Append(i);
      //    }
      // }

      // mesh->GeneralRefinement(refs, 1);
      // refs.DeleteAll();
      
      mesh->Finalize(true);
   }

   const int mesh_NE = mesh->GetNE();
   if (mpi.Root())
   {
      cout << "Number of zones in the serial mesh: " << mesh_NE << endl;
   }

   // mesh->GetBoundingBox(bb_min, bb_max, max(param.mesh.order_v, 1));

   // Parallel partitioning of the mesh.
   ParMesh *pmesh = nullptr;
   const int num_tasks = mpi.WorldSize(); int unit = 1;
   int *nxyz = new int[dim];
   switch (param.mesh.partition_type)
   {
      case 0:
         for (int d = 0; d < dim; d++) { nxyz[d] = unit; }
         break;
      case 11:
      case 111:
         unit = static_cast<int>(floor(pow(num_tasks, 1.0 / dim) + 1e-2));
         for (int d = 0; d < dim; d++) { nxyz[d] = unit; }
         break;
      case 21: // 2D
         unit = static_cast<int>(floor(pow(num_tasks / 2, 1.0 / 2) + 1e-2));
         nxyz[0] = 2 * unit; nxyz[1] = unit;
         break;
      case 31: // 2D
         unit = static_cast<int>(floor(pow(num_tasks / 3, 1.0 / 2) + 1e-2));
         nxyz[0] = 3 * unit; nxyz[1] = unit;
         break;
      case 32: // 2D
         unit = static_cast<int>(floor(pow(2 * num_tasks / 3, 1.0 / 2) + 1e-2));
         nxyz[0] = 3 * unit / 2; nxyz[1] = unit;
         break;
      case 49: // 2D
         unit = static_cast<int>(floor(pow(9 * num_tasks / 4, 1.0 / 2) + 1e-2));
         nxyz[0] = 4 * unit / 9; nxyz[1] = unit;
         break;
      case 51: // 2D
         unit = static_cast<int>(floor(pow(num_tasks / 5, 1.0 / 2) + 1e-2));
         nxyz[0] = 5 * unit; nxyz[1] = unit;
         break;
      case 211: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 2, 1.0 / 3) + 1e-2));
         nxyz[0] = 2 * unit; nxyz[1] = unit; nxyz[2] = unit;
         break;
      case 221: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 4, 1.0 / 3) + 1e-2));
         nxyz[0] = 2 * unit; nxyz[1] = 2 * unit; nxyz[2] = unit;
         break;
      case 311: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 3, 1.0 / 3) + 1e-2));
         nxyz[0] = 3 * unit; nxyz[1] = unit; nxyz[2] = unit;
         break;
      case 321: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 6, 1.0 / 3) + 1e-2));
         nxyz[0] = 3 * unit; nxyz[1] = 2 * unit; nxyz[2] = unit;
         break;
      case 322: // 3D.
         unit = static_cast<int>(floor(pow(2 * num_tasks / 3, 1.0 / 3) + 1e-2));
         nxyz[0] = 3 * unit / 2; nxyz[1] = unit; nxyz[2] = unit;
         break;
      case 432: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 3, 1.0 / 3) + 1e-2));
         nxyz[0] = 2 * unit; nxyz[1] = 3 * unit / 2; nxyz[2] = unit;
         break;
      case 511: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 5, 1.0 / 3) + 1e-2));
         nxyz[0] = 5 * unit; nxyz[1] = unit; nxyz[2] = unit;
         break;
      case 521: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 10, 1.0 / 3) + 1e-2));
         nxyz[0] = 5 * unit; nxyz[1] = 2 * unit; nxyz[2] = unit;
         break;
      case 522: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 20, 1.0 / 3) + 1e-2));
         nxyz[0] = 5 * unit; nxyz[1] = 2 * unit; nxyz[2] = 2 * unit;
         break;
      case 911: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 9, 1.0 / 3) + 1e-2));
         nxyz[0] = 9 * unit; nxyz[1] = unit; nxyz[2] = unit;
         break;
      case 921: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 18, 1.0 / 3) + 1e-2));
         nxyz[0] = 9 * unit; nxyz[1] = 2 * unit; nxyz[2] = unit;
         break;
      case 922: // 3D.
         unit = static_cast<int>(floor(pow(num_tasks / 36, 1.0 / 3) + 1e-2));
         nxyz[0] = 9 * unit; nxyz[1] = 2 * unit; nxyz[2] = 2 * unit;
         break;
      default:
         if (myid == 0)
         {
            cout << "Unknown partition type: " << param.mesh.partition_type << '\n';
         }
         delete mesh;
         MPI_Finalize();
         return 3;
   }
   int product = 1;
   for (int d = 0; d < dim; d++) { product *= nxyz[d]; }
   const bool cartesian_partitioning = (cxyz.Size()>0)?true:false;
   if (product == num_tasks || cartesian_partitioning)
   {
      if (cartesian_partitioning)
      {
         int cproduct = 1;
         for (int d = 0; d < dim; d++) { cproduct *= cxyz[d]; }
         MFEM_VERIFY(!cartesian_partitioning || cxyz.Size() == dim,
                     "Expected " << mesh->SpaceDimension() << " integers with the "
                     "option --cartesian-partitioning.");
         MFEM_VERIFY(!cartesian_partitioning || num_tasks == cproduct,
                     "Expected cartesian partitioning product to match number of ranks.");
      }
      int *partitioning = cartesian_partitioning ?
                          mesh->CartesianPartitioning(cxyz):
                          mesh->CartesianPartitioning(nxyz);
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);
      delete [] partitioning;
   }
   else
   {
      if (myid == 0)
      {
         cout << "Non-Cartesian partitioning through METIS will be used.\n";
#ifndef MFEM_USE_METIS
         cout << "MFEM was built without METIS. "
              << "Adjust the number of tasks to use a Cartesian split." << endl;
#endif
      }
#ifndef MFEM_USE_METIS
      return 1;
#endif
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   }
   delete [] nxyz;
   delete mesh;

   // Refine the mesh further in parallel to increase the resolution.
   for (int lev = 0; lev < param.mesh.rp_levels; lev++) { pmesh->UniformRefinement(); }

   // pmesh->Rebalance();

   int NE = pmesh->GetNE(), ne_min, ne_max;
   MPI_Reduce(&NE, &ne_min, 1, MPI_INT, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&NE, &ne_max, 1, MPI_INT, MPI_MAX, 0, pmesh->GetComm());
   if (myid == 0)
   { cout << "Zones min/max: " << ne_min << " " << ne_max << endl; }

   pmesh->GetBoundingBox(bb_min, bb_max, max(param.mesh.order_v, 1));

   if(dim == 2)
   {
      bb_center[0] = (bb_min[0]+bb_max[0])*0.5;
      bb_center[1] = (bb_min[1]+bb_max[1])*0.5;

      bb_length[0] = (bb_max[0]-bb_min[0]);
      bb_length[1] = (bb_max[1]-bb_min[1]);

      // if (myid == 0)
      // { 
      //    cout << "min " << bb_min[0] << " " << bb_min[1] << endl; 
      //    cout << "max " << bb_max[0] << " " << bb_max[1] << endl;
      //    cout << "center " << bb_center[0] << " " << bb_center[1] << endl;
      //    cout << "length " << bb_length[0] << " " << bb_length[1] << endl;
      // }
   }
   else
   {
      bb_center[0] = (bb_min[0]+bb_max[0])*0.5;
      bb_center[1] = (bb_min[1]+bb_max[1])*0.5;
      bb_center[2] = (bb_min[2]+bb_max[2])*0.5;

      bb_length[0] = (bb_max[0]-bb_min[0]);
      bb_length[1] = (bb_max[1]-bb_min[1]);
      bb_length[2] = (bb_max[2]-bb_min[2]);
   }
   
   // Define the parallel finite element spaces. We use:
   // - H1 (Gauss-Lobatto, continuous) for position and velocity.
   // - L2 (Bernstein, discontinuous) for specific internal energy.
   // L2_FECollection L2FEC(param.mesh.order_e, dim, BasisType::Positive); // Bernstein polynomials.
   // L2_FECollection L2FEC(param.mesh.order_e, dim, BasisType::GaussLegendre); // Open type.
   L2_FECollection L2FEC(param.mesh.order_e, dim, BasisType::GaussLobatto); // Closed type.

   L2_FECollection l2_fec(param.mesh.order_e, pmesh->Dimension(), BasisType::Positive); // Non-positive basis drives oscillation while interpolation wtihin DG type elements. 
   ParFiniteElementSpace l2_fes(pmesh, &l2_fec);

   // switch (param.mesh.l2_basis)
   // {
   //    case 1: L2_FECollection L2FEC(param.mesh.order_e, dim, BasisType::GaussLobatto); break;
   //    case 2: L2_FECollection L2FEC(param.mesh.order_e, dim, BasisType::GaussLegendre); break;
   //    case 3: L2_FECollection L2FEC(param.mesh.order_e, dim, BasisType::Positive); break;
   //    default:
   //       if (myid == 0)
   //       {
   //          cout << "Unknown l2 basis type: " << param.mesh.l2_basis << '\n';
   //       }
   //       delete pmesh;
   //       MPI_Finalize();
   //       return 3;
   // }

   H1_FECollection H1FEC(param.mesh.order_v, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, pmesh->Dimension());
   ParFiniteElementSpace L2FESpace_stress(pmesh, &L2FEC, 3*(dim-1)); // three varibles for 2D, six varibles for 3D

   // Boundary conditions: all tests use v.n = 0 on the boundary, and we assume
   // that the boundaries are straight.
   // Remove square brackets and spaces
   param.bc.bc_ids.erase(std::remove(param.bc.bc_ids.begin(), param.bc.bc_ids.end(), '['), param.bc.bc_ids.end());
   param.bc.bc_ids.erase(std::remove(param.bc.bc_ids.begin(), param.bc.bc_ids.end(), ']'), param.bc.bc_ids.end());
   param.bc.bc_ids.erase(std::remove(param.bc.bc_ids.begin(), param.bc.bc_ids.end(), ' '), param.bc.bc_ids.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss(param.bc.bc_ids);
   std::vector<int> bc_id;
    
   // Temporary variable to store each token
   std::string token;

   // std::cout <<"check point1"<<std::endl;
   // Tokenize the string and convert tokens to integers
   while (getline(ss, token, ',')) 
   {bc_id.push_back(std::stoi(token)); // Convert string to int and add to vector
   }

   // std::cout <<"check point2"<<std::endl;

   if(pmesh->bdr_attributes.Max() != bc_id.size())
   {
      if (myid == 0){cout << "The number of boundaries are not consistent with the given mesh. \nBC indicator from mesh is " << pmesh->bdr_attributes.Max() << " but input is " << bc_id.size() << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }

   // Boundary velocity of x component
   param.bc.bc_vxs.erase(std::remove(param.bc.bc_vxs.begin(), param.bc.bc_vxs.end(), '['), param.bc.bc_vxs.end());
   param.bc.bc_vxs.erase(std::remove(param.bc.bc_vxs.begin(), param.bc.bc_vxs.end(), ']'), param.bc.bc_vxs.end());
   param.bc.bc_vxs.erase(std::remove(param.bc.bc_vxs.begin(), param.bc.bc_vxs.end(), ' '), param.bc.bc_vxs.end());

   // Create a stringstream to tokenize the string
   std::stringstream vxs(param.bc.bc_vxs);
   std::vector<double> bc_vx;
    
   // Tokenize the string and convert tokens to integers
   while (getline(vxs, token, ',')) 
   {bc_vx.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // Boundary velocity of y component
   param.bc.bc_vys.erase(std::remove(param.bc.bc_vys.begin(), param.bc.bc_vys.end(), '['), param.bc.bc_vys.end());
   param.bc.bc_vys.erase(std::remove(param.bc.bc_vys.begin(), param.bc.bc_vys.end(), ']'), param.bc.bc_vys.end());
   param.bc.bc_vys.erase(std::remove(param.bc.bc_vys.begin(), param.bc.bc_vys.end(), ' '), param.bc.bc_vys.end());

   // Create a stringstream to tokenize the string
   std::stringstream vys(param.bc.bc_vys);
   std::vector<double> bc_vy;

   // Tokenize the string and convert tokens to integers
   while (getline(vys, token, ',')) 
   {bc_vy.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // Boundary velocity of z component
   param.bc.bc_vzs.erase(std::remove(param.bc.bc_vzs.begin(), param.bc.bc_vzs.end(), '['), param.bc.bc_vzs.end());
   param.bc.bc_vzs.erase(std::remove(param.bc.bc_vzs.begin(), param.bc.bc_vzs.end(), ']'), param.bc.bc_vzs.end());
   param.bc.bc_vzs.erase(std::remove(param.bc.bc_vzs.begin(), param.bc.bc_vzs.end(), ' '), param.bc.bc_vzs.end());

   // Create a stringstream to tokenize the string
   std::stringstream vzs(param.bc.bc_vzs);
   std::vector<double> bc_vz;

   // Tokenize the string and convert tokens to integers
   while (getline(vzs, token, ',')) 
   {bc_vz.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   if(param.bc.bc_unit =="cm/yr")
   {v_unit = v_unit/100.0;}
   else if(param.bc.bc_unit =="mm/yr")
   {v_unit = v_unit/1000.0;}
   else if(param.bc.bc_unit =="m/s")
   {v_unit = v_unit*1.0;}
   else if(param.bc.bc_unit =="cm/s")
   {v_unit = v_unit*0.01;}

   // Dirichlet type boundary condition (i.e., fixing velocity component at boundaries)
   Array<int> ess_tdofs, ess_vdofs;
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max()), dofs_marker, dofs_list;
      for (int i = 0; i < bc_id.size(); ++i) 
      {
        ess_bdr = 0;
        if(bc_id[i] > 0)
        {
            if(dim == 2)
            {
               switch (bc_id[i])
               {
                  // case 1 : x compoent is constained
                  // case 2 : y compoent is constained
                  // case 3 : all compoents are constained
                  case 1: ess_bdr[i] = 1; H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,0); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); break;
                  case 2: ess_bdr[i] = 1; H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,1); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); break;
                  case 3: ess_bdr[i] = 1; H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); break;
                  default:
                     if (myid == 0)
                     {
                        cout << "Unknown boundary type: " << bc_id[i] << '\n';
                     }
                     delete pmesh;
                     MPI_Finalize();
                     return 3;
               }
            }
            else
            {
               switch (bc_id[i])
               {
                  // case 1 : x compoent is constained
                  // case 2 : y compoent is constained
                  // case 3 : z compoent is constained
                  // case 4 : all compoents are constained
                  // case 5 : x and y compoents are constained
                  // case 6 : x and z compoents are constained
                  // case 7 : y and z compoents are constained
                  case 1: ess_bdr[i] = 1; H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,0); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); break;
                  case 2: ess_bdr[i] = 1; H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,1); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); break;
                  case 3: ess_bdr[i] = 1; H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,2); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,2); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); break;
                  case 4: ess_bdr[i] = 1; H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); break;
                  case 5:
                     ess_bdr[i] = 1; 
                     H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,0); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); 
                     H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,1); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); 
                     break;
                  case 6:
                     ess_bdr[i] = 1; 
                     H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,0); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); 
                     H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,2); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,2); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); 
                     break;
                  case 7: 
                     ess_bdr[i] = 1; 
                     H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,1); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); 
                     H1FESpace.GetEssentialTrueDofs(ess_bdr, dofs_list,2); ess_tdofs.Append(dofs_list); H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,2); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list); ess_vdofs.Append(dofs_list); 
                     break; 

                  default:
                     if (myid == 0)
                     {
                        cout << "Unknown boundary type: " << bc_id[i] << '\n';
                     }
                     delete pmesh;
                     MPI_Finalize();
                     return 3;
               }
            }
        }
      }
   }

   Vector bc_id_pa(pmesh->bdr_attributes.Max());
   for (int i = 0; i < bc_id.size(); ++i){bc_id_pa[i]=bc_id[i];}

   // Define the explicit ODE solver used for time integration.
   ODESolver *ode_solver = NULL;
   switch (param.solver.ode_solver_type)
   {
      case 1: ode_solver = new ForwardEulerSolver; break;
      case 2: ode_solver = new RK2Solver(0.5); break;
      case 3: ode_solver = new RK3SSPSolver; break;
      case 4: ode_solver = new RK4Solver; break;
      case 6: ode_solver = new RK6Solver; break;
      case 7: ode_solver = new RK2AvgSolver; break;
      default:
         if (myid == 0)
         {
            cout << "Unknown ODE solver type: " << param.solver.ode_solver_type << '\n';
         }
         delete pmesh;
         MPI_Finalize();
         return 3;
   }

   // 4. Define the ODE solver for submesh used for time integration. Several implicit
   //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
   //    explicit Runge-Kutta methods are available.
   int ode_solver_type = 12;
   ODESolver *ode_solver_sub;
   ODESolver *ode_solver_sub2;
   switch (ode_solver_type)
   {
      // Implicit L-stable methods
      case 1:  ode_solver_sub = new BackwardEulerSolver; break;
      case 2:  ode_solver_sub = new SDIRK23Solver(2); break;
      case 3:  ode_solver_sub = new SDIRK33Solver; break;
      // Explicit methods
      case 11: ode_solver_sub = new ForwardEulerSolver; break;
      case 12: ode_solver_sub = new RK2Solver(0.5); break; // midpoint method
      case 13: ode_solver_sub = new RK3SSPSolver; break;
      case 14: ode_solver_sub = new RK4Solver; break;
      case 15: ode_solver_sub = new GeneralizedAlphaSolver(0.5); break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver_sub = new ImplicitMidpointSolver; break;
      case 23: ode_solver_sub = new SDIRK23Solver; break;
      case 24: ode_solver_sub = new SDIRK34Solver; break;
      default:
         if (myid == 0)
         {
            cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         }
         delete pmesh;
         MPI_Finalize();
         return 3;
   }

   switch (ode_solver_type)
   {
      // Implicit L-stable methods
      case 1:  ode_solver_sub2 = new BackwardEulerSolver; break;
      case 2:  ode_solver_sub2 = new SDIRK23Solver(2); break;
      case 3:  ode_solver_sub2 = new SDIRK33Solver; break;
      // Explicit methods
      case 11: ode_solver_sub2 = new ForwardEulerSolver; break;
      case 12: ode_solver_sub2 = new RK2Solver(0.5); break; // midpoint method
      case 13: ode_solver_sub2 = new RK3SSPSolver; break;
      case 14: ode_solver_sub2 = new RK4Solver; break;
      case 15: ode_solver_sub2 = new GeneralizedAlphaSolver(0.5); break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver_sub2 = new ImplicitMidpointSolver; break;
      case 23: ode_solver_sub2 = new SDIRK23Solver; break;
      case 24: ode_solver_sub2 = new SDIRK34Solver; break;
      default:
         if (myid == 0)
         {
            cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         }
         delete pmesh;
         MPI_Finalize();
         return 3;
   }

   const HYPRE_Int glob_size_l2 = L2FESpace.GlobalTrueVSize();
   const HYPRE_Int glob_size_h1 = H1FESpace.GlobalTrueVSize();
   if (mpi.Root())
   {
      cout << "Number of kinematic (position, velocity) dofs: "
           << glob_size_h1 << endl;
      cout << "Number of specific internal energy dofs: "
           << glob_size_l2 << endl;
   }

   // The monolithic BlockVector stores unknown fields as:
   // - 0 -> position
   // - 1 -> velocity
   // - 2 -> specific internal energy
   // - 3 -> stress
   const int Vsize_l2 = L2FESpace.GetVSize();
   const int Vsize_h1 = H1FESpace.GetVSize();
   Array<int> offset(5); // when you change this number, you should chnage block offset in solver.cpp too
   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   offset[2] = offset[1] + Vsize_h1;
   offset[3] = offset[2] + Vsize_l2;
   offset[4] = offset[3] + Vsize_l2*3*(dim-1);
   // offset[5] = offset[4] + Vsize_h1;
   BlockVector S(offset, Device::GetMemoryType());

   // Define GridFunction objects for the position, velocity and specific
   // internal energy. There is no function for the density, as we can always
   // compute the density values given the current mesh position, using the
   // property of pointwise mass conservation.
   ParGridFunction x_gf, v_gf, e_gf, s_gf;
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   v_gf.MakeRef(&H1FESpace, S, offset[1]);
   e_gf.MakeRef(&L2FESpace, S, offset[2]);
   s_gf.MakeRef(&L2FESpace_stress, S, offset[3]);
   pmesh->SetNodalGridFunction(&x_gf);
   // Sync the data location of x_gf with its base, S
   x_gf.SyncAliasMemory(S);

   // Create a "sub mesh" from the boundary elements with attribute 3 (top boundary) for Surface process
   Array<int> bdr_attrs(1);
   bdr_attrs[0] = 4;
   ParSubMesh submesh(ParSubMesh::CreateFromBoundary(*pmesh, bdr_attrs));
   ParFiniteElementSpace sub_fespace0(&submesh, &H1FEC, pmesh->Dimension()); // nodes of submesh
   ParFiniteElementSpace sub_fespace1(&submesh, &H1FEC); // topography

   // Solve a Poisson problem on the boundary. This just follows ex0p.
   Array<int> boundary_dofs;
   sub_fespace1.GetBoundaryTrueDofs(boundary_dofs);
   ParGridFunction x_top(&sub_fespace0); 
   // ParGridFunction x_top_old(&sub_fespace0);
   ParGridFunction topo(&sub_fespace1);
   submesh.SetNodalGridFunction(&x_top);
   // submesh.SetNodalGridFunction(&x_top_old);
   for (int i = 0; i < topo.Size(); i++){topo[i] = x_top[i+topo.Size()];}

   Vector topo_t, topo_t_old;
   topo.GetTrueDofs(topo_t); topo_t_old=topo;

   // Create a "sub mesh" from the boundary elements with attribute 2 (bottom boundary) for flattening 
   Array<int> bdr_attrs_b(1);
   bdr_attrs_b[0] = 3;
   ParSubMesh submesh_bottom(ParSubMesh::CreateFromBoundary(*pmesh, bdr_attrs_b));
   ParFiniteElementSpace sub_fespace2(&submesh_bottom, &H1FEC, pmesh->Dimension()); // nodes of submesh
   ParFiniteElementSpace sub_fespace3(&submesh_bottom, &H1FEC); // topography
   
   // Solve a Poisson problem on the boundary. This just follows ex0p.
   Array<int> boundary_dofs_bot;
   sub_fespace2.GetBoundaryTrueDofs(boundary_dofs_bot);
   ParGridFunction x_bottom(&sub_fespace2);
   // ParGridFunction x_bottom_old(&sub_fespace2);
   ParGridFunction bottom(&sub_fespace3);
   submesh_bottom.SetNodalGridFunction(&x_bottom);
   // submesh_bottom.SetNodalGridFunction(&x_bottom_old);
   for (int i = 0; i < bottom.Size(); i++){bottom[i] = x_bottom[i+bottom.Size()];}

   Vector bottom_t, bottom_t_old;
   bottom.GetTrueDofs(bottom_t); bottom_t_old=bottom;


   // Create a "sub mesh" from the boundary elements with attribute 0
   Array<int> bdr_attrs_x0(1);
   bdr_attrs_x0[0] = 1;
   ParSubMesh submesh_x0(ParSubMesh::CreateFromBoundary(*pmesh, bdr_attrs_x0));
   ParFiniteElementSpace sub_fespace4(&submesh_x0, &H1FEC, pmesh->Dimension()); // right sidewall 
   ParGridFunction x0_side(&sub_fespace4);
   submesh_x0.SetNodalGridFunction(&x0_side);

   // Create a "sub mesh" from the boundary elements with attribute 0
   Array<int> bdr_attrs_x1(1);
   bdr_attrs_x1[0] = 2;
   ParSubMesh submesh_x1(ParSubMesh::CreateFromBoundary(*pmesh, bdr_attrs_x1));
   ParFiniteElementSpace sub_fespace5(&submesh_x1, &H1FEC, pmesh->Dimension()); // right sidewall 
   ParGridFunction x1_side(&sub_fespace5);
   submesh_x1.SetNodalGridFunction(&x1_side);

   if(dim == 3)
   {
      // Create a "sub mesh" from the boundary elements with attribute 0
      Array<int> bdr_attrs_y0(1);
      bdr_attrs_y0[0] = 5;
      ParSubMesh submesh_y0(ParSubMesh::CreateFromBoundary(*pmesh, bdr_attrs_y0));
      ParFiniteElementSpace sub_fespace4(&submesh_y0, &H1FEC, pmesh->Dimension()); // right sidewall 
      ParGridFunction y0_side(&sub_fespace4);
      submesh_y0.SetNodalGridFunction(&y0_side);

      // Create a "sub mesh" from the boundary elements with attribute 0
      Array<int> bdr_attrs_y1(1);
      bdr_attrs_y1[0] = 6;
      ParSubMesh submesh_y1(ParSubMesh::CreateFromBoundary(*pmesh, bdr_attrs_y1));
      ParFiniteElementSpace sub_fespace5(&submesh_y1, &H1FEC, pmesh->Dimension()); // right sidewall 
      ParGridFunction y1_side(&sub_fespace5);
      submesh_y1.SetNodalGridFunction(&y1_side);

   }

   // 
   // ParGridFunction x_bottom(&sub_fespace2);
   // ParGridFunction bottom(&sub_fespace3);
   // submesh_bottom.SetNodalGridFunction(&x_bottom);
   // for (int i = 0; i < bottom.Size(); i++){bottom[i] = x_bottom[i+bottom.Size()];}

   // 9. Initialize the conduction operator for surface diffusion
   ConductionOperator oper_sub(  sub_fespace1, param.bc.surf_alpha, param.bc.surf_diff, topo_t   );
   ConductionOperator oper_sub2( sub_fespace3, param.bc.base_alpha, param.bc.base_diff, bottom_t );

   // xyz coordinates in L2 space
   ParFiniteElementSpace L2FESpace_xyz(pmesh, &l2_fec, dim); //
   ParGridFunction xyz_gf_l2(&L2FESpace_xyz);
   VectorFunctionCoefficient xyz_coeff(pmesh->Dimension(), xyz0);
   xyz_gf_l2.ProjectCoefficient(xyz_coeff);

   int nSize = 1, nAspr = 1, nSkew = 1;
   if (dim == 3)
   {
      nAspr = 2;
      nSkew = 3;
   }

   // Total number of geometric parameters; for now we skip orientation.
   const int nTotalParams = nSize + nAspr + nSkew;
   // Define a GridFunction for all geometric parameters associated with the
   // mesh.
   ParFiniteElementSpace L2FESpace_geometric(pmesh, &L2FEC, nTotalParams); // must order byNodes
   ParGridFunction quality(&L2FESpace_geometric);
   // Vector quality; quality.SetSize(e_gf.Size()*nTotalParams);
   
   DenseMatrix jacobian(dim);
   Vector geomParams(nTotalParams);
   Array<int> vdofs;
   Vector allVals;
   // Compute the geometric parameter at the dofs of each element.
   for (int e = 0; e < pmesh->GetNE(); e++)
   {
      const FiniteElement *fe = L2FESpace_geometric.GetFE(e);
      const IntegrationRule &ir = fe->GetNodes();
      L2FESpace_geometric.GetElementVDofs(e, vdofs);
      allVals.SetSize(vdofs.Size());
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         pmesh->GetElementJacobian(e, jacobian, &ip);
         double sizeVal;
         Vector asprVals, skewVals, oriVals;
         pmesh->GetGeometricParametersFromJacobian(jacobian, sizeVal,
                                                  asprVals, skewVals, oriVals);
         allVals(q + 0) = sizeVal;
         for (int n = 0; n < nAspr; n++)
         {
            if(asprVals(n) > 1.0){allVals(q + (n+1)*ir.GetNPoints()) = asprVals(n);}
            else{allVals(q + (n+1)*ir.GetNPoints()) = 1/asprVals(n);}
            
         }
         for (int n = 0; n < nSkew; n++)
         {
            allVals(q + (n+1+nAspr)*ir.GetNPoints()) = skewVals(n);
         }
      }
      quality.SetSubVector(vdofs, allVals);
   }
   ParGridFunction vol_ini_gf(&L2FESpace);
   ParGridFunction skew_ini_gf(&L2FESpace);
   for (int i = 0; i < vol_ini_gf.Size(); i++){vol_ini_gf[i] = quality[i]; skew_ini_gf[i] = quality[i + e_gf.Size()];}
   

   // Initialize the velocity.
   v_gf = 0.0;
   // PlasticCoefficient p_coeff(dim, xyz_gf_l2, weak_location, param.mat.weak_rad, param.mat.ini_pls);
   VectorFunctionCoefficient v_coeff(pmesh->Dimension(), v0);
   v_gf.ProjectCoefficient(v_coeff);

   double max_vbc_val = param.control.max_vbc_val;
   for (int i = 0; i < bc_id.size(); ++i) 
   // for (int i = bc_id.size() -1; i > -1; --i) 
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max()), dofs_marker, dofs_list1, dofs_list2, dofs_list3;
      if(bc_id[i] > 0)
      {
         
         if(dim == 2)
            max_vbc_val = std::max(max_vbc_val,  sqrt(pow(v_unit*bc_vx[i], 2) + pow(v_unit*bc_vy[i], 2)));
         else
            max_vbc_val = std::max(max_vbc_val,  sqrt(pow(v_unit*bc_vx[i], 2) + pow(v_unit*bc_vy[i], 2) + pow(v_unit*bc_vz[i], 2)));

         ess_bdr = 0;
         if(dim == 2)
         {
            switch (bc_id[i])
            {
               // case 1 : x compoent is constained
               // case 2 : y compoent is constained
               // case 3 : all compoents are constained
               case 1: ess_bdr[i] = 1; H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list1); for (int j = 0; j < dofs_list1.Size(); j++){v_gf(dofs_list1[j]) = v_unit*bc_vx[i];} break;
               case 2: ess_bdr[i] = 1; H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list2); for (int j = 0; j < dofs_list2.Size(); j++){v_gf(dofs_list2[j]) = v_unit*bc_vy[i];} break;
               case 3: ess_bdr[i] = 1; 
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list1); for (int j = 0; j < dofs_list1.Size(); j++){v_gf(dofs_list1[j]) = v_unit*bc_vx[i];}
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list2); for (int j = 0; j < dofs_list2.Size(); j++){v_gf(dofs_list2[j]) = v_unit*bc_vy[i];} break;
               
               default:
                  if (myid == 0)
                  {
                        cout << "Unknown boundary type: " << bc_id[i] << '\n';
                  }
                  delete pmesh;
                  MPI_Finalize();
                  return 3;
            }
         }
         else
         {
            switch (bc_id[i])
            {
               // case 1 : x compoent is constained
               // case 2 : y compoent is constained
               // case 3 : z compoent is constained
               // case 4 : all compoents are constained
               // case 5 : x and y compoents are constained
               // case 6 : x and z compoents are constained
               // case 7 : y and z compoents are constained
               case 1: ess_bdr[i] = 1; H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list1); for (int j = 0; j < dofs_list1.Size(); j++){v_gf(dofs_list1[j]) = v_unit*bc_vx[i];} break;
               case 2: ess_bdr[i] = 1; H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list2); for (int j = 0; j < dofs_list2.Size(); j++){v_gf(dofs_list2[j]) = v_unit*bc_vy[i];} break;
               case 3: ess_bdr[i] = 1; H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,2); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list3); for (int j = 0; j < dofs_list3.Size(); j++){v_gf(dofs_list3[j]) = v_unit*bc_vz[i];} break;
               case 4: ess_bdr[i] = 1; 
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list1); for (int j = 0; j < dofs_list1.Size(); j++){v_gf(dofs_list1[j]) = v_unit*bc_vx[i];}
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list2); for (int j = 0; j < dofs_list2.Size(); j++){v_gf(dofs_list2[j]) = v_unit*bc_vy[i];}
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,2); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list3); for (int j = 0; j < dofs_list3.Size(); j++){v_gf(dofs_list3[j]) = v_unit*bc_vz[i];} break;
               case 5: ess_bdr[i] = 1; 
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list1); for (int j = 0; j < dofs_list1.Size(); j++){v_gf(dofs_list1[j]) = v_unit*bc_vx[i];}
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list2); for (int j = 0; j < dofs_list2.Size(); j++){v_gf(dofs_list2[j]) = v_unit*bc_vy[i];} break;
               case 6: ess_bdr[i] = 1;
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,0); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list1); for (int j = 0; j < dofs_list1.Size(); j++){v_gf(dofs_list1[j]) = v_unit*bc_vx[i];}
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,2); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list3); for (int j = 0; j < dofs_list3.Size(); j++){v_gf(dofs_list3[j]) = v_unit*bc_vz[i];} break;
               case 7: ess_bdr[i] = 1;
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,1); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list2); for (int j = 0; j < dofs_list2.Size(); j++){v_gf(dofs_list2[j]) = v_unit*bc_vy[i];}
                       H1FESpace.GetEssentialVDofs(ess_bdr, dofs_marker,2); FiniteElementSpace::MarkerToList(dofs_marker, dofs_list3); for (int j = 0; j < dofs_list3.Size(); j++){v_gf(dofs_list3[j]) = v_unit*bc_vz[i];} break;

               default:
                  if (myid == 0)
                  {
                     cout << "Unknown boundary type: " << bc_id[i] << '\n';
                  }
                  delete pmesh;
                  MPI_Finalize();
                  return 3;
            }
         }
      }
   }

   

   // for (int i = 0; i < ess_vdofs.Size(); i++)
   // {
   //    v_gf(ess_vdofs[i]) = 0.0;
   // }

   // For the Sedov test, we use a delta function at the origin.
   // Vector dir(pmesh->Dimension());
   // L2_FECollection h1_fec(order_v, pmesh->Dimension());
   // ParFiniteElementSpace h1_fes(pmesh, &h1_fec);
   // ParGridFunction h1_v(&h1_fes);
   // dir=0.0;
   // dir(1)=1.0;
   // // VectorDeltaCoefficient v2_coeff(pmesh->Dimension());
   // // v2_coeff.SetScale(blast_energy2);
   // // v2_coeff.SetDeltaCenter(cent);
   // // v2_coeff.SetDirection(dir);
   // VectorDeltaCoefficient v2_coeff(dir, blast_position2[0], blast_position2[1], blast_energy2);
  
   // h1_v.ProjectCoefficient(v2_coeff);
   // // v_gf.ProjectGridFunction(h1_v);

   // Sync the data location of v_gf with its base, S
   v_gf.SyncAliasMemory(S);

   // String (Material) extraction
   param.mat.rho.erase(std::remove(param.mat.rho.begin(), param.mat.rho.end(), '['), param.mat.rho.end());
   param.mat.rho.erase(std::remove(param.mat.rho.begin(), param.mat.rho.end(), ']'), param.mat.rho.end());
   param.mat.rho.erase(std::remove(param.mat.rho.begin(), param.mat.rho.end(), ' '), param.mat.rho.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss2(param.mat.rho);
   std::vector<double> rho_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss2, token, ',')) 
   {rho_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   param.mat.lambda.erase(std::remove(param.mat.lambda.begin(), param.mat.lambda.end(), '['), param.mat.lambda.end());
   param.mat.lambda.erase(std::remove(param.mat.lambda.begin(), param.mat.lambda.end(), ']'), param.mat.lambda.end());
   param.mat.lambda.erase(std::remove(param.mat.lambda.begin(), param.mat.lambda.end(), ' '), param.mat.lambda.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss3(param.mat.lambda);
   std::vector<double> lambda_vec;

   // Tokenize the string and convert tokens to integers
   while (getline(ss3, token, ',')) 
   {
      lambda_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.mu.erase(std::remove(param.mat.mu.begin(), param.mat.mu.end(), '['), param.mat.mu.end());
   param.mat.mu.erase(std::remove(param.mat.mu.begin(), param.mat.mu.end(), ']'), param.mat.mu.end());
   param.mat.mu.erase(std::remove(param.mat.mu.begin(), param.mat.mu.end(), ' '), param.mat.mu.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss4(param.mat.mu);
   std::vector<double> mu_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss4, token, ',')) 
   {mu_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.tension_cutoff.erase(std::remove(param.mat.tension_cutoff.begin(), param.mat.tension_cutoff.end(), '['), param.mat.tension_cutoff.end());
   param.mat.tension_cutoff.erase(std::remove(param.mat.tension_cutoff.begin(), param.mat.tension_cutoff.end(), ']'), param.mat.tension_cutoff.end());
   param.mat.tension_cutoff.erase(std::remove(param.mat.tension_cutoff.begin(), param.mat.tension_cutoff.end(), ' '), param.mat.tension_cutoff.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss5(param.mat.tension_cutoff);
   std::vector<double> tension_cutoff_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss5, token, ',')) 
   {tension_cutoff_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }
   
   // String (Material) extraction
   param.mat.cohesion0.erase(std::remove(param.mat.cohesion0.begin(), param.mat.cohesion0.end(), '['), param.mat.cohesion0.end());
   param.mat.cohesion0.erase(std::remove(param.mat.cohesion0.begin(), param.mat.cohesion0.end(), ']'), param.mat.cohesion0.end());
   param.mat.cohesion0.erase(std::remove(param.mat.cohesion0.begin(), param.mat.cohesion0.end(), ' '), param.mat.cohesion0.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss6(param.mat.cohesion0);
   std::vector<double> cohesion0_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss6, token, ',')) 
   {cohesion0_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.cohesion1.erase(std::remove(param.mat.cohesion1.begin(), param.mat.cohesion1.end(), '['), param.mat.cohesion1.end());
   param.mat.cohesion1.erase(std::remove(param.mat.cohesion1.begin(), param.mat.cohesion1.end(), ']'), param.mat.cohesion1.end());
   param.mat.cohesion1.erase(std::remove(param.mat.cohesion1.begin(), param.mat.cohesion1.end(), ' '), param.mat.cohesion1.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss7(param.mat.cohesion1);
   std::vector<double> cohesion1_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss7, token, ',')) 
   {cohesion1_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.friction_angle0.erase(std::remove(param.mat.friction_angle0.begin(), param.mat.friction_angle0.end(), '['), param.mat.friction_angle0.end());
   param.mat.friction_angle0.erase(std::remove(param.mat.friction_angle0.begin(), param.mat.friction_angle0.end(), ']'), param.mat.friction_angle0.end());
   param.mat.friction_angle0.erase(std::remove(param.mat.friction_angle0.begin(), param.mat.friction_angle0.end(), ' '), param.mat.friction_angle0.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss8(param.mat.friction_angle0);
   std::vector<double> friction_angle0_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss8, token, ',')) 
   {friction_angle0_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.friction_angle1.erase(std::remove(param.mat.friction_angle1.begin(), param.mat.friction_angle1.end(), '['), param.mat.friction_angle1.end());
   param.mat.friction_angle1.erase(std::remove(param.mat.friction_angle1.begin(), param.mat.friction_angle1.end(), ']'), param.mat.friction_angle1.end());
   param.mat.friction_angle1.erase(std::remove(param.mat.friction_angle1.begin(), param.mat.friction_angle1.end(), ' '), param.mat.friction_angle1.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss81(param.mat.friction_angle1);
   std::vector<double> friction_angle1_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss81, token, ',')) 
   {friction_angle1_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.dilation_angle0.erase(std::remove(param.mat.dilation_angle0.begin(), param.mat.dilation_angle0.end(), '['), param.mat.dilation_angle0.end());
   param.mat.dilation_angle0.erase(std::remove(param.mat.dilation_angle0.begin(), param.mat.dilation_angle0.end(), ']'), param.mat.dilation_angle0.end());
   param.mat.dilation_angle0.erase(std::remove(param.mat.dilation_angle0.begin(), param.mat.dilation_angle0.end(), ' '), param.mat.dilation_angle0.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss9(param.mat.dilation_angle0);
   std::vector<double> dilation_angle0_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss9, token, ',')) 
   {dilation_angle0_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.dilation_angle1.erase(std::remove(param.mat.dilation_angle1.begin(), param.mat.dilation_angle1.end(), '['), param.mat.dilation_angle1.end());
   param.mat.dilation_angle1.erase(std::remove(param.mat.dilation_angle1.begin(), param.mat.dilation_angle1.end(), ']'), param.mat.dilation_angle1.end());
   param.mat.dilation_angle1.erase(std::remove(param.mat.dilation_angle1.begin(), param.mat.dilation_angle1.end(), ' '), param.mat.dilation_angle1.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss91(param.mat.dilation_angle1);
   std::vector<double> dilation_angle1_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss91, token, ',')) 
   {dilation_angle1_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.pls0.erase(std::remove(param.mat.pls0.begin(), param.mat.pls0.end(), '['), param.mat.pls0.end());
   param.mat.pls0.erase(std::remove(param.mat.pls0.begin(), param.mat.pls0.end(), ']'), param.mat.pls0.end());
   param.mat.pls0.erase(std::remove(param.mat.pls0.begin(), param.mat.pls0.end(), ' '), param.mat.pls0.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss10(param.mat.pls0);
   std::vector<double> pls0_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss10, token, ',')) 
   {pls0_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.pls1.erase(std::remove(param.mat.pls1.begin(), param.mat.pls1.end(), '['), param.mat.pls1.end());
   param.mat.pls1.erase(std::remove(param.mat.pls1.begin(), param.mat.pls1.end(), ']'), param.mat.pls1.end());
   param.mat.pls1.erase(std::remove(param.mat.pls1.begin(), param.mat.pls1.end(), ' '), param.mat.pls1.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss11(param.mat.pls1);
   std::vector<double> pls1_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss11, token, ',')) 
   {pls1_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // String (Material) extraction
   param.mat.plastic_viscosity.erase(std::remove(param.mat.plastic_viscosity.begin(), param.mat.plastic_viscosity.end(), '['), param.mat.plastic_viscosity.end());
   param.mat.plastic_viscosity.erase(std::remove(param.mat.plastic_viscosity.begin(), param.mat.plastic_viscosity.end(), ']'), param.mat.plastic_viscosity.end());
   param.mat.plastic_viscosity.erase(std::remove(param.mat.plastic_viscosity.begin(), param.mat.plastic_viscosity.end(), ' '), param.mat.plastic_viscosity.end());

   // Create a stringstream to tokenize the string
   std::stringstream ss12(param.mat.plastic_viscosity);
   std::vector<double> plastic_viscosity_vec;
    
   // Tokenize the string and convert tokens to integers
   while (getline(ss12, token, ',')) 
   {plastic_viscosity_vec.push_back(std::stod(token)); // Convert string to int and add to vector
   }

   // Initialize density and specific internal energy values. We interpolate in
   // a non-positive basis to get the correct values at the dofs. Then we do an
   // L2 projection to the positive basis in which we actually compute. The goal
   // is to get a high-order representation of the initial condition. Note that
   // this density is a temporary function and it will not be updated during the
   // time evolution.
   num_materials =pmesh->attributes.Max();
   Vector z_rho(pmesh->attributes.Max());
   Vector s_rho(pmesh->attributes.Max());

   double pseudo_speed =  max_vbc_val * param.control.mscale;
   double pseudo_speed_sqrd =  pseudo_speed * pseudo_speed;

   // EC: Why consider this case separately? Whatever the size is, shouldn't it be the same with attributes.Max()?
   if(rho_vec.size() == 1) { 
      z_rho = rho_vec[0]; 
      s_rho = (lambda_vec[0] + 2*mu_vec[0]) / pseudo_speed_sqrd;
   }
   else if(rho_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0)
         cout << "The number of rho are not consistent with material ID in the given mesh." << endl;        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else 
   {
      for (int i = 0; i < pmesh->attributes.Max(); i++) {
         z_rho[i] = rho_vec[i]; 
         s_rho[i] = (lambda_vec[i] + 2*mu_vec[i]) / pseudo_speed_sqrd;
      }
   }
   
   // z_rho = 2700.0;
   // s_rho = 2700.0 * param.control.mscale;

   // GridFunctionCoefficient u_coeff(&u_alpha_gf);
   
   ParGridFunction rho0_gf(&L2FESpace);
   ParGridFunction fictitious_rho0_gf(&L2FESpace);

   PWConstCoefficient rho0_coeff(z_rho);
   PWConstCoefficient scale_rho0_coeff(s_rho);
   ParGridFunction l2_rho0_gf(&l2_fes),l2_e(&l2_fes);
   l2_rho0_gf.ProjectCoefficient(rho0_coeff);
   rho0_gf.ProjectGridFunction(l2_rho0_gf);
   l2_rho0_gf.ProjectCoefficient(scale_rho0_coeff);
   fictitious_rho0_gf.ProjectGridFunction(l2_rho0_gf);

   // rho_ini_gf.ProjectGridFunction(l2_rho0_gf);

   /*
   ParGridFunction rho0_gf(&L2FESpace);
   FunctionCoefficient rho0_coeff(rho0);
   // FunctionCoefficient super_rho0_coeff(rho0);
   L2_FECollection l2_fec(order_e, pmesh->Dimension());
   ParFiniteElementSpace l2_fes(pmesh, &l2_fec);
   ParGridFunction l2_rho0_gf(&l2_fes), l2_e(&l2_fes);
   l2_rho0_gf.ProjectCoefficient(rho0_coeff);
   rho0_gf.ProjectGridFunction(l2_rho0_gf);
   */
   
   if (param.sim.problem == 1)
   {
      // For the Sedov test, we use a delta function at the origin.
      DeltaCoefficient e_coeff(blast_position[0], blast_position[1],
                               blast_position[2], blast_energy);
      l2_e.ProjectCoefficient(e_coeff);
   }
   else
   {
      FunctionCoefficient e_coeff(e0);
      l2_e.ProjectCoefficient(e_coeff);
   }
   e_gf.ProjectGridFunction(l2_e);
   // e_gf = 0;
   // Sync the data location of e_gf with its base, S
   e_gf.SyncAliasMemory(S);

   // Piecewise constant elastic stiffness over the Lagrangian mesh.
   // Lambda and Mu is Lame's first and second constants
   Vector lambda(pmesh->attributes.Max());
   Vector mu(pmesh->attributes.Max());
   // lambda = param.mat.lambda;
   // mu = param.mat.mu;
   if(lambda_vec.size() == 1) {lambda =lambda_vec[0];}
   else if(lambda_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of lambda are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {lambda[i] = lambda_vec[i];}}

   if(mu_vec.size() == 1) {mu =mu_vec[0];}
   else if(mu_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of mu are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {mu[i] = mu_vec[i];}}

   PWConstCoefficient lambda_func(lambda);
   PWConstCoefficient mu_func(mu);
   
   // Project PWConstCoefficient to grid function
   L2_FECollection lambda_fec(param.mesh.order_e, pmesh->Dimension());
   ParFiniteElementSpace lambda_fes(pmesh, &lambda_fec);
   ParGridFunction lambda0_gf(&lambda_fes);
   lambda0_gf.ProjectCoefficient(lambda_func);
   
   L2_FECollection mu_fec(param.mesh.order_e, pmesh->Dimension());
   ParFiniteElementSpace mu_fes(pmesh, &mu_fec);
   ParGridFunction mu0_gf(&mu_fes);
   mu0_gf.ProjectCoefficient(mu_func);

   // Piecewise constant for material index
   Vector mat(pmesh->attributes.Max());
   for (int i = 0; i < pmesh->attributes.Max(); i++)
   {
      mat[i] = i;
   }
   PWConstCoefficient mat_func(mat);

   // Project PWConstCoefficient to grid function
   L2_FECollection mat_fec(param.mesh.order_e, pmesh->Dimension());
   ParFiniteElementSpace mat_fes(pmesh, &mat_fec);
   ParGridFunction mat_gf(&mat_fes);
   mat_gf.ProjectCoefficient(mat_func);

   // Composition
   ParFiniteElementSpace L2FESpace_mat(pmesh, &L2FEC, num_materials); // material composition 
   ParGridFunction comp_gf(&L2FESpace_mat); ParGridFunction comp_ref_gf(&L2FESpace_mat);
   CompoCoefficient comp_coeff(num_materials, mat_gf);
   comp_gf.ProjectCoefficient(comp_coeff); // Initialize the composition with material indicators

   // for (int i = 0; i < num_materials; i++)
   // {
   //    for (int j = 0; j < e_gf.Size(); j++)
   //    {
   //       if(mat_gf[j] == i){comp_gf[j+e_gf.Size()*i] = 1.0;}
   //    }
   // }

   comp_ref_gf = comp_gf;

   // Material properties of Plasticity
   Vector tension_cutoff(pmesh->attributes.Max());
   Vector cohesion0(pmesh->attributes.Max());
   Vector cohesion1(pmesh->attributes.Max());
   Vector friction_angle0(pmesh->attributes.Max());
   Vector friction_angle1(pmesh->attributes.Max());
   Vector dilation_angle0(pmesh->attributes.Max());
   Vector dilation_angle1(pmesh->attributes.Max());
   Vector plastic_viscosity(pmesh->attributes.Max());
   Vector pls0(pmesh->attributes.Max());
   Vector pls1(pmesh->attributes.Max());

   if(tension_cutoff_vec.size() == 1) {tension_cutoff =tension_cutoff_vec[0];}
   else if(tension_cutoff_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of tension_cutoff are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {tension_cutoff[i] = tension_cutoff_vec[i];}}

   if(cohesion0_vec.size() == 1) {cohesion0 =cohesion0_vec[0];}
   else if(cohesion0_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of cohesion0 are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {cohesion0[i] = cohesion0_vec[i];}}

   if(cohesion1_vec.size() == 1) {cohesion1 =cohesion1_vec[0];}
   else if(cohesion1_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of cohesion1 are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {cohesion1[i] = cohesion1_vec[i];}}

   if(friction_angle0_vec.size() == 1) {friction_angle0 =friction_angle0_vec[0];}
   else if(friction_angle0_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of friction_angle0 are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {friction_angle0[i] = friction_angle0_vec[i];}}

   if(friction_angle1_vec.size() == 1) {friction_angle1 =friction_angle1_vec[0];}
   else if(friction_angle1_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of friction_angle1 are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {friction_angle1[i] = friction_angle1_vec[i];}}

   if(dilation_angle0_vec.size() == 1) {dilation_angle0 =dilation_angle0_vec[0];}
   else if(dilation_angle0_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of dilation_angle0 are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {dilation_angle0[i] = dilation_angle0_vec[i];}}

   if(dilation_angle1_vec.size() == 1) {dilation_angle1 =dilation_angle1_vec[0];}
   else if(dilation_angle1_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of dilation_angle1 are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {dilation_angle1[i] = dilation_angle1_vec[i];}}

   if(pls0_vec.size() == 1) {pls0 =pls0_vec[0];}
   else if(pls0_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of pls0 are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {pls0[i] = pls0_vec[i];}}

   if(pls1_vec.size() == 1) {pls1 =pls1_vec[0];}
   else if(pls1_vec.size() != pmesh->attributes.Max())
   {
      if (myid == 0){cout << "The number of pls1 are not consistent with material ID in the given mesh." << endl; }        
      delete pmesh;
      MPI_Finalize();
      return 3;
   }
   else {for (int i = 0; i < pmesh->attributes.Max(); i++) {pls1[i] = pls1_vec[i];}}

   if(param.mat.viscoplastic)
   {
      if(plastic_viscosity_vec.size() == 1) {plastic_viscosity =plastic_viscosity_vec[0];}
      else if(plastic_viscosity_vec.size() != pmesh->attributes.Max())
      {
         if (myid == 0){cout << "The number of plastic_viscosity are not consistent with material ID in the given mesh." << endl; }        
         delete pmesh;
         MPI_Finalize();
         return 3;
      }
      else {for (int i = 0; i < pmesh->attributes.Max(); i++) {plastic_viscosity[i] = plastic_viscosity_vec[i];}}
   }
   else
   {
      if (myid == 0){cout << "viscoplasticity is not activate." << endl; }        
      plastic_viscosity = 1.0e+300;
   }
   
   // lithostatic pressure
   s_gf=0.0;
   ATMCoefficient ATM_coeff(dim, xyz_gf_l2, rho0_gf, param.control.gravity, param.control.thickness);
   s_gf.ProjectCoefficient(ATM_coeff);

   if(param.control.lithostatic)
   {
      LithostaticCoefficient Lithostatic_coeff(dim, xyz_gf_l2, rho0_gf, param.control.gravity, param.control.thickness);
      s_gf.ProjectCoefficient(Lithostatic_coeff);
   }   

   s_gf.SyncAliasMemory(S);

   // for(int i = 0; i < s_gf.Size()/3; i++)
   // {
   //    cout << "sxx " << s_gf[i+0*s_gf.Size()/3] << ", syy " << s_gf[i+1*s_gf.Size()/3] << ", sxy " << s_gf[i+2*s_gf.Size()/3] << endl;;
   // }

   ParGridFunction s_old_gf(&L2FESpace_stress);
   s_old_gf=s_gf;

   // Copy initial cooriante to x_gf
   ParGridFunction x_ini_gf(&H1FESpace); 
   x_ini_gf = x_gf; 
   // x_ini_gf.SyncAliasMemory(S);
   ParGridFunction x_old_gf(&H1FESpace);
   x_old_gf = 0.0;

   // Plastic strain (J2 strain invariant)
   ParGridFunction p_gf(&L2FESpace);
   ParGridFunction p_gf_old(&L2FESpace);
   p_gf=0.0; p_gf_old = 0.0;
   Vector weak_location(dim);
   if(dim == 2){weak_location[0] = param.mat.weak_x; weak_location[1] = param.mat.weak_y;}
   else if(dim ==3){weak_location[0] = param.mat.weak_x; weak_location[1] = param.mat.weak_y; weak_location[2] = param.mat.weak_z;}

   Vector ini_weakzone(pmesh->attributes.Max());
   for (int i = 0; i < pmesh->attributes.Max(); i++) {ini_weakzone[i] = 0.0;}
   ini_weakzone[pmesh->attributes.Max()-1] = 0.5;
   PWConstCoefficient weak_func(ini_weakzone);

   PlasticCoefficient p_coeff(dim, xyz_gf_l2, weak_location, param.mat.weak_rad, param.mat.ini_pls);
   // p_gf.ProjectCoefficient(p_coeff);
   // // interpolation using non-basis function
   ParGridFunction l2_p_gf(&l2_fes);
   l2_p_gf.ProjectCoefficient(p_coeff);
   p_gf.ProjectGridFunction(l2_p_gf);

   // p_gf.ProjectCoefficient(weak_func);

   p_gf_old = p_gf;
   // Non-initial plastic strain
   ParGridFunction ini_p_gf(&L2FESpace);
   ParGridFunction ini_p_old_gf(&L2FESpace);
   ParGridFunction n_p_gf(&L2FESpace);
   ini_p_gf=p_gf; ini_p_old_gf=p_gf;
   n_p_gf=0.0;

   
   ParGridFunction u_gf(&H1FESpace);  // Displacment
   u_gf = 0.0;


   // // The numerical sandbox case
   // if (param.tmop.tmop)
   // {
   //    ParGridFunction x_mod_gf(&H1FESpace); 
   //    // Store source mesh positions.
   //    ParMesh *pmesh_copy =  new ParMesh(*pmesh);
   //    x_mod_gf = x_gf;
            
   //    if(myid == 0){cout << "First Remeshing " << endl;}
   //    HR_adaptivity(pmesh_copy, x_mod_gf, ess_tdofs, myid, param.tmop.mesh_poly_deg, param.mesh.rs_levels, param.mesh.rp_levels, param.tmop.jitter, param.tmop.metric_id, param.tmop.target_id,\
   //                   param.tmop.lim_const, param.tmop.adapt_lim_const, param.tmop.quad_type, param.tmop.quad_order, param.tmop.solver_type, param.tmop.solver_iter, param.tmop.solver_rtol, \
   //                   param.tmop.solver_art_type, param.tmop.lin_solver, param.tmop.max_lin_iter, param.tmop.move_bnd, param.tmop.combomet, param.tmop.bal_expl_combo, param.tmop.hradaptivity, \
   //                   param.tmop.h_metric_id, param.tmop.normalization, param.tmop.verbosity_level, param.tmop.fdscheme, param.tmop.adapt_eval, param.tmop.exactaction, param.solver.p_assembly, \
   //                   param.tmop.n_hr_iter, param.tmop.n_h_iter, param.tmop.mesh_node_ordering, param.tmop.barrier_type, param.tmop.worst_case_type);

      
   //    x_gf = *pmesh_copy->GetNodes();  x_gf *= param.tmop.ale; x_gf.Add(1.0 - param.tmop.ale, x_old_gf);
   //    pmesh->NewNodes(x_gf, false); 
   //    delete pmesh_copy;

   //    xyz_gf_l2.ProjectCoefficient(xyz_coeff);

   //    if(param.control.lithostatic)
   //    {
   //       LithostaticCoefficient Lithostatic_coeff(dim, xyz_gf_l2, rho0_gf, param.control.gravity, param.control.thickness);
   //       s_gf.ProjectCoefficient(Lithostatic_coeff);
   //    }

   //    CompoCoefficient comp_coeff(dim, xyz_gf_l2, param.control.thickness);
   //    comp_gf.ProjectCoefficient(comp_coeff);
   //    rho0_gf = 0.0; lambda0_gf = 0.0; mu0_gf = 0.0;

   //    for(int j = 0; j < rho0_gf.Size(); j++ )
   //    {
   //       for (int i = 0; i < pmesh->attributes.Max(); i++)
   //       {
   //          rho0_gf[j] = rho0_gf[j] + z_rho[i]*comp_gf[j+rho0_gf.Size()*i];
   //          lambda0_gf[j] = lambda0_gf[j] + lambda[i]*comp_gf[j+rho0_gf.Size()*i];
   //          mu0_gf[j] = mu0_gf[j] + mu[i]*comp_gf[j+rho0_gf.Size()*i];
   //       }
   //    }
   // }

   // // Flattening node
   // ParLinearForm flattening(&H1FESpace); 
   // Array<int> nbc_bdr(pmesh->bdr_attributes.Max());   
   // nbc_bdr = 0; nbc_bdr[2] = 1; // bottome boundary
   // VectorArrayCoefficient bottom_node(dim);
   // for (int i = 0; i < dim-1; i++)
   // {
   //   bottom_node.Set(i, new ConstantCoefficient(0.0));
   // }

   // Vector bottom_node_id(pmesh->bdr_attributes.Max());
   // bottom_node_id = 0.0;
   // bottom_node_id(2) = 1.0;
   // bottom_node.Set(dim-1, new PWConstCoefficient(bottom_node_id));
   // flattening.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(bottom_node), nbc_bdr);
   // flattening.Assemble();

   // L2_FECollection mat_fec(0, pmesh->Dimension());
   // ParFiniteElementSpace mat_fes(pmesh, &mat_fec);
   // ParGridFunction mat_gf(&mat_fes);
   // FunctionCoefficient mat_coeff(gamma_func);
   // mat_gf.ProjectCoefficient(mat_coeff);

   int source = 0; bool visc = false, vorticity = false;
   // switch (param.sim.problem)
   // {
   //    case 0: if (pmesh->Dimension() == 2) { source = 1; } visc = false; break;
   //    case 1: visc = true; break;
   //    case 2: visc = true; break;
   //    case 3: visc = true; S.HostRead(); break;
   //    case 4: visc = false; break;
   //    case 5: visc = true; break;
   //    case 6: visc = true; break;
   //    case 7: source = 2; visc = true; vorticity = true;  break;
   //    default: MFEM_ABORT("Wrong problem specification!");
   // }
   if (param.solver.impose_visc) { visc = true; }

   
   geodynamics::LagrangianGeoOperator geo(S.Size(),
                                          H1FESpace, L2FESpace, L2FESpace_stress, ess_tdofs,
                                          // rho0_coeff, scale_rho0_coeff, rho0_gf,
                                          rho0_gf, fictitious_rho0_gf,
                                          mat_gf, source, param.solver.cfl,
                                          visc, vorticity, param.solver.p_assembly,
                                          param.solver.cg_tol, param.solver.cg_max_iter, 
                                          param.solver.ftz_tol,
                                          param.mesh.order_q, lambda0_gf, mu0_gf, param.control.mscale, param.control.gravity, param.control.thickness,
                                          param.bc.winkler_foundation, param.bc.winkler_rho, param.control.dyn_damping, param.control.dyn_factor, bc_id_pa, max_vbc_val);
    

   socketstream vis_rho, vis_v, vis_e;
   char vishost[] = "localhost";
   int  visport   = 19916;

   // ParGridFunction rho_gf;
   // if (param.sim.visualization || param.sim.visit || param.sim.paraview) { geo.ComputeDensity(rho_gf); }
   // if (mass_bal) { geo.ComputeDensity(rho_gf); }
   const double energy_init = geo.InternalEnergy(e_gf) +
                              geo.KineticEnergy(v_gf);

   if (param.sim.visualization)
   {
      // Make sure all MPI ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh->GetComm());
      vis_rho.precision(8);
      vis_v.precision(8);
      vis_e.precision(8);
      int Wx = 0, Wy = 0; // window position
      const int Ww = 350, Wh = 350; // window size
      int offx = Ww+10; // window offsets
      if (param.sim.problem != 0 && param.sim.problem != 4)
      {
         geodynamics::VisualizeField(vis_rho, vishost, visport, rho0_gf,
                                       "Density", Wx, Wy, Ww, Wh);
      }
      Wx += offx;
      geodynamics::VisualizeField(vis_v, vishost, visport, v_gf,
                                    "Velocity", Wx, Wy, Ww, Wh);
      Wx += offx;
      geodynamics::VisualizeField(vis_e, vishost, visport, e_gf,
                                    "Specific Internal Energy", Wx, Wy, Ww, Wh);
   }

   // Save data for VisIt visualization.
   VisItDataCollection visit_dc(param.sim.basename, pmesh);
   if (param.sim.visit)
   {
      visit_dc.RegisterField("Density",  &rho0_gf);
      visit_dc.RegisterField("Displacement", &u_gf);
      visit_dc.RegisterField("Velocity", &v_gf);
      visit_dc.RegisterField("Specific Internal Energy", &e_gf);
      visit_dc.RegisterField("Stress", &s_gf);
      visit_dc.RegisterField("Plastic Strain", &p_gf);
      visit_dc.RegisterField("Non-inital Plastic Strain", &n_p_gf);
      // visit_dc.RegisterField("Geometric Parameters", &quality);
      visit_dc.RegisterField("Composition", &comp_gf);
      visit_dc.RegisterField("Lambda", &lambda0_gf);
      visit_dc.RegisterField("Mu", &mu0_gf);
      visit_dc.SetCycle(0);
      visit_dc.SetTime(0.0);
      visit_dc.Save();
   }

   ParaViewDataCollection *pd = NULL;
   if (param.sim.paraview)
   {
      pd = new ParaViewDataCollection(param.sim.basename, pmesh);
      // pd->SetPrefixPath("ParaView");
      pd->RegisterField("Density",  &rho0_gf);
      pd->RegisterField("Displacement", &u_gf);
      pd->RegisterField("Velocity", &v_gf);
      pd->RegisterField("Specific Internal Energy", &e_gf);
      pd->RegisterField("Stress", &s_gf);
      pd->RegisterField("Plastic Strain", &p_gf);
      pd->RegisterField("inital Plastic Strain", &ini_p_gf);
      pd->RegisterField("Non-inital Plastic Strain", &n_p_gf);
      pd->RegisterField("Geometric Parameters", &quality);
      pd->RegisterField("Composition", &comp_gf);
      pd->RegisterField("Lambda", &lambda0_gf);
      pd->RegisterField("Mu", &mu0_gf);
      pd->SetLevelsOfDetail(param.mesh.order_v);
      pd->SetDataFormat(VTKFormat::BINARY);
      // pd->SetDataFormat(VTKFormat::ASCII);
      pd->SetHighOrderOutput(true);
      pd->SetCycle(0);
      pd->SetTime(0.0);
      pd->Save();
   }

   // Perform time-integration (looping over the time iterations, ti, with a
   // time-step dt). The object oper is of type LagrangianGeoOperator that
   // defines the Mult() method that used by the time integrators.
   ode_solver->Init(geo);
   geo.ResetTimeStepEstimate();
   double t = 0.0, dt = 0.0, t_old, dt_old = 0.0, h_min = 1.0;
   // dt = geo.GetTimeStepEstimate(S, dt); // To provide dt before the estimate, initializing is necessary
   // h_min = geo.GetLengthEstimate(S, dt); // To provide dt before the estimate, initializing is necessary
   dt = geo.GetTimeStepEstimate(S); // To provide dt before the estimate, initializing is necessary
   h_min = geo.GetLengthEstimate(S); // To provide dt before the estimate, initializing is necessar
   double ini_h_min = h_min;
   dt = init_dt;
   bool last_step = false;
   int steps = 0;
   BlockVector S_old(S);
   long mem=0, mmax=0, msum=0;
   int checks = 0;

   // 10. Perform time-integration (looping over the time iterations, ti, with a
   //     time-step dt).
   ode_solver_sub->Init(oper_sub); ode_solver_sub2->Init(oper_sub2);

   //   const double internal_energy = geo.InternalEnergy(e_gf);
   //   const double kinetic_energy = geo.KineticEnergy(v_gf);
   //   if (mpi.Root())
   //   {
   //      cout << std::fixed;
   //      cout << "step " << std::setw(5) << 0
   //            << ",\tt = " << std::setw(5) << std::setprecision(4) << t
   //            << ",\tdt = " << std::setw(5) << std::setprecision(6) << dt
   //            << ",\t|IE| = " << std::setprecision(10) << std::scientific
   //            << internal_energy
   //            << ",\t|KE| = " << std::setprecision(10) << std::scientific
   //            << kinetic_energy
   //            << ",\t|E| = " << std::setprecision(10) << std::scientific
   //            << kinetic_energy+internal_energy;
   //      cout << std::fixed;
   //      if (mem_usage)
   //      {
   //         cout << ", mem: " << mmax << "/" << msum << " MB";
   //      }
   //      cout << endl;
   //   }

   if (mpi.Root()) 
   {
      std::cout<<""<<std::endl;
      std::cout<<"simulation starts"<<std::endl;
   }

   // Estimate element errors using the Zienkiewicz-Zhu error estimator.
   // Vector errors(pmesh->GetNE()); errors=0.0;
   // if(param.tmop.amr)
   // {
   //    // 11. As in Example 6p, we also need a refiner. This time the refinement
   //    //     strategy is based on a fixed threshold that is applied locally to each
   //    //     element. The global threshold is turned off by setting the total error
   //    //     fraction to zero. We also enforce a maximum refinement ratio between
   //    //     adjacent elements.
   //    ThresholdRefiner refiner(&errors);
   //    refiner.SetTotalErrorFraction(0.0); // use purely local threshold
   //    refiner.SetLocalErrorGoal(1e-4);
   //    refiner.PreferConformingRefinement();
   //    refiner.SetNCLimit(2);

   //    // 12. A derefiner selects groups of elements that can be coarsened to form
   //    //     a larger element. A conservative enough threshold needs to be set to
   //    //     prevent derefining elements that would immediately be refined again.
   //    ThresholdDerefiner derefiner(&errors);
   //    derefiner.SetThreshold(1e-4 * 0.25);
   //    derefiner.SetNCLimit(2);
   // }


   for (int ti = 1; !last_step; ti++)
   {
      if (t + dt >= param.sim.t_final)
      {
         dt = param.sim.t_final - t;
         last_step = true;
      }
      if (steps == param.sim.max_tsteps) { last_step = true; }
      S_old = S;
      t_old = t;
      double year = t/86400/365.25;
      p_gf_old = p_gf; ini_p_old_gf = ini_p_gf; x_old_gf = x_gf;
      geo.ResetTimeStepEstimate();
      // S is the vector of dofs, t is the current time, and dt is the time step
      // to advance.

      // // Adjusting nodes on side walls based on original mesh. 
      // ParSubMesh::Transfer(x_ini_gf, x0_side);
      // ParSubMesh::Transfer(x_ini_gf, x1_side);

      if(param.control.pseudo_transient)
      {
         for (int i = 0; i < param.control.transient_num; i++)
         {
            x_gf = x_old_gf; // back to orignal mesh to fix mesh during pseudo transient loop
            s_gf = s_old_gf; // 
            // e_gf = e_old_gf; // 

            ode_solver->Step(S, t, dt);
         }
         t = t - dt*(param.control.transient_num-1.0);
      }
      else
      {
         ode_solver->Step(S, t, dt);
      }

      // // Keep the domain size for ALE
      // ParSubMesh::Transfer(x0_side, x_gf);
      // ParSubMesh::Transfer(x1_side, x_gf);

      if(param.bc.surf_proc)
      {
         ParSubMesh::Transfer(x_gf, x_top); // update current mesh to submesh
         for (int i = 0; i < topo.Size(); i++){topo[i] = x_top[i+topo.Size()];}
         topo.GetTrueDofs(topo_t); 
         // x_top_old=x_top; 
         topo_t_old=topo; 
         ode_solver_sub->Step(topo_t, t, dt); t=t-dt;
         topo.SetFromTrueDofs(topo_t);
         for (int i = 0; i < topo.Size(); i++){x_top[i+topo.Size()] = topo[i];}
         submesh.NewNodes(x_top, false);
         ParSubMesh::Transfer(x_top, x_gf); // update adjusted nodes on top boundary 
      }

      if(param.bc.winkler_foundation & param.bc.base_proc)
      {
         ParSubMesh::Transfer(x_gf, x_bottom); // update current mesh to submesh
         for (int i = 0; i < bottom.Size(); i++){bottom[i] = x_bottom[i+bottom.Size()];}
         bottom.GetTrueDofs(bottom_t); 
         // x_bottom_old=x_bottom; 
         bottom_t_old=bottom; 
         ode_solver_sub2->Step(bottom_t, t, dt); t=t-dt;
         bottom.SetFromTrueDofs(bottom_t);
         if(param.bc.winkler_flat)
         {
            for (int i = 0; i < bottom.Size(); i++){x_bottom[i+bottom.Size()] = 0.0;}
         }
         else
         {
            for (int i = 0; i < bottom.Size(); i++){x_bottom[i+bottom.Size()] = bottom[i];}
         }
         submesh_bottom.NewNodes(x_bottom, false);
         ParSubMesh::Transfer(x_bottom, x_gf); // update adjusted nodes on top boundary 
      }

      if(param.mat.plastic)
      {
         /*
         // temporay stress and plastic strain array
         ParFiniteElementSpace L2FESpace_temp(pmesh, &L2FEC, 3*(dim-1) + 1); 
         ParGridFunction temp_gf(&L2FESpace_temp); ParGridFunction temp2_gf(&L2FESpace); 
         temp_gf = 1.0;
          
         if(dim == 2)
         {
            // simplified model is on test. 
            // ReturnMapping2D_simple_Coefficient Return_coeff(dim, h_min, dt_old, param.mat.viscoplastic, comp_gf, s_gf, s_old_gf, p_gf, mat_gf, z_rho, lambda, mu, tension_cutoff, cohesion0, cohesion1, pls0, pls1, friction_angle0, friction_angle1, dilation_angle0, dilation_angle1, plastic_viscosity);
            ReturnMapping2DCoefficient Return_coeff(dim, h_min, dt_old, param.mat.viscoplastic, comp_gf, s_gf, s_old_gf, p_gf, mat_gf, z_rho, lambda, mu, tension_cutoff, cohesion0, cohesion1, pls0, pls1, friction_angle0, friction_angle1, dilation_angle0, dilation_angle1, plastic_viscosity);
            temp_gf.ProjectCoefficient(Return_coeff);
         }
         else
         {
            ReturnMapping3DCoefficient Return_coeff(dim, h_min, dt_old, param.mat.viscoplastic, comp_gf, s_gf, s_old_gf, p_gf, mat_gf, z_rho, lambda, mu, tension_cutoff, cohesion0, cohesion1, pls0, pls1, friction_angle0, friction_angle1, dilation_angle0, dilation_angle1, plastic_viscosity);
            temp_gf.ProjectCoefficient(Return_coeff);
         }
         StressMappingCoefficient Stress_coeff(dim, temp_gf);
         s_gf.ProjectCoefficient(Stress_coeff);
         PlasticityMappingCoefficient Plasticity_coeff(dim, temp_gf);
         temp2_gf.ProjectCoefficient(Plasticity_coeff);
         p_gf.Add(1.0, temp2_gf);
         */

         if(dim == 2){Returnmapping2d (comp_gf, s_gf, s_old_gf, p_gf, mat_gf, dim, h_min, z_rho, lambda, mu, tension_cutoff, cohesion0, cohesion1, pls0, pls1, friction_angle0, friction_angle1, dilation_angle0, dilation_angle1, plastic_viscosity, param.mat.viscoplastic, dt_old);}
         else{Returnmapping3d (comp_gf, s_gf, s_old_gf, p_gf, mat_gf, dim, h_min, z_rho, lambda, mu, tension_cutoff, cohesion0, cohesion1, pls0, pls1, friction_angle0, friction_angle1, dilation_angle0, dilation_angle1, plastic_viscosity, param.mat.viscoplastic, dt_old);}   
         n_p_gf  = ini_p_gf;
         n_p_gf -= p_gf;
         n_p_gf.Neg();
      }

      steps++;
      dt_old = dt;

      // Adaptive time step control.
      // const double dt_est = geo.GetTimeStepEstimate(S, dt);
      // double dt_est = geo.GetTimeStepEstimate(S, dt);
      // h_min = geo.GetLengthEstimate(S, dt);

      double dt_est = geo.GetTimeStepEstimate(S);
      h_min = geo.GetLengthEstimate(S);
      double cond_num = 0.0;
      // cond_num = ini_h_min/h_min;
      double global_min_vol = 1.0e5;

      if (param.tmop.tmop)
      {

         if(ti == 1 || (ti % param.tmop.remesh_steps) == 0 || cond_num > param.tmop.tmop_cond_num || global_min_vol < 1e3)
         {
            if(myid == 0)
            {
               if ((ti % param.tmop.remesh_steps) == 0){cout << "*** calling remeshing due to constant remeshing step " << param.tmop.remesh_steps << endl;}
               else if (cond_num > param.tmop.tmop_cond_num){cout << "*** calling remeshing due to relative aspect ratio is greater than " << param.tmop.tmop_cond_num << endl;}
               else if (ti == 1){cout << "*** Initial remehsing *** " << endl;}
               else if (global_min_vol < 1e3){cout << "*** calling remeshing due to small jacobian " << global_min_vol << endl;}
            }

            ti = ti -1;
            if (param.sim.visit)
            {
                  visit_dc.SetCycle(ti);
                  visit_dc.SetTime(t*0.995);
                  if(param.sim.year){visit_dc.SetTime(year-1);}
                  visit_dc.Save();
            }

            if (param.sim.paraview)
            {
                  pd->SetCycle(ti);
                  pd->SetTime(t*0.995);
                  if(param.sim.year){pd->SetTime(year-1);}
                  pd->Save();
            }

            ti = ti+1;

            // mass balance
            CompMassCoefficient CompBalance(num_materials, comp_ref_gf, vol_ini_gf, quality);
            comp_gf.ProjectCoefficient(CompBalance); // Initialize the composition with material indicators
            ParGridFunction x_mod_gf(&H1FESpace); ParGridFunction x_mod2_gf(&H1FESpace);
            // Store source mesh positions.
            ParMesh *pmesh_copy =  new ParMesh(*pmesh);
            ParMesh *pmesh_copy_old =  new ParMesh(*pmesh);
            ParMesh *pmesh_old  =  new ParMesh(*pmesh);

            x_old_gf = *pmesh->GetNodes();
            x_mod_gf = x_ini_gf;

            pmesh->GetBoundingBox(bb_min2, bb_max2, max(param.mesh.order_v, 1));
            if(dim == 2)
            {
               bb_center2[0] = (bb_min2[0]+bb_max2[0])*0.5;
               bb_center2[1] = (bb_min2[1]+bb_max2[1])*0.5;
               bb_length2[0] = (bb_max2[0]-bb_min2[0]); // x width
               bb_length2[1] = (bb_max2[1]-bb_min2[1]); // y height
               stretching_factor[0] = bb_length2[0]/bb_length[0];
               stretching_factor[1] = bb_length2[1]/bb_length[1];

                  if (myid == 0)
                  { 
                     cout << "streching factor x: " << stretching_factor[0] << ", streching factor y:" << stretching_factor[1] << endl; 
                  }
            }
            else
            {
               bb_center2[0] = (bb_min2[0]+bb_max2[0])*0.5;
               bb_center2[1] = (bb_min2[1]+bb_max2[1])*0.5;
               bb_center2[2] = (bb_min2[2]+bb_max2[2])*0.5;
               bb_length2[0] = (bb_max2[0]-bb_min2[0]); // x width
               bb_length2[1] = (bb_max2[1]-bb_min2[1]); // y width
               bb_length2[2] = (bb_max2[2]-bb_min2[2]); // z height
               stretching_factor[0] = bb_length2[0]/bb_length[0];
               stretching_factor[1] = bb_length2[1]/bb_length[1];
               stretching_factor[2] = bb_length2[2]/bb_length[2];
            }

            for( int i = 0; i < x_mod_gf.Size()/dim; i++ )
            {
               if(dim == 2)
               {
                  x_mod_gf[i] = (x_mod_gf[i] - bb_center[0])*stretching_factor[0] + bb_center2[0];
                  x_mod_gf[i + x_mod_gf.Size()/dim] = (x_mod_gf[i + x_mod_gf.Size()/dim] - bb_center[1])*stretching_factor[1] + bb_center2[1];
               }  
               else
               {
                  x_mod_gf[i] = (x_mod_gf[i] - bb_center[0])*stretching_factor[0] + bb_center2[0];
                  x_mod_gf[i + x_mod_gf.Size()/dim] = (x_mod_gf[i + x_mod_gf.Size()/dim] - bb_center[1])*stretching_factor[1] + bb_center2[1];
                  x_mod_gf[i + 2*x_mod_gf.Size()/dim] = (x_mod_gf[i + 2*x_mod_gf.Size()/dim] - bb_center[2])*stretching_factor[2] + bb_center2[2];
               }
            }

            x_mod2_gf = x_mod_gf; // store streched initial mesh

            /*
            // Calculating thickness of side walls of current mesh and streched mesh to adujust node postion of stresched mesh.
            ParSubMesh::Transfer(x_gf, x0_side);
            ParSubMesh::Transfer(x_gf, x1_side); 

            Vector x0(x0_side.Size());  
            Vector x1(x1_side.Size());  
            double local_min_x0{0.0}, local_max_x0{0.0};
            double local_min_x1{0.0}, local_max_x1{0.0};
            double global_min_x0, global_max_x0;
            double global_min_x1, global_max_x1;
            for (int i = 0; i < x0_side.Size()/dim; i++)
            {
               if(dim ==2)
               {
                  local_min_x0 = std::min(local_min_x0, x0_side[i+1*x0_side.Size()/dim]);
                  local_max_x0 = std::max(local_max_x0, x0_side[i+1*x0_side.Size()/dim]);
               }
               else
               {
                  local_min_x0 = std::min(local_min_x0, x0_side[i+2*x0_side.Size()/dim]);
                  local_max_x0 = std::max(local_max_x0, x0_side[i+2*x0_side.Size()/dim]);
               }
            }

            for (int i = 0; i < x1_side.Size()/dim; i++)
            {
               if(dim ==2)
               {
                  local_min_x1 = std::min(local_min_x1, x1_side[i+1*x1_side.Size()/dim]);
                  local_max_x1 = std::max(local_max_x1, x1_side[i+1*x1_side.Size()/dim]);
               }
               else
               {
                  local_min_x1 = std::min(local_min_x1, x1_side[i+2*x1_side.Size()/dim]);
                  local_max_x1 = std::max(local_max_x1, x1_side[i+2*x1_side.Size()/dim]);
               }
            }

            MPI_Reduce(&local_min_x0, &global_min_x0, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_min_x0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local_max_x0, &global_max_x0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_max_x0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local_min_x1, &global_min_x1, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_min_x1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local_max_x1, &global_max_x1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_max_x1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


            double x0_thickness = global_max_x0 - global_min_x0;
            double x1_thickness = global_max_x1 - global_min_x1;

            ParSubMesh::Transfer(x_mod2_gf, x0_side);
            ParSubMesh::Transfer(x_mod2_gf, x1_side);

            for (int i = 0; i < x0_side.Size()/dim; i++)
            {
               if(dim ==2)
               {
                  local_min_x0 = std::min(local_min_x0, x0_side[i+1*x0_side.Size()/dim]);
                  local_max_x0 = std::max(local_max_x0, x0_side[i+1*x0_side.Size()/dim]);
               }
               else
               {
                  local_min_x0 = std::min(local_min_x0, x0_side[i+2*x0_side.Size()/dim]);
                  local_max_x0 = std::max(local_max_x0, x0_side[i+2*x0_side.Size()/dim]);
               }
            }

            for (int i = 0; i < x1_side.Size()/dim; i++)
            {
               if(dim ==2)
               {
                  local_min_x1 = std::min(local_min_x1, x1_side[i+1*x1_side.Size()/dim]);
                  local_max_x1 = std::max(local_max_x1, x1_side[i+1*x1_side.Size()/dim]);
               }
               else
               {
                  local_min_x1 = std::min(local_min_x1, x1_side[i+2*x1_side.Size()/dim]);
                  local_max_x1 = std::max(local_max_x1, x1_side[i+2*x1_side.Size()/dim]);
               }
            }

            MPI_Reduce(&local_min_x0, &global_min_x0, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_min_x0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local_max_x0, &global_max_x0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_max_x0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local_min_x1, &global_min_x1, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_min_x1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local_max_x1, &global_max_x1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_max_x1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            double x0_ini_thickness = global_max_x0 - global_min_x0;
            double x1_ini_thickness = global_max_x1 - global_min_x1;
            */
           
            double global_max_top = bb_max2[1];
            double global_min_bot = bb_min2[1];

            {
               // Projecting top and bottom boundaries on flat surface
               ParSubMesh::Transfer(x_mod_gf, x_top); // update current mesh to submesh
               for (int i = 0; i < topo.Size(); i++){x_top[i+topo.Size()] = global_max_top;}
               ParSubMesh::Transfer(x_top, x_mod_gf);


               ParSubMesh::Transfer(x_mod_gf, x_bottom); // update current mesh to submesh
               if(param.bc.winkler_foundation)
               {
                  for (int i = 0; i < bottom.Size(); i++){x_bottom[i+bottom.Size()] = global_min_bot;}
                  ParSubMesh::Transfer(x_bottom, x_mod_gf);
               }
            }
            
            {
               x_gf = x_old_gf;
               // Projecting top and bottom boundaries on flat surface
               ParSubMesh::Transfer(x_gf, x_top); // update current mesh to submesh
               for (int i = 0; i < topo.Size(); i++){x_top[i+topo.Size()] = global_max_top;}
               ParSubMesh::Transfer(x_top, x_gf);


               ParSubMesh::Transfer(x_gf, x_bottom); // update current mesh to submesh
               if(param.bc.winkler_foundation)
               {
                  for (int i = 0; i < bottom.Size(); i++){x_bottom[i+bottom.Size()] = global_min_bot;}
                  ParSubMesh::Transfer(x_bottom, x_gf);
               }
            }

            pmesh_copy->NewNodes(x_mod_gf, false); // Deformined mesh for H1 interpolation 
            pmesh_copy_old->NewNodes(x_gf, false); // 

        
            Vector vxyz;
            int point_ordering;
            vxyz = *pmesh_copy->GetNodes(); // from target mesh
            point_ordering = pmesh_copy->GetNodes()->FESpace()->GetOrdering();
            FindPointsGSLIB finder(MPI_COMM_WORLD);
            finder.Setup(*pmesh_copy_old); // source mesh
            Vector interp_vals(x_gf.Size());               
            finder.Interpolate(vxyz, x_old_gf, interp_vals, point_ordering);
            x_mod_gf = interp_vals; 
            
            x_gf = x_old_gf; // back to original gridfunction

            /*
            // Adjusting nodes on side walls based on original mesh. 
            ParSubMesh::Transfer(x_mod2_gf, x0_side);
            ParSubMesh::Transfer(x_mod2_gf, x1_side);
            if(dim == 2)
            {
               for (int i = 0; i < x0_side.Size()/dim; i++){x0_side[i+x0_side.Size()/dim] = x0_thickness*x0_side[i+1*x0_side.Size()/dim]/x0_ini_thickness;}
               for (int i = 0; i < x1_side.Size()/dim; i++){x1_side[i+x1_side.Size()/dim] = x1_thickness*x1_side[i+1*x1_side.Size()/dim]/x1_ini_thickness;}
            }
            else
            {
               for (int i = 0; i < x0_side.Size()/dim; i++){x0_side[i+x0_side.Size()/dim] = x0_thickness*x0_side[i+2*x0_side.Size()/dim]/x0_ini_thickness;}
               for (int i = 0; i < x1_side.Size()/dim; i++){x1_side[i+x1_side.Size()/dim] = x1_thickness*x1_side[i+2*x1_side.Size()/dim]/x1_ini_thickness;}
            }
            ParSubMesh::Transfer(x0_side, x_gf);
            ParSubMesh::Transfer(x1_side, x_gf); 
            */

            // transfer interpolate coord to submesh
            ParSubMesh::Transfer(x_mod_gf, x_top); 
            ParSubMesh::Transfer(x_mod_gf, x_bottom);
            // transfer interpolate coord in submesh to original mesh  
            ParSubMesh::Transfer(x_top, x_gf);
            ParSubMesh::Transfer(x_bottom, x_gf);
            pmesh_copy->NewNodes(x_gf, false);
            
            if(myid == 0){cout << "First Remeshing " << endl;}
            HR_adaptivity(pmesh_copy, x_mod_gf, ess_tdofs, myid, param.tmop.mesh_poly_deg, param.mesh.rs_levels, param.mesh.rp_levels, param.tmop.jitter, param.tmop.metric_id, param.tmop.target_id,\
                           param.tmop.lim_const, param.tmop.adapt_lim_const, param.tmop.quad_type, param.tmop.quad_order, param.tmop.solver_type, param.tmop.solver_iter, param.tmop.solver_rtol, \
                           param.tmop.solver_art_type, param.tmop.lin_solver, param.tmop.max_lin_iter, param.tmop.move_bnd, param.tmop.combomet, param.tmop.bal_expl_combo, param.tmop.hradaptivity, \
                           param.tmop.h_metric_id, param.tmop.normalization, param.tmop.verbosity_level, param.tmop.fdscheme, param.tmop.adapt_eval, param.tmop.exactaction, param.solver.p_assembly, \
                           param.tmop.n_hr_iter, param.tmop.n_h_iter, param.tmop.mesh_node_ordering, param.tmop.barrier_type, param.tmop.worst_case_type);

            mesh_changed = true;

            if(param.tmop.move_bnd)
            {
               Vector vxyz;
               int point_ordering;
               // vxyz = *pmesh_copy_old->GetNodes(); // from target mesh
               vxyz = *pmesh_copy->GetNodes(); // from target mesh
               // point_ordering = pmesh_copy_old->GetNodes()->FESpace()->GetOrdering();
               point_ordering = pmesh_copy->GetNodes()->FESpace()->GetOrdering();
               FindPointsGSLIB finder(MPI_COMM_WORLD);
               // finder.Setup(*pmesh_copy); // source mesh
               finder.Setup(*pmesh_copy_old); // source mesh
               Vector interp_vals(x_gf.Size());               
               finder.Interpolate(vxyz, x_old_gf, interp_vals, point_ordering);
               x_mod_gf = interp_vals; 
               x_gf = x_old_gf; // back to original gridfunction

               // 
               ParSubMesh::Transfer(x_mod_gf, x_top); 
               ParSubMesh::Transfer(x_mod_gf, x_bottom);
               ParSubMesh::Transfer(x_top, x_gf);
               ParSubMesh::Transfer(x_bottom, x_gf);
               pmesh_copy->NewNodes(x_gf, false); 

               // fixed boundary TMOP
               double lim_const  = 0.0;
               param.tmop.move_bnd = false;
               if(myid == 0){cout << "Second Remeshing " << endl;}
               HR_adaptivity(pmesh_copy, x_mod_gf, ess_tdofs, myid, param.tmop.mesh_poly_deg, param.mesh.rs_levels, param.mesh.rp_levels, param.tmop.jitter, param.tmop.metric_id, param.tmop.target_id,\
                              param.tmop.lim_const, param.tmop.adapt_lim_const, param.tmop.quad_type, param.tmop.quad_order, param.tmop.solver_type, param.tmop.solver_iter, param.tmop.solver_rtol, \
                              param.tmop.solver_art_type, param.tmop.lin_solver, param.tmop.max_lin_iter, param.tmop.move_bnd, param.tmop.combomet, param.tmop.bal_expl_combo, param.tmop.hradaptivity, \
                              param.tmop.h_metric_id, param.tmop.normalization, param.tmop.verbosity_level, param.tmop.fdscheme, param.tmop.adapt_eval, param.tmop.exactaction, param.solver.p_assembly, \
                              param.tmop.n_hr_iter, param.tmop.n_h_iter, param.tmop.mesh_node_ordering, param.tmop.barrier_type, param.tmop.worst_case_type);
               param.tmop.move_bnd = true;
            }
            
            x_gf = *pmesh_copy->GetNodes();  x_gf *= param.tmop.ale; x_gf.Add(1.0 - param.tmop.ale, x_old_gf);
            pmesh->NewNodes(x_gf, false); 
            pmesh_copy->NewNodes(x_gf, false); // Deformined mesh for H1 interpolation
            {
               ParGridFunction U(&H1FESpace); 
               U =0.0; U = x_old_gf;
               ParGridFunction S1(&L2FESpace); ParGridFunction S2(&L2FESpace); ParGridFunction S3(&L2FESpace);
               ParGridFunction S4(&L2FESpace); ParGridFunction S5(&L2FESpace); ParGridFunction S6(&L2FESpace);
               S1 =0.0; S2 =0.0; S3 =0.0; 
               S4 =0.0; S5 =0.0; S6 =0.0;
               ParGridFunction comps(&L2FESpace); comps =0.0;
               ParGridFunction rmass(&L2FESpace); rmass =1.0; // to prevent mass leaking or adding after remeshing
               double all_comp = 0.0;

               // if (param.control.mass_bal && ti > 1) { geo.ComputeDensity(rho_gf); rho0_gf = rho_gf;}

               if(dim == 2){for(int i = 0; i < S1.Size(); i++ ){S1[i] = s_gf[i+S1.Size()*0];S2[i] = s_gf[i+S1.Size()*1];S3[i] = s_gf[i+S1.Size()*2];} }
               else{for(int i = 0; i < S1.Size(); i++ ){S1[i] = s_gf[i+S1.Size()*0];S2[i] = s_gf[i+S1.Size()*1];S3[i] = s_gf[i+S1.Size()*2]; S4[i] = s_gf[i+S1.Size()*3];S5[i] = s_gf[i+S1.Size()*4];S6[i] = s_gf[i+S1.Size()*5];}}
               
               if(myid==0){std::cout << "remapping for L2" << std::endl;}
               
               {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, rmass, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}

               for (int i = 0; i < pmesh->attributes.Max(); i++)
               {
                  for(int j = 0; j < comps.Size(); j++ ){comps[j] = comp_gf[j+comps.Size()*i];}
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, comps, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
                  for(int j = 0; j < comps.Size(); j++ ){comp_gf[j+comps.Size()*i] = comps[j]/rmass[j];}
                  comps =0.0;
               }
               
               {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, e_gf, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
               {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, p_gf, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
               {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, ini_p_gf, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
               {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, rho0_gf, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
               {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, fictitious_rho0_gf, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
               
               if(dim == 2)
               {
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, S1, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, S2, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, S3, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
               }
               else
               {
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, S1, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, S2, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, S3, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, S4, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, S5, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
                  {ParMesh *pmesh_old1 =  new ParMesh(*pmesh_old); Remapping(pmesh_old1, U, x_gf, S6, param.mesh.order_v, param.mesh.order_e, param.solver.p_assembly,param.mesh.local_refinement); delete pmesh_old1; U = x_old_gf;}
               }
               lambda0_gf = 0.0; mu0_gf = 0.0;
               for(int j = 0; j < comps.Size(); j++ )
               {
                  all_comp = 0.0;
                  for (int i = 0; i < pmesh->attributes.Max(); i++){all_comp = all_comp + comp_gf[j+comps.Size()*i];}
                  e_gf[j] = e_gf[j]/rmass[j]; p_gf[j] = p_gf[j]/rmass[j]; ini_p_gf[j] = ini_p_gf[j]/rmass[j];
                  rho0_gf[j] = rho0_gf[j]/rmass[j]; 
                  fictitious_rho0_gf[j] = fictitious_rho0_gf[j]/rmass[j]; //
                  for (int i = 0; i < pmesh->attributes.Max(); i++)
                  {
                    comp_gf[j+comps.Size()*i] = comp_gf[j+comps.Size()*i]/all_comp;
                  //   comp_gf[j+comps.Size()*i] = comp_gf[j+comps.Size()*i]/rmass[j];
                  //   rho_gf[j] = rho_gf[j] + z_rho[i]*comp_gf[j+comps.Size()*i];
                    lambda0_gf[j] = lambda0_gf[j] + lambda[i]*comp_gf[j+comps.Size()*i];
                    mu0_gf[j] = mu0_gf[j] + mu[i]*comp_gf[j+comps.Size()*i];
                  }

                  if(dim == 2){s_gf[j+S1.Size()*0]=S1[j]/rmass[j];s_gf[j+S1.Size()*1]=S2[j]/rmass[j];s_gf[j+S1.Size()*2]=S3[j]/rmass[j];}
                  else{s_gf[j+S1.Size()*0]=S1[j]/rmass[j];s_gf[j+S1.Size()*1]=S2[j]/rmass[j];s_gf[j+S1.Size()*2]=S3[j]/rmass[j];s_gf[j+S1.Size()*3]=S4[j]/rmass[j];s_gf[j+S1.Size()*4]=S5[j]/rmass[j];s_gf[j+S1.Size()*5]=S6[j]/rmass[j];}
               }
               
               if(myid==0){std::cout << "remapping for H1" << std::endl;}
               Vector vxyz;
               int point_ordering;
               vxyz = *pmesh_copy->GetNodes(); // from target mesh
               point_ordering = pmesh_copy->GetNodes()->FESpace()->GetOrdering();

               // // Find and Interpolate FE function values on the desired points.
               // // Vector interp_vals(nodes_cnt*tar_ncomp);
               FindPointsGSLIB finder(MPI_COMM_WORLD);
               finder.Setup(*pmesh_old);
               Vector interp_vals(v_gf.Size());
               finder.Interpolate(vxyz, v_gf, interp_vals, point_ordering); for(int i = 0; i < interp_vals.Size(); i++ ){if(interp_vals[i] != 0.0){v_gf[i] = interp_vals[i];}}
               finder.Interpolate(vxyz, u_gf, interp_vals, point_ordering); for(int i = 0; i < interp_vals.Size(); i++ ){if(interp_vals[i] != 0.0){u_gf[i] = interp_vals[i];}}
            }
            delete pmesh_old; 
            delete pmesh_copy;
            delete pmesh_copy_old;

            if (ti == 1)
            {
               x_ini_gf = *pmesh->GetNodes(); // copy optimized initial mesh
            }

            {
               // Compute the geometric parameter at the dofs of each element.
               for (int e = 0; e < pmesh->GetNE(); e++)
               {
                  const FiniteElement *fe = L2FESpace_geometric.GetFE(e);
                  const IntegrationRule &ir = fe->GetNodes();
                  L2FESpace_geometric.GetElementVDofs(e, vdofs);
                  allVals.SetSize(vdofs.Size());
                  for (int q = 0; q < ir.GetNPoints(); q++)
                  {
                     const IntegrationPoint &ip = ir.IntPoint(q);
                     pmesh->GetElementJacobian(e, jacobian, &ip);
                     double sizeVal;
                     Vector asprVals, skewVals, oriVals;
                     pmesh->GetGeometricParametersFromJacobian(jacobian, sizeVal,
                                                            asprVals, skewVals, oriVals);
                     allVals(q + 0) = sizeVal;
                     for (int n = 0; n < nAspr; n++)
                     {
                        allVals(q + (n+1)*ir.GetNPoints()) = asprVals(n);
                     }
                     for (int n = 0; n < nSkew; n++)
                     {
                        allVals(q + (n+1+nAspr)*ir.GetNPoints()) = skewVals(n);
                     }
                  }
                  quality.SetSubVector(vdofs, allVals);
               }

               Vector vol_vec(e_gf.Size());
               Vector skew_vec(e_gf.Size());
               vol_vec = 1.0; skew_vec = 0.0;
               for(int i = 0; i < e_gf.Size(); i++)
               {  
                  vol_vec[i] = quality[i]; 
                  if(quality[e_gf.Size() + i] < 1.0)
                  {
                     quality[e_gf.Size() + i] = 1/quality[e_gf.Size() + i];  
                  }
                  skew_vec[i] = quality[e_gf.Size() + i];
                  skew_ini_gf[i] = quality[e_gf.Size() + i]; // reinitalize conidtion number after remeshing.
               }

               skew_vec.Add(-1.0, skew_ini_gf);
               // skew_vec.Abs();

               double local_min_vol = vol_vec.Min();
               MPI_Reduce(&local_min_vol, &global_min_vol, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
               MPI_Bcast(&global_min_vol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

               if(global_min_vol < 0){MFEM_ABORT("Negative Jacobian (volume) occurs!");}
               
               // double local_max_skew = skew_vec.Max();
               double local_max_skew = std::max( skew_vec.Max(), -skew_vec.Min() );
               double global_max_skew;
               
               MPI_Reduce(&local_max_skew, &global_max_skew, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
               MPI_Bcast(&global_max_skew, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
               
               cond_num = global_max_skew;

               // if(myid == 0)
               // {
               //    cout << "*** after remeshing :  " << cond_num << endl;
               // }
         }

            // if(param.control.winkler_foundation & param.control.winkler_flat)
            // {
             
            //    ParSubMesh::Transfer(x_gf, x_bottom); // update current mesh to submesh
            //    for (int i = 0; i < bottom.Size(); i++){bottom[i] = x_bottom[i+bottom.Size()];}
            //    // double local_sum_bot = bottom.Sum();
            //    // int local_bot_size = bottom.Size();
            //    // double global_sum_bot;
            //    // int global_sum_size;
            //    // MPI_Reduce(&local_sum_bot, &global_sum_bot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            //    // MPI_Reduce(&local_bot_size, &global_sum_size, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            //    // global_sum_bot = global_sum_bot/global_sum_size; // average height of bottom surface
            //    // MPI_Bcast(&global_sum_bot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            //    for (int i = 0; i < bottom.Size(); i++){x_bottom[i+bottom.Size()] = 0.0;}
            //    // for (int i = 0; i < bottom.Size(); i++){x_bottom[i+bottom.Size()] = global_sum_bot;}
            //    submesh_bottom.NewNodes(x_bottom, false);
            //    ParSubMesh::Transfer(x_bottom, x_gf); // update adjusted nodes on top boundary 
            // }

            /*
            for (int i = 0; i < vol_ini_gf.Size(); i++){vol_ini_gf[i] = quality[i];}
            comp_ref_gf = comp_gf;

            Vector skew_vec(e_gf.Size());
            skew_vec = 1.0;
            for(int i = 0; i < e_gf.Size(); i++){skew_vec[i] = quality[e_gf.Size() + i];}

            double local_max_skew = skew_vec.Max();
            double local_min_skew = skew_vec.Min();
            double global_max_skew;
            double global_min_skew;

            MPI_Reduce(&local_max_skew, &global_max_skew, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local_min_skew, &global_min_skew, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_max_skew, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_min_skew, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            cond_num = std::max(global_max_skew, 1.0/global_min_skew);
            */

            // h_min = geo.GetLengthEstimate(S);
            // cond_num = ini_h_min/h_min;
            // if(myid==0){std::cout << "skew number : " << cond_num << std::endl;}
            
            // if(dim ==2) {delete S1; delete S2; delete S3;}
            // else {delete S1; delete S2; delete S3; delete S4; delete S5; delete S6;}
            // if (mesh_changed & param.tmop.amr)
            // {
            //    // update state and operator
            //    TMOPUpdate(S, S_old, offset, x_gf, v_gf, e_gf, s_gf, x_ini_gf, p_gf, n_p_gf, ini_p_gf, u_gf, rho0_gf, lambda0_gf, mu0_gf, mat_gf, dim, param.tmop.amr);
            //    geo.TMOPUpdate(S, true);
            //    pmesh->Rebalance();
            //    TMOPUpdate(S, S_old, offset, x_gf, v_gf, e_gf, s_gf, x_ini_gf, p_gf, n_p_gf, ini_p_gf, u_gf, rho0_gf, lambda0_gf, mu0_gf, mat_gf, dim, param.tmop.amr);
            //    geo.TMOPUpdate(S, false);
            //    ode_solver->Init(geo);
            // }

            if (param.sim.visit)
            {
                  visit_dc.SetCycle(ti);
                  visit_dc.SetTime(t);
                  if(param.sim.year){visit_dc.SetTime(year);}
                  visit_dc.Save();
            }

            if (param.sim.paraview)
            {
                  pd->SetCycle(ti);
                  pd->SetTime(t);
                  if(param.sim.year){pd->SetTime(year);}
                  pd->Save();
            }

            // ti = ti+1;

          }
      }

      // if(mesh_changed){mesh_changed=false; dt_est=dt_est*1e-5;}

      // const double dt_est = geo.GetTimeStepEstimate(S, dt, mpi.Root());
      // const double dt_est = geo.GetTimeStepEstimate(S);

      // if (dt < std::numeric_limits<double>::epsilon())
      // { 
      //    if (param.sim.visit)
      //    {
      //          visit_dc.SetCycle(ti);
      //          visit_dc.SetTime(t);
      //          visit_dc.Save();
      //    }

      //    if (param.sim.paraview)
      //    {
      //          pd->SetCycle(ti);
      //          pd->SetTime(t);
      //          pd->Save();
      //    }

      //    MFEM_ABORT("The time step crashed!"); 
      // }

      if(mesh_changed)
      {
      mesh_changed = false;
      }
      else
      {
         if (dt_est < dt)
         {
         // Repeat (solve again) with a decreased time step - decrease of the
         // time estimate suggests appearance of oscillations.
         // dt *= 0.50;
         dt  = dt_est; 
         // if (dt < std::numeric_limits<double>::epsilon())
         if (dt < 1.0E-38)
         { 
            if (param.sim.visit)
            {
                  visit_dc.SetCycle(ti);
                  visit_dc.SetTime(t);
                  visit_dc.Save();
            }

            if (param.sim.paraview)
            {
                  pd->SetCycle(ti);
                  pd->SetTime(t);
                  pd->Save();
            }

            MFEM_ABORT("The time step crashed!"); 
         }
            t = t_old;
            S = S_old;
            p_gf = p_gf_old; ini_p_gf = ini_p_old_gf;
            // if(surface_diff){x_top=x_top_old; topo=topo_t_old;}
            geo.ResetQuadratureData();
            // if (mpi.Root()) { cout << "Repeating step " << ti << ", dt " << dt/86400/365.25 << std::setprecision(6) << std::scientific << " yr" << endl; }
            if (steps < param.sim.max_tsteps) { last_step = false; }
            ti--; continue;
         }
         else if (dt_est > 1.25 * dt) { dt *= 1.02; }
      }

      // Ensure the sub-vectors x_gf, v_gf, and e_gf know the location of the
      // data in S. This operation simply updates the Memory validity flags of
      // the sub-vectors to match those of S.

      // // Estimate element errors using the Zienkiewicz-Zhu error estimator.
      // if(param.tmop.amr)
      // {
      //    if (myid == 0) { std::cout << "Estimating Error ... " << flush; }
      //    Vector errors(pmesh->GetNE());
      //    // geo.GetErrorEstimates(e_gf, errors);
      //    if (myid == 0) { std::cout << "done." << std::endl; }

      //    double local_max_err = errors.Max();
      //    double global_max_err;
      //    MPI_Allreduce(&local_max_err, &global_max_err, 1,MPI_DOUBLE, MPI_MAX, pmesh.GetComm());

      //    // Refine the elements whose error is larger than a fraction of the
      //    // maximum element error.
      //    const double frac = 0.7;
      //    double threshold = frac * global_max_err;
      //    if (Mpi::Root()) { cout << "Refining ..." << endl; }
      //    pmesh->RefineByError(errors, threshold);

      //    // update state and operator
      //    TMOPUpdate(S, S_old, offset, x_gf, v_gf, e_gf, s_gf, x_ini_gf, p_gf, n_p_gf, ini_p_gf, u_gf, rho0_gf, lambda0_gf, mu0_gf, mat_gf, flattening, dim, param.tmop.amr);
      //    geo.TMOPUpdate(S, true);
      //    pmesh->Rebalance();
      //    TMOPUpdate(S, S_old, offset, x_gf, v_gf, e_gf, s_gf, x_ini_gf, p_gf, n_p_gf, ini_p_gf, u_gf, rho0_gf, lambda0_gf, mu0_gf, mat_gf, flattening, dim, param.tmop.amr);
      //    geo.TMOPUpdate(S, false);
      //    ode_solver->Init(geo);
      // }


      // if(param.control.winkler_foundation & param.control.winkler_flat) {for( int i = 0; i < x_gf.Size(); i++ ){if(flattening[i] > 0.0){x_gf[i] = x_ini_gf[i];}}}
      
      x_gf.SyncAliasMemory(S);
      v_gf.SyncAliasMemory(S);
      e_gf.SyncAliasMemory(S);
      s_gf.SyncAliasMemory(S);

      s_old_gf = s_gf; // storing old Caushy stress

      
      // Adding stress increment to total stress and storing spin rate
      // Make sure that the mesh corresponds to the new solution state. This is
      // needed, because some time integrators use different S-type vectors
      // and the oper object might have redirected the mesh positions to those.

      u_gf.Add(dt, v_gf);
      
      
      {
         // Compute the geometric parameter at the dofs of each element.
         for (int e = 0; e < pmesh->GetNE(); e++)
         {
            const FiniteElement *fe = L2FESpace_geometric.GetFE(e);
            const IntegrationRule &ir = fe->GetNodes();
            L2FESpace_geometric.GetElementVDofs(e, vdofs);
            allVals.SetSize(vdofs.Size());
            for (int q = 0; q < ir.GetNPoints(); q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               pmesh->GetElementJacobian(e, jacobian, &ip);
               double sizeVal;
               Vector asprVals, skewVals, oriVals;
               pmesh->GetGeometricParametersFromJacobian(jacobian, sizeVal,
                                                      asprVals, skewVals, oriVals);
               allVals(q + 0) = sizeVal;
               for (int n = 0; n < nAspr; n++)
               {
                  allVals(q + (n+1)*ir.GetNPoints()) = asprVals(n);
               }
               for (int n = 0; n < nSkew; n++)
               {
                  allVals(q + (n+1+nAspr)*ir.GetNPoints()) = skewVals(n);
               }
            }
            quality.SetSubVector(vdofs, allVals);
         }

         Vector vol_vec(e_gf.Size());
         Vector skew_vec(e_gf.Size());
         vol_vec = 1.0; skew_vec = 0.0;
         for(int i = 0; i < e_gf.Size(); i++)
         {  
            vol_vec[i] = quality[i]; 
            if(quality[e_gf.Size() + i] < 1.0)
            {
               quality[e_gf.Size() + i] = 1/quality[e_gf.Size() + i];  
            }
            skew_vec[i] = quality[e_gf.Size() + i];
         }

         skew_vec.Add(-1.0, skew_ini_gf);
         // skew_vec.Abs();

         double local_min_vol = vol_vec.Min();
         MPI_Reduce(&local_min_vol, &global_min_vol, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
         MPI_Bcast(&global_min_vol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

         if(global_min_vol < 0){MFEM_ABORT("Negative Jacobian (volume) occurs!");}
         
         // double local_max_skew = skew_vec.Max();
         double local_max_skew = std::max( skew_vec.Max(), -skew_vec.Min() );
         double global_max_skew;
         
         MPI_Reduce(&local_max_skew, &global_max_skew, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
         MPI_Bcast(&global_max_skew, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         
         cond_num = global_max_skew;

         // if(myid == 0)
         // {
         //    cout << "*** check condition number :  " << ti << ","<< cond_num << endl;
         // }
      }
      
      
      if (param.tmop.tmop)
      {
         if(param.control.mass_bal && ti > 1)
         {
            geo.TMOPUpdate(S, true); // update mass matrix and density to keep same. 
         }
         else
         {  
            geo.TMOPUpdate(S, false); // update mass matrix and density to keep same. 
         }
      }
      else
      {
         // initialize density (test)
         geo.TMOPUpdate(S, false); // update mass matrix and density to keep same. 
      }
      
      if (last_step || (ti % param.sim.vis_steps) == 0)
      {
         double lnorm = e_gf * e_gf, norm;
         MPI_Allreduce(&lnorm, &norm, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
         if (param.sim.mem_usage)
         {
            mem = GetMaxRssMB();
            MPI_Reduce(&mem, &mmax, 1, MPI_LONG, MPI_MAX, 0, pmesh->GetComm());
            MPI_Reduce(&mem, &msum, 1, MPI_LONG, MPI_SUM, 0, pmesh->GetComm());
         }
         const double internal_energy = geo.InternalEnergy(e_gf);
         const double kinetic_energy = geo.KineticEnergy(v_gf);
         Vector vel_mag(v_gf.Size()/dim);
         for (int i = 0; i < v_gf.Size()/dim; i++)
         {
            if(dim == 2){vel_mag[i] = sqrt(pow(v_gf[i], 2) + pow(v_gf[i+v_gf.Size()/dim], 2));}
            else{vel_mag[i] = sqrt(pow(v_gf[i], 2) + pow(v_gf[i+v_gf.Size()/dim], 2) + pow(v_gf[i+2*v_gf.Size()/dim], 2));}
         }

         double local_max_vel = vel_mag.Max();
         double global_max_vel;

         MPI_Reduce(&local_max_vel, &global_max_vel, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
         MPI_Bcast(&global_max_vel, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

         if(param.sim.year)
         {
            if (mpi.Root())
            {
            const double sqrt_norm = sqrt(norm);

            cout << std::fixed;
            cout << "step " << std::setw(5) << ti
                 << ",\tt = " << std::setw(5) << std::setprecision(4) << t/86400/365.25
                 << ",\tdt (yr) = " << std::setw(5) << std::setprecision(6) << std::scientific << dt/86400/365.25
                 << ",\t|e| = " << std::setw(5) << std::setprecision(3) << std::scientific
                 << sqrt_norm
                 << ", max_vel (cm/yr) = " << std::setw(5) << std::setprecision(3) << std::scientific
                 << global_max_vel*86400*365*100
                 << ", relative max_skew = " << std::setw(5) << std::setprecision(3) << std::scientific
                 << cond_num
                 << ", h_min = " << std::setw(5) << std::setprecision(3) << std::scientific
                 << h_min;
            //  << ",\t|IE| = " << std::setprecision(10) << std::scientific
            //  << internal_energy
            //   << ",\t|KE| = " << std::setprecision(10) << std::scientific
            //  << kinetic_energy
            //   << ",\t|E| = " << std::setprecision(10) << std::scientific
            //  << kinetic_energy+internal_energy;
            cout << std::fixed;
            if (param.sim.mem_usage)
               {
                  cout << ", mem: " << mmax << "/" << msum << " MB";
               }
            cout << endl;
            }
         }
         else
         {
            if (mpi.Root())
            {
            const double sqrt_norm = sqrt(norm);

            cout << std::fixed;
            cout << "step " << std::setw(5) << ti
                 << ",\tt = " << std::setw(5) << std::setprecision(4) << t
                 << ",\tdt (sec) = " << std::setw(5) << std::setprecision(6) << std::scientific << dt
                 << ",\t|e| = " << std::setw(5) << std::setprecision(3) << std::scientific
                 << sqrt_norm
                 << ", max_vel (m/sec) = " << std::setw(5) << std::setprecision(3) << std::scientific
                 << global_max_vel*1
                 << ", relative max_skew = " << std::setw(5) << std::setprecision(3) << std::scientific
                 << cond_num
                 << ", h_min = " << std::setw(5) << std::setprecision(3) << std::scientific
                 << h_min;
            //  << ",\t|IE| = " << std::setprecision(10) << std::scientific
            //  << internal_energy
            //   << ",\t|KE| = " << std::setprecision(10) << std::scientific
            //  << kinetic_energy
            //   << ",\t|E| = " << std::setprecision(10) << std::scientific
            //  << kinetic_energy+internal_energy;
            cout << std::fixed;
            if (param.sim.mem_usage)
               {
                  cout << ", mem: " << mmax << "/" << msum << " MB";
               }
            cout << endl;
            }
         }
         
         // Make sure all ranks have sent their 'v' solution before initiating
         // another set of GLVis connections (one from each rank):
         MPI_Barrier(pmesh->GetComm());

         // if (param.sim.visualization || param.sim.visit || param.sim.gfprint || param.sim.paraview) { geo.ComputeDensity(rho_gf); }
         // if (param.control.mass_bal) { geo.ComputeDensity(rho_gf); }
         // geo.ComputeDensity(rho_gf);
         if (param.sim.visualization)
         {
            int Wx = 0, Wy = 0; // window position
            int Ww = 350, Wh = 350; // window size
            int offx = Ww+10; // window offsets
            if (param.sim.problem != 0 && param.sim.problem != 4)
            {
               geodynamics::VisualizeField(vis_rho, vishost, visport, rho0_gf,
                                             "Density", Wx, Wy, Ww, Wh);
            }
            Wx += offx;
            geodynamics::VisualizeField(vis_v, vishost, visport,
                                          v_gf, "Velocity", Wx, Wy, Ww, Wh);
            Wx += offx;
            geodynamics::VisualizeField(vis_e, vishost, visport, e_gf,
                                          "Specific Internal Energy",
                                          Wx, Wy, Ww,Wh);
            Wx += offx;
         }


         if (param.sim.visit)
         {
            visit_dc.SetCycle(ti);
            visit_dc.SetTime(t);
            if(param.sim.year){visit_dc.SetTime(year);}
            visit_dc.Save();
         }

         if (param.sim.paraview)
         {
            pd->SetCycle(ti);
            pd->SetTime(t);
            if(param.sim.year){pd->SetTime(year);}
            pd->Save();
         }

         if (param.sim.gfprint)
         {
            std::ostringstream mesh_name, rho_name, v_name, e_name;
            mesh_name << param.sim.basename << "_" << ti << "_mesh";
            rho_name  << param.sim.basename << "_" << ti << "_rho";
            v_name << param.sim.basename << "_" << ti << "_v";
            e_name << param.sim.basename << "_" << ti << "_e";

            std::ofstream mesh_ofs(mesh_name.str().c_str());
            mesh_ofs.precision(8);
            pmesh->PrintAsOne(mesh_ofs);
            mesh_ofs.close();

            std::ofstream rho_ofs(rho_name.str().c_str());
            rho_ofs.precision(8);
            rho0_gf.SaveAsOne(rho_ofs);
            rho_ofs.close();

            std::ofstream v_ofs(v_name.str().c_str());
            v_ofs.precision(8);
            v_gf.SaveAsOne(v_ofs);
            v_ofs.close();

            std::ofstream e_ofs(e_name.str().c_str());
            e_ofs.precision(8);
            e_gf.SaveAsOne(e_ofs);
            e_ofs.close();
         }
      }

      // Problems checks
      if (param.sim.check)
      {
         double lnorm = e_gf * e_gf, norm;
         MPI_Allreduce(&lnorm, &norm, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
         const double e_norm = sqrt(norm);
         MFEM_VERIFY(param.mesh.rs_levels==0 && param.mesh.rp_levels==0, "check: rs, rp");
         MFEM_VERIFY(param.mesh.order_v==2, "check: order_v");
         MFEM_VERIFY(param.mesh.order_e==1, "check: order_e");
         MFEM_VERIFY(param.solver.ode_solver_type==4, "check: ode_solver_type");
         MFEM_VERIFY(param.sim.t_final == 0.6, "check: t_final");
         MFEM_VERIFY(param.solver.cfl==0.5, "check: cfl");
         MFEM_VERIFY(param.mesh.mesh_file.compare("default") == 0, "check: mesh_file");
         MFEM_VERIFY(dim==2 || dim==3, "check: dimension");
         Checks(ti, e_norm, checks);
      }
   }
   MFEM_VERIFY(!param.sim.check || checks == 2, "Check error!");

   switch (param.solver.ode_solver_type)
   {
      case 2: steps *= 2; break;
      case 3: steps *= 3; break;
      case 4: steps *= 4; break;
      case 6: steps *= 6; break;
      case 7: steps *= 2;
   }

   geo.PrintTimingData(mpi.Root(), steps, param.sim.fom);

   if (param.sim.mem_usage)
   {
      mem = GetMaxRssMB();
      MPI_Reduce(&mem, &mmax, 1, MPI_LONG, MPI_MAX, 0, pmesh->GetComm());
      MPI_Reduce(&mem, &msum, 1, MPI_LONG, MPI_SUM, 0, pmesh->GetComm());
   }

   const double energy_final = geo.InternalEnergy(e_gf) +
                               geo.KineticEnergy(v_gf);
   if (mpi.Root())
   {
      cout << endl;
      cout << "Energy  diff: " << std::scientific << std::setprecision(2)
           << fabs(energy_init - energy_final) << endl;
      if (param.sim.mem_usage)
      {
         cout << "Maximum memory resident set size: "
              << mmax << "/" << msum << " MB" << endl;
      }
   }

   // Print the error.
   // For problems 0 and 4 the exact velocity is constant in time.
   if (param.sim.problem == 0 || param.sim.problem == 4)
   {
      const double error_max = v_gf.ComputeMaxError(v_coeff),
                   error_l1  = v_gf.ComputeL1Error(v_coeff),
                   error_l2  = v_gf.ComputeL2Error(v_coeff);
      if (mpi.Root())
      {
         cout << "L_inf  error: " << error_max << endl
              << "L_1    error: " << error_l1 << endl
              << "L_2    error: " << error_l2 << endl;
      }
   }

   if (param.sim.visualization)
   {
      vis_v.close();
      vis_e.close();
   }

   // Free the used memory.
   delete ode_solver;
   delete pmesh;
   delete ode_solver_sub;
   delete ode_solver_sub2;

   return 0;
}

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
         {{5, 6.546538624534384e+00}, { 27, 7.588576357792927e+00}},
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
         {{5, 2.041491591302486e+01}, { 59, 3.443180411803796e+01}},
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