            __                __               __ 
           / /   ____ _____ _/ /_  ____  _____/ /_
          / /   / __ `/ __ `/ __ \/ __ \/ ___/ __/
         / /___/ /_/ / /_/ / / / / /_/ (__  ) /_  
        /_____/\__,_/\__, /_/ /_/\____/____/\__/  
                    /____/                        
        Lagrangian High-order Solver for Tectonics

## Purpose

**Laghost** (LAGrangian High-Order Solver for Tectnoics) solves the
time-dependent momentum balance of geological media in a moving
Lagrangian frame using unstructured high-order finite element spatial
discretization and explicit high-order time-stepping.
[MFEM](http://mfem.org) is a modular parallel C++ library to enable high-performance scalable finite element discretization. 
LAGHOST extends the capabilities of the [Laghos](https://github.com/CEED/Laghos) (Lagrangian High-Order Solver) one of mini-apps of MFEM, which solves 
the time-dependent Euler equations of compressible gas dynamics in a moving Lagrangian frame 
using high-order finite element spatial discretization and explicit time-stepping (Runge-Kutta method).

> Veselin A. Dobrev, Tzanio V. Kolev, and Robert N. Riebenn <br>
> [High-order curvilinear finite element methods for Lagrangian hydrodynamics](https://doi.org/10.1137/120864672) <br>
> *SIAM Journal on Scientific Computing*, (34) 2012, pp. B606–B641.

> Robert W. Anderson, Veselin A. Dobrev, Tzanio V. Kolev, Robert N. Rieben, and Vladimir Z. <br>
> [High-Order Multi-Material ALE Hydrodynamics](https://doi.org/10.1137/17M1116453) <br>
> *Computational Methods in Science and Engineering*, (40) 2018.

<!-- Laghos captures the basic structure of many compressible shock hydrocodes,
including the [BLAST code](http://llnl.gov/casc/blast) at [Lawrence Livermore
National Laboratory](http://llnl.gov). The miniapp is built on top of a general
discretization library, [MFEM](http://mfem.org), thus separating the pointwise
physics from finite element and meshing concerns.

The Laghos miniapp is part of the [CEED software suite](http://ceed.exascaleproject.org/software),
a collection of software benchmarks, miniapps, libraries and APIs for
efficient exascale discretizations based on high-order finite element
and spectral element methods. See http://github.com/ceed for more
information and source code availability.

The CEED research is supported by the [Exascale Computing Project](https://exascaleproject.org/exascale-computing-project)
(17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
organizations (Office of Science and the National Nuclear Security
Administration) responsible for the planning and preparation of a
[capable exascale ecosystem](https://exascaleproject.org/what-is-exascale),
including software, applications, hardware, advanced system engineering and early
testbed platforms, in support of the nation’s exascale computing imperative. -->

## Characteristics

The problem that Laghost is solving is formulated as a big (block) system of
ordinary differential equations (ODEs) for the unknown (high-order) velocity,
internal energy, stress and mesh nodes (position). The left-hand side of this system of
ODEs is controlled by *mass matrices* (one for velocity and one for energy and stress),
while the right-hand side is constructed from a *force matrix*.

Laghost supports two options for deriving and solving the ODE system, namely the
*full assembly* and the *partial assembly* methods. Partial assembly is the main
algorithm of interest for high orders. For low orders (e.g. 2nd order in 3D),
both algorithms are of interest.

The full assembly option relies on constructing and utilizing global mass and
force matrices stored in compressed sparse row (CSR) format.  In contrast, the
[partial assembly](http://ceed.exascaleproject.org/ceed-code) option defines
only the local action of those matrices, which is then used to perform all
necessary operations. As the local action is defined by utilizing the tensor
structure of the finite element spaces, the amount of data storage, memory
transfers, and FLOPs are lower (especially for higher orders).

The mother code, Laghos, implementation includes support for hardware devices, such
as GPUs, and programming models, such as CUDA, OCCA, RAJA and OpenMP,
based on [MFEM](http://mfem.org), version 4.1 or later. These device
backends are selectable at runtime, see the `-d/--device` command-line
option. So, Laghost share those capability, however, they are not tested enough yet.

Other computational motives in Laghost include the following:

- Support for unstructured meshes, in 2D and 3D, with quadrilateral and
  hexahedral elements (triangular and tetrahedral elements can also be used, but
  with the less efficient full assembly option). Serial and parallel mesh
  refinement options can be set via a command-line flag.
- Explicit time-stepping loop with a specialized Runge-Kutta method of order 2 
  that ensures exact energy conservation on fully discrete level (RK2Avg).
- Continuous and discontinuous high-order finite element discretization spaces
  of runtime-specified order.
- Moving (high-order) meshes.
- Separation between the assembly and the quadrature point-based computations.
- Point-wise definition of mesh size, time-step estimate and artificial
  viscosity coefficient.
- Constant-in-time velocity mass operator that is inverted iteratively on
  each time step. This is an example of an operator that is prepared once (fully
  or partially assembled), but is applied many times. The application cost is
  dominant for this operator.
- Time-dependent force matrix that is prepared every time step (fully or
  partially assembled) and is applied just twice per "assembly". Both the
  preparation and the application costs are important for this operator.
- Domain-decomposed MPI parallelism.
- Optional in-situ visualization with [GLVis](http:/glvis.org) and data output
  for visualization and data analysis with [VisIt](http://visit.llnl.gov) and [ParaView](https://www.paraview.org/).
- Rock rhelogies : Compressible elastic medium, Mohr-Coulomb rate-independnt and rate-independent plasticity, 
  plastic softening based on accumulated plastic strain for cohesion, friction coefficient, and dilation coefficient.
- Mass scaling for *mass matrices* to achieve year-length time step size.
- Dynamic relaxation (a.k.a. Cundall's damping).
- Enabling the application of a Winkler foundation or spring boundary condition for the bottom boundary.
- Multi-material tracking based on composition field
- Remeshing and improving the quality of high-order finite element meshes based on the TMOP (Target-Matrix Optimization Paradigm)
- Remapping high-order continuous (velocity and mesh nodes) and discontinous variables (energy, stress, composition, plastic strain) 
  from source mesh (before remeshing) to new mesh (after remeshing) 
  using [GSLIB](https://mfem.org/howto/findpts/) and [Remhos](https://github.com/CEED/Remhos).
- Input file system (default.cfg) based on boost library (1.42 or newer version).

## Main Code Structure

- The file `laghost.cpp` contains the main driver with the time integration loop.
- In each time step, the ODE system of interest is constructed and solved by
  the class `LagrangianGeoOperator`, defined in `laghost.cpp`
  and implemented in files `laghost_solver.hpp` and `laghost_solver.cpp`.
- In `LagrangianGeoOperator::RK2AvgSolver::Step`, `UpdateMesh`, `SolveVelocity`, `SolveEnergy`, and `SolveStress`
  are sequentially called.
- All quadrature-based computations are performed in the function
  `LagrangianGeoOperator::UpdateQuadratureData` in `laghost_solver.cpp`.
- In `UpdateQuadratureData`, total stress and stress increment based on objective stress rate (Jaumann stress rate)
  are calculated to construct work matrix `F_ij` (force x length; i and j for continous and discontinous space).
- In `SolveVelocity`, a vector,`rhs`, is assembled by multiplying the work matrix `F_ij` and the unity vector of the discontinuous space. 
  Then, taking the negative sign on the `rhs` vector and adding damping force, which is stored in a new vector based on the current force vector, the `rhs`.
- Depending on the chosen option (`-pa` for partial assembly or `-fa` for full
  assembly), the function `LagrangianGeoOperator::Mult` uses the corresponding
  method to construct and solve the final ODE system.
- The full assembly computations for all mass matrices are performed by the MFEM
  library, e.g., classes `MassIntegrator` and `VectorMassIntegrator`.  Full
  assembly of the ODE's right hand side is performed by utilizing the class
  `ForceIntegrator` defined in `laghost_assembly.hpp`.
- The partial assembly computations are performed by the classes
  `ForcePAOperator` and `MassPAOperator` defined in `laghost_assembly.hpp`.
- When partial assembly is used, the main computational kernels are the
  `Mult*` functions of the classes `MassPAOperator` and `ForcePAOperator`
  implemented in file `laghost_assembly.cpp`. These functions have specific
  versions for quadrilateral and hexahedral elements.
- The orders of the velocity and position (continuous kinematic space)
  and the internal energy, stress, composition and plastic strain 
  (discontinuous thermodynamic space) are given by the `-ok` and `-ot` input parameters, respectively.

## Building

Laghost has the following external dependencies:

-  Working MPI compiler
-  hypre, used for parallel linear algebra<br>
   https://github.com/hypre-space/hypre
-  METIS, used for parallel domain decomposition<br>
   https://github.com/KarypisLab/METIS
-  boost-program-options, used for input file system<br>
   https://www.boost.org/
-  GSLIb for  
-  MFEM, core library for arbitrary-order finite elements<br>
   https://github.com/mfem/mfem

### Build boost:

```sh
apt install libboost-program-options-dev
```

or

```sh
~> tar -zxvf boost_1_84_0.tar.gz
~> cd boost_1_84_0/src/
~/boost_1_84_0/src> ./bootstrap.sh
~/boost_1_84_0/src> ./b2 --with-program_options -q
~/boost_1_84_0/src> cd ..
```

### hypre and METIS 

Put them on the same level as the `Laghost` directory: e.g.,
```sh
~> ls
Laghost/  hypre-2.11.2.tar.gz  metis-4.0.3.tar.gz
```

#### Build hypre: e.g.,

```sh
git clone https://github.com/hypre-space/hypre
cd hypre/src
./configure --disable-fortran
make -j
```

or

```sh
~> tar -zxvf hypre-2.11.2.tar.gz
~> cd hypre-2.11.2/src/
~/hypre-2.11.2/src> ./configure --disable-fortran
~/hypre-2.11.2/src> make -j
~/hypre-2.11.2/src> cd ../..
```
For large runs (problem size above 2 billion unknowns), add the
`--enable-bigint` option to the above `configure` line.

#### Build METIS:

```sh
git clone https://github.com/KarypisLab/GKlib
cd GKlib
make config prefix=./
make install
cd ..
git clone https://github.com/KarypisLab/METIS
cd METIS
make config cc=mpicc gklib_path=../GKlib/build/Linux-x86_64 prefix=./
make install
```

This build is optional, as MFEM can be build without METIS by specifying
`MFEM_USE_METIS = NO` below.

<!--```sh
~> tar -zxvf metis-4.0.3.tar.gz
~> cd metis-4.0.3
~/metis-4.0.3> make
~/metis-4.0.3> cd ..
~> ln -s metis-4.0.3 metis-4.0
```-->

### Build GSLIB:

```sh
~> git clone https://github.com/CEED/GSLIB.git
~> cd GSLIB
~/GSLIB> make CC=mpicc
~/GSLIB> cd ..
~> ln -s GSLIB gslib
```

### Build MFEM

Clone and build the parallel version of MFEM:
```sh
~> git clone https://github.com/mfem/mfem.git ./mfem
~> cd mfem/
~/mfem> git checkout master
~/mfem> cp ../Laghost/mfem_modification/vector* ./linalg/
~/mfem> make parallel -j MFEM_USE_GSLIB=YES
~/mfem> cd ..
```

Clone and build the cuda version of MFEM:
```sh
~> git clone https://github.com/mfem/mfem.git ./mfem
~> cd mfem/
~/mfem> git checkout master
~/mfem> cp ../Laghost/mfem_modification/vector* ./linalg/
~/mfem> make pcuda -j MFEM_USE_GSLIB=YES
~/mfem> cd ..
```

The above uses the `master` branch of MFEM.
See the [MFEM building page](http://mfem.org/building/) for additional details.

<!-- (Optional) Clone and build GLVis:
```sh
~> git clone https://github.com/GLVis/glvis.git ./glvis
~> cd glvis/
~/glvis> make
~/glvis> cd ..
```
The easiest way to visualize Laghost results is to have GLVis running in a
separate terminal. Then the `-vis` option in Laghos will stream results directly
to the GLVis socket.
 -->
 
### Build Laghost

```sh
~> cd Laghost/
~/Laghost> make -j
```

<!-- This can be followed by `make test` and `make install` to check and install the
build respectively. See `make help` for additional options.

See also the `make setup` target that can be used to automated the
download and building of hypre, METIS and MFEM. -->

<!--## Versions

In addition to the main MPI-based CPU implementation in https://github.com/CEED/Laghost,
the following versions of Laghost have been developed

<!-- - **SERIAL** version in the [serial/](./serial/README.md) directory.
- **AMR** version in the [amr/](./amr/README.md) directory.
  This version supports dynamic adaptive mesh refinement.
 -->
-->

## Running

**TBA**

## Contact

You can reach the Laghost team by emailing slee29@memphis.edu or sungho91123@gmail.com or by leaving a
comment in the [issue tracker](https://github.com/GeoFLAC/Laghost/issues).

## Copyright

The following copyright applies to each file in the CEED software suite,
unless otherwise stated in the file:

> Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at the
> Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights reserved.

See files LICENSE and NOTICE for details.
