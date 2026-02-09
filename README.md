            __                __               __ 
           / /   ____ _____ _/ /_  ____  _____/ /_
          / /   / __ `/ __ `/ __ \/ __ \/ ___/ __/
         / /___/ /_/ / /_/ / / / / /_/ (__  ) /_  
        /_____/\__,_/\__, /_/ /_/\____/____/\__/  
                    /____/                        
        Lagrangian High-order Solver for Tectonics

[![Unit Tests](https://github.com/GeoFLAC/Laghost/workflows/Unit%20Tests/badge.svg)](https://github.com/GeoFLAC/Laghost/actions?query=workflow%3A%22Unit+Tests%22)
[![CI](https://github.com/GeoFLAC/Laghost/workflows/Laghost%20CI/badge.svg)](https://github.com/GeoFLAC/Laghost/actions?query=workflow%3A%22Laghost+CI%22)
[![CI Fast](https://github.com/GeoFLAC/Laghost/workflows/Laghost%20CI%20-%20Fast/badge.svg)](https://github.com/GeoFLAC/Laghost/actions?query=workflow%3A%22Laghost+CI+-+Fast%22)
[![Comprehensive CI](https://github.com/GeoFLAC/Laghost/workflows/Comprehensive%20CI/badge.svg)](https://github.com/GeoFLAC/Laghost/actions?query=workflow%3A%22Comprehensive+CI%22)

## Purpose

**Laghost** (LAGrangian High-Order Solver for Tectnoics) solves the
time-dependent momentum balance of geological media in a moving
Lagrangian frame using unstructured high-order finite element spatial
discretization and explicit high-order time-stepping.

Laghost extends the capabilities of the [Laghos](https://github.com/CEED/Laghos) (Lagrangian High-Order Solver) one of mini-apps of [MFEM](http://mfem.org), a modular parallel C++ library to enable high-performance scalable finite element discretization. Laghos solves 
the time-dependent Euler equations of compressible gas dynamics in a moving Lagrangian frame 
using high-order finite element spatial discretization and explicit time-stepping (Runge-Kutta method).
Laghost inherits most of these features.

> Veselin A. Dobrev, Tzanio V. Kolev, and Robert N. Riebenn <br>
> [High-order curvilinear finite element methods for Lagrangian hydrodynamics](https://doi.org/10.1137/120864672) <br>
> *SIAM Journal on Scientific Computing*, (34) 2012, pp. B606–B641.

> Robert W. Anderson, Veselin A. Dobrev, Tzanio V. Kolev, Robert N. Rieben, and Vladimir Z. <br>
> [High-Order Multi-Material ALE Hydrodynamics](https://doi.org/10.1137/17M1116453) <br>
> *Computational Methods in Science and Engineering*, (40) 2018.

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

Like the parent code, Laghos, Laghost can support, in principle, hardware devices, such
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
- Data output for visualization and data analysis with [VisIt](http://visit.llnl.gov) and [ParaView](https://www.paraview.org/).
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

The parallel build of MFEM has the following external dependencies:

-  MPI compiler
-  hypre: https://github.com/hypre-space/hypre
-  METIS: See below

Laghost has these additional dependencies:
-  GSLIb, used for remeshing. See below
-  MFEM, core library for arbitrary-order finite elements<br>
   https://github.com/mfem/mfem
-  boost-program-options, used for input file system<br>
   https://www.boost.org/

The MFEM library has a serial and an MPI-based parallel version, which largely
share the same code base. The only prerequisite for building the serial version
of MFEM is a (modern) C++ compiler, such as g++. The parallel version of MFEM
requires an MPI C++ compiler, hypre and METIS.

### Clone MFEM

```sh
git clone https://github.com/mfem/mfem.git
```

#### Build hypre

hypre is expected to be on the same level as the `Laghost` directory: e.g.,

**Option 1: From Git Repository**
```sh
$ ls
mfem
$ git clone https://github.com/hypre-space/hypre
$ ls
hypre  mfem
$ cd hypre/src
$ ./configure --disable-fortran
$ make -j
```

**Option 2: From Tarball**
```sh
$ ls
mfem
$ wget https://github.com/hypre-space/hypre/archive/v2.28.0.tar.gz
$ tar -xzf v2.28.0.tar.gz
$ ls
hypre-2.28.0  mfem  v2.28.0.tar.gz
$ cd hypre-2.28.0/src
$ ./configure --disable-fortran
$ make -j
$ cd ../..
$ ln -sf hypre-2.28.0 hypre
$ ls
hypre  hypre-2.28.0  mfem  v2.28.0.tar.gz
```

After building, ensure the `hypre` symbolic link points to the source directory so MFEM can find it.

#### Build METIS

From [mfem INSTALL document](https://github.com/mfem/mfem/blob/master/INSTALL):

- METIS (a family of multilevel partitioning algorithms)
  https://github.com/mfem/tpls

  Note: We recommend using the MFEM third-party libraries mirror at 
  https://github.com/mfem/tpls/raw/gh-pages/metis-4.0.3.tar.gz because the original METIS
  webpage, http://glaros.dtc.umn.edu/gkhome/metis/metis/overview, is often down
  and we don't support yet the new repo https://github.com/KarypisLab/METIS.

- Follow https://mfem.org/building/#parallel-mpi-version-of-mfem
  ```sh
  $ ls
  hypre  mfem
  $ git clone https://github.com/mfem/tpls.git mfem-tpls
  $ ls
  hypre  mfem  mfem-tpls
  $ cd mfem-tpls
  $ tar zxvf metis-4.0.3.tar.gz
  $ cd metis-4.0.3
  $ make OPTFLAGS=-Wno-error=implicit-function-declaration
  $ cd ../..
  $ ln -s mfem-tpls/metis-4.0.3 metis-4.0
  ```
- This build is optional but recommended.

### Build GSLIB

From [mfem INSTALL document](https://github.com/mfem/mfem/blob/master/INSTALL):

>  GSLIB (optional), used when MFEM_USE_GSLIB = YES. The gslib library must be
>  built prior to the MFEM build, as follows: download gslib-1.0.9, untar it at
>  the same level as MFEM and create a symbolic link: "ln -s gslib-1.0.9 gslib".
>  Build gslib in parallel or in serial based on the desired MFEM build: "make
>  clean; make CC=mpicc" or "make clean; make CC=gcc MPI=0". Build MFEM with
>  MFEM_USE_GSLIB=YES.
  
- URL: https://github.com/Nek5000/gslib.git

Follow the above instruction. The whole process might be as follows:

```sh
$ ls
hypre  metis-4.0  mfem  mfem-tpls 
$ git clone https://github.com/Nek5000/gslib.git
$ cd gslib
$ make CC=mpicc
```

### Build MFEM

Build the parallel version of MFEM:

```sh
$ ls
gslib-1.0.9  gslib  hypre  metis-4.0  mfem  mfem-tpls
$ cd mfem
$ make parallel -j MFEM_USE_GSLIB=YES MFEM_USE_METIS_5=NO
```

To build the cuda version of MFEM:

```sh
$ make pcuda -j MFEM_USE_GSLIB=YES MFEM_USE_METIS_5=NO
```

The above uses the `master` branch of MFEM.
See the [MFEM building page](http://mfem.org/building/) for additional details.

### Install boost:

```sh
apt install libboost-program-options-dev
```

Or download a release package and install it locally: e.g.,

```sh
$ tar xzvf boost_1_88_0.tar.gz
$ cd boost_1_88_0
$ ./bootstrap.sh
$ ./b2 --with-program_options -q
```

### Clone Laghost

```sh
$ git clone https://github.com/GeoFLAC/Laghost.git
$ ls
Laghost  gslib-1.0.9  gslib  hypre  metis-4.0  mfem  mfem-tpls
```

### Build Laghost

```sh
$ cd Laghost/
$ make -j
```

If `libboost-program-options.so` is locally installed, specify its location as follows:

```sh
make -j PROGRAMOPTIONS_LIBDIR=../boost_1_88_0/stage/lib
```

<!-- This can be followed by `make test` and `make install` to check and install the
build respectively. See `make help` for additional options.

See also the `make setup` target that can be used to automated the
download and building of hypre, METIS and MFEM. -->

<!--## Versions

The main MPI-based CPU implementation is the primary version of Laghost.

<!-- Note: The serial/ and amr/ directories have been removed as they were obsolete -->

<!-- - **SERIAL** version in the [serial/](./serial/README.md) directory.
- **AMR** version in the [amr/](./amr/README.md) directory.
  This version supports dynamic adaptive mesh refinement.
 -->

## Running Laghost

```sh
./laghost 
```
Parameters in `defaults.cfg` will be used.

```sh
mpirun -np 8 ./laghost -i ./input_parameters.cfg
```
to use a user-provided input file, `input_parameters.cfg` and run laghost on 8 cores.

For other available command-line options, 

```sh
./laghost -h
```

## Visualizing Laghost output

Use ParaView to load `results/Laghost/Laghost.pvd`

## Contact

Leave a comment or ask a question in the [issue tracker](https://github.com/GeoFLAC/Laghost/issues).

## Copyright

The following copyright applies to each file in the CEED software suite,
unless otherwise stated in the file:

> Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at the
> Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights reserved.

See files LICENSE and NOTICE for details.
