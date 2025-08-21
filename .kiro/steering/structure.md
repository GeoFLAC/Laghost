# Project Structure

## Source Code Organization

### `/src/main/`
- **`laghost.cpp`**: Main driver with time integration loop

### `/src/core/`
Core solver components:
- **`laghost_solver.{hpp,cpp}`**: `LagrangianGeoOperator` class - core solver logic
- **`laghost_assembly.{hpp,cpp}`**: Force/mass matrix assembly (full/partial)
- **`laghost_constants.hpp`**: Physical and mathematical constants

### `/src/physics/`
Physics modules:
- **`laghost_rheology.{hpp,cpp}`**: Material models and constitutive relations
- **`laghost_function.{hpp,cpp}`**: Problem-specific functions and coefficients
- **`laghost_tmop.{hpp,cpp}`**: Mesh optimization (TMOP) functionality

### `/src/io/`
Input/Output handling:
- **`laghost_parameters.hpp`**: Parameter structures and definitions
- **`laghost_input.{hpp,cpp}`**: Configuration file parsing

### `/src/remhos/`
Remapping functionality:
- **`laghost_remhos.{hpp,cpp}`**: Main remapping interface
- **`remhos_*.{hpp,cpp}`**: Specialized remapping algorithms (FCT, HO, LO, etc.)

### `/src/utils/`
Utility components:
- **`array2d.hpp`**: 2D array utilities
- **`mesh-optimizer.hpp`**: Mesh optimization utilities

## Key Directories

### `/external/`
External dependencies and modifications:
- **`external/common/`**: Shared MFEM miniapp utilities
  
  *(Note: external/mfem/ directory removed as obsolete)*

### `/tests/`
Testing framework:
- **`tests/unit/`**: Unit tests for individual components
- **`tests/benchmarks/`**: Scientific validation benchmarks
- **`tests/integration/`**: Integration and workflow tests

### `/config/`
Configuration files:
- **`config/defaults/`**: Default configuration files

### `/data/`
Data files:
- **`data/meshes/`**: Mesh files for testing and examples

### Alternative Implementations (REMOVED)
*(Note: serial/ and amr/ directories have been removed as they were obsolete)*

## Build Configuration
- **`makefile`**: Main build configuration with new source paths
- **`.travis.yml`**: CI configuration

## Code Organization Patterns

### Finite Element Spaces
- **H1 (continuous)**: Position and velocity fields
- **L2 (discontinuous)**: Energy, stress, material properties
- **Mixed formulation**: Kinematic (H1) + thermodynamic (L2) spaces

### Assembly Methods
- **Full Assembly (`-fa`)**: Traditional sparse matrix approach
- **Partial Assembly (`-pa`)**: Matrix-free, tensor-product optimized

### Solver Structure
1. **Time Integration**: Explicit Runge-Kutta methods (RK2Avg default)
2. **Spatial Discretization**: High-order finite elements
3. **Linear Solvers**: CG with appropriate preconditioners
4. **Quadrature**: Tensor-product rules for efficiency

### Naming Conventions
- **Files**: `laghost_<module>.{hpp,cpp}` pattern
- **Classes**: CamelCase (e.g., `LagrangianGeoOperator`)
- **Functions**: CamelCase for methods, snake_case for some utilities
- **Variables**: Mixed conventions, often descriptive names