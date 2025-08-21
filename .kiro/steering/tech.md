# Technology Stack

## Core Dependencies
- **MFEM**: Modular parallel C++ finite element library (v4.5+)
- **HYPRE**: Parallel linear algebra and solvers (v2.11.2+)
- **METIS**: Graph partitioning for domain decomposition (v4.0.3)
- **Boost**: Program options library (v1.42+)
- **GSLIB**: High-order interpolation library
- **MPI**: Parallel communication

## Build System
- **Make-based**: Uses MFEM's configuration system
- **Compiler**: Inherits from MFEM (typically GCC/Intel with MPI wrappers)
- **GPU Support**: CUDA backend available via MFEM

## Common Build Commands

### Setup Dependencies
```bash
# Build all dependencies automatically
make setup

# Build with CUDA support
make setup MFEM_BUILD=pcuda

# Manual dependency build order
make hypre
make metis  
make mfem
```

### Build Laghost
```bash
# Standard parallel build
make -j 4

# Debug build (modify makefile CCC flags to include -g)
make clean && make -j 4

# Check build status
make status
```

### Testing & Validation
```bash
# Quick test
make test

# Full test suite
make tests

# Verification checks
make checks

# Performance benchmarks (various problem sizes/ranks)
make 1  # 1 MPI rank
make 4  # 4 MPI ranks
```

## Configuration
- **Input files**: `.cfg` format using Boost program options
- **Default config**: `defaults.cfg`
- **Runtime options**: Command-line parameters override config file settings

## Compilation Flags
- Inherits MFEM's compiler and linker flags
- Links against: `$(MFEM_LIBS) $(MFEM_EXT_LIBS) -lboost_program_options`
- Debug builds: Add `-O0 -g` to CXXFLAGS