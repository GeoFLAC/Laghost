# Proposed Laghost Directory Structure

## Current Issues
- 25+ source files in root directory
- Mixed naming conventions (laghost_* vs remhos_*)
- Unclear module boundaries
- Difficult to navigate and understand dependencies

## Proposed Structure

```
laghost/
├── README.md
├── LICENSE
├── CHANGELOG
├── makefile
├── defaults.cfg
├── .gitignore
├── .travis.yml
│
├── src/                          # Main source code
│   ├── main/
│   │   └── laghost.cpp          # Main driver
│   ├── core/                    # Core solver components
│   │   ├── solver.{hpp,cpp}     # Main solver (laghost_solver)
│   │   ├── assembly.{hpp,cpp}   # Assembly operations (laghost_assembly)
│   │   └── constants.hpp        # Physical/mathematical constants
│   ├── physics/                 # Physics modules
│   │   ├── rheology.{hpp,cpp}   # Material models (laghost_rheology)
│   │   ├── functions.{hpp,cpp}  # Problem functions (laghost_function)
│   │   └── tmop.{hpp,cpp}       # Mesh optimization (laghost_tmop)
│   ├── io/                      # Input/Output
│   │   ├── input.{hpp,cpp}      # Config parsing (laghost_input)
│   │   └── parameters.hpp       # Parameter structures (laghost_parameters)
│   ├── remhos/                  # Remapping functionality
│   │   ├── remhos.{hpp,cpp}     # Main remhos (laghost_remhos)
│   │   ├── fct.{hpp,cpp}        # FCT methods (remhos_fct)
│   │   ├── ho.{hpp,cpp}         # High-order (remhos_ho)
│   │   ├── lo.{hpp,cpp}         # Low-order (remhos_lo)
│   │   ├── mono.{hpp,cpp}       # Monotonic (remhos_mono)
│   │   ├── sync.{hpp,cpp}       # Synchronization (remhos_sync)
│   │   └── tools.{hpp,cpp}      # Utilities (remhos_tools)
│   └── utils/                   # Utilities
│       ├── array2d.hpp          # 2D array utilities
│       └── mesh_optimizer.hpp   # Mesh optimization utilities
│
├── include/                     # Public headers (if needed for library use)
│   └── laghost.hpp             # Main public header
│
├── external/                    # External dependencies/modifications
│   ├── mfem/                   # MFEM modifications
│   │   ├── vector.hpp
│   │   └── vector.cpp
│   └── common/                 # Shared MFEM miniapp utilities
│       ├── mfem-common.hpp
│       ├── fem_extras.{hpp,cpp}
│       ├── mesh_extras.{hpp,cpp}
│       ├── pfem_extras.{hpp,cpp}
│       ├── dist_solver.{hpp,cpp}
│       ├── makefile
│       └── CMakeLists.txt
│
├── variants/                   # Alternative implementations (REMOVED)
│   # Note: serial/ and amr/ directories have been removed as they were obsolete
│
├── data/                      # Mesh files and test data
│   ├── meshes/               # Mesh files
│   │   ├── 1d/
│   │   ├── 2d/
│   │   └── 3d/
│   ├── reference/            # Reference solutions
│   └── images/               # Visualization images
│
├── tests/                    # All testing
│   ├── unit/                # Unit tests
│   │   ├── test_*.cpp
│   │   ├── unit_tests.cpp
│   │   ├── Makefile
│   │   ├── CMakeLists.txt
│   │   └── README.md
│   ├── integration/         # Integration tests
│   │   ├── verification/    # Verification problems
│   │   ├── performance/     # Performance benchmarks
│   │   └── regression/      # Regression tests
│   └── benchmarks/          # Scientific benchmarks
│       ├── elastic_column/
│       ├── plastic_oedometer_test/
│       └── viscoplastic_simple_shear_test/
│
├── docs/                     # Documentation
│   ├── user_guide/
│   ├── developer_guide/
│   ├── api/                 # API documentation
│   └── examples/            # Usage examples
│
├── scripts/                  # Build and utility scripts
│   ├── build/               # Build scripts
│   ├── setup/               # Setup scripts
│   └── ci/                  # CI/CD scripts
│
├── config/                   # Configuration files
│   ├── defaults/            # Default configurations
│   ├── examples/            # Example configurations
│   └── templates/           # Configuration templates
│
└── .kiro/                    # Kiro IDE configuration
    └── steering/
        ├── product.md
        ├── tech.md
        └── structure.md
```

## Benefits of Proposed Structure

### 1. **Clear Module Organization**
- **src/core/**: Essential solver components
- **src/physics/**: Physics-specific modules
- **src/io/**: Input/output handling
- **src/remhos/**: Remapping functionality (clearly separated)
- **src/utils/**: Utility functions

### 2. **Improved Maintainability**
- Related files grouped together
- Clear dependency hierarchy
- Easier to locate specific functionality
- Better separation of concerns

### 3. **Enhanced Testing**
- **tests/unit/**: Unit tests for individual components
- **tests/integration/**: Full workflow testing
- **tests/benchmarks/**: Scientific validation cases
- Clear separation of test types

### 4. **Better Documentation**
- Centralized documentation in **docs/**
- API documentation separate from user guides
- Examples and tutorials organized

### 5. **Cleaner Root Directory**
- Only essential files at top level
- Configuration files in dedicated directory
- Scripts organized by purpose

### 6. **External Dependencies**
- **external/**: Clear separation of external code
- **variants/**: Alternative implementations isolated
- No mixing of core and external code

## Migration Strategy

### Phase 1: Core Reorganization
1. Create new directory structure
2. Move core source files to **src/core/**
3. Move physics modules to **src/physics/**
4. Update makefiles and includes

### Phase 2: Module Separation
1. Move remhos files to **src/remhos/**
2. Organize I/O components in **src/io/**
3. Create utilities directory **src/utils/**

### Phase 3: Testing Reorganization
1. Separate unit and integration tests
2. Move benchmarks to **tests/benchmarks/**
3. Update test build system

### Phase 4: Documentation and Scripts
1. Create **docs/** structure
2. Move scripts to **scripts/**
3. Organize configuration files

## Implementation Considerations

### Build System Updates
- Update makefile include paths
- Modify CMakeLists.txt for new structure
- Update CI/CD scripts

### Header Dependencies
- Update #include paths in source files
- Create convenience headers if needed
- Maintain backward compatibility during transition

### Testing
- Ensure all tests pass after reorganization
- Update test file paths
- Verify CI/CD pipeline works

### Documentation
- Update README files
- Create migration guide
- Update developer documentation