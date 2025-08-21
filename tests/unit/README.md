# Laghost Unit Tests

This directory contains comprehensive unit tests for the Laghost computational geodynamics solver.

## Overview

The test suite covers all core components of Laghost:

- **Parameters**: Configuration structures and constants validation
- **Input**: Configuration file parsing and command-line argument handling  
- **Functions**: Mathematical functions, coefficients, and problem-specific implementations
- **Rheology**: Material models and constitutive relations (elastic, plastic, viscoplastic)
- **Assembly**: Finite element assembly operations (full and partial assembly)
- **Solver**: Main solver components, time integration, and energy computations

## Requirements

### Dependencies
- **Google Test**: C++ testing framework
- **MFEM**: Finite element library (v4.5+)
- **HYPRE**: Parallel linear algebra (v2.11.2+)
- **METIS**: Graph partitioning (v4.0.3, optional)
- **Boost**: Program options library (v1.42+)
- **MPI**: Parallel communication

### Installation (Ubuntu/Debian)
```bash
# Install Google Test
sudo apt-get install libgtest-dev cmake
cd /usr/src/gtest
sudo cmake CMakeLists.txt
sudo make
sudo cp *.a /usr/lib

# Install Boost
sudo apt-get install libboost-program-options-dev
```

## Building and Running Tests

### Method 1: Using Make (Recommended)
```bash
# Build and run tests
make test

# Just build
make all

# Clean
make clean
```

### Method 2: Using CMake
```bash
# Create build directory
mkdir build && cd build

# Configure
cmake ..

# Build
make

# Run tests
./run_tests

# Or use CTest
ctest
```

## Test Structure

### Core Test Files

- **`unit_tests.cpp`**: Main test runner with MPI initialization
- **`test_parameters.cpp`**: Parameter structures and constants
- **`test_input.cpp`**: Configuration parsing and validation
- **`test_function.cpp`**: Mathematical functions and coefficients
- **`test_rheology.cpp`**: Material models and return mapping algorithms
- **`test_assembly.cpp`**: Finite element assembly operations
- **`test_solver.cpp`**: Main solver components and time integration

### Test Categories

#### 1. Parameter Tests
- Default initialization validation
- Boundary condition parsing
- Material parameter parsing
- Constants verification

#### 2. Input Tests  
- Configuration file parsing
- Command-line argument processing
- String parsing for arrays and boundary conditions
- Error handling for invalid inputs

#### 3. Function Tests
- Basic mathematical functions (e0, p0, rho0, etc.)
- Vector functions (v0, xyz0)
- Coefficient evaluations (PlasticCoefficient, LithostaticCoefficient, etc.)
- Principal stress calculations

#### 4. Rheology Tests
- Elastic behavior validation
- Plastic yielding detection
- Viscoplastic behavior
- Material property validation
- Return mapping algorithms (2D and 3D)

#### 5. Assembly Tests
- Quadrature data initialization
- Integrator creation and usage
- Partial assembly operators (Force, Stress, Mass)
- Operator dimension consistency
- Matrix-vector operations

#### 6. Solver Tests
- LagrangianGeoOperator initialization
- Time step estimation
- Energy computations (kinetic and internal)
- Time integration schemes (RK2Avg)
- Coefficient evaluations

## Running Specific Tests

```bash
# Run all tests
./run_tests

# Run tests with filter
./run_tests --gtest_filter="ParametersTest.*"

# Run with verbose output
./run_tests --gtest_verbose

# List all tests
./run_tests --gtest_list_tests
```

## Test Configuration

### Environment Variables
- `MFEM_DIR`: Path to MFEM installation (default: `../../mfem`)
- `GTEST_DIR`: Path to Google Test source (default: `/usr/src/googletest`)

### Parallel Testing
Tests are designed to run with MPI but can also run in serial mode. The test framework automatically initializes MPI and creates appropriate parallel mesh and finite element spaces.

## Adding New Tests

### Test File Template
```cpp
#include <gtest/gtest.h>
#include "../your_header.hpp"

class YourTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize test data
    }
    
    void TearDown() override {
        // Clean up
    }
    
    // Test data members
};

TEST_F(YourTest, TestName) {
    // Test implementation
    EXPECT_EQ(expected, actual);
}
```

### Guidelines
1. **Isolation**: Each test should be independent
2. **Naming**: Use descriptive test and fixture names
3. **Setup**: Initialize required MFEM objects in SetUp()
4. **Cleanup**: Delete allocated objects in TearDown()
5. **Assertions**: Use appropriate EXPECT_* macros
6. **Documentation**: Comment complex test logic

## Continuous Integration

Tests are designed to run in CI environments:
- Automatic MPI initialization
- Graceful handling of missing dependencies
- Proper cleanup of resources
- Exit codes for success/failure

## Troubleshooting

### Common Issues

1. **MFEM not found**: Set `MFEM_DIR` environment variable
2. **Google Test not found**: Install libgtest-dev package
3. **MPI errors**: Ensure MPI is properly installed
4. **Boost errors**: Install libboost-program-options-dev

### Debug Mode
```bash
# Build with debug symbols
make clean
CXXFLAGS="-g -O0" make test

# Run with debugger
gdb ./run_tests
```

### Memory Checking
```bash
# Run with Valgrind
valgrind --leak-check=full ./run_tests
```

## Contributing

When adding new features to Laghost:
1. Add corresponding unit tests
2. Ensure tests pass in both serial and parallel modes
3. Update this README if adding new test categories
4. Follow existing naming conventions and patterns