# Laghost Testing Guide

## Overview

This document provides comprehensive guidance for testing the Laghost computational geodynamics solver. The testing framework includes both unit tests and integration tests to ensure code reliability and correctness.

## Test Architecture

### Test Types

1. **Unit Tests**: Test individual components in isolation
   - Parameter structures and validation
   - Mathematical functions and coefficients
   - Material models and rheology
   - Finite element assembly operations
   - Solver components

2. **Integration Tests**: Test complete workflows
   - Full simulation runs with known solutions
   - Multi-physics coupling validation
   - Performance benchmarks

### Test Framework

- **Google Test**: C++ unit testing framework
- **MPI Support**: Parallel testing capabilities
- **MFEM Integration**: Finite element library testing
- **Automated CI**: Continuous integration support

## Quick Start

### Prerequisites
```bash
# Install dependencies (Ubuntu/Debian)
sudo apt-get install libgtest-dev cmake libboost-program-options-dev

# Build Google Test
cd /usr/src/gtest
sudo cmake CMakeLists.txt && sudo make
sudo cp *.a /usr/lib
```

### Running Tests
```bash
# From project root - run all tests
make test-all

# Unit tests only
make unit-tests

# Integration tests only  
make tests

# Custom test runs
cd test
./run_tests.sh --help
```

## Unit Test Details

### Test Coverage

| Component | File | Coverage |
|-----------|------|----------|
| Parameters | `test_parameters.cpp` | Structure initialization, constants validation |
| Input | `test_input.cpp` | Config parsing, command-line args, string parsing |
| Functions | `test_function.cpp` | Math functions, coefficients, problem setup |
| Rheology | `test_rheology.cpp` | Material models, plasticity, return mapping |
| Assembly | `test_assembly.cpp` | FE assembly, operators, quadrature data |
| Solver | `test_solver.cpp` | Time integration, energy computation, main solver |

### Key Test Cases

#### Parameter Validation
```cpp
TEST_F(ParametersTest, DefaultSimParametersInitialization) {
    Sim sim = {};
    EXPECT_EQ(sim.problem, 0);
    EXPECT_EQ(sim.dim, 0);
    EXPECT_FALSE(sim.visualization);
}
```

#### Rheology Testing
```cpp
TEST_F(RheologyTest, PlasticYielding2D) {
    // Set high stress to trigger yielding
    stress_2d[0] = -200.0e6;  // High compression
    
    Returnmapping2d(stress_2d, stress_old, strain_inc, plastic_strain, ...);
    
    // Verify plastic strain increases
    EXPECT_GT(plastic_strain[0], initial_plastic_strain);
}
```

#### Assembly Operations
```cpp
TEST_F(AssemblyTest, MassPAOperatorMult) {
    MassPAOperator mass_pa(*H1FESpace, int_rule, rho_coeff);
    
    Vector x(mass_pa.Width()), y(mass_pa.Height());
    x = 1.0;
    mass_pa.Mult(x, y);
    
    // Mass matrix should be positive definite
    EXPECT_GT(y.Norml2(), 0.0);
}
```

## Integration Test Details

### Verification Problems

The integration tests use established benchmarks with known analytical or reference solutions:

1. **Taylor-Green Vortex** (`-p 0`): Smooth velocity field evolution
2. **Sedov Blast** (`-p 1`): Shock wave propagation  
3. **Noh Problem** (`-p 2`): Spherical implosion
4. **Triple Point** (`-p 3`): Multi-material shock interaction
5. **Gresho Vortex** (`-p 4`): Vorticity preservation

### Reference Solutions

| Problem | Dimension | Final Time | Expected Energy | Tolerance |
|---------|-----------|------------|-----------------|-----------|
| Taylor-Green 2D | 2 | 0.75 | 4.9695537349e+01 | 1e-8 |
| Taylor-Green 3D | 3 | 0.75 | 3.3909635545e+03 | 1e-8 |
| Sedov 2D | 2 | 0.8 | 4.6303396053e+01 | 1e-8 |
| Sedov 3D | 3 | 0.6 | 1.3408616722e+02 | 1e-8 |

### Running Verification Tests
```bash
# Run specific verification problem
mpirun -np 4 ./laghost -p 1 -dim 2 -rs 3 -tf 0.8 -pa

# Automated verification suite
make checks

# Performance benchmarks
make 1  # 1 MPI rank
make 4  # 4 MPI ranks
```

## Test Development

### Adding New Unit Tests

1. **Create Test File**: Follow naming convention `test_<component>.cpp`
2. **Include Headers**: Add necessary Laghost and MFEM headers
3. **Test Fixture**: Create class inheriting from `::testing::Test`
4. **Setup/Teardown**: Initialize/cleanup in `SetUp()`/`TearDown()`
5. **Test Cases**: Use `TEST_F(FixtureName, TestName)` macro

#### Template
```cpp
#include <gtest/gtest.h>
#include "../laghost_component.hpp"

class ComponentTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize test objects
        mesh = new Mesh(...);
        // ... other setup
    }
    
    void TearDown() override {
        delete mesh;
        // ... cleanup
    }
    
    Mesh* mesh;
    // ... other test data
};

TEST_F(ComponentTest, BasicFunctionality) {
    // Test implementation
    EXPECT_EQ(expected, actual);
    EXPECT_NEAR(expected_double, actual_double, tolerance);
    EXPECT_TRUE(condition);
}
```

### Test Guidelines

#### Best Practices
1. **Independence**: Tests should not depend on each other
2. **Determinism**: Tests should produce consistent results
3. **Speed**: Keep tests fast (< 1 second each)
4. **Clarity**: Use descriptive names and comments
5. **Coverage**: Test both success and failure cases

#### Common Patterns
```cpp
// Numerical tolerance testing
EXPECT_NEAR(computed_value, expected_value, 1e-10);

// Exception testing
EXPECT_THROW(risky_function(), std::exception);
EXPECT_NO_THROW(safe_function());

// Vector/array testing
for (int i = 0; i < vector.Size(); i++) {
    EXPECT_DOUBLE_EQ(vector[i], expected[i]);
}

// Positive definiteness (for mass matrices)
EXPECT_GT(matrix_norm, 0.0);
```

### Adding Integration Tests

1. **Problem Setup**: Define new problem in `laghost_function.cpp`
2. **Reference Solution**: Compute or obtain analytical solution
3. **Test Case**: Add to verification table in makefile
4. **Tolerance**: Set appropriate numerical tolerance

## Debugging Tests

### Common Issues

#### Build Failures
```bash
# Check MFEM configuration
make status

# Verify dependencies
pkg-config --libs mfem
ldconfig -p | grep gtest
```

#### Runtime Failures
```bash
# Run single test with debugging
gdb ./run_tests
(gdb) run --gtest_filter="SpecificTest.*"

# Memory checking
valgrind --leak-check=full ./run_tests

# MPI debugging
mpirun -np 1 gdb ./run_tests
```

#### Numerical Issues
```bash
# Increase tolerance for floating-point comparisons
EXPECT_NEAR(a, b, 1e-12);  // Instead of EXPECT_EQ

# Check for NaN/Inf values
EXPECT_TRUE(std::isfinite(result));
EXPECT_FALSE(std::isnan(result));
```

### Test Output Analysis

#### Successful Run
```
[==========] Running 45 tests from 6 test suites.
[----------] Global test environment set-up.
[----------] 8 tests from ParametersTest
[ RUN      ] ParametersTest.DefaultSimParametersInitialization
[       OK ] ParametersTest.DefaultSimParametersInitialization (0 ms)
...
[==========] 45 tests from 6 test suites ran. (1234 ms total)
[  PASSED  ] 45 tests.
```

#### Failed Test
```
[ RUN      ] RheologyTest.PlasticYielding2D
test_rheology.cpp:123: Failure
Expected: (plastic_strain[0]) > (initial_plastic_strain), actual: 0 vs 0
[  FAILED  ] RheologyTest.PlasticYielding2D (5 ms)
```

## Continuous Integration

### GitHub Actions
```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v2
    - name: Install dependencies
      run: |
        sudo apt-get update
        sudo apt-get install -y libgtest-dev cmake libboost-program-options-dev
    - name: Build MFEM
      run: make setup
    - name: Build Laghost
      run: make -j 4
    - name: Run tests
      run: make test-all
```

### Local CI Simulation
```bash
# Clean build and test
make distclean
make setup
make -j 4
make test-all
```

## Performance Testing

### Benchmarking
```bash
# Time complete test suite
time make test-all

# Profile specific tests
perf record ./run_tests --gtest_filter="SolverTest.*"
perf report

# Memory usage
/usr/bin/time -v ./run_tests
```

### Scaling Tests
```bash
# Test parallel scaling
for np in 1 2 4 8; do
    echo "Testing with $np processes"
    mpirun -np $np ./run_tests
done
```

## Maintenance

### Regular Tasks
1. **Update Reference Solutions**: When algorithms change
2. **Add New Tests**: For new features
3. **Review Coverage**: Ensure comprehensive testing
4. **Performance Monitoring**: Track test execution time
5. **Dependency Updates**: Keep test framework current

### Release Testing
```bash
# Full test suite before release
make distclean
make setup
make -j 4
make test-all
make checks

# Multi-platform testing
# - Different compilers (GCC, Clang, Intel)
# - Different MPI implementations (OpenMPI, MPICH)
# - Different architectures (x86_64, ARM)
```