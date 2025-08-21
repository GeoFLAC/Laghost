# Laghost Continuous Integration

This directory contains GitHub Actions workflows that automatically test Laghost on every commit and pull request.

## Workflows

### 1. Unit Tests (`unit-tests.yml`)
- **Trigger**: Pull requests affecting source code, tests, or build system
- **Purpose**: Fast feedback for code changes
- **Tests**: Unit tests only
- **Duration**: ~10-15 minutes (with caching)

### 2. Standard CI (`ci.yml`) 
- **Trigger**: Push to master/main branches
- **Purpose**: Comprehensive validation for main branch
- **Tests**: Unit tests + integration tests + functionality checks
- **Duration**: ~20-30 minutes

### 3. Fast CI (`ci-fast.yml`)
- **Trigger**: Manual or on-demand
- **Purpose**: Optimized testing with job separation
- **Tests**: Parallel unit and integration test jobs
- **Duration**: ~15-25 minutes

### 4. Comprehensive CI (`comprehensive-ci.yml`)
- **Trigger**: Push to master + weekly schedule (Sunday 2 AM UTC)
- **Purpose**: Complete validation including memory checks
- **Tests**: Full test suite + memory analysis + multi-MPI testing
- **Duration**: ~30-45 minutes

## Dependencies Built

All workflows automatically build the required dependencies:

- **HYPRE** v2.28.0 - Parallel linear solvers
- **METIS** v4.0.3 - Graph partitioning  
- **GSLIB** - High-order interpolation
- **MFEM** (latest) - Finite element library
- **Google Test** - C++ testing framework

## Caching Strategy

Dependencies are cached to speed up builds:
- Cache key includes OS and workflow file hash
- Typical cache hit reduces build time by 60-80%
- Cache automatically invalidates when workflows change

## Test Coverage

### Unit Tests (55 tests across 5 suites):
- ✅ Parameters (9 tests)
- ✅ Input/Output (8 tests) 
- ✅ Functions (11 tests)
- ✅ Assembly (16 tests)
- ✅ Solver (11 tests)

### Integration Tests:
- ✅ Basic MFEM integration
- ✅ MPI parallel execution
- ✅ Command line argument parsing
- ✅ Configuration file loading
- ✅ Short simulation runs

### Functionality Tests:
- ✅ Help command execution
- ✅ Configuration file parsing
- ✅ Multi-MPI process execution
- ✅ Memory leak detection (basic)

## Status Badges

Add these to your main README.md:

```markdown
[![Unit Tests](https://github.com/GeoFLAC/Laghost/workflows/Unit%20Tests/badge.svg)](https://github.com/GeoFLAC/Laghost/actions?query=workflow%3A%22Unit+Tests%22)
[![CI](https://github.com/GeoFLAC/Laghost/workflows/Laghost%20CI/badge.svg)](https://github.com/GeoFLAC/Laghost/actions?query=workflow%3A%22Laghost+CI%22)
[![Comprehensive CI](https://github.com/GeoFLAC/Laghost/workflows/Comprehensive%20CI/badge.svg)](https://github.com/GeoFLAC/Laghost/actions?query=workflow%3A%22Comprehensive+CI%22)
```

## Local Testing

To run the same tests locally:

```bash
# Unit tests only
make unit-tests

# Integration test  
make test

# All tests
make test-all
```

## Troubleshooting

### Common Issues:

1. **Cache invalidation**: If builds fail unexpectedly, clear cache by updating cache key version numbers in workflows

2. **Dependency build failures**: Usually due to network issues or dependency updates. Check the build logs in the "Build MFEM dependencies" step

3. **Test failures**: Check unit test output in the "Run unit tests" step. Most failures are due to configuration or dependency issues

4. **Memory issues**: The comprehensive workflow includes basic memory checking. For detailed analysis, run valgrind locally

### Updating Workflows:

1. Modify workflow files in `.github/workflows/`
2. Test changes on a branch first
3. Update cache key versions if changing dependencies
4. Update this documentation if adding new workflows

## Performance Notes

- Typical fresh build time: 25-35 minutes
- With cache hit: 8-15 minutes  
- Unit tests execution: 1-3 minutes
- Integration tests: 2-5 minutes

The workflows are optimized for the GitHub Actions free tier limits (2000 minutes/month for public repositories).
