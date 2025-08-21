#!/bin/bash

# Laghost Unit Test Runner Script

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Default values
MFEM_DIR=${MFEM_DIR:-"../../mfem"}
BUILD_TESTS=true
RUN_TESTS=true
PARALLEL=false
NP=1
VERBOSE=false
FILTER=""

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
show_usage() {
    echo "Laghost Unit Test Runner"
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help              Show this help message"
    echo "  -b, --build-only        Only build tests, don't run them"
    echo "  -r, --run-only          Only run tests, don't build them"
    echo "  -p, --parallel          Run tests in parallel with MPI"
    echo "  -n, --np NUM            Number of MPI processes (default: 1)"
    echo "  -v, --verbose           Verbose test output"
    echo "  -f, --filter PATTERN    Run only tests matching pattern"
    echo "  --mfem-dir DIR          Path to MFEM directory"
    echo ""
    echo "Examples:"
    echo "  $0                      Build and run all tests"
    echo "  $0 -p -n 4              Run tests with 4 MPI processes"
    echo "  $0 -f \"ParametersTest.*\" Run only parameter tests"
    echo "  $0 -v                   Run with verbose output"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_usage
            exit 0
            ;;
        -b|--build-only)
            RUN_TESTS=false
            shift
            ;;
        -r|--run-only)
            BUILD_TESTS=false
            shift
            ;;
        -p|--parallel)
            PARALLEL=true
            shift
            ;;
        -n|--np)
            NP="$2"
            PARALLEL=true
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -f|--filter)
            FILTER="$2"
            shift 2
            ;;
        --mfem-dir)
            MFEM_DIR="$2"
            shift 2
            ;;
        *)
            print_error "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Check if we're in the test directory
if [[ ! -f "Makefile" ]] || [[ ! -f "unit_tests.cpp" ]]; then
    print_error "This script must be run from the test directory"
    exit 1
fi

# Check MFEM directory
if [[ ! -d "$MFEM_DIR" ]]; then
    print_error "MFEM directory not found: $MFEM_DIR"
    print_error "Please set MFEM_DIR environment variable or use --mfem-dir option"
    exit 1
fi

export MFEM_DIR

print_status "Laghost Unit Test Runner"
print_status "MFEM Directory: $MFEM_DIR"
print_status "Parallel: $PARALLEL (np=$NP)"
print_status "Verbose: $VERBOSE"
if [[ -n "$FILTER" ]]; then
    print_status "Filter: $FILTER"
fi

# Build tests
if [[ "$BUILD_TESTS" == true ]]; then
    print_status "Building tests..."
    
    # Clean previous build
    make clean > /dev/null 2>&1 || true
    
    # Build
    if make all; then
        print_status "Build successful"
    else
        print_error "Build failed"
        exit 1
    fi
fi

# Run tests
if [[ "$RUN_TESTS" == true ]]; then
    print_status "Running tests..."
    
    # Check if test executable exists
    if [[ ! -f "run_tests" ]]; then
        print_error "Test executable not found. Please build first."
        exit 1
    fi
    
    # Prepare test command
    TEST_CMD="./run_tests"
    
    # Add Google Test options
    if [[ "$VERBOSE" == true ]]; then
        TEST_CMD="$TEST_CMD --gtest_verbose"
    fi
    
    if [[ -n "$FILTER" ]]; then
        TEST_CMD="$TEST_CMD --gtest_filter=\"$FILTER\""
    fi
    
    # Run with or without MPI
    if [[ "$PARALLEL" == true ]]; then
        if command -v mpirun > /dev/null 2>&1; then
            print_status "Running with MPI ($NP processes)..."
            eval "mpirun -np $NP $TEST_CMD"
        else
            print_warning "mpirun not found, running in serial mode"
            eval "$TEST_CMD"
        fi
    else
        print_status "Running in serial mode..."
        eval "$TEST_CMD"
    fi
    
    # Check test results
    if [[ $? -eq 0 ]]; then
        print_status "All tests passed!"
    else
        print_error "Some tests failed"
        exit 1
    fi
fi

print_status "Done"