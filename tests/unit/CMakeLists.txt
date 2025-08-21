cmake_minimum_required(VERSION 3.10)
project(LaghostTests)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(MPI REQUIRED)
find_package(PkgConfig REQUIRED)

# Find Google Test
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Find Boost
find_package(Boost REQUIRED COMPONENTS program_options)

# Find MFEM
set(MFEM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../mfem" CACHE PATH "Path to MFEM installation")
find_library(MFEM_LIBRARY mfem PATHS ${MFEM_DIR}/lib ${MFEM_DIR})
find_path(MFEM_INCLUDE_DIR mfem.hpp PATHS ${MFEM_DIR}/include ${MFEM_DIR})

if(NOT MFEM_LIBRARY OR NOT MFEM_INCLUDE_DIR)
    message(FATAL_ERROR "MFEM not found. Please set MFEM_DIR to the MFEM installation directory.")
endif()

# Find HYPRE (required by MFEM)
find_library(HYPRE_LIBRARY HYPRE PATHS ${MFEM_DIR}/../hypre-*/src/hypre/lib)
find_path(HYPRE_INCLUDE_DIR HYPRE.h PATHS ${MFEM_DIR}/../hypre-*/src/hypre/include)

# Find METIS (optional, used by MFEM)
find_library(METIS_LIBRARY metis PATHS ${MFEM_DIR}/../metis-*/libmetis)

# Include directories
include_directories(${MFEM_INCLUDE_DIR})
include_directories(${HYPRE_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall")
add_definitions(-DMFEM_USE_MPI)

# Source files
set(LAGHOST_SOURCES
    ../laghost_input.cpp
    ../laghost_function.cpp
    ../laghost_rheology.cpp
    ../laghost_assembly.cpp
    ../laghost_solver.cpp
    ../laghost_tmop.cpp
    ../laghost_remhos.cpp
)

set(TEST_SOURCES
    unit_tests.cpp
    test_parameters.cpp
    test_input.cpp
    test_function.cpp
    test_rheology.cpp
    test_assembly.cpp
    test_solver.cpp
)

# Create test executable
add_executable(run_tests ${TEST_SOURCES} ${LAGHOST_SOURCES})

# Link libraries
target_link_libraries(run_tests
    ${GTEST_LIBRARIES}
    ${MFEM_LIBRARY}
    ${HYPRE_LIBRARY}
    ${MPI_LIBRARIES}
    ${Boost_LIBRARIES}
    pthread
)

if(METIS_LIBRARY)
    target_link_libraries(run_tests ${METIS_LIBRARY})
endif()

# Add MPI compile flags
if(MPI_COMPILE_FLAGS)
    set_target_properties(run_tests PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(run_tests PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

# Enable testing
enable_testing()

# Add test
add_test(NAME LaghostUnitTests COMMAND run_tests)

# Custom target for running tests
add_custom_target(test_laghost
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/run_tests
    DEPENDS run_tests
    COMMENT "Running Laghost unit tests"
)

# Installation
install(TARGETS run_tests DESTINATION bin)

# Print configuration summary
message(STATUS "Laghost Tests Configuration Summary:")
message(STATUS "  MFEM_DIR: ${MFEM_DIR}")
message(STATUS "  MFEM_LIBRARY: ${MFEM_LIBRARY}")
message(STATUS "  MFEM_INCLUDE_DIR: ${MFEM_INCLUDE_DIR}")
message(STATUS "  HYPRE_LIBRARY: ${HYPRE_LIBRARY}")
message(STATUS "  HYPRE_INCLUDE_DIR: ${HYPRE_INCLUDE_DIR}")
if(METIS_LIBRARY)
    message(STATUS "  METIS_LIBRARY: ${METIS_LIBRARY}")
else()
    message(STATUS "  METIS_LIBRARY: Not found (optional)")
endif()
message(STATUS "  Boost_LIBRARIES: ${Boost_LIBRARIES}")
message(STATUS "  GTEST_LIBRARIES: ${GTEST_LIBRARIES}")