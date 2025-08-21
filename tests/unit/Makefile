# Laghost Unit Tests Makefile

# Use the MFEM build directory
MFEM_DIR ?= ../../../mfem
CONFIG_MK = $(MFEM_DIR)/config/config.mk

# Include MFEM configuration
-include $(CONFIG_MK)

# Google Test configuration
GTEST_DIR ?= /usr/src/googletest
GTEST_INCLUDE = -I$(GTEST_DIR)/googletest/include
GTEST_LIB = -lgtest -lgtest_main -lpthread

# Boost configuration (for input parsing tests)
BOOST_LIBS = -lboost_program_options

# Compiler and flags
CXX = $(MFEM_CXX)
CPPFLAGS = $(MFEM_CPPFLAGS) $(GTEST_INCLUDE) -I../../src
CXXFLAGS = $(MFEM_CXXFLAGS) -std=c++17
LAGHOST_FLAGS = $(CPPFLAGS) $(CXXFLAGS) $(MFEM_INCFLAGS)

# Extra include dir for MFEM headers
EXTRA_INC_DIR = $(or $(wildcard $(MFEM_DIR)/include/mfem),$(MFEM_DIR))
CCC = $(strip $(CXX) $(LAGHOST_FLAGS) $(if $(EXTRA_INC_DIR),-I$(EXTRA_INC_DIR)))

# Libraries
LAGHOST_LIBS = $(MFEM_LIBS) $(MFEM_EXT_LIBS)
LIBS = $(strip $(LAGHOST_LIBS) $(LDFLAGS) $(BOOST_LIBS) $(GTEST_LIB))

# Source files
TEST_SOURCES = unit_tests.cpp test_parameters.cpp test_input.cpp test_function.cpp \
               test_assembly.cpp test_solver.cpp

# Object files
TEST_OBJECTS = $(TEST_SOURCES:.cpp=.o)

# Laghost source files (needed for linking)
LAGHOST_SOURCES = ../../src/io/laghost_input.cpp ../../src/physics/laghost_function.cpp \
                  ../../src/physics/laghost_rheology.cpp ../../src/core/laghost_assembly.cpp \
                  ../../src/core/laghost_solver.cpp ../../src/physics/laghost_tmop.cpp
LAGHOST_OBJECTS = $(LAGHOST_SOURCES:.cpp=.o)

# Main target
.PHONY: all test clean

all: run_tests

# Build and run tests
test: run_tests
	./run_tests

# Test executable
run_tests: $(TEST_OBJECTS) $(LAGHOST_OBJECTS) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_OBJECTS) $(LAGHOST_OBJECTS) $(LIBS)

# Compile test files
%.o: %.cpp
	$(CCC) -c $< -o $@

# Compile Laghost source files
../../src/%.o: ../../src/%.cpp
	$(CCC) -c $< -o $@

# Clean targets
clean:
	rm -f $(TEST_OBJECTS) $(LAGHOST_OBJECTS) run_tests

distclean: clean
	rm -f *~

# Dependencies
$(TEST_OBJECTS): $(wildcard ../../src/*/*.hpp)
$(LAGHOST_OBJECTS): $(wildcard ../../src/*/*.hpp)

# Help target
help:
	@echo "Laghost Unit Tests Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all       - Build test executable"
	@echo "  test      - Build and run tests"
	@echo "  clean     - Remove object files and executable"
	@echo "  distclean - Remove all generated files"
	@echo "  help      - Show this help message"
	@echo ""
	@echo "Requirements:"
	@echo "  - Google Test library installed"
	@echo "  - MFEM library built and configured"
	@echo "  - Boost program_options library"
	@echo ""
	@echo "Environment Variables:"
	@echo "  MFEM_DIR   - Path to MFEM directory (default: ../../mfem)"
	@echo "  GTEST_DIR  - Path to Google Test source (default: /usr/src/googletest)"

# Generate an error message if the MFEM library is not configured
$(CONFIG_MK):
	$(error The MFEM library is not configured. Please build MFEM first.)