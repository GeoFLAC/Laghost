#include <gtest/gtest.h>
#include <mpi.h>
#include "mfem.hpp"

// Include all Laghost headers
#include "../../src/io/laghost_parameters.hpp"
#include "../../src/core/laghost_constants.hpp"
#include "../../src/io/laghost_input.hpp"
#include "../../src/physics/laghost_function.hpp"
#include "../../src/physics/laghost_rheology.hpp"
#include "../../src/core/laghost_assembly.hpp"
#include "../../src/core/laghost_solver.hpp"

using namespace mfem;
using namespace mfem::geodynamics;

// Global test fixture for MPI initialization
class LaghostTestEnvironment : public ::testing::Environment {
public:
    void SetUp() override {
        int provided;
        MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SINGLE, &provided);
        Hypre::Init();
    }
    
    void TearDown() override {
        Hypre::Finalize();
        MPI_Finalize();
    }
};

// Base test fixture with common setup
class LaghostTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for all tests
        dim = 2;
        mesh = nullptr;
        pmesh = nullptr;
    }
    
    void TearDown() override {
        delete pmesh;
        delete mesh;
    }
    
    // Helper methods
    void CreateSimpleMesh(int nx = 2, int ny = 2) {
        mesh = new Mesh(Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL, true));
        pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
        delete mesh;
        mesh = nullptr;
    }
    
    int dim;
    Mesh* mesh;
    ParMesh* pmesh;
};

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new LaghostTestEnvironment);
    return RUN_ALL_TESTS();
}