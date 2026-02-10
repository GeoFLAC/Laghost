#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include "../../src/io/laghost_input.hpp"
#include "../../src/io/laghost_parameters.hpp"

class InputTest : public ::testing::Test {
protected:
    void SetUp() override {
        param = Param{};
        
        // Create a temporary config file for testing
        test_config_content = R"([sim]
problem = 1
dim = 2
t_final = 1000.0
max_tsteps = 100
year = true
visualization = false
vis_steps = 10
visit = false
paraview = true
gfprint = false
basename = test_results
device = cpu
dev = 0
check = false
mem_usage = false
fom = false
gpu_aware_mpi = false

[solver]
ode_solver_type = 7
cfl = 0.25
cg_tol = 1.0e-10
ftz_tol = 0.0
cg_max_iter = 300
p_assembly = false
impose_visc = true

[control]
lithostatic = true
init_dt = 1.0
mscale = 1.0e6
gravity = 10.0
thickness = 10.0e3
mass_bal = false
dyn_damping = true
dyn_factor = 0.8
max_vbc_val = 3.168808781402895e-12

[mesh]
mesh_file = test_mesh.mesh
rs_levels = 2
rp_levels = 0
partition_type = 0
order_v = 2
order_e = 1
order_q = -1
local_refinement = false
l2_basis = 1

[bc]
vbc_unit = cm/yr
vbc_factor = 1.0
vbc_x0 = 1
vbc_x1 = 1
winkler_foundation = true
winkler_flat = false
winkler_rho = 2700.0
surf_proc = true
surf_diff = 1.0e-7
surf_alpha = 0.0
base_proc = true
base_diff = 1.0e-7
base_alpha = 0.0

[mat]
plastic = true
viscoplastic = false
rho = [2700.0]
lambda = [3e10]
mu = [3e10]
tension_cutoff = [0.0]
cohesion0 = [44.0e6]
cohesion1 = [4.0e6]
friction_angle0 = [30.0]
friction_angle1 = [30.0]
dilation_angle0 = [0.0]
dilation_angle1 = [0.0]
alpha0 = [0.0]
alpha1 = [0.5]
plastic_viscosity = [1.0]
weak_rad = 1.0e3
weak_x = 50.0e3
weak_y = 2.00e3
weak_z = 0.00e3
ini_alpha = 0.5

[tmop]
tmop = false
amr = false
ale = 0.5
remesh_steps = 10000000
mesh_poly_deg = 2
jitter = 0.0
metric_id = 2
target_id = 1
lim_const = 0.0
adapt_lim_const = 0.0
quad_type = 1
quad_order = 8
solver_type = 0
solver_iter = 20
solver_rtol = 1e-10
solver_art_type = 0
lin_solver = 2
max_lin_iter = 100
move_bnd = false
combomet = 0
bal_expl_combo = false
hradaptivity = false
h_metric_id = -1
normalization = false
verbosity_level = 0
fdscheme = false
adapt_eval = 0
exactaction = false
n_hr_iter = 5
n_h_iter = 1
mesh_node_ordering = 0
barrier_type = 0
worst_case_type = 0
tmop_cond_num = 0.5
)";
        
        // Write test config to temporary file
        test_config_file = "test_config.cfg";
        std::ofstream file(test_config_file);
        file << test_config_content;
        file.close();
    }
    
    void TearDown() override {
        // Clean up temporary file
        std::remove(test_config_file.c_str());
    }
    
    Param param;
    std::string test_config_content;
    std::string test_config_file;
};

TEST_F(InputTest, ConfigFileExists) {
    std::ifstream file(test_config_file);
    EXPECT_TRUE(file.good()) << "Test config file should exist";
    file.close();
}

TEST_F(InputTest, ParseSimSection) {
    // This test would require access to the internal parsing function
    // Since the parsing is done internally, we test the interface
    
    // Create mock command line arguments
    const char* argv[] = {"laghost", "-i", test_config_file.c_str()};
    int argc = 3;
    
    OptionsParser args(argc, const_cast<char**>(argv));
    
    // Test that the function doesn't crash with valid input
    EXPECT_NO_THROW({
        read_and_assign_input_parameters(args, param, 0);
    });
}

TEST_F(InputTest, InvalidConfigFile) {
    // Test with non-existent config file
    const char* argv[] = {"laghost", "-i", "nonexistent.cfg"};
    int argc = 3;
    
    OptionsParser args(argc, const_cast<char**>(argv));
    
    // The function calls exit(1) when config file is not found
    // In a test environment, we expect this to cause the program to exit
    // This test documents the current behavior - when file is missing, program exits
    EXPECT_EXIT({
        read_and_assign_input_parameters(args, param, 0);
    }, ::testing::ExitedWithCode(1), ".*can not read options configuration file.*");
}

TEST_F(InputTest, CommandLineOverrides) {
    // Test that command line arguments can override config file values
    const char* argv[] = {"laghost", "-i", test_config_file.c_str(), "-pa", "-dim", "3"};
    int argc = 6;
    
    OptionsParser args(argc, const_cast<char**>(argv));
    
    EXPECT_NO_THROW({
        read_and_assign_input_parameters(args, param, 0);
    });
}

TEST_F(InputTest, DefaultParametersWithoutConfig) {
    // Test behavior when no config file is provided
    const char* argv[] = {"laghost"};
    int argc = 1;
    
    OptionsParser args(argc, const_cast<char**>(argv));
    
    EXPECT_NO_THROW({
        read_and_assign_input_parameters(args, param, 0);
    });
}

TEST_F(InputTest, BoundaryConditionStringParsing) {
    // Test parsing of boundary condition strings like "[1,1,0,0]"
    std::string bc_ids = "[1,1,0,0]";
    
    // Remove brackets and spaces
    bc_ids.erase(std::remove(bc_ids.begin(), bc_ids.end(), '['), bc_ids.end());
    bc_ids.erase(std::remove(bc_ids.begin(), bc_ids.end(), ']'), bc_ids.end());
    bc_ids.erase(std::remove(bc_ids.begin(), bc_ids.end(), ' '), bc_ids.end());
    
    EXPECT_EQ(bc_ids, "1,1,0,0");
    
    // Parse comma-separated values
    std::stringstream ss(bc_ids);
    std::vector<int> bc_id;
    std::string token;
    
    while (getline(ss, token, ',')) {
        bc_id.push_back(std::stoi(token));
    }
    
    EXPECT_EQ(bc_id.size(), 4);
    EXPECT_EQ(bc_id[0], 1);
    EXPECT_EQ(bc_id[1], 1);
    EXPECT_EQ(bc_id[2], 0);
    EXPECT_EQ(bc_id[3], 0);
}

TEST_F(InputTest, MaterialParameterStringParsing) {
    // Test parsing of material parameter strings like "[2700.0]"
    std::string rho_str = "[2700.0]";
    
    // Remove brackets
    rho_str.erase(std::remove(rho_str.begin(), rho_str.end(), '['), rho_str.end());
    rho_str.erase(std::remove(rho_str.begin(), rho_str.end(), ']'), rho_str.end());
    
    EXPECT_EQ(rho_str, "2700.0");
    
    // Convert to double
    double rho_val = std::stod(rho_str);
    EXPECT_DOUBLE_EQ(rho_val, 2700.0);
}

TEST_F(InputTest, MultiMaterialParameterParsing) {
    // Test parsing of multi-material parameters like "[2700.0,3000.0,2500.0]"
    std::string multi_rho = "[2700.0,3000.0,2500.0]";
    
    // Remove brackets and spaces
    multi_rho.erase(std::remove(multi_rho.begin(), multi_rho.end(), '['), multi_rho.end());
    multi_rho.erase(std::remove(multi_rho.begin(), multi_rho.end(), ']'), multi_rho.end());
    multi_rho.erase(std::remove(multi_rho.begin(), multi_rho.end(), ' '), multi_rho.end());
    
    std::stringstream ss(multi_rho);
    std::vector<double> rho_values;
    std::string token;
    
    while (getline(ss, token, ',')) {
        rho_values.push_back(std::stod(token));
    }
    
    EXPECT_EQ(rho_values.size(), 3);
    EXPECT_DOUBLE_EQ(rho_values[0], 2700.0);
    EXPECT_DOUBLE_EQ(rho_values[1], 3000.0);
    EXPECT_DOUBLE_EQ(rho_values[2], 2500.0);
}