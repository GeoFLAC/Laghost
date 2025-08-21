#ifndef LAGHOST_INPUT_HPP
#define LAGHOST_INPUT_HPP

#include "mfem.hpp"
#include "laghost_parameters.hpp"
using namespace mfem;

void read_and_assign_input_parameters( mfem::OptionsParser& args, Param& param, const int &myid );
static void get_input_parameters( const char*, Param& );

#endif
