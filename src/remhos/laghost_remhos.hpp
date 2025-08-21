#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;

void Remapping(ParMesh *, ParGridFunction &, ParGridFunction &, ParGridFunction &, int &, int &, bool &, bool &);
// void NCRemapping(ParNCMesh *, ParGridFunction &, ParGridFunction &, ParGridFunction &, int &, int &, bool &);

// void Remapping_stress(ParMesh *, ParGridFunction &, ParGridFunction &, ParGridFunction &, int &, int &, bool &);