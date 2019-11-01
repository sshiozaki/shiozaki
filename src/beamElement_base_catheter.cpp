#include "beamElement_base.h"
using namespace std;

void beamElement_base::read_beam_geometry_catheter(TextParser &tp,const int number)
{
  string base="/Catheter"+to_string(number);
  readGeometryFromFile(tp,base);
  set_initial_t(t,x);
  set_initial_t(t0,x0);
  set_initial_t(t_ref,x_ref);

  setBoundaryCondition(tp,base);
  set_physicalCondition(tp,base);
}
