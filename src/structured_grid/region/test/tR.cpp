#include <iostream>
#include <fstream>
#include <iomanip>

// closing DSO objects
#include "VerboseObject_objs.hh"

using std::cout;
using std::endl;
#include <ParmParse.H>
#include <VisMF.H>

#include <RegionManager.H>

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  RegionManager rm;
  std::cout << rm << std::endl;

  Array<Real> plo(BL_SPACEDIM,0);
  Array<Real> dx(BL_SPACEDIM,1);
  const Array<const Region*> regions = rm.RegionPtrArray();
  for (int i=0; i<regions.size(); ++i) {
    std::cout << regions[i]->approximate_bounds(plo,dx) << '\n';
  }


  return 0;
}
