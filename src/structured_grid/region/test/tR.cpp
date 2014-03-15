#include <iostream>
#include <fstream>
#include <iomanip>
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
  const Array<const Region*> regions = rm.RegionPtrArray();
  for (int i=0; i<regions.size(); ++i) {
    std::cout << *regions[i] << std::endl;
  }
  return 0;
}
