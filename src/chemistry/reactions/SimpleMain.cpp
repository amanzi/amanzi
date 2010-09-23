#include "Geochemistry.hpp"
#include "Glenn.hpp"

#include <iostream>
#include <vector>

using namespace std;

int main (int argc, char **args) {

  std::string filename;

  std::vector<double> total;

  // create geochemistry object
  Geochemistry g;
#if 0
  // set up simple 2-species carbonate system (H,HCO3-)
  createCarbonateSystem(&total,&g);
#else
  filename = "reaction.dat";
  readChemistryFromFile(filename,&g);
  filename = "target_total.dat";
  readTargetTotalFromFile(filename,g.get_ncomp(),&total) ;
#endif
  // solve for free-ion concentrations
  g.speciate(total);

  cout << "Done!\n";

}
