#include "Geochemistry.hpp"
#include "Glenn.hpp"

#include <iostream>
#include <vector>

using namespace std;

int main (int argc, char **args) {

  std::string filename;

  int verbose = 1;

  std::vector<double> total;

  // create geochemistry object
  Beaker g;
#if 1
  // set up simple 2-species carbonate system (H,HCO3-)
  createCarbonateSystem(&total,&g);
#else
  filename = "reaction.dat";
  readChemistryFromFile(filename,&g);
  filename = "target_total.dat";
  readTargetTotalFromFile(filename,g.ncomp(),&total) ;
#endif
  // solve for free-ion concentrations
  g.verbose(verbose);
  g.react(total,1.,1.,1.,1.);
  g.print_results();
  g.speciate(total);
  g.print_results();

  cout << "Done!\n";

}
