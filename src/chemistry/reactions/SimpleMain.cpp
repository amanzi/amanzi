#include "Geochemistry.hpp"
#include "Glenn.hpp"

#include <iostream>
#include <vector>

using namespace std;

int main (int argc, char **args) {
  int verbose = 1;
  std::vector<double> total;

  // create geochemistry object
  Geochemistry g;
  // set up simple 2-species carbonate system (H,HCO3-)
  createCarbonateSystem(&total,&g);
  // solve for free-ion concentrations
  g.verbose(verbose);
  g.speciate(total);
  g.print_results();

  cout << "Done!\n";

}
