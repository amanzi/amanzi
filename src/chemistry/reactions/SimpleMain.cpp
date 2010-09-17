#include "Geochemistry.hpp"
#include "Glenn.hpp"

#include <iostream>
#include <vector>

using namespace std;

int main (int argc, char **args) {

  std::vector<double> total;

  // create geochemistry object
  Geochemistry g;
  // set up simple 2-species carbonate system (H,HCO3-)
  createCarbonateSystem(&total,&g);
  // solve for free-ion concentrations
  g.speciate(total);

  cout << "Done!\n";

}
