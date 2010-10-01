#include "Geochemistry.hpp"
#include "Beaker.hpp"
#include "Glenn.hpp"

#include <iostream>
#include <vector>

using namespace std;

int main (int argc, char **args) {

  std::string filename;

  int verbose = 1;

  std::vector<double> total;

  // create geochemistry object
  Beaker beaker;
  Beaker::BeakerParameters parameters = beaker.GetDefaultParameters();
  parameters.porosity = 0.25;
  parameters.saturation = 1.;
  parameters.water_density = 997.205133945901; // kg / m^3
  parameters.volume = 0.25; // m^3

#if 1
  // set up simple 2-species carbonate system (H,HCO3-)
  createCarbonateSystem(&total,&beaker);
#else
  filename = "reaction.dat";
  readChemistryFromFile(filename,&beaker);
  filename = "target_total.dat";
  readTargetTotalFromFile(filename,beaker.ncomp(),&total) ;
  // convert totals from molality [mol/kg water] -> molarity [mol/L water]
  for (int i = 0; i < (int)total.size(); i++)
    total[i] *= water_density/1000.;
#endif
  // solve for free-ion concentrations
  beaker.verbosity(verbose);

  // to test react with unitary time step
//  beaker.react(total, parameters, 1.0);
//  beaker.print_results();

  // to test speciation
//  beaker.speciate(total, parameters.water_density);
//  beaker.print_results();

  // to test a time stepping loop with kinetic reactions
  Glenn g(&beaker);

  double final_time = 0.25 * 365. * 24. * 3600; // 1/4 year
  double time_step = final_time / 100;

  g.solve(total, final_time, time_step, parameters);

  cout << "Done!\n";

}
