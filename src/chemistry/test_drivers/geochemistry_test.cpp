/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "SimpleThermoDatabase.hpp"
#include "Beaker.hpp"
#include "ActivityModelFactory.hpp"
#include "Verbosity.hpp"
#include "ChemistryException.hpp"

#include "geochemistry_test.hpp"

const std::string kCrunch("crunch");
const std::string kPflotran("pflotran");


int main(int argc, char **argv) {
  Verbosity verbosity = kTerse;
  int test = 0;
  std::string model("");
  int error = EXIT_SUCCESS;

  Beaker *chem = NULL;
  Beaker::BeakerComponents components;
  components.free_ion.clear();
  components.minerals.clear();
  components.ion_exchange_sites.clear();
  components.total.clear();
  components.total_sorbed.clear();

  std::string thermo_database_file("");
  std::string activity_model_name("");
  double delta_time = 0;  // seconds
  int num_time_steps = 0;
  int output_interval = 0;

  error = CommandLineOptions(argc, argv, &verbosity, &test, &model);

  if (error == EXIT_SUCCESS) {
    switch (test) {
      case 0: {
        // this is expected to fail with error checking
        if (verbosity == kTerse) {
          std::cout << "Running component size sanity check." << std::endl;
        }
        thermo_database_file = "input/carbonate.bgd";
        activity_model_name = ActivityModelFactory::unit;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(1.0e-3);  // HCO3-
        components.total.push_back(0.0);
        components.free_ion.push_back(0.0);
        components.minerals.push_back(0.0);
        components.ion_exchange_sites.push_back(0.0);
        break;
      }
      case 1: {
        // set up simple 2-species carbonate system (H,HCO3-),
        // unit activity coefficients
        if (verbosity == kTerse) {
          std::cout << "Running simple carbonate example, "
                    << "unit activity coefficients." << std::endl;
        }
        thermo_database_file = "input/carbonate.bgd";
        activity_model_name = ActivityModelFactory::unit;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(1.0e-3);  // HCO3-
        break;
      }
      case 2: {
        // set up simple 2-species carbonate system (H,HCO3-),
        // debye-huckel activity coefficients
        if (verbosity == kTerse) {
          std::cout << "Running simple carbonate example, debye-huckel." 
                    << std::endl;
        }
        thermo_database_file = "input/carbonate.bgd";
        activity_model_name = ActivityModelFactory::debye_huckel;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(1.0e-3);  // HCO3-
        break;
      }
      case 3: {
        // larger carbonate system, 3 components, 9 secondary, 
        // unit activity coefficients
        if (verbosity == kTerse) {
          std::cout << "Running large carbonate speciation example, "
                    << "unit activity coefficients." << std::endl;
        }
        thermo_database_file = "input/ca-carbonate.bgd";
        activity_model_name = ActivityModelFactory::unit;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(3.0e-3);  // HCO3-
        components.total.push_back(1.0e-3);  // Ca++
        break;
      }
      case 4: {
        // larger carbonate system, 3 components, 9 secondary,
        // debye-huckel activity coefficients
        if (verbosity == kTerse) {
          std::cout << "Running large carbonate speciation example, "
                    << "debye-huckel activity coefficients." << std::endl;
        }
        thermo_database_file = "input/ca-carbonate.bgd";
        activity_model_name = ActivityModelFactory::debye_huckel;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(3.0e-3);  // HCO3-
        components.total.push_back(1.0e-3);  // Ca++
        break;
      }
      case 5: {
        calcite_kinetics(verbosity,
                         &thermo_database_file,
                         &activity_model_name,
                         &components,
                         &delta_time,
                         &num_time_steps,
                         &output_interval);
        break;
      }
      case 6: {
        // Na+ Ca++ ion exchange
        if (verbosity == kTerse) {
          std::cout << "Running Na-Ca ion exchange problem." << std::endl;
        }
        thermo_database_file = "input/na-ca-ion-exchange.bgd";
        activity_model_name = ActivityModelFactory::debye_huckel;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(3.0e-3);  // HCO3-
        components.total.push_back(1.0e-3);  // Ca++
        components.total.push_back(1.0e-3);  // Na+
        components.total.push_back(1.0e-3);  // Cl-
        components.ion_exchange_sites.push_back(15.0);  // X-, units? equivalents per 100 grams solid?
        break;
      }
      case 7: {
        // fbasin, initial condition, no minerals
        fbasin_initial_speciation(verbosity,
                                  &thermo_database_file,
                                  &activity_model_name,
                                  &components);
        break;
      }
      case 8: {
        // fbasin, initial condition, mineral saturation state, no kinetics
        fbasin_initial_speciation(verbosity,
                                  &thermo_database_file,
                                  &activity_model_name,
                                  &components);
        fbasin_minerals(verbosity,
                        &thermo_database_file,
                        &activity_model_name,
                        &components);
        break;
      }
      case 9: {
        // fbasin, initial condition, mineral kinetics
        fbasin_initial_all(verbosity,
                           &thermo_database_file,
                           &activity_model_name,
                           &components);
        break;
      }
      case 10: {
        // calcite kinetics, long time steps
        calcite_kinetics_large_time_steps(verbosity,
                                          &thermo_database_file,
                                          &activity_model_name,
                                          &components,
                                          &delta_time,
                                          &num_time_steps,
                                          &output_interval);
        break;
      }
      case 11: {
        // surface complexation
        surface_complexation(verbosity,
                             &thermo_database_file,
                             &activity_model_name,
                             &components,
                             &delta_time,
                             &num_time_steps,
                             &output_interval);
        break;
      }
      default: {
        std::cout << "Invalid test number specified on command line. "
                  << "try using the \'-h\' option." << std::endl;
        break;
      }
    }
  }

  try {
    if (thermo_database_file.size() != 0) {
      chem = new SimpleThermoDatabase();
      chem->verbosity(verbosity);

      Beaker::BeakerParameters parameters = chem->GetDefaultParameters();
      parameters.thermo_database_file = thermo_database_file;
      parameters.activity_model_name = activity_model_name;
      parameters.porosity = 0.5;  // -
      parameters.saturation = 1.0;  // -
      parameters.volume = 1.0;  // m^3
      ModelSpecificParameters(model, &parameters);

      if (components.free_ion.size() == 0) {
        // initialize free-ion concentrations, these are actual
        // concentrations, so the value must be positive or ln(free_ion)
        // will return a nan!
        components.free_ion.resize(components.total.size(), 1.0e-9);
      }

      chem->Setup(components, parameters);

      if (verbosity >= kVerbose) {
        chem->Display();
        chem->DisplayComponents(components);
      }

      // solve for free-ion concentrations
      chem->Speciate(components, parameters);
      if (verbosity >= kTerse) {
        chem->DisplayResults();
      }

      if (num_time_steps != 0) {
        std::cout << "-- Test Beaker Reaction Stepping -------------------------------------" << std::endl;
        chem->DisplayTotalColumnHeaders();
        chem->DisplayTotalColumns(0.0, components.total);
        //parameters.max_iterations = 2;
        for (int time_step = 0; time_step < num_time_steps; time_step++) {
          chem->ReactionStep(&components, parameters, delta_time);
          if ((time_step+1) % output_interval == 0) {
            chem->DisplayTotalColumns((time_step+1) * delta_time, 
                                      components.total);
          }
        }
        std::cout << "---- Final Speciation" << std::endl;
        chem->Speciate(components, parameters);
        if (verbosity >= kTerse) {
          chem->DisplayResults();
        }
      }
    }
  }
  catch (const ChemistryException& geochem_error) {
    std::cout << geochem_error.what() << std::endl;
    error = EXIT_FAILURE;
  }
  catch (const std::runtime_error& rt_error) {
    std::cout << rt_error.what() << std::endl;
    error = EXIT_FAILURE;
  }
  catch (const std::logic_error& lg_error) {
    std::cout << lg_error.what() << std::endl;
    error = EXIT_FAILURE;
  }

  // cleanup memory
  delete chem;

  std::cout << "Done!\n";
  return error;
}  // end main()


void ModelSpecificParameters(const std::string model, 
                             Beaker::BeakerParameters* parameters)
{
  if (model == kCrunch) {
    parameters->water_density = 1000.0; // kg / m^3
  } else if (model == kPflotran) {
    parameters->water_density = 997.16; // kg / m^3
    // where did this number come from? 
    // default parameters->water_density = 997.205133945901; // kg / m^3
  } else {
    // bad model name, how did we get here....
  }
}  // end ModelSpecificParameters()

void PrintDoubleVector(const std::vector<double> &total)
{
  std::cout << "[ ";
  std::vector<double>::const_iterator i;
  for (i = total.begin(); i != total.end(); i++) {
    std::cout << std::scientific << std::setprecision(10) << *i << ", ";
  }
  std::cout << " ]" << std::endl;
}  // end PrintDoubleVector()

int CommandLineOptions(int argc, char **argv, 
                       Verbosity* verbosity, int* test, std::string* model)
{
  int error = -2;
  int option;
  extern char *optarg;

  while ((option = getopt(argc, argv, "m:ht:v:?")) != -1) {
    switch (option) {
      case 'm': {
        model->assign(optarg);
        error = EXIT_SUCCESS;
        break;
      }
      case 't': {
        /* specify the test that should be run */
        *test = std::atoi(optarg);
        error = EXIT_SUCCESS;
        break;
      }
      case 'v': {
        *verbosity = static_cast<Verbosity>(std::atoi(optarg));
        break;
      }
      case '?': case 'h': {  /* help mode */
        /* print some help stuff and exit without doing anything */
        std::cout << argv[0] << " command line options:" << std::endl;
        std::cout << "    -m string " << std::endl;
        std::cout << "         adjusts parameters to match validation models. Options:" << std::endl;
        std::cout << "             " << kCrunch << std::endl;
        std::cout << "             " << kPflotran << std::endl;
        std::cout << std::endl;
        std::cout << "    -t integer " << std::endl;
        std::cout << "         run a test case. valid test numbers are: " 
                  << std::endl;
        std::cout << "             0: error test" << std::endl;
        std::cout << "             1: simple carbonate speciation, unit activity coeff" << std::endl;
        std::cout << "             2: simple carbonate speciation, debye-huckel" << std::endl;
        std::cout << "             3: larger carbonate speciation, unit activity coeff" << std::endl;
        std::cout << "             4: larger carbonate speciation, debye-huckel" << std::endl;
        std::cout << "             5: calcite kinetics, TST rate law" << std::endl;
        std::cout << "             6: Na+ / Ca++ ion exchange" << std::endl;
        std::cout << "             7: fbasin initial condition" << std::endl;
        std::cout << "             8: fbasin initial condition with minerals" << std::endl;
        std::cout << "             9: fbasin initial condition with mineral kinetics" << std::endl;
        std::cout << "            10: calcite kinetics with large time steps" << std::endl;
        std::cout << "            11: surface complexation" << std::endl;
        std::cout << std::endl;
        std::cout << "    -v integer" << std::endl;
        std::cout << "         verbose output:" << std::endl;
        std::cout << "             0: silent" << std::endl;
        std::cout << "             1: terse" << std::endl;
        std::cout << "             2: verbose" << std::endl;
        std::cout << "             3: debug" << std::endl;
        std::cout << "             4: debug beaker" << std::endl;
        std::cout << "             5: debug input file" << std::endl;
        std::cout << "             6: debug mineral kinetics" << std::endl;
        std::cout << "             7: debug ion exchange" << std::endl;
        error = -1;
        break;
      }
      default: {
        /* no options */
        break;
      }
    }
  }

  if (error != -1 && *test == 0) {
    std::cout << "No test number specified on command line. Try \""
              <<  argv[0] << " -h \" for help." << std::endl;
  }

  if (model->c_str() != kCrunch && model->c_str() != kPflotran) {
    std::cout << "Invalid model name \'" << *model << "\'." << std::endl;
    std::cout << "Run \"" <<  argv[0] << " -h \" for help." << std::endl;
  }

  if (*verbosity >= kVerbose) {
    std::cout << "Command Line Options: " << std::endl;
    std::cout << "\tModel: " << *model << std::endl;
    std::cout << "\tTest number: " << *test << std::endl;
    std::cout << "\tVerbosity: " << *verbosity << std::endl;
  }
  std::cout << std::endl << std::endl;

  return error;
}  // end commandLineOptions()


void calcite_kinetics(const Verbosity& verbosity,
                      std::string* thermo_database_file,
                      std::string* activity_model_name,
                      Beaker::BeakerComponents* components,
                      double* delta_time,
                      int* num_time_steps,
                      int* output_interval)
{
        // calcite TST kinetics
        if (verbosity == kTerse) {
          std::cout << "Running calcite kinetics tst problem." << std::endl;
        }
        *thermo_database_file = "input/calcite.bgd";
        *activity_model_name = ActivityModelFactory::debye_huckel;
        components->total.push_back(-2.0e-5);  // H+
        components->total.push_back(1.0e-3);  // HCO3-
        components->total.push_back(1.0e-4);  // Ca++
        components->minerals.push_back(0.2);  // calcite
        *delta_time = 0.5;
        *num_time_steps = 500;
        *output_interval = 1;
}  // end calcite_kinetics()


void calcite_kinetics_large_time_steps(const Verbosity& verbosity,
                                       std::string* thermo_database_file,
                                       std::string* activity_model_name,
                                       Beaker::BeakerComponents* components,
                                       double* delta_time,
                                       int* num_time_steps,
                                       int* output_interval)
{
        // calcite TST kinetics
        if (verbosity == kTerse) {
          std::cout << "Running calcite kinetics tst problem." << std::endl;
        }
        *thermo_database_file = "input/calcite.bgd";
        *activity_model_name = ActivityModelFactory::debye_huckel;
        components->total.push_back(1.0e-3);  // H+
        components->total.push_back(1.0e-3);  // HCO3-
        components->total.push_back(0.0);  // Ca++
        components->minerals.push_back(0.2);  // calcite
        *delta_time = 30.0 * 24.0 * 3600.0;
        *num_time_steps = 12;
        *output_interval = 1;
}  // end calcite_kinetics_large_time_steps()


void surface_complexation(const Verbosity& verbosity,
                      std::string* thermo_database_file,
                      std::string* activity_model_name,
                      Beaker::BeakerComponents* components,
                      double* delta_time,
                      int* num_time_steps,
                      int* output_interval)
{
        // calcite TST kinetics
        if (verbosity == kTerse) {
          std::cout << "Running surface complexation problem." << std::endl;
        }
        *thermo_database_file = "input/surface-complexation.bgd";
        *activity_model_name = ActivityModelFactory::debye_huckel;
        components->total.push_back(1.0e-3);  // H+
        components->total.push_back(3.0e-3);  // Cl-
        components->total.push_back(1.0e-3);  // UO2++
        components->total_sorbed.push_back(1.0e-4);  // H+
        components->total_sorbed.push_back(1.0e-4);  // Cl-
        components->total_sorbed.push_back(1.0e-4);  // UO2++
        *delta_time = 1.0;
        *num_time_steps = 10;
        *output_interval = 1;
}  // end surface_complexation()

void fbasin_initial_all(const Verbosity& verbosity,
                        std::string* thermo_database_file,
                        std::string* activity_model_name,
                        Beaker::BeakerComponents* components)
{
  fbasin_initial_speciation(verbosity, thermo_database_file,
                            activity_model_name, components);
  fbasin_minerals(verbosity, thermo_database_file,
                  activity_model_name, components);

  *thermo_database_file = "input/fbasin-initial-mineral-kinetics.bgd";
}  // end fbasin_initial_full

void fbasin_initial_speciation(const Verbosity& verbosity,
                               std::string* thermo_database_file,
                               std::string* activity_model_name,
                               Beaker::BeakerComponents* components)
{
        if (verbosity == kTerse) {
          std::cout << "Running fbasin speciation problem, \'initial\' conditions." << std::endl;
        }
        *thermo_database_file = "input/fbasin-initial.bgd";
        *activity_model_name = ActivityModelFactory::debye_huckel;
        components->total.push_back(1.0000E-05);  // Na+
        components->total.push_back(1.0000E-05);  // Ca++
        components->total.push_back(8.4757E-08);  // Fe++
        components->total.push_back(1.9211E-04);  // K+
        components->total.push_back(6.6779E-09);  // Al+++
        components->total.push_back(1.2683E-05);  // H+
        components->total.push_back(1.0000E-05);  // N2(aq)
        components->total.push_back(1.0000E-05);  // NO3-
        components->total.push_back(2.0999E-04);  // HCO3-
        components->total.push_back(1.0000E-05);  // Cl-
        components->total.push_back(1.0000E-06);  // SO4--
        components->total.push_back(1.0000E-06);  // HPO4--
        components->total.push_back(1.0000E-06);  // F-
        components->total.push_back(1.8703E-04);  // SiO2(aq)
        components->total.push_back(1.0000E-15);  // UO2++
        components->total.push_back(2.5279E-04);  // O2(aq)
        components->total.push_back(1.0000E-15);  // Tracer

}  // end fbasin_initial_speciation()

void fbasin_minerals(const Verbosity& verbosity,
                               std::string* thermo_database_file,
                               std::string* activity_model_name,
                               Beaker::BeakerComponents* components)
{
        if (verbosity == kTerse) {
          std::cout << "Adding fbasin kinetic minerals." << std::endl;
        }
        *thermo_database_file = "input/fbasin-initial-minerals.bgd";
        components->minerals.push_back(0.0);  // Gibbsite
        components->minerals.push_back(0.21);  // Quartz
        components->minerals.push_back(0.15);  // K-Feldspar
        components->minerals.push_back(0.0);  // Jurbanite
        components->minerals.push_back(0.1);  // Ferrihydrite
        components->minerals.push_back(0.15);  // Kaolinite
        components->minerals.push_back(0.0);  // Schoepite
        components->minerals.push_back(0.0);  // (UO2)3(PO4)2.4H2O
        components->minerals.push_back(0.0);  // Soddyite
        components->minerals.push_back(0.0);  // Calcite
        components->minerals.push_back(0.0);  // Chalcedony

}  // end fbasin_minerals()
