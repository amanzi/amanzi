/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "ThermoDatabase.hpp"
#include "Beaker.hpp"
#include "ActivityModelFactory.hpp"
#include "Verbosity.hpp"

const std::string kCrunch("crunch");
const std::string kPflotran("pflotran");

int CommandLineOptions(int argc, char **argv, Verbosity* verbosity, int* test, std::string* model);
void ModelSpecificParameters(const std::string model, Beaker::BeakerParameters* parameters);
void PrintDoubleVector(const std::vector<double> &total);

int main(int argc, char **argv) {
  Verbosity verbosity = kTerse;
  int test = 0;
  std::string model("");
  int error = EXIT_SUCCESS;

  Beaker *chem = NULL;
  std::vector<double> total;
  total.clear();

  std::string thermo_database_file("");
  std::string mineral_kinetics_file("");
  std::string activity_model_name("");
  
  error = CommandLineOptions(argc, argv, &verbosity, &test, &model);

  if (error == EXIT_SUCCESS) {
    switch (test) {
      case 1: {
        // set up simple 2-species carbonate system (H,HCO3-) unit activity coefficients
        if (verbosity == kTerse) {
          std::cout << "Running simple carbonate example, unit activity coefficients." << std::endl;
        }
        thermo_database_file = "input/carbonate.bgd";
        activity_model_name = ActivityModelFactory::unit;
        total.push_back(1.0e-3);  // H+
        total.push_back(1.0e-3);  // HCO3-
        break;
      }
      case 2: {
        // set up simple 2-species carbonate system (H,HCO3-), debye-huckel activity coefficients
        if (verbosity == kTerse) {
          std::cout << "Running simple carbonate example, debye-huckel." << std::endl;
        }
        thermo_database_file = "input/carbonate.bgd";
        activity_model_name = ActivityModelFactory::debye_huckel;
        total.push_back(1.0e-3);  // H+
        total.push_back(1.0e-3);  // HCO3-
        break;
      }
      case 3: {
        // larger carbonate system, 3 components, 9 secondary, unit activity coefficients
        if (verbosity == kTerse) {
          std::cout << "Running large carbonate speciation example, unit activity coefficients." << std::endl;
        }
        thermo_database_file = "input/ca-carbonate.bgd";
        activity_model_name = ActivityModelFactory::unit;
        total.push_back(1.0e-3);  // H+
        total.push_back(3.0e-3);  // HCO3-
        total.push_back(1.0e-3);  // Ca++
        break;
      }
      case 4: {
        // larger carbonate system, 3 components, 9 secondary, debye-huckel activity coefficients
        if (verbosity == kTerse) {
          std::cout << "Running large carbonate speciation example, debye-huckel activity coefficients." << std::endl;
        }
        thermo_database_file = "input/ca-carbonate.bgd";
        activity_model_name = ActivityModelFactory::debye_huckel;
        total.push_back(1.0e-3);  // H+
        total.push_back(3.0e-3);  // HCO3-
        total.push_back(1.0e-3);  // Ca++
        break;
      }
      case 5: {
        // calcite TST kinetics
        if (verbosity == kTerse) {
          std::cout << "Running calcite kinetics tst problem." << std::endl;
        }
        thermo_database_file = "input/calcite.bgd";
        activity_model_name = ActivityModelFactory::debye_huckel;
        mineral_kinetics_file = "input/calcite-kinetics-tst.ain";
        total.push_back(1.0e-3);  // H+
        total.push_back(3.0e-3);  // HCO3-
        total.push_back(1.0e-3);  // Ca++
        break;
      }
      case 6: {
        // Na+ Ca++ ion exchange
        if (verbosity == kTerse) {
          std::cout << "Running Na-Ca ion exchange problem." << std::endl;
        }
        thermo_database_file = "input/na-ca-ion-exchange.bgd";
        activity_model_name = ActivityModelFactory::debye_huckel;
        total.push_back(1.0e-3);  // H+
        total.push_back(3.0e-3);  // HCO3-
        total.push_back(1.0e-3);  // Ca++
        total.push_back(1.0e-3);  // Na+
        total.push_back(1.0e-3);  // Cl-
        break;
      }
      default: {
        std::cout << "Invalid test number specified on command line. try using the \'-h\' option." << std::endl;
        break;
      }
    }
  }

  if (thermo_database_file.size() != 0) {
    chem = new ThermoDatabase();
    chem->verbosity(verbosity);
    Beaker::BeakerParameters parameters = chem->GetDefaultParameters();
    parameters.thermo_database_file = thermo_database_file;
    parameters.mineral_kinetics_file = mineral_kinetics_file;
    parameters.activity_model_name = activity_model_name;
    parameters.porosity = 0.5;  // -
    parameters.saturation = 1.0;  // - 
    parameters.volume = 1.0;  // m^3
    ModelSpecificParameters(model, &parameters);
    chem->setup(total, parameters);
    if (verbosity >= kVerbose) {
      chem->Display();
    }

    // solve for free-ion concentrations
    chem->speciate(total, parameters.water_density);
    if (verbosity >= kTerse) {
      chem->DisplayResults();
    }
  
    if (mineral_kinetics_file.size() != 0) {
      std::cout << "-- Test Beaker Reaction Stepping -------------------------------------" << std::endl;
      chem->DisplayTotalColumnHeaders();
      chem->DisplayTotalColumns(0.0, total);
      double delta_time = 3660.0;  // seconds
      int num_time_steps = 12;
      for (int time_step = 0; time_step <= num_time_steps; time_step++) {
        chem->ReactionStep(total, parameters, delta_time);        
        chem->DisplayTotalColumns(time_step+1 * delta_time, total);
      }
    }
  }
  // cleanup memory
  if (chem != NULL) {
    delete chem;
  }

  std::cout << "Done!\n";
}  // end main()


void ModelSpecificParameters(const std::string model, Beaker::BeakerParameters* parameters)
{
  if (model == kCrunch) {
    parameters->water_density = 1000.0; // kg / m^3    
  } else if (model == kPflotran) {
    parameters->water_density = 997.16; // kg / m^3
    // where did this number come from? parameters->water_density = 997.205133945901; // kg / m^3
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

int CommandLineOptions(int argc, char **argv, Verbosity* verbosity, int* test,
                       std::string* model)
{
  int error = -2;
  int option;
  extern char *optarg;

  while ((option = getopt(argc, argv, "m:ht:v:?")) != EOF) {
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
        std::cout << "         run a test case. valid test numbers are: " << std::endl;
        std::cout << "             1: simple carbonate speciation, unit activity coeff" << std::endl;
        std::cout << "             2: simple carbonate speciation, debye-huckel" << std::endl;
        std::cout << "             3: larger carbonate speciation, unit activity coeff" << std::endl;
        std::cout << "             4: larger carbonate speciation, debye-huckel" << std::endl;
        std::cout << "             5: calcite kinetics, TST rate law" << std::endl;
        std::cout << "             6: Na+ / Ca++ ion exchange" << std::endl;
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
    std::cout << "Invalid model name \'" << model << "\'." << std::endl;
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

