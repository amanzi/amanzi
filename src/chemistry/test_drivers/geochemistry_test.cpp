/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>

#include <iostream>
#include <vector>
#include <string>

#include "SimpleCarbonate.hpp"
#include "LargeCarbonate.hpp"
#include "Beaker.hpp"
#include "ActivityModelFactory.hpp"
#include "Verbosity.hpp"

int CommandLineOptions(int argc, char **argv, Verbosity& verbosity, int& test);

int main(int argc, char **argv) {
  Verbosity verbosity = kTerse;
  int test = 0;
  int error = EXIT_SUCCESS;

  Beaker *chem = NULL;
  std::string activity_model_name;
  std::vector<KineticRate*> mineral_rates;
  mineral_rates.clear();
  std::string mineral_kinetics_file;
  mineral_kinetics_file.clear();

  error = CommandLineOptions(argc, argv, verbosity, test);

  if (error == EXIT_SUCCESS) {
    switch (test) {
      case 1: {
        // set up simple 2-species carbonate system (H,HCO3-) unit activity coefficients
        if (verbosity >= kTerse) {
          std::cout << "Running simple carbonate example, unit activity coefficients." << std::endl;
        }
        chem = new SimpleCarbonate();
        activity_model_name = ActivityModelFactory::unit;
        break;
      }
      case 2: {
        // set up simple 2-species carbonate system (H,HCO3-), debye-huckel activity coefficients
        if (verbosity >= kTerse) {
          std::cout << "Running simple carbonate example, debye-huckel." << std::endl;
        }
        chem = new SimpleCarbonate();
        activity_model_name = ActivityModelFactory::debye_huckel;
        break;
      }
      case 3: {
        // larger carbonate system, 3 components, 9 secondary, unit activity coefficients
        if (verbosity >= kTerse) {
          std::cout << "Running large carbonate speciation example, unit activity coefficients." << std::endl;
        }
        chem = new LargeCarbonate();
        activity_model_name = ActivityModelFactory::unit;
        break;
      }
      case 4: {
        // larger carbonate system, 3 components, 9 secondary, debye-huckel activity coefficients
        if (verbosity >= kTerse) {
          std::cout << "Running large carbonate speciation example, debye-huckel activity coefficients." << std::endl;
        }
        chem = new LargeCarbonate();
        activity_model_name = ActivityModelFactory::debye_huckel;
        break;
      }
      case 5: {
        // calcite TST kinetics
        if (verbosity >= kTerse) {
          std::cout << "Running calcite kinetics tst problem." << std::endl;
        }
        chem = new LargeCarbonate();
        activity_model_name = ActivityModelFactory::debye_huckel;
        mineral_kinetics_file = "calcite-kinetics-tst.txt";
      }
      default: {
        std::cout << "Invalid test number specified on command line. try using the \'-h\' option." << std::endl;
        break;
      }
    }
  }

  if (chem != NULL) {
    std::vector<double> total;
    chem->verbosity(verbosity);
    chem->SetupActivityModel(activity_model_name);
    chem->setup(total, mineral_kinetics_file);
    if (verbosity >= kVerbose) {
      chem->display();
    }

    // solve for free-ion concentrations
    double water_density = 997.205133945901; // kg / m^3
    chem->speciate(total, water_density);
    if (verbosity >= kTerse) {
      chem->print_results();
    }
  }



  // cleanup memory
  if (chem != NULL) {
    delete chem;
  }


  std::cout << "Done!\n";
}  // end main()


int CommandLineOptions(int argc, char **argv, Verbosity& verbosity, int& test)
{
  int error = -2;
  int option;
  extern char *optarg;

  while ((option = getopt(argc, argv, "ht:v:?")) != EOF) {
    switch (option) {
      case 't': {
        /* specify the test that should be run */
        test = std::atoi(optarg);
        error = EXIT_SUCCESS;
        break;
      }
      case 'v': {
        verbosity = static_cast<Verbosity>(std::atoi(optarg));
        break;
      }
      case '?': case 'h': {  /* help mode */
        /* print some help stuff and exit without doing anything */
        std::cout << argv[0] << " command line options:" << std::endl;
        std::cout << "    -t integer " << std::endl
                  << "         run a test case. valid test numbers are: " << std::endl
                  << "             1: simple carbonate speciation, unit activity coeff" << std::endl
                  << "             2: simple carbonate speciation, debye-huckel" << std::endl
                  << "             3: larger carbonate speciation, unit activity coeff" << std::endl
                  << "             4: larger carbonate speciation, debye-huckel" << std::endl
                  << "             5: calcite kinetics, TST rate law" << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "    -v integer" << std::endl;
        std::cout << "         verbose output:" << std::endl;
        std::cout << "             0: silent" << std::endl;
        std::cout << "             1: terse" << std::endl;
        std::cout << "             2: verbose" << std::endl;
        std::cout << "             3: debug" << std::endl;
        std::cout << "             4: debug beaker" << std::endl;
        std::cout << "             5: debug mineral kinetics" << std::endl;
        error = -1;
        break;
      }
      default: {
        /* no options */
        break;
      }
    }
  }

  if (error != -1 && test == 0) {
    std::cout << "No test number specified on command line. Try \""
              <<  argv[0] << " -h \" for help." << std::endl;
  }

  if (verbosity >= kVerbose) {
    std::cout << "Command Line Options: " << std::endl;
    std::cout << "\tTest number: " << test << std::endl;
    std::cout << "\tVerbosity: " << verbosity << std::endl;
  }
  std::cout << std::endl << std::endl;

  return error;
}  // end commandLineOptions()

