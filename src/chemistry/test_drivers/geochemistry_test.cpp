#include <cstdlib>
#include <iostream>
#include <vector>


#include "SimpleCarbonate.hpp"
//#include "LargeCarbonate.hpp"
#include "Geochemistry.hpp"

int commandLineOptions(int argc, char **argv, int& verbose, int& test);


int main (int argc, char **argv) {

  int verbose = 1;
  int test = 0;
  int error = EXIT_SUCCESS;

  Geochemistry *chem = NULL;
  error = commandLineOptions(argc, argv, verbose, test);

  if (error == EXIT_SUCCESS) {
    switch (test){
    case 1:
      // set up simple 2-species carbonate system (H,HCO3-)
      if (verbose > 0) {
	std::cout << "Running simple carbonate example." << std::endl;
      }
      chem = new SimpleCarbonate();
      break;
    case 2:
      // larger carbonate system, 3 components, 9 secondary
      if (verbose > 0) {
	std::cout << "Running large carbonate speciation example." << std::endl;
      }
      //    chem = new ComplexCarbonate();
      break;
    default:
      break;
    }
  }

  if (chem != NULL) {
    std::vector<double> total;
    chem->verbose(verbose);
    chem->setup(&total);
    // solve for free-ion concentrations
    chem->speciate(total);
    if (verbose > 1) {
      chem->print_results();
    }
  }

  if (chem != NULL) {
    delete chem;
  }
  std::cout << "Done!\n";

}


int commandLineOptions(int argc, char **argv, int& verbose, int& test)
{
  int error = -2;
  int option;
  extern char *optarg;

  while ((option = getopt(argc, argv, "ht:v:?")) != EOF) {
    switch (option) {
    case 't':
      /* specify the test that should be run */
      test = std::atoi(optarg);
      error = EXIT_SUCCESS;
      break;
    case 'v':
      verbose = std::atoi(optarg);
      break;
    case '?': case 'h':  /* help mode */
      /* print some help stuff and exit without doing anything */
      std::cout << argv[0] << " command line options:" << std::endl;
      std::cout << "    -t integer " << std::endl
		<< "         run a test case. valid test numbers are: " << std::endl
		<< "             1: simple carbonate speciation" << std::endl
		<< "             2: larger carbonate speciation" << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "    -v integer" << std::endl;
      std::cout << "         verbose output:" << std::endl;
      std::cout << "             0: silent" << std::endl;
      std::cout << "             1: terse" << std::endl;
      std::cout << "             2: verbose" << std::endl;
      std::cout << "             3: debug" << std::endl;
      error = -1;
      break;
    default:
      /* no options */
      break;
    }
  }

  if (error != -1 && test == 0) {
    std::cout << "No test number specified on command line. Try \"" 
	      <<  argv[0] << " -h \" for help." << std::endl;
  }

  if (verbose > 1) {
    std::cout << "Command Line Options: " << std::endl;
    std::cout << "\tTest number: " << test << std::endl;
    std::cout << "\tVerbosity: " << verbose << std::endl;
  }
  std::cout << std::endl << std::endl;

  return error; /* end commandLineOptions() */
}

