/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "geochemistry_test.hh"

#include <unistd.h>

#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>

#include "simple_thermo_database.hh"
#include "beaker.hh"
#include "activity_model_factory.hh"
#include "verbosity.hh"
#include "chemistry_exception.hh"

namespace ac = amanzi::chemistry;

const std::string kCrunch("crunch");
const std::string kPflotran("pflotran");


int main(int argc, char** argv) {
  ac::Verbosity verbosity = ac::kTerse;
  int test = 0;
  std::string model("");
  int error = EXIT_SUCCESS;

  ac::Beaker* chem = NULL;
  ac::Beaker::BeakerComponents components;
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
        if (verbosity == ac::kTerse) {
          std::cout << "Running component size sanity check." << std::endl;
        }
        thermo_database_file = "input/carbonate.bgd";
        activity_model_name = ac::ActivityModelFactory::unit;
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
        if (verbosity == ac::kTerse) {
          std::cout << "Running simple carbonate example, "
                    << "unit activity coefficients." << std::endl;
        }
        thermo_database_file = "input/carbonate.bgd";
        activity_model_name = ac::ActivityModelFactory::unit;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(1.0e-3);  // HCO3-
        break;
      }
      case 2: {
        // set up simple 2-species carbonate system (H,HCO3-),
        // debye-huckel activity coefficients
        if (verbosity == ac::kTerse) {
          std::cout << "Running simple carbonate example, debye-huckel."
                    << std::endl;
        }
        thermo_database_file = "input/carbonate.bgd";
        activity_model_name = ac::ActivityModelFactory::debye_huckel;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(1.0e-3);  // HCO3-
        break;
      }
      case 3: {
        // larger carbonate system, 3 components, 9 secondary,
        // unit activity coefficients
        if (verbosity == ac::kTerse) {
          std::cout << "Running large carbonate speciation example, "
                    << "unit activity coefficients." << std::endl;
        }
        thermo_database_file = "input/ca-carbonate.bgd";
        activity_model_name = ac::ActivityModelFactory::unit;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(3.0e-3);  // HCO3-
        components.total.push_back(1.0e-3);  // Ca++
        break;
      }
      case 4: {
        // larger carbonate system, 3 components, 9 secondary,
        // debye-huckel activity coefficients
        if (verbosity == ac::kTerse) {
          std::cout << "Running large carbonate speciation example, "
                    << "debye-huckel activity coefficients." << std::endl;
        }
        thermo_database_file = "input/ca-carbonate.bgd";
        activity_model_name = ac::ActivityModelFactory::debye_huckel;
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
        if (verbosity == ac::kTerse) {
          std::cout << "Running Na-Ca ion exchange problem." << std::endl;
        }
        thermo_database_file = "input/na-ca-ion-exchange.bgd";
        activity_model_name = ac::ActivityModelFactory::debye_huckel;
        components.total.push_back(1.0e-3);  // H+
        components.total.push_back(3.0e-3);  // HCO3-
        components.total.push_back(1.0e-3);  // Ca++
        components.total.push_back(1.0e-3);  // Na+
        components.total.push_back(1.0e-3);  // Cl-
        components.ion_exchange_sites.push_back(15.0);  // X-, units? equivalents per 100 grams solid?
        break;
      }
      case 7: {
        // fbasin 'initial' condition
        fbasin_initial(verbosity,
                       &thermo_database_file,
                       &activity_model_name,
                       &components,
                       &delta_time,
                       &num_time_steps,
                       &output_interval);
        break;
      }
      case 8: {
        // fbasin 'infiltration' condition
        fbasin_infiltration(verbosity,
                            &thermo_database_file,
                            &activity_model_name,
                            &components,
                            &delta_time,
                            &num_time_steps,
                            &output_interval);
        break;
      }
      case 9: {
        // fbasin 'source' condition
        fbasin_source(verbosity,
                      &thermo_database_file,
                      &activity_model_name,
                      &components,
                      &delta_time,
                      &num_time_steps,
                      &output_interval);
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
      case 12: {
        // empty test slot
        break;
      }
      case 13: {
        // uo2_5_component 'initial' condition
        uo2_5_component_initial(verbosity,
                                &thermo_database_file,
                                &activity_model_name,
                                &components,
                                &delta_time,
                                &num_time_steps,
                                &output_interval);
        break;
      }
      case 14: {
        // uo2_5_component 'outlet' condition
        uo2_5_component_outlet(verbosity,
                               &thermo_database_file,
                               &activity_model_name,
                               &components,
                               &delta_time,
                               &num_time_steps,
                               &output_interval);
        break;
      }
      case 15: {
        // uo2_5_component 'source' condition
        uo2_5_component_source(verbosity,
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
      chem = new ac::SimpleThermoDatabase();
      chem->verbosity(verbosity);

      ac::Beaker::BeakerParameters parameters = chem->GetDefaultParameters();
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

      if (verbosity >= ac::kVerbose) {
        chem->Display();
        chem->DisplayComponents(components);
      }

      // solve for free-ion concentrations
      chem->Speciate(components, parameters);
      chem->UpdateComponents(&components);
      if (verbosity >= ac::kTerse) {
        chem->DisplayResults();
      }

      if (num_time_steps != 0) {
        std::cout << "-- Test Beaker Reaction Stepping -------------------------------------" << std::endl;
        chem->DisplayTotalColumnHeaders();
        chem->DisplayTotalColumns(0.0, components.total);
        // parameters.max_iterations = 2;
        for (int time_step = 0; time_step < num_time_steps; time_step++) {
          chem->ReactionStep(&components, parameters, delta_time);
          if ((time_step + 1) % output_interval == 0) {
            chem->DisplayTotalColumns((time_step + 1) * delta_time,
                                      components.total);
          }
          if (verbosity >= ac::kDebugNewtonSolver) {
            ac::Beaker::SolverStatus status = chem->status();
            std::cout << "Timestep: " << time_step << std::endl;
            std::cout << "    number of rhs evaluations: " << status.num_rhs_evaluations << std::endl;
            std::cout << "    number of jacobian evaluations: " << status.num_jacobian_evaluations << std::endl;
            std::cout << "    number of newton iterations: " << status.num_newton_iterations << std::endl;
            std::cout << "    solution converged: " << status.converged << std::endl;
          }
        }
        std::cout << "---- Final Speciation" << std::endl;
        chem->Speciate(components, parameters);
        if (verbosity >= ac::kTerse) {
          chem->DisplayResults();
        }
      }
    }
  } catch (const ac::ChemistryException& geochem_error) {
    std::cout << geochem_error.what() << std::endl;
    error = EXIT_FAILURE;
  } catch (const std::runtime_error& rt_error) {
    std::cout << rt_error.what() << std::endl;
    error = EXIT_FAILURE;
  } catch (const std::logic_error& lg_error) {
    std::cout << lg_error.what() << std::endl;
    error = EXIT_FAILURE;
  }

  // cleanup memory
  delete chem;

  std::cout << "Done!\n";
  return error;
}  // end main()


void ModelSpecificParameters(const std::string model,
                             ac::Beaker::BeakerParameters* parameters) {
  if (model == kCrunch) {
    parameters->water_density = 1000.0;  // kg / m^3
  } else if (model == kPflotran) {
    parameters->water_density = 997.16;  // kg / m^3
    // where did this number come from?
    // default parameters->water_density = 997.205133945901;  // kg / m^3
  } else {
    // bad model name, how did we get here....
  }
}  // end ModelSpecificParameters()

void PrintDoubleVector(const std::vector<double> &total) {
  std::cout << "[ ";
  std::vector<double>::const_iterator i;
  for (i = total.begin(); i != total.end(); i++) {
    std::cout << std::scientific << std::setprecision(10) << *i << ", ";
  }
  std::cout << " ]" << std::endl;
}  // end PrintDoubleVector()

int CommandLineOptions(int argc, char** argv,
                       ac::Verbosity* verbosity, int* test, std::string* model) {
  int error = -2;
  int option;
  extern char* optarg;

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
        *verbosity = static_cast<ac::Verbosity>(std::atoi(optarg));
        break;
      }
      case '?':
      case 'h': {  /* help mode */
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
        std::cout << "             7: fbasin 17 component initial condition" << std::endl;
        std::cout << "             8: fbasin 17 component infiltration condition" << std::endl;
        std::cout << "             9: fbasin 17 component source condition" << std::endl;
        std::cout << "            10: calcite kinetics with large time steps" << std::endl;
        std::cout << "            11: surface complexation" << std::endl;
        std::cout << "            12: empty test slot" << std::endl;
        std::cout << "            13: UO2 5 component initial condition" << std::endl;
        std::cout << "            14: UO2 5 component source condition" << std::endl;
        std::cout << "            15: UO2 5 component outlet condition" << std::endl;
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
        std::cout << "             8: debug newton solver" << std::endl;
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

  if (*verbosity >= ac::kVerbose) {
    std::cout << "Command Line Options: " << std::endl;
    std::cout << "\tModel: " << *model << std::endl;
    std::cout << "\tTest number: " << *test << std::endl;
    std::cout << "\tVerbosity: " << *verbosity << std::endl;
  }
  std::cout << std::endl << std::endl;

  return error;
}  // end commandLineOptions()


void calcite_kinetics(const ac::Verbosity& verbosity,
                      std::string* thermo_database_file,
                      std::string* activity_model_name,
                      ac::Beaker::BeakerComponents* components,
                      double* delta_time,
                      int* num_time_steps,
                      int* output_interval) {
  // calcite TST kinetics
  if (verbosity == ac::kTerse) {
    std::cout << "Running calcite kinetics tst problem." << std::endl;
  }
  *thermo_database_file = "input/calcite.bgd";
  *activity_model_name = ac::ActivityModelFactory::debye_huckel;
  components->total.push_back(-2.0e-5);  // H+
  components->total.push_back(1.0e-3);  // HCO3-
  components->total.push_back(1.0e-4);  // Ca++
  components->minerals.push_back(0.2);  // calcite
  *delta_time = 0.5;
  *num_time_steps = 500;
  *output_interval = 1;
}  // end calcite_kinetics()


void calcite_kinetics_large_time_steps(const ac::Verbosity& verbosity,
                                       std::string* thermo_database_file,
                                       std::string* activity_model_name,
                                       ac::Beaker::BeakerComponents* components,
                                       double* delta_time,
                                       int* num_time_steps,
                                       int* output_interval) {
  // calcite TST kinetics
  if (verbosity == ac::kTerse) {
    std::cout << "Running calcite kinetics tst problem." << std::endl;
  }
  *thermo_database_file = "input/calcite.bgd";
  *activity_model_name = ac::ActivityModelFactory::debye_huckel;
  components->total.push_back(1.0e-3);  // H+
  components->total.push_back(1.0e-3);  // HCO3-
  components->total.push_back(0.0);  // Ca++
  components->minerals.push_back(0.2);  // calcite
  *delta_time = 30.0 * 24.0 * 3600.0;
  *num_time_steps = 12;
  *output_interval = 1;
}  // end calcite_kinetics_large_time_steps()


void surface_complexation(const ac::Verbosity& verbosity,
                          std::string* thermo_database_file,
                          std::string* activity_model_name,
                          ac::Beaker::BeakerComponents* components,
                          double* delta_time,
                          int* num_time_steps,
                          int* output_interval) {
  // calcite TST kinetics
  if (verbosity == ac::kTerse) {
    std::cout << "Running surface complexation problem." << std::endl;
  }
  *thermo_database_file = "input/surface-complexation.bgd";
  *activity_model_name = ac::ActivityModelFactory::debye_huckel;
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

/*******************************************************************************
 **
 ** fbasin UO2 17 component problem
 **
 ******************************************************************************/
void fbasin_initial(const ac::Verbosity& verbosity,
                    std::string* thermo_database_file,
                    std::string* activity_model_name,
                    ac::Beaker::BeakerComponents* components,
                    double* delta_time,
                    int* num_time_steps,
                    int* output_interval) {
  if (verbosity == ac::kTerse) {
    std::cout << "Running fbasin problem, \'initial\' conditions." << std::endl;
  }
  *thermo_database_file = "input/fbasin-17.bgd";
  *activity_model_name = ac::ActivityModelFactory::debye_huckel;
  *delta_time = 30.0 * 24.0 * 3600.0;
  *num_time_steps = 12;
  *output_interval = 1;

  fbasin_aqueous_initial(components);
  fbasin_free_ions(components);
  fbasin_minerals(components);
  fbasin_sorbed(components);
}  // end fbasin_initial

void fbasin_infiltration(const ac::Verbosity& verbosity,
                         std::string* thermo_database_file,
                         std::string* activity_model_name,
                         ac::Beaker::BeakerComponents* components,
                         double* delta_time,
                         int* num_time_steps,
                         int* output_interval) {
  if (verbosity == ac::kTerse) {
    std::cout << "Running fbasin problem, \'infiltration\' conditions." << std::endl;
  }
  *thermo_database_file = "input/fbasin-17.bgd";
  *activity_model_name = ac::ActivityModelFactory::debye_huckel;
  *delta_time = 30.0 * 24.0 * 3600.0;
  *num_time_steps = 12;
  *output_interval = 1;

  fbasin_aqueous_infiltration(components);
  fbasin_free_ions(components);
  fbasin_minerals(components);
  fbasin_sorbed(components);
}  // end fbasin_infiltration

void fbasin_source(const ac::Verbosity& verbosity,
                   std::string* thermo_database_file,
                   std::string* activity_model_name,
                   ac::Beaker::BeakerComponents* components,
                   double* delta_time,
                   int* num_time_steps,
                   int* output_interval) {
  if (verbosity == ac::kTerse) {
    std::cout << "Running fbasin problem, \'source\' conditions." << std::endl;
  }
  *thermo_database_file = "input/fbasin-17.bgd";
  *activity_model_name = ac::ActivityModelFactory::debye_huckel;
  *delta_time = 30.0 * 24.0 * 3600.0;
  *num_time_steps = 12;
  *output_interval = 1;

  fbasin_aqueous_source(components);
  fbasin_free_ions(components);
  fbasin_minerals(components);
  fbasin_sorbed(components);
}  // end fbasin_source

void fbasin_aqueous_initial(ac::Beaker::BeakerComponents* components) {
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
}  // end fbasin_aqueous_initial

void fbasin_aqueous_infiltration(ac::Beaker::BeakerComponents* components) {
  // constraint: infiltration
  components->total.push_back(1.3132E-04);  // Na+
  components->total.push_back(1.0000E-05);  // Ca++
  components->total.push_back(1.0000E-12);  // Fe++
  components->total.push_back(1.0000E-06);  // K+
  components->total.push_back(1.0000E-12);  // Al+++
  components->total.push_back(1.0716E-05);  // H+
  components->total.push_back(1.0000E-05);  // N2(aq)
  components->total.push_back(1.0000E-05);  // NO3-
  components->total.push_back(6.0081E-05);  // HCO3-
  components->total.push_back(1.0000E-05);  // Cl-
  components->total.push_back(1.0000E-06);  // SO4--
  components->total.push_back(1.0000E-06);  // HPO4--
  components->total.push_back(7.8954E-05);  // F-
  components->total.push_back(1.0000E-05);  // SiO2(aq)
  components->total.push_back(1.0000E-15);  // UO2++
  components->total.push_back(2.5277E-04);  // O2(aq)
  components->total.push_back(1.0000E-15);  // Tracer
}  // end fbasin_aqueous_infiltration

void fbasin_aqueous_source(ac::Beaker::BeakerComponents* components) {
  // constraint: source
  components->total.push_back(3.4363E-02);  // Na+
  components->total.push_back(1.2475E-05);  // Ca++
  components->total.push_back(3.0440E-05);  // Fe++
  components->total.push_back(1.7136E-05);  // K+
  components->total.push_back(2.8909E-05);  // Al+++
  components->total.push_back(3.6351E-03);  // H+
  components->total.push_back(1.3305E-03);  // N2(aq)
  components->total.push_back(3.4572E-02);  // NO3-
  components->total.push_back(2.1830E-03);  // HCO3-
  components->total.push_back(3.3848E-05);  // Cl-
  components->total.push_back(6.2463E-04);  // SO4--
  components->total.push_back(7.1028E-05);  // HPO4--
  components->total.push_back(7.8954E-05);  // F-
  components->total.push_back(2.5280E-04);  // SiO2(aq)
  components->total.push_back(3.5414E-05);  // UO2++
  components->total.push_back(2.6038E-04);  // O2(aq)
  components->total.push_back(3.5414E-05);  // Tracer
}  // end fbasin_aqueous_source

void fbasin_free_ions(ac::Beaker::BeakerComponents* components) {
  // free ion concentrations (better initial guess)
  components->free_ion.push_back(9.9969E-06);  // Na+
  components->free_ion.push_back(9.9746E-06);  // Ca++
  components->free_ion.push_back(2.2405E-18);  // Fe++
  components->free_ion.push_back(1.8874E-04);  // K+
  components->free_ion.push_back(5.2970E-16);  // Al+++
  components->free_ion.push_back(3.2759E-08);  // H+
  components->free_ion.push_back(1.0000E-05);  // N2(aq)
  components->free_ion.push_back(1.0000E-05);  // NO3-
  components->free_ion.push_back(1.9282E-04);  // HCO3-
  components->free_ion.push_back(9.9999E-06);  // Cl-
  components->free_ion.push_back(9.9860E-07);  // SO4--
  components->free_ion.push_back(9.9886E-07);  // HPO4--
  components->free_ion.push_back(1.0000E-06);  // F-
  components->free_ion.push_back(1.8703E-04);  // SiO2(aq)
  components->free_ion.push_back(1.7609E-20);  // UO2++
  components->free_ion.push_back(2.5277E-04);  // O2(aq)
  components->free_ion.push_back(1.0000E-15);  // Tracer
}  // end fbasin_free_ions()

void fbasin_minerals(ac::Beaker::BeakerComponents* components) {
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

void fbasin_sorbed(ac::Beaker::BeakerComponents* components) {
  components->total_sorbed.resize(components->total.size(), 0.0);
}  // end fbasin_sorbed


/*******************************************************************************
 **
 ** fbasin UO2 5 component problem
 **
 ******************************************************************************/
void uo2_5_component_initial(const ac::Verbosity& verbosity,
                             std::string* thermo_database_file,
                             std::string* activity_model_name,
                             ac::Beaker::BeakerComponents* components,
                             double* delta_time,
                             int* num_time_steps,
                             int* output_interval) {
  if (verbosity == ac::kTerse) {
    std::cout << "Running uo2_5_component problem, \'initial\' conditions." << std::endl;
  }
  *thermo_database_file = "input/uo2-5-component.bgd";
  *activity_model_name = ac::ActivityModelFactory::debye_huckel;
  *delta_time = 30.0 * 24.0 * 3600.0;
  *num_time_steps = 12;
  *output_interval = 1;

  uo2_5_component_aqueous_initial(components);
  uo2_5_component_free_ions(components);
  uo2_5_component_minerals(components);
  uo2_5_component_sorbed(components);
}  // end uo2_5_component_initial

void uo2_5_component_outlet(const ac::Verbosity& verbosity,
                            std::string* thermo_database_file,
                            std::string* activity_model_name,
                            ac::Beaker::BeakerComponents* components,
                            double* delta_time,
                            int* num_time_steps,
                            int* output_interval) {
  if (verbosity == ac::kTerse) {
    std::cout << "Running uo2_5_component problem, \'outlet\' conditions." << std::endl;
  }
  *thermo_database_file = "input/uo2-5-component.bgd";
  *activity_model_name = ac::ActivityModelFactory::debye_huckel;
  *delta_time = 30.0 * 24.0 * 3600.0;
  *num_time_steps = 12;
  *output_interval = 1;

  uo2_5_component_aqueous_outlet(components);
  uo2_5_component_free_ions(components);
  uo2_5_component_minerals(components);
  uo2_5_component_sorbed(components);
}  // end uo2_5_component_outlet

void uo2_5_component_source(const ac::Verbosity& verbosity,
                            std::string* thermo_database_file,
                            std::string* activity_model_name,
                            ac::Beaker::BeakerComponents* components,
                            double* delta_time,
                            int* num_time_steps,
                            int* output_interval) {
  if (verbosity == ac::kTerse) {
    std::cout << "Running uo2_5_component problem, \'source\' conditions." << std::endl;
  }
  *thermo_database_file = "input/uo2-5-component.bgd";
  *activity_model_name = ac::ActivityModelFactory::debye_huckel;
  *delta_time = 30.0 * 24.0 * 3600.0;
  *num_time_steps = 12;
  *output_interval = 1;

  uo2_5_component_aqueous_source(components);
  uo2_5_component_free_ions(components);
  uo2_5_component_minerals(components);
  uo2_5_component_sorbed(components);
}  // end uo2_5_component_source

void uo2_5_component_aqueous_initial(ac::Beaker::BeakerComponents* components) {
  components->total.push_back(6.5874E-09);  // Al+++
  components->total.push_back(-3.1408E-07);  // H+
  components->total.push_back(1.0000E-06);  // HPO4--
  components->total.push_back(1.8703E-04);  // SiO2(aq)
  components->total.push_back(1.0000E-15);  // UO2++
}  // end uo2_5_component_aqueous_initial

void uo2_5_component_aqueous_outlet(ac::Beaker::BeakerComponents* components) {
  // constraint: outlet
  components->total.push_back(1.0000E-12);  // Al+++
  components->total.push_back(-1.1407E-09);  // H+
  components->total.push_back(1.0000E-06);  // HPO4--
  components->total.push_back(1.0000E-05);  // SiO2(aq)
  components->total.push_back(1.0000E-15);  // UO2++
}  // end uo2_5_component_aqueous_outlet

void uo2_5_component_aqueous_source(ac::Beaker::BeakerComponents* components) {
  // constraint: source
  components->total.push_back(2.8909E-05);  // Al+++
  components->total.push_back(1.2786E-03);  // H+
  components->total.push_back(7.1028E-05);  // HPO4--
  components->total.push_back(2.5280E-04);  // SiO2(aq)
  components->total.push_back(3.5414E-05);  // UO2++
}  // end uo2_5_component_aqueous_source

void uo2_5_component_free_ions(ac::Beaker::BeakerComponents* components) {
  // free ion concentrations (better initial guess)
  components->free_ion.push_back(5.2970E-16);  // Al+++
  components->free_ion.push_back(3.2759E-08);  // H+
  components->free_ion.push_back(9.9886E-07);  // HPO4--
  components->free_ion.push_back(1.8703E-04);  // SiO2(aq)
  components->free_ion.push_back(1.7609E-20);  // UO2++
}  // end uo2_5_component_free_ions()

void uo2_5_component_minerals(ac::Beaker::BeakerComponents* components) {
  components->minerals.push_back(0.15);  // Kaolinite
  components->minerals.push_back(0.21);  // Quartz
  components->minerals.push_back(0.0);  // (UO2)3(PO4)2.4H2O
}  // end uo2_5_component_minerals()

void uo2_5_component_sorbed(ac::Beaker::BeakerComponents* components) {
  components->total_sorbed.resize(components->total.size(), 0.0);
}  // end uo2_5_component_sorbed
