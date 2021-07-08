/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
**
** test program for openmp threading of multiple beaker objects
**
*/
#include "beaker-threads.hh"

#include <unistd.h>
#include <omp.h>

#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>

#include "Beaker.hh"
#include "chemistry_output.hh"
#include "chemistry_exception.hh"
#include "SimpleThermoDatabase.hh"

namespace ac = Amanzi::AmanziChemistry;

int main(int argc, char** argv) {
  ac::OutputOptions output_options;
  output_options.use_stdout = true;
  output_options.file_name = "chemistry-unit-test-results.txt";
  ac::chem_out = new ac::ChemistryOutput();
  ac::chem_out->Initialize(output_options);

  int test = 0;
  int error = EXIT_SUCCESS;

  std::cout << "Serial setup:" << std::endl;

  error = CommandLineOptions(argc, argv);

  // hard code some basic info we need for chemistry
  ac::Beaker::BeakerParameters parameters;
  parameters.thermo_database_file = "input/fbasin-17.bgd";
  parameters.activity_model_name = "debye-huckel";
  parameters.porosity = 0.5;  // -
  parameters.saturation = 1.0;  // -
  parameters.volume = 1.0;  // m^3
  parameters.water_density = 997.16;
  parameters.tolerance = 1.0e-12;
  parameters.max_iterations = 250;

  // setup and initialize a beaker object for each thread
  int num_threads = omp_get_max_threads();
  std::cout << "  number of threads: " << num_threads << std::endl;
  std::vector<ac::Beaker*> mixing_cells(num_threads);
  std::vector<ac::Beaker::BeakerState> cell_state(num_threads);

  for (int thread = 0; thread < num_threads; thread++) {
    // apply fbasin 'source' condition to all the state
    fbasin_source(&(cell_state[thread]));

    // create and setup the beakers
    mixing_cells[thread] = new ac::SimpleThermoDatabase();
    mixing_cells[thread]->verbosity(verbosity);
    mixing_cells[thread]->Setup(cell_state[thread], parameters);

    {
      std::cout << "Thread: " << thread << std::endl;
      mixing_cells[thread]->Display();
      mixing_cells[thread]->DisplayComponents(cell_state[thread]);
    }

    // solve for initial free-ion concentrations
    mixing_cells[thread]->Speciate(&(cell_state[thread]), parameters);
    mixing_cells[thread]->CopyBeakerToState(&cell_state[thread]);
    mixing_cells[thread]->DisplayResults();
  }

  // create a bunch of computational work
  double delta_time = 30.0 * 24.0 * 3600.0;  // seconds
  int num_time_steps = 1;
  int num_grid_blocks = 1000;

  std::cout << "Parallel work region:" << std::endl;
  for (int time_step = 0; time_step < num_time_steps; time_step++) {
#pragma omp parallel for schedule(dynamic, 1)
    for (int grid_block = 0; grid_block < num_grid_blocks; grid_block++) {
      int thread = omp_get_thread_num();
      std::cout << "grid block: " << grid_block << "   on thread: " << thread << std::endl;

      try {
        mixing_cells[thread]->ReactionStep(&cell_state[thread], parameters, delta_time);
        mixing_cells[thread]->Speciate(&(cell_state[thread]), parameters);
      } catch (const ac::ChemistryException& geochem_error) {
        std::cout << geochem_error.what() << std::endl;
      } catch (const std::runtime_error& rt_error) {
        std::cout << rt_error.what() << std::endl;
      } catch (const std::logic_error& lg_error) {
        std::cout << lg_error.what() << std::endl;
      }
    }
  }
  std::cout << "Serial cleanup:" << std::endl;
  // cleanup memory
  for (int thread = 0; thread < num_threads; thread++) {
    delete mixing_cells[thread];
  }

  std::cout << "Done!\n";
  return error;
}  // end main()


int CommandLineOptions(int argc, char** argv) {
  int error = EXIT_SUCCESS;
  int option, verbosity;
  extern char* optarg;

  while ((option = getopt(argc, argv, "hv:?")) != -1) {
    switch (option) {
      case 'v': {
        verbosity = std::atoi(optarg);
        break;
      }
      case '?':
      case 'h': {  /* help mode */
        /* print some help stuff and exit without doing anything */
        std::cout << argv[0] << " command line options:" << std::endl;
        std::cout << std::endl;
        std::cout << "    -v integer" << std::endl;
        std::cout << "         verbose output:" << std::endl;
        std::cout << "             0: silent" << std::endl;
        std::cout << "             1: terse" << std::endl;
        std::cout << "             2: verbose" << std::endl;
        std::cout << "             3: debug" << std::endl;
        std::cout << "             4: debug beaker" << std::endl;
        std::cout << "             5......" << std::endl;
        error = -1;
        break;
      }
      default: {
        /* no options */
        break;
      }
    }
  }

  return error;
}


void fbasin_source(ac::Beaker::BeakerState* state) {
  fbasin_aqueous_source(state);
  fbasin_free_ions(state);
  fbasin_minerals(state);
  fbasin_sorbed(state);
}


void fbasin_aqueous_source(ac::Beaker::BeakerState* state) {
  // constraint: source
  state->total.push_back(3.4363E-02);  // Na+
  state->total.push_back(1.2475E-05);  // Ca++
  state->total.push_back(3.0440E-05);  // Fe++
  state->total.push_back(1.7136E-05);  // K+
  state->total.push_back(2.8909E-05);  // Al+++
  state->total.push_back(3.6351E-03);  // H+
  state->total.push_back(1.3305E-03);  // N2(aq)
  state->total.push_back(3.4572E-02);  // NO3-
  state->total.push_back(2.1830E-03);  // HCO3-
  state->total.push_back(3.3848E-05);  // Cl-
  state->total.push_back(6.2463E-04);  // SO4--
  state->total.push_back(7.1028E-05);  // HPO4--
  state->total.push_back(7.8954E-05);  // F-
  state->total.push_back(2.5280E-04);  // SiO2(aq)
  state->total.push_back(3.5414E-05);  // UO2++
  state->total.push_back(2.6038E-04);  // O2(aq)
  state->total.push_back(3.5414E-05);  // Tracer
}  // end fbasin_aqueous_source

void fbasin_free_ions(ac::Beaker::BeakerState* state) {
  // free ion concentrations (better initial guess)
  state->free_ion.push_back(9.9969E-06);  // Na+
  state->free_ion.push_back(9.9746E-06);  // Ca++
  state->free_ion.push_back(2.2405E-18);  // Fe++
  state->free_ion.push_back(1.8874E-04);  // K+
  state->free_ion.push_back(5.2970E-16);  // Al+++
  state->free_ion.push_back(3.2759E-08);  // H+
  state->free_ion.push_back(1.0000E-05);  // N2(aq)
  state->free_ion.push_back(1.0000E-05);  // NO3-
  state->free_ion.push_back(1.9282E-04);  // HCO3-
  state->free_ion.push_back(9.9999E-06);  // Cl-
  state->free_ion.push_back(9.9860E-07);  // SO4--
  state->free_ion.push_back(9.9886E-07);  // HPO4--
  state->free_ion.push_back(1.0000E-06);  // F-
  state->free_ion.push_back(1.8703E-04);  // SiO2(aq)
  state->free_ion.push_back(1.7609E-20);  // UO2++
  state->free_ion.push_back(2.5277E-04);  // O2(aq)
  state->free_ion.push_back(1.0000E-15);  // Tracer
}  // end fbasin_free_ions()

void fbasin_minerals(ac::Beaker::BeakerState* state) {
  state->mineral_volume_fraction.push_back(0.0);  // Gibbsite
  state->mineral_volume_fraction.push_back(0.21);  // Quartz
  state->mineral_volume_fraction.push_back(0.15);  // K-Feldspar
  state->mineral_volume_fraction.push_back(0.0);  // Jurbanite
  state->mineral_volume_fraction.push_back(0.1);  // Ferrihydrite
  state->mineral_volume_fraction.push_back(0.15);  // Kaolinite
  state->mineral_volume_fraction.push_back(0.0);  // Schoepite
  state->mineral_volume_fraction.push_back(0.0);  // (UO2)3(PO4)2.4H2O
  state->mineral_volume_fraction.push_back(0.0);  // Soddyite
  state->mineral_volume_fraction.push_back(0.0);  // Calcite
  state->mineral_volume_fraction.push_back(0.0);  // Chalcedony
}  // end fbasin_minerals()

void fbasin_sorbed(ac::Beaker::BeakerState* state) {
  state->total_sorbed.resize(state->total.size(), 0.0);
}  // end fbasin_sorbed
