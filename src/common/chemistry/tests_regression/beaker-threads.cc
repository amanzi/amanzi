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

#include "simple_thermo_database.hh"
#include "beaker.hh"
#include "chemistry_output.hh"
#include "chemistry_exception.hh"

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
  std::vector<ac::Beaker::BeakerComponents> cell_components(num_threads);

  for (int thread = 0; thread < num_threads; thread++) {
    // apply fbasin 'source' condition to all the components
    fbasin_source(&(cell_components[thread]));

    // create and setup the beakers
    mixing_cells[thread] = new ac::SimpleThermoDatabase();
    mixing_cells[thread]->verbosity(verbosity);
    mixing_cells[thread]->Setup(cell_components[thread], parameters);

    {
      std::cout << "Thread: " << thread << std::endl;
      mixing_cells[thread]->Display();
      mixing_cells[thread]->DisplayComponents(cell_components[thread]);
    }

    // solve for initial free-ion concentrations
    mixing_cells[thread]->Speciate(&(cell_components[thread]), parameters);
    mixing_cells[thread]->CopyBeakerToComponents(&cell_components[thread]);
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
        mixing_cells[thread]->ReactionStep(&cell_components[thread], parameters, delta_time);
        mixing_cells[thread]->Speciate(&(cell_components[thread]), parameters);
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


void fbasin_source(ac::Beaker::BeakerComponents* components) {
  fbasin_aqueous_source(components);
  fbasin_free_ions(components);
  fbasin_minerals(components);
  fbasin_sorbed(components);
}


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
  components->mineral_volume_fraction.push_back(0.0);  // Gibbsite
  components->mineral_volume_fraction.push_back(0.21);  // Quartz
  components->mineral_volume_fraction.push_back(0.15);  // K-Feldspar
  components->mineral_volume_fraction.push_back(0.0);  // Jurbanite
  components->mineral_volume_fraction.push_back(0.1);  // Ferrihydrite
  components->mineral_volume_fraction.push_back(0.15);  // Kaolinite
  components->mineral_volume_fraction.push_back(0.0);  // Schoepite
  components->mineral_volume_fraction.push_back(0.0);  // (UO2)3(PO4)2.4H2O
  components->mineral_volume_fraction.push_back(0.0);  // Soddyite
  components->mineral_volume_fraction.push_back(0.0);  // Calcite
  components->mineral_volume_fraction.push_back(0.0);  // Chalcedony
}  // end fbasin_minerals()

void fbasin_sorbed(ac::Beaker::BeakerComponents* components) {
  components->total_sorbed.resize(components->total.size(), 0.0);
}  // end fbasin_sorbed
