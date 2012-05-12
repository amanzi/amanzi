/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "batch_chem.hh"

#ifdef WINDOWS
#include "xgetopt.hh"
#else
#include <unistd.h>
#endif

#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

#include "simple_thermo_database.hh"
#include "beaker.hh"
#include "activity_model_factory.hh"
#include "chemistry_verbosity.hh"
#include "chemistry_output.hh"
#include "chemistry_containers.hh"
#include "chemistry_utilities.hh"
#include "chemistry_exception.hh"
#include "string_tokenizer.hh"

// create a global ChemistryOutput* pointer in the amanzi::chemisry
// namespace that can be used by an other chemistry object
namespace amanzi {
namespace chemistry {
ChemistryOutput* chem_out;
}  // end namespace chemistry
}  // end namespace amanzi

namespace ac = amanzi::chemistry;

const std::string kCrunch("crunch");
const std::string kPflotran("pflotran");

/* TODO: might be worth switching over to reading the component values
   into a map rather than a vector, then order of components in the
   cfg file wouldn't matter, but we need to request an name-id map
   from the beaker.  */

int main(int argc, char** argv) {
  ac::OutputOptions output_options;
  output_options.use_stdout = true;
  output_options.file_name = "chemistry-unit-test-results.txt";
  output_options.verbosity_levels.push_back(ac::strings::kVerbosityError);
  output_options.verbosity_levels.push_back(ac::strings::kVerbosityWarning);
  output_options.verbosity_levels.push_back(ac::strings::kVerbosityVerbose);

  ac::chem_out = new ac::ChemistryOutput();
  ac::chem_out->Initialize(output_options);
  std::stringstream message;

  bool debug_batch_driver(false);
  std::string verbosity_name("");
  std::string input_file_name("");
  std::string template_file_name("");
  int error = EXIT_SUCCESS;

  error = CommandLineOptions(argc, argv,
                             &verbosity_name,
                             &input_file_name,
                             &template_file_name,
                             &debug_batch_driver);

  if (!template_file_name.empty()) {
    WriteTemplateFile(template_file_name);
    exit(EXIT_SUCCESS);
  }

  ac::Beaker::BeakerComponents components;
  components.free_ion.clear();
  components.minerals.clear();
  components.ion_exchange_sites.clear();
  components.total.clear();
  components.total_sorbed.clear();

  SimulationParameters simulation_params;
  simulation_params.mineral_ssa.clear();
  simulation_params.site_density.clear();
  simulation_params.cation_exchange_capacity = -1.0;
  if (!input_file_name.empty()) {
    ReadInputFile(input_file_name, &simulation_params, &components);
  }

  if (components.total.size() == 0) {
    message.str("");
    message << "Must have a non-zero number of total component values." << std::endl;
    ac::chem_out->Write(ac::kError, message);
    abort();
  }

  ac::VerbosityMap verbosity_map = ac::CreateVerbosityMap();
  simulation_params.verbosity = verbosity_map[simulation_params.verbosity_name];
  // if verbosity was specified on the command line, override the
  // value specified in the input file
  if (!verbosity_name.empty()) {
    simulation_params.verbosity = verbosity_map[verbosity_name];
  }

  if (debug_batch_driver) {
    PrintInput(simulation_params, components);
  }
  ac::Beaker* chem = NULL;

  try {
    if (simulation_params.database_file.size() != 0) {
      chem = new ac::SimpleThermoDatabase();
      chem->verbosity(simulation_params.verbosity);

      ac::Beaker::BeakerParameters parameters = chem->GetDefaultParameters();
      parameters.thermo_database_file = simulation_params.database_file;
      parameters.activity_model_name = simulation_params.activity_model;
      parameters.porosity = simulation_params.porosity;  // -
      parameters.saturation = simulation_params.saturation;  // -
      parameters.volume = simulation_params.volume;  // m^3
      ModelSpecificParameters(simulation_params.comparison_model, &parameters);
      OverrideParameters(simulation_params, &parameters);

      if (components.free_ion.size() == 0) {
        // initialize free-ion concentrations, these are actual
        // concentrations, so the value must be positive or ln(free_ion)
        // will return a nan!
        components.free_ion.resize(components.total.size(), 1.0e-9);
      }

      chem->Setup(components, parameters);

      if (simulation_params.verbosity >= ac::kVerbose) {
        chem->Display();
        chem->DisplayComponents(components);
      }

      // solve for free-ion concentrations
      chem->Speciate(components, parameters);
      chem->UpdateComponents(&components);
      if (simulation_params.verbosity >= ac::kTerse) {
        chem->DisplayResults();
      }

      if (simulation_params.num_time_steps != 0) {
        message.str("");
        message << "-- Test Beaker Reaction Stepping -------------------------------------" << std::endl;
        ac::chem_out->Write(ac::kVerbose, message);
        chem->DisplayTotalColumnHeaders();
        chem->DisplayTotalColumns(0.0, components);
        // parameters.max_iterations = 2;
        for (int time_step = 0; time_step < simulation_params.num_time_steps;
             time_step++) {
          chem->ReactionStep(&components, parameters, simulation_params.delta_time);
          if ((time_step + 1) % simulation_params.output_interval == 0) {
            chem->DisplayTotalColumns((time_step + 1) * simulation_params.delta_time,
                                      components);
          }
          if (simulation_params.verbosity >= ac::kDebugNonlinearSolver) {
            message.str("");
            ac::Beaker::SolverStatus status = chem->status();
            message << "Timestep: " << time_step << std::endl;
            message << "    number of rhs evaluations: " << status.num_rhs_evaluations << std::endl;
            message << "    number of jacobian evaluations: " << status.num_jacobian_evaluations << std::endl;
            message << "    number of newton iterations: " << status.num_newton_iterations << std::endl;
            message << "    solution converged: " << status.converged << std::endl;
            ac::chem_out->Write(ac::kVerbose, message);
          }
        }
        ac::chem_out->Write(ac::kVerbose, "---- Final Speciation\n");
        chem->Speciate(components, parameters);
        if (simulation_params.verbosity >= ac::kTerse) {
          chem->DisplayResults();
        }
      }
    } else {
      ac::chem_out->Write(ac::kVerbose, "No database file specified in input file.\n");
    }
  } catch (const ac::ChemistryException& geochem_error) {
    ac::chem_out->Write(ac::kError, geochem_error.what());
    error = EXIT_FAILURE;
  } catch (const std::runtime_error& rt_error) {
    ac::chem_out->Write(ac::kError, rt_error.what());
    error = EXIT_FAILURE;
  } catch (const std::logic_error& lg_error) {
    ac::chem_out->Write(ac::kError, lg_error.what());
    error = EXIT_FAILURE;
  }

  if (!error) {
    ac::chem_out->Write(ac::kVerbose, "Success!\n");
  } else {
    ac::chem_out->Write(ac::kVerbose, "Failed!\n");
  }

  // cleanup memory
  delete chem;
  delete ac::chem_out;

  return error;
}  // end main()


void ModelSpecificParameters(const std::string model,
                             ac::Beaker::BeakerParameters* parameters) {
  if (model == kCrunch) {
    parameters->water_density = 997.075;  // kg / m^3
  } else if (model == kPflotran) {
    parameters->water_density = 997.16;  // kg / m^3
    // where did this number come from?
    // default parameters->water_density = 997.205133945901;  // kg / m^3
  } else {
    // bad model name, how did we get here....
  }
}  // end ModelSpecificParameters()

void OverrideParameters(const SimulationParameters& simulation_params,
                        ac::Beaker::BeakerParameters* parameters) {
  if (simulation_params.mineral_ssa.size() > 0) {
    parameters->override_database = true;
    parameters->mineral_specific_surface_area.assign(
        simulation_params.mineral_ssa.begin(),
        simulation_params.mineral_ssa.end());
  }

  if (simulation_params.site_density.size() > 0) {
    parameters->override_database = true;
    parameters->sorption_site_density.assign(
        simulation_params.site_density.begin(),
        simulation_params.site_density.end());
  }

  if (simulation_params.cation_exchange_capacity > 0.0) {
    parameters->override_database = true;
    parameters->cation_exchange_capacity = 
        simulation_params.cation_exchange_capacity;
  }
}  // end OverrideParameters()

/*******************************************************************************
 **
 **  Commandline
 **
 *******************************************************************************/
int CommandLineOptions(int argc, char** argv,
                       std::string* verbosity_name,
                       std::string* input_file_name,
                       std::string* template_file_name,
                       bool* debug_batch_driver)
{
  int error = -2;
  int option;
  extern char* optarg;

  while ((option = getopt(argc, argv, "di:ht:v:?")) != -1) {
    switch (option) {
      case 'd': {
        *debug_batch_driver = true;
        break;
      }
      case 'i': {
        /* input file name */
        input_file_name->assign(optarg);
        error = EXIT_SUCCESS;
        break;
      }
      case 't': {
        /* template file name */
        template_file_name->assign(optarg);
        error = EXIT_SUCCESS;
        break;
      }
      case 'v': {
        verbosity_name->assign(optarg);
        break;
      }
      case '?':
      case 'h': {  /* help mode */
        /* print some help stuff and exit without doing anything */
        std::cout << argv[0] << " command line options:" << std::endl;
        std::cout << "    -d" << std::endl;
        std::cout << "         debugging flag for batch driver" << std::endl;
        std::cout << "    -i string " << std::endl;
        std::cout << "         input file name" << std::endl;
        std::cout << std::endl;
        std::cout << "    -t string" << std::endl;
        std::cout << "         write a template input file" << std::endl;
        std::cout << std::endl;
        std::cout << "    -v string" << std::endl;
        std::cout << "         override verbosity from input file:" << std::endl;
        std::cout << "            silent" << std::endl;
        std::cout << "            terse" << std::endl;
        std::cout << "            verbose" << std::endl;
        std::cout << "            debug" << std::endl;
        std::cout << "            debug_beaker" << std::endl;
        std::cout << "            debug_database" << std::endl;
        std::cout << "            debug_mineral_kinetics" << std::endl;
        std::cout << "            debug_ion_exchange" << std::endl;
        std::cout << "            debug_newton_solver" << std::endl;
        error = -1;
        break;
      }
      default: {
        /* no options */
        break;
      }
    }
  }

  if (!input_file_name->c_str() && !template_file_name->c_str()) {
    std::cout << "An input or template file name must be specified." << std::endl;
    std::cout << "Run \"" <<  argv[0] << " -h \" for help." << std::endl;
  }

  if (*debug_batch_driver) {
    std::stringstream message;
    message << "- Command Line Options -----------------------------------------------" << std::endl;
    message << "\tdebug batch driver: " << *debug_batch_driver << std::endl;
    message << "\tinput file name: " << *input_file_name << std::endl;
    message << "\ttemplate file name: " << *template_file_name << std::endl;
    message << "\tverbosity name: " << *verbosity_name << std::endl;
    message << "----------------------------------------------- Command Line Options -" << std::endl;
    message << std::endl << std::endl;
    ac::chem_out->Write(ac::kDebugDriver, message);
  }
  return error;
}  // end commandLineOptions()


/*******************************************************************************
 **
 **  Input file parser
 **
 *******************************************************************************/
void ReadInputFile(const std::string& file_name,
                   SimulationParameters* simulation_params,
                   ac::Beaker::BeakerComponents* components)
{
  std::stringstream message;
  std::ifstream input_file(file_name.c_str());
  if (!input_file) {
    message.str("");
    message << "batch_chem: \n";
    message << "input file \'" << file_name
              << "\' could not be opened." << std::endl;
    ac::chem_out->Write(ac::kError, message);
    abort();
  }

  enum LineType {
    kCommentLine,
    kSection,
    kParameter
  } line_type;

  enum SectionType {
    kSectionSimulation,
    kSectionTotal,
    kSectionMineral,
    kSectionSorbed,
    kSectionFreeIon,
    kSectionIonExchange,
    kSectionSiteDensity,
    kSectionSpecificSurfaceArea,
    kSectionCationExchangeCapacity
  } current_section;

  int count = 0;
  const int max_lines = 500;
  while (!input_file.eof() && count < max_lines) {
    count++;
    std::string raw_line;
    getline(input_file, raw_line);
    //std::cout << raw_line << std::endl;
    if ((raw_line.size() > 0) && (raw_line.at(raw_line.size() - 1) == '\r')) {
      // getline only searches for \n line ends. windows files use \r\n
      // check for a hanging \r and remove it if it is there
      raw_line.resize(raw_line.size() - 1);
    }

    char first = '\0';
    if (raw_line.length() > 0) first = raw_line[0];
    if (first == '#' || first == '\0') {
      line_type = kCommentLine;
    } else if (first == '[') {
      line_type = kSection;
    } else {
      line_type = kParameter;
    }

    if (line_type == kSection) {
      size_t first = raw_line.find_first_not_of('[');
      size_t last = raw_line.find_last_of(']');
      last--;
      std::string section_name = raw_line.substr(first, last);
      if (section_name.compare(kSimulationSection) == 0) {
        current_section = kSectionSimulation;
      } else if (section_name.compare(kTotalSection) == 0) {
        current_section = kSectionTotal;
      } else if (section_name.compare(kMineralSection) == 0) {
        current_section = kSectionMineral;
      } else if (section_name.compare(kIonExchangeSection) == 0) {
        current_section = kSectionIonExchange;
      } else if (section_name.compare(kSorbedSection) == 0) {
        current_section = kSectionSorbed;
      } else if (section_name.compare(kFreeIonSection) == 0) {
        current_section = kSectionFreeIon;
      } else if (section_name.compare(kSiteDensitySection) == 0) {
        current_section = kSectionSiteDensity;
      } else if (section_name.compare(kSpecificSurfaceAreaSection) == 0) {
        current_section = kSectionSpecificSurfaceArea;
      } else if (section_name.compare(kCationExchangeCapacitySection) == 0) {
        current_section = kSectionCationExchangeCapacity;
      } else {
        message.str("");
        message << "batch_chem::ReadInputFile(): ";
        message << "unknown section found on line " << count << ":";
        message << "\'" << raw_line << "\'"<< std::endl;
        ac::chem_out->Write(ac::kDebugInputFile, message);
      }
    } else if (line_type == kParameter) {
      // assume parameter line, but it may be empty (just spaces or missing an = )...
      if (current_section == kSectionSimulation) {
        ParseSimulationParameter(raw_line, simulation_params);
      } else if (current_section == kSectionTotal) {
        ParseComponentValue(raw_line, &(components->total));
      } else if (current_section == kSectionMineral) {
        ParseComponentValue(raw_line, &(components->minerals));
      } else if (current_section == kSectionIonExchange) {
        ParseComponentValue(raw_line, &(components->ion_exchange_sites));
      } else if (current_section == kSectionSorbed) {
        ParseComponentValue(raw_line, &(components->total_sorbed));
      } else if (current_section == kSectionFreeIon) {
        ParseComponentValue(raw_line, &(components->free_ion));
      } else if (current_section == kSectionSiteDensity) {
        ParseComponentValue(raw_line, &(simulation_params->site_density));
      } else if (current_section == kSectionSpecificSurfaceArea) {
        ParseComponentValue(raw_line, &(simulation_params->mineral_ssa));
      } else if (current_section == kSectionCationExchangeCapacity) {
        ParseComponentValue(raw_line, &(simulation_params->cation_exchange_capacity));
      }
    }
  }

  input_file.close();
}  // end ReadInputFile()


void ParseSimulationParameter(const std::string& raw_line,
                              SimulationParameters* params)
{
  std::string equal("=:");
  std::string spaces(" \t");
  ac::StringTokenizer param(raw_line, equal);
  //std::cout << "\'" << raw_line << "\'" << std::endl;
  // if param.size() == 0 then we have a blank line
  if (param.size() != 0) {
    ac::StringTokenizer param_value(param.at(1), spaces);
    std::string value("");
    if (param_value.size() > 0) {
      value.assign(param_value.at(0));
    }
    if (param.at(0).find(kDescriptionParam) != std::string::npos) {
      // the description probably has spaces in it, so we want to use
      // the raw parameter value from param.at(1) rather than the
      // version in value, which has been tokenized by spaces!
      params->description.assign(param.at(1));
    } else if (param.at(0).find(kVerbosityParam) != std::string::npos) {
      params->verbosity_name.assign(value);
    } else if (param.at(0).find(kComparisonModelParam) != std::string::npos) {
      params->comparison_model.assign(value);
    } else if (param.at(0).find(kDatabaseTypeParam) != std::string::npos) {
      params->database_type.assign(value);
    } else if (param.at(0).find(kDatabaseFileParam) != std::string::npos) {
      params->database_file.assign(value);
    } else if (param.at(0).find(kActivityModelParam) != std::string::npos) {
      params->activity_model.assign(value);
    } else if (param.at(0).find(kPorosityParam) != std::string::npos) {
      params->porosity = std::atof(value.c_str());
    } else if (param.at(0).find(kSaturationParam) != std::string::npos) {
      params->saturation = std::atof(value.c_str());
    } else if (param.at(0).find(kVolumeParam) != std::string::npos) {
      params->volume = std::atof(value.c_str());
    } else if (param.at(0).find(kDeltaTimeParam) != std::string::npos) {
      params->delta_time = std::atof(value.c_str());
    } else if (param.at(0).find(kNumTimeStepsParam) != std::string::npos) {
      params->num_time_steps = std::atoi(value.c_str());
    } else if (param.at(0).find(kOutputIntervalParam) != std::string::npos) {
      params->output_interval = std::atoi(value.c_str());
    }
  }

}  // end ParseSimulationParameter()


void ParseComponentValue(const std::string& raw_line,
                         std::vector<double>* component)
{
  // for now we assume that the order of the component is the
  // same as the order in the database file
  std::string equal("=:");
  std::string spaces(" \t");
  ac::StringTokenizer param(raw_line, equal);
  //std::cout << "\'" << raw_line << "\'" << std::endl;
  // if param.size() == 0 then we have a blank line
  if (param.size() != 0) {
    ac::StringTokenizer param_value(param.at(1), spaces);
    double value;
    if (param_value.size() > 0) {
      value = std::atof(param_value.at(0).c_str());
    }
    component->push_back(value);
  }

}  // end ParseComponentValue();

void ParseComponentValue(const std::string& raw_line,
                         double* component)
{
  // this is intended for a single value, not a c-style array!
  std::string equal("=:");
  std::string spaces(" \t");
  ac::StringTokenizer param(raw_line, equal);
  //std::cout << "\'" << raw_line << "\'" << std::endl;
  // if param.size() == 0 then we have a blank line
  if (param.size() != 0) {
    ac::StringTokenizer param_value(param.at(1), spaces);
    double value;
    if (param_value.size() > 0) {
      value = std::atof(param_value.at(0).c_str());
    }
    *component = value;
  }

}  // end ParseComponentValue();


/*******************************************************************************
 **
 **  Output related functions
 **
 *******************************************************************************/
void WriteTemplateFile(const std::string& file_name)
{
  std::ofstream template_file(file_name.c_str());
  if (!template_file) {
    std::stringstream message;
    message << "batch_chem: \n";
    message << "template file \'" << file_name
              << "\' could not be opened." << std::endl;
    ac::chem_out->Write(ac::kError, message);
    abort();
  }
  template_file << "[" << kSimulationSection << "]" << std::endl;
  template_file << kDescriptionParam << " = " << std::endl;
  template_file << kVerbosityParam << " = verbose" << std::endl;
  template_file << kComparisonModelParam << " = pflotran" << std::endl;
  template_file << std::endl;
  template_file << kDatabaseTypeParam << " = simple" << std::endl;
  template_file << kDatabaseFileParam << " = " << std::endl;
  template_file << kActivityModelParam << " = debye-huckel" << std::endl;
  template_file << kPorosityParam << " = " << std::endl;
  template_file << kSaturationParam << " = " << std::endl;
  template_file << kVolumeParam << " = " << std::endl;
  template_file << kDeltaTimeParam << " = " << std::endl;
  template_file << kNumTimeStepsParam << " = " << std::endl;
  template_file << kOutputIntervalParam << " = " << std::endl;
  template_file << std::endl;
  template_file << "# all component values must be in the same order as the database file" << std::endl;
  template_file << "[" << kTotalSection << "]" << std::endl;
  template_file << std::endl;
  template_file << "[" << kMineralSection << "]" << std::endl;
  template_file << std::endl;
  template_file << "[" << kSorbedSection << "]" << std::endl;
  template_file << std::endl;
  template_file << "[" << kFreeIonSection << "]" << std::endl;
  template_file << std::endl;
  template_file << "[" << kIonExchangeSection << "]" << std::endl;
  template_file << std::endl;
  template_file.close();
}  // end WriteTemplateFile()


void PrintInput(const SimulationParameters& params,
                const amanzi::chemistry::Beaker::BeakerComponents& components)
{
  ac::chem_out->Write(ac::kVerbose, "- Input File ---------------------------------------------------------\n");
  PrintSimulationParameters(params);
  PrintComponents(components);
  ac::chem_out->Write(ac::kVerbose, "--------------------------------------------------------- Input File -\n");
}  // end PrintInput()


void PrintSimulationParameters(const SimulationParameters& params)
{
  std::stringstream message;
  message << "-- Simulation parameters:" << std::endl;
  message << "\tdescription: " << params.description << std::endl;
  message << "\tverbosity name: " << params.verbosity_name << std::endl;
  message << "\tverbosity enum: " << params.verbosity << std::endl;
  message << "\tcomparison model: " << params.comparison_model << std::endl;
  message << "\tdatabase type: " << params.database_type << std::endl;
  message << "\tdatabase file: " << params.database_file << std::endl;
  message << "\tactivity model: " << params.activity_model << std::endl;
  message << "\tporosity: " << params.porosity << std::endl;
  message << "\tsaturation: " << params.saturation << std::endl;
  message << "\tvolume: " << params.volume << std::endl;
  message << "\tdelta time: " << params.delta_time << std::endl;
  message << "\tnum time steps: " << params.num_time_steps << std::endl;
  message << "\toutput interval: " << params.output_interval << std::endl;
  message << "-- Database override parameters:" << std::endl;
  ac::chem_out->Write(ac::kVerbose, message);
  ac::utilities::PrintVector("  Site Density", params.site_density);
  ac::utilities::PrintVector("  Specific Surface Area", params.mineral_ssa);
}


void PrintComponents(const ac::Beaker::BeakerComponents& components)
{
  ac::chem_out->Write(ac::kVerbose, "-- Input components: \n");
  ac::utilities::PrintVector("  Totals", components.total);
  ac::utilities::PrintVector("  Minerals", components.minerals);
  ac::utilities::PrintVector("  Total sorbed", components.total_sorbed);
  ac::utilities::PrintVector("  Free Ion", components.free_ion);
  ac::utilities::PrintVector("  Ion Exchange", components.ion_exchange_sites);

}  // end PrintComponents()

