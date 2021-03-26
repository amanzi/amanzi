#include <unistd.h>

//#define ABORT_ON_FLOATING_POINT_EXCEPTIONS
#ifdef __APPLE__
  #include <xmmintrin.h>
#endif

#include <cstdlib>
#include <cctype>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

// TPLs
#include "VerboseObject_objs.hh"
#include "VerboseObject.hh"
#include "boost/algorithm/string.hpp"

// Chemistry
#include "simple_thermo_database.hh"
#include "beaker.hh"
#include "activity_model_factory.hh"
#include "chemistry_utilities.hh"
#include "chemistry_exception.hh"
#include "string_tokenizer.hh"

#include "batch_chem.hh"

namespace ac = Amanzi::AmanziChemistry;

const std::string kCrunch("crunch");
const std::string kPflotran("pflotran");

/* TODO: might be worth switching over to reading the component values
   into a map rather than a vector, then order of state in the
   cfg file wouldn't matter, but we need to request an name-id map
   from the beaker.  */

int main(int argc, char** argv) {
#ifdef ABORT_ON_FLOATING_POINT_EXCEPTIONS
#ifdef __APPLE__
  // Make floating point exceptions abort the program. runtime error
  // message isn't helpful, but running in gdb will stop at the
  // correct line. This may code may not be apple specific....
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_DENORM);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_DIV_ZERO);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_OVERFLOW);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_UNDERFLOW);
#endif
#endif

  std::stringstream message;

  bool debug_batch_driver(false);
  std::string verbosity_name("");
  std::string input_file_name("");
  std::string template_file_name("");
  int error = EXIT_SUCCESS;

  // if verbosity was specified on the command line, add the level
  Amanzi::VerboseObject::global_hide_line_prefix = false;  // two default value
  Amanzi::VerboseObject::global_default_level = Teuchos::VERB_MEDIUM;

  Teuchos::ParameterList plist;
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Chemistry PK", plist));

  error = CommandLineOptions(argc, argv,
                             &verbosity_name,
                             &input_file_name,
                             &template_file_name,
                             &debug_batch_driver,
                             vo);

  if (!template_file_name.empty()) {
    WriteTemplateFile(template_file_name, vo);
    exit(EXIT_SUCCESS);
  }

  ac::Beaker::BeakerState state;

  SimulationParameters simulation_params;
  if (!input_file_name.empty()) {
    ReadInputFile(input_file_name, &simulation_params, &state, vo);
  }

  if (state.total.size() == 0) {
    message << "Must have a non-zero number of total component values.\n";
    vo->WriteWarning(Teuchos::VERB_LOW, message);
    abort();
  }

  if (debug_batch_driver) {
    PrintInput(simulation_params, state, vo);
  }

  double time_units_conversion = 1.0;
  char time_units = 's';
  std::fstream text_output;
  if (simulation_params.text_output.size() > 0) {
    SetupTextOutput(simulation_params, input_file_name,
                    &text_output, &time_units, &time_units_conversion);
  }

  ac::Beaker* chem = NULL;

  try {
    if (simulation_params.database_file.size() != 0) {
      chem = new ac::SimpleThermoDatabase(vo);

      ac::Beaker::BeakerParameters parameters = chem->GetDefaultParameters();
      parameters.thermo_database_file = simulation_params.database_file;
      parameters.activity_model_name = simulation_params.activity_model;
      parameters.max_iterations = simulation_params.max_iterations;
      parameters.tolerance = simulation_params.tolerance;
      parameters.porosity = simulation_params.porosity;  // -
      parameters.saturation = simulation_params.saturation;  // -
      parameters.volume = simulation_params.volume;  // m^3
      ModelSpecificParameters(simulation_params.comparison_model, &parameters);
      if (state.free_ion.size() != state.total.size()) {
        state.free_ion.resize(state.total.size(), 1.0e-9);
      }
      chem->Setup(state, parameters);

      if (vo->getVerbLevel() >= Teuchos::VERB_HIGH) {
        chem->Display();
        chem->DisplayComponents(state);
      }

      // solve for free-ion concentrations
      chem->Speciate(&state, parameters);
      chem->CopyBeakerToState(&state);
      if (vo->getVerbLevel() >= Teuchos::VERB_EXTREME) {
        chem->DisplayResults();
      }
      bool using_sorption = false;
      if (state.total_sorbed.size() > 0) {
        using_sorption = true;
      }
      if (simulation_params.num_time_steps != 0) {
        message.str("");
        message << "-- Test Beaker Reaction Stepping -------------------------------------" << std::endl;
        vo->Write(Teuchos::VERB_HIGH, message.str());

        // write out the headers info and initial conditions
        chem->DisplayTotalColumnHeaders(simulation_params.display_free_columns);
        chem->DisplayTotalColumns(0.0, state,
                                  simulation_params.display_free_columns);
        std::vector<std::string> names;
        chem->GetPrimaryNames(&names);
        WriteTextOutputHeader(&text_output, time_units, names, using_sorption);
        WriteTextOutput(&text_output, 0.0, state);

        // parameters.max_iterations = 2;
        for (int time_step = 0; time_step < simulation_params.num_time_steps;
             time_step++) {
          chem->ReactionStep(&state, parameters, simulation_params.delta_time);
          if ((time_step + 1) % simulation_params.output_interval == 0) {
            double time = (time_step + 1) * simulation_params.delta_time;
            chem->DisplayTotalColumns(time, state, 
                                      simulation_params.display_free_columns);
            WriteTextOutput(&text_output, time * time_units_conversion, state);
          }
          if (vo->getVerbLevel() >= Teuchos::VERB_HIGH) {
            message.str("");
            ac::Beaker::SolverStatus status = chem->status();
            message << "Timestep: " << time_step << std::endl;
            message << "    number of rhs evaluations: " << status.num_rhs_evaluations << std::endl;
            message << "    number of jacobian evaluations: " << status.num_jacobian_evaluations << std::endl;
            message << "    number of newton iterations: " << status.num_newton_iterations << std::endl;
            message << "    solution converged: " << status.converged << std::endl;
            vo->Write(Teuchos::VERB_HIGH, message.str());
          }
        }
        vo->Write(Teuchos::VERB_HIGH, "---- Final Speciation\n");
        chem->Speciate(&state, parameters);
        if (vo->getVerbLevel() >= Teuchos::VERB_EXTREME) {
          chem->DisplayResults();
        }
      }
    } else {
      vo->Write(Teuchos::VERB_HIGH, "No database file specified in input file.\n");
    }
  } catch (const ac::ChemistryException& geochem_error) {
    vo->WriteWarning(Teuchos::VERB_LOW, geochem_error.what());
    error = EXIT_FAILURE;
  } catch (const std::runtime_error& rt_error) {
    vo->WriteWarning(Teuchos::VERB_LOW, rt_error.what());
    error = EXIT_FAILURE;
  } catch (const std::logic_error& lg_error) {
    vo->WriteWarning(Teuchos::VERB_LOW, lg_error.what());
    error = EXIT_FAILURE;
  }

  if (!error) {
    vo->Write(Teuchos::VERB_HIGH, "Success!\n");
  } else {
    vo->Write(Teuchos::VERB_HIGH, "Failed!\n");
  }

  text_output.close();
  // cleanup memory
  delete chem;

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


/*******************************************************************************
 **
 **  Commandline
 **
 *******************************************************************************/
int CommandLineOptions(int argc, char** argv,
                       std::string* verbosity_name,
                       std::string* input_file_name,
                       std::string* template_file_name,
                       bool* debug_batch_driver,
                       const Teuchos::RCP<Amanzi::VerboseObject>& vo)
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
        std::cout << "         additional verbosity level:" << std::endl;
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
    vo->Write(Teuchos::VERB_EXTREME, message.str());
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
                   ac::Beaker::BeakerState* state,
                   const Teuchos::RCP<Amanzi::VerboseObject>& vo)
{
  std::stringstream message;
  std::ifstream input_file(file_name.c_str());
  if (!input_file) {
    message.str("");
    message << "batch_chem: \n";
    message << "input file \'" << file_name
              << "\' could not be opened." << std::endl;
    vo->WriteWarning(Teuchos::VERB_LOW, message);
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
    kSectionIsotherms
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

    char sym_first = '\0';
    if (raw_line.length() > 0) sym_first = raw_line[0];
    if (sym_first == '#' || sym_first == '\0') {
      line_type = kCommentLine;
    } else if (sym_first == '[') {
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
      } else if (section_name.compare(kIsothermSection) == 0) {
        current_section = kSectionIsotherms;
      } else {
        message.str("");
        message << "batch_chem::ReadInputFile(): ";
        message << "unknown section found on line " << count << ":";
        message << "\'" << raw_line << "\'"<< std::endl;
        vo->Write(Teuchos::VERB_LOW, message.str());
      }
    } else if (line_type == kParameter) {
      // assume parameter line, but it may be empty (just spaces or missing an = )...
      if (current_section == kSectionSimulation) {
        ParseSimulationParameter(raw_line, simulation_params);
      } else if (current_section == kSectionTotal) {
        ParseComponentValue(raw_line, &(state->total));
      } else if (current_section == kSectionSorbed) {
        ParseComponentValue(raw_line, &(state->total_sorbed));
      } else if (current_section == kSectionFreeIon) {
        ParseComponentValue(raw_line, &(state->free_ion));
      } else if (current_section == kSectionMineral) {
        ParseComponentValue(raw_line, &(state->mineral_volume_fraction));
      } else if (current_section == kSectionSpecificSurfaceArea) {
        ParseComponentValue(raw_line, &(state->mineral_specific_surface_area));
      } else if (current_section == kSectionIonExchange) {
        ParseComponentValue(raw_line, &(state->ion_exchange_sites));
      } else if (current_section == kSectionSiteDensity) {
        ParseComponentValue(raw_line, &(state->surface_site_density));
      } else if (current_section == kSectionIsotherms) {
        // TODO: need to figure out the format of this data...
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
    ac::StringTokenizer values(param.at(1), ",");
    std::string value("");
    if (values.size() == 1) {
      value.assign(values.at(0));
      ac::utilities::RemoveLeadingAndTrailingWhitespace(&value);
    }
    //std::cout << "Parsing ----->  '" << param.at(0) << "'" << std::endl;
    if (param.at(0).find(kDescriptionParam) != std::string::npos) {
      // the description probably has spaces in it, so we want to use
      // the raw parameter value from param.at(1) rather than the
      // version in value, which has been tokenized by spaces!
      params->description.assign(param.at(1));
    } else if (param.at(0).find(kTextOutputParam) != std::string::npos) {
      params->text_output.assign(value) ;
    } else if (param.at(0).find(kTextTimeUnitsParam) != std::string::npos) {
      params->text_time_units.assign(value) ;
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
    } else if (param.at(0).find(kToleranceParam) != std::string::npos) {
      params->tolerance = std::atof(value.c_str());
    } else if (param.at(0).find(kMaxIterationsParam) != std::string::npos) {
      params->max_iterations = std::atoi(value.c_str());
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
    double value = 0.;
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
void WriteTemplateFile(const std::string& file_name,
                       const Teuchos::RCP<Amanzi::VerboseObject>& vo)
{
  std::ofstream template_file(file_name.c_str());
  if (!template_file) {
    std::stringstream message;
    message << "batch_chem: \n";
    message << "template file \'" << file_name
              << "\' could not be opened." << std::endl;
    vo->WriteWarning(Teuchos::VERB_LOW, message);
    abort();
  }
  template_file << "[" << kSimulationSection << "]" << std::endl;
  template_file << kDescriptionParam << " = " << std::endl;
  template_file << "# verbosity can be a comma seperated list." << std::endl;
  template_file << kComparisonModelParam << " = pflotran" << std::endl;
  template_file << kTextOutputParam << " = true" << std::endl;
  template_file << kTextTimeUnitsParam << " = days" << std::endl;
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
  template_file << "[" << kIonExchangeSection << "] # CEC" << std::endl;
  template_file << std::endl;
  template_file << "[" << kIsothermSection << "]" << std::endl;
  template_file << std::endl;
  template_file.close();
}  // end WriteTemplateFile()

void SetupTextOutput(const SimulationParameters& simulation_params,
                     const std::string& input_file_name,
                     std::fstream* text_output, char* time_units,
                     double* time_units_conversion) {
  // are we writting to observations to a text file?
  if (boost::iequals(simulation_params.text_output, "true") ||
      boost::iequals(simulation_params.text_output, "yes") ||
      boost::iequals(simulation_params.text_output, "on")) {
    // generate the output file name:
    size_t position = input_file_name.find_last_of('.');
    std::string text_output_name = input_file_name.substr(0, position) + ".txt";

    text_output->open(text_output_name.c_str(), std::fstream::out);

    // do we want to change the time units for the output?
    if (simulation_params.text_time_units.size() > 0) {
      *time_units = std::tolower(simulation_params.text_time_units.at(0));
      switch (*time_units) {
        case 's':
          break;
        case 'm':
          *time_units_conversion = 60.0;
          break;
        case 'h':
          *time_units_conversion = 60.0 * 60.0;
          break;
        case 'd':
          *time_units_conversion = 60.0 * 60.0 * 24.0;
          break;
        case 'y':
          *time_units_conversion = 60.0 * 60.0 * 24.0 * 365.25;
          break;
        default:
          break;
      }
    }
    *time_units_conversion = 1.0 / (*time_units_conversion);
  }
}


void WriteTextOutputHeader(std::fstream* text_output, const char time_units,
                           const std::vector<std::string>& names,
                           const bool using_sorption) {
  if (text_output->is_open()) {
    *text_output << "# Time(" << time_units << ")";
    for (std::vector<std::string>::const_iterator name = names.begin();
         name != names.end(); ++name) {
      *text_output <<  " , " << *name;
    }
    if (using_sorption) {
      for (std::vector<std::string>::const_iterator name = names.begin();
           name != names.end(); ++name) {
        *text_output <<  " , " << *name << "_sorbed";
      }
    }
    *text_output << std::endl;
  }
}


void WriteTextOutput(std::fstream* text_output, const double time, 
                     const Amanzi::AmanziChemistry::Beaker::BeakerState& state) {
  if (text_output->is_open()) {
    std::string seperator(" , ");
    *text_output << std::scientific << std::setprecision(6) << std::setw(15) << time;
    for (int i = 0; i < state.total.size(); ++i) {
      *text_output << seperator << state.total.at(i);
    }
    for (int i = 0; i < state.total_sorbed.size(); ++i) {
      *text_output << seperator << state.total_sorbed.at(i);
    }
    *text_output << std::endl;
  }
}


void PrintInput(const SimulationParameters& params,
                const Amanzi::AmanziChemistry::Beaker::BeakerState& state,
                const Teuchos::RCP<Amanzi::VerboseObject>& vo)
{
  vo->Write(Teuchos::VERB_HIGH, "- Input File ---------------------------------------------------------\n");
  PrintSimulationParameters(params, vo);
  ac::Display(state, "-- Input state: \n", vo);
  vo->Write(Teuchos::VERB_HIGH, "--------------------------------------------------------- Input File -\n");
}  // end PrintInput()


void PrintSimulationParameters(const SimulationParameters& params,
                               const Teuchos::RCP<Amanzi::VerboseObject>& vo)
{
  std::stringstream message;
  message << "-- Simulation parameters:" << std::endl;
  message << "\tdescription: " << params.description << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
  message.str("");
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
  message << "\tmax iterations: " << params.max_iterations << std::endl;
  message << "\ttolerance: " << params.tolerance << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}



