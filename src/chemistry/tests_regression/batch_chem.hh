/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_BATCH_CHEM_HH_
#define AMANZI_CHEMISTRY_BATCH_CHEM_HH_

#include <string>
#include <vector>

#include "beaker.hh"
#include "chemistry_verbosity.hh"

struct SimulationParameters {
  std::string description;
  std::vector<std::string> verbosity_names;
  amanzi::chemistry::Verbosity verbosity;
  std::string text_output;
  std::string text_time_units;
  std::string comparison_model;
  std::string database_type;
  std::string database_file;
  std::string activity_model;
  double porosity;  // [-]
  double saturation;  // [-]
  double volume;  // [m^3]
  double delta_time;  // [s]
  int num_time_steps;  // [-]
  int output_interval;  // [steps]
  bool display_free_columns;
  std::vector<double> mineral_ssa;  // specific surface area []
  std::vector<double> site_density;  // sorption site density []
  std::vector<double> isotherm_kd;
  std::vector<double> isotherm_langmuir_b;
  std::vector<double> isotherm_freundlich_n;

  SimulationParameters()
      : description(""),
        verbosity_names(),
        verbosity(amanzi::chemistry::kVerbose), 
        text_output(""),
        text_time_units("s"),
        comparison_model("pflotran"),
        database_type("simple"),
        database_file(""),
        activity_model(""),
        porosity(1.0),
        saturation(1.0),
        volume(1.0),
        delta_time(1.0),
        num_time_steps(0),
        output_interval(1),
        display_free_columns(false),
        mineral_ssa(),
        site_density(),
        isotherm_kd(),
        isotherm_langmuir_b(),
        isotherm_freundlich_n() {
  }
};

static const std::string kSimulationSection("simulation parameters");
static const std::string kDescriptionParam("description");
static const std::string kVerbosityParam("verbosity");
static const std::string kTextOutputParam("text_output");
static const std::string kTextTimeUnitsParam("text_time_units");
static const std::string kComparisonModelParam("comparison_model");
static const std::string kDatabaseTypeParam("database_type");
static const std::string kDatabaseFileParam("database_file");
static const std::string kActivityModelParam("activity_model");
static const std::string kPorosityParam("porosity");
static const std::string kSaturationParam("saturation");
static const std::string kVolumeParam("volume");
static const std::string kDeltaTimeParam("delta_time");
static const std::string kNumTimeStepsParam("num_time_steps");
static const std::string kOutputIntervalParam("output_interval");

static const std::string kTotalSection("total");
static const std::string kMineralSection("mineral");
static const std::string kSorbedSection("total_sorbed");
static const std::string kFreeIonSection("free_ion");
static const std::string kIonExchangeSection("ion_exchange");

static const std::string kSiteDensitySection("site_density");
static const std::string kSpecificSurfaceAreaSection("specific_surface_area");
static const std::string kIsothermSection("sorption_isotherms");

int CommandLineOptions(int argc, char** argv,
                       std::string* verbosity_name,
                       std::string* input_file_name,
                       std::string* template_file_name,
                       bool* debug_batch_driver);

void WriteTemplateFile(const std::string& file_name);

void SetupChemistryOutput(void);

void SetupTextOutput(const SimulationParameters& simulation_params,
                     const std::string& input_file_name,
                     std::fstream* text_output, char* time_units,
                     double* time_units_conversion);

void WriteTextOutputHeader(std::fstream* text_output,
                           const char time_units,
                           const std::vector<std::string>& names);

void WriteTextOutput(std::fstream* text_output,
                     const double time,
                     const amanzi::chemistry::Beaker::BeakerComponents& components);

void ReadInputFile(const std::string& file_name,
                   SimulationParameters* simulation_params,
                   amanzi::chemistry::Beaker::BeakerComponents* components);

void ParseSimulationParameter(const std::string& raw_line,
                              SimulationParameters* params);

void ParseComponentValue(const std::string& raw_line,
                         std::vector<double>* component);
void ParseComponentValue(const std::string& raw_line,
                         double* component);

void ModelSpecificParameters(const std::string model,
                             amanzi::chemistry::Beaker::BeakerParameters* parameters);
void CopySimulationParametersToBeakeParameters(
    const SimulationParameters& simulation_params,
    amanzi::chemistry::Beaker::BeakerParameters* parameters);

void PrintInput(const SimulationParameters& params,
                const amanzi::chemistry::Beaker::BeakerComponents& components);
void PrintSimulationParameters(const SimulationParameters& params);
void PrintComponents(const amanzi::chemistry::Beaker::BeakerComponents& components);


#endif  /* AMANZI_CHEMISTRY_BATCH_CHEM_HH_ */
