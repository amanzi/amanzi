#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include <xercesc/dom/DOM.hpp>
#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziNewInput {

#define AMANZI_INPUT_VERSION_MAJOR 2
#define AMANZI_INPUT_VERSION_MINOR 1
#define AMANZI_INPUT_VERSION_MICRO 1

// validate sets
const Teuchos::Array<std::string> verbosityStrings = Teuchos::tuple<std::string>("None", "Low", "Medium", "High", "Extreme");
const Teuchos::Array<std::string> meshfileStrings = Teuchos::tuple<std::string>("exodus ii", "exodus II", "Exodus II", "Exodus ii", "H5M", "h5m");
  
//Teuchos::ParameterList translate (const std::string& xmlfilename, const std::string& xmlSchemafile);
Teuchos::ParameterList translate (const std::string& xmlfilename);

Teuchos::ParameterList get_verbosity(xercesc::DOMDocument* xmlDoc);
Teuchos::ParameterList get_constants(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
std::string get_amanzi_version(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
void get_sim_type(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList* def_list);
Teuchos::ParameterList get_model_description(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
Teuchos::ParameterList get_Mesh(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
Teuchos::ParameterList get_execution_controls(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList* def_list);
Teuchos::ParameterList get_phases(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
Teuchos::ParameterList get_regions(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList* def_list);
Teuchos::ParameterList get_materials(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
Teuchos::ParameterList get_initial_conditions(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
Teuchos::ParameterList get_boundary_conditions(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
Teuchos::ParameterList get_sources(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
Teuchos::ParameterList get_output(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);
Teuchos::ParameterList get_echo_translated(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list);

Teuchos::ParameterList get_file_info(Teuchos::ParameterList propertyList, xercesc::DOMElement* propElement, std::string propName, std::string sectionName);
void get_gslib_info(Teuchos::ParameterList& propertyList,const Teuchos::ParameterList& defList, const xercesc::DOMElement* propElement, const std::string& propName, const std::string& sectionName, const std::string& dfile_DEF);

double get_time_value(std::string time_value, Teuchos::ParameterList def_list);
double convert_time_value(char* time_value);
double get_double_constant(std::string pos_name, Teuchos::ParameterList def_list);
int get_int_constant(std::string pos_name, Teuchos::ParameterList def_list);
Teuchos::Array<int> make_int_list(char* char_array);
Teuchos::Array<std::string> make_regions_list(char* char_array);
bool compare_region_names(Teuchos::Array<std::string> regions, Teuchos::ParameterList def_list);
std::string trim_string(char* tmp);
Teuchos::Array<double> make_coordinates(char* char_array, Teuchos::ParameterList def_list);
Teuchos::ParameterList make_chemistry(Teuchos::ParameterList def_list);

void write_BDG_file(Teuchos::ParameterList solute_list, Teuchos::ParameterList def_list);

void throw_error_str_ustr(std::string section, std::string element_type, std::string sim_type);
void throw_error_illformed(std::string section, std::string element_type, std::string ill_formed);
void throw_error_illformed(std::string section, std::string element_type, std::string ill_formed, std::string options);
void throw_error_missattr(std::string section, std::string att_elem_type, std::string missing, std::string elem_name);
void throw_error_missattr(std::string section, std::string att_elem_type, std::string missing, std::string elem_name, std::string options);
void throw_warning_skip(std::string element);

static bool isUnstr_ ;
static int dimension_;
static Teuchos::Array<std::string> regionNames_string_;
static Teuchos::RCP<VerboseObject> voI_;



}
}
