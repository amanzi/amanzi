#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include <xercesc/dom/DOM.hpp>

namespace Amanzi {
namespace AmanziNewInput {

#define AMANZI_INPUT_VERSION_MAJOR 2
#define AMANZI_INPUT_VERSION_MINOR 0
#define AMANZI_INPUT_VERSION_MICRO 0



Teuchos::ParameterList translate (const std::string& xmlfilename, const std::string& xmlSchemafile);

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
double get_time_value(std::string time_value, Teuchos::ParameterList def_list);
double get_double_constant(std::string pos_name, Teuchos::ParameterList def_list);
int get_int_constant(std::string pos_name, Teuchos::ParameterList def_list);
Teuchos::Array<std::string> make_regions_list(char* char_array);
bool compare_region_names(Teuchos::Array<std::string> regions, Teuchos::ParameterList def_list);
Teuchos::Array<double> make_coordinates(char* char_array, Teuchos::ParameterList def_list);
Teuchos::ParameterList make_chemistry(Teuchos::ParameterList def_list);
void write_BDG_file(Teuchos::ParameterList solute_list, Teuchos::ParameterList def_list);

static bool isUnstr_ ;
static int dimension_;

}
}
