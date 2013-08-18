#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include <xercesc/dom/DOM.hpp>

namespace Amanzi {
namespace AmanziNewInput {

#define AMANZI_INPUT_VERSION_MAJOR 1
#define AMANZI_INPUT_VERSION_MINOR 1
#define AMANZI_INPUT_VERSION_MICRO 0



Teuchos::ParameterList translate (const std::string& xmlfilename);

Teuchos::ParameterList get_model_description(xercesc::DOMDocument* xmlDoc);
Teuchos::ParameterList get_Mesh(xercesc::DOMDocument* xmlDoc);
Teuchos::ParameterList get_execution_controls(xercesc::DOMDocument* xmlDoc);
Teuchos::ParameterList get_phases(xercesc::DOMDocument* xmlDoc);
Teuchos::ParameterList get_regions(xercesc::DOMDocument* xmlDoc);
Teuchos::ParameterList get_materials(xercesc::DOMDocument* xmlDoc);
Teuchos::ParameterList get_initial_conditions(xercesc::DOMDocument* xmlDoc);
Teuchos::ParameterList get_boundary_conditions(xercesc::DOMDocument* xmlDoc);
Teuchos::ParameterList get_output(xercesc::DOMDocument* xmlDoc);

static bool isUnstr_ ;
static int dimension_;

}
}
