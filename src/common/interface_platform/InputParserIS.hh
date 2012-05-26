#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace Amanzi {
namespace AmanziInput {

  Teuchos::ParameterList translate (Teuchos::ParameterList* plist, int numproc);

Teuchos::ParameterList get_Time_Macro (const std::string& macro_name, Teuchos::ParameterList* plist );
Teuchos::Array<int> get_Cycle_Macro ( const std::string& macro_name, Teuchos::ParameterList* plist );
Teuchos::Array<std::string> get_Variable_Macro ( const std::string& macro_name, Teuchos::ParameterList* plist );
void init_global_info( Teuchos::ParameterList* plist );

Teuchos::ParameterList create_Checkpoint_Data_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_Visualization_Data_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_Observation_Data_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList get_Regions_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList get_Mesh_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList translate_Mesh_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList get_Domain_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_MPC_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_Transport_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_Flow_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_WRM_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_DPC_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_SS_FlowBC_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_State_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_Verbosity_List ( const std::string& vlevel );
Teuchos::ParameterList CreateChemistryList ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_Preconditioners_List ( Teuchos::ParameterList* plist );
Teuchos::ParameterList create_Solvers_List ( Teuchos::ParameterList* plist );

static std::string phase_name;
static std::string phase_comp_name;
static Teuchos::Array<std::string> comp_names;
static std::map<std::string, int> comp_names_map;
static Teuchos::Array<std::string> mineral_names_;
static Teuchos::Array<std::string> sorption_site_names_;
static std::string verbosity_level("low");
static int numproc_;
  
}
}
