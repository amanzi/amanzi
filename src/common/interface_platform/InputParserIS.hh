#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace Amanzi {
namespace AmanziInput {

#define AMANZI_OLD_INPUT_VERSION_MAJOR 1
#define AMANZI_OLD_INPUT_VERSION_MINOR 2
#define AMANZI_OLD_INPUT_VERSION_MICRO 2

Teuchos::ParameterList Translate(Teuchos::ParameterList* plist, int numproc);

Teuchos::ParameterList CreateCheckpointDataList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateWalkaboutDataList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateVisualizationDataList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateObservationDataList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateMPC_List(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateTransportList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateTransportSrcList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateFlowList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateWRM_List(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateFlowSrcList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateSS_FlowBC_List(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateStateList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateVerbosityList(const std::string& vlevel);
Teuchos::ParameterList CreateChemistryList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreatePreconditionersList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateDPC_List(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateBILU_List(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateHypreAMG_List(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateSolversList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateNonlinearSolversList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateTimePeriodControlList(Teuchos::ParameterList* plist);
Teuchos::ParameterList CreateFlowOperatorList(const std::string& disc_method);
Teuchos::ParameterList CreatePartitionList(Teuchos::ParameterList* plist);

Teuchos::ParameterList TranslateMeshList(Teuchos::ParameterList* plist);
Teuchos::Array<std::string> TranslateForms (Teuchos::Array<std::string>& forms);

Teuchos::ParameterList get_Time_Macro(const std::string& macro_name, Teuchos::ParameterList* plist);
Teuchos::ParameterList get_Cycle_Macro(const std::string& macro_name, Teuchos::ParameterList* plist);
Teuchos::ParameterList get_Domain_List(Teuchos::ParameterList* plist);
Teuchos::ParameterList get_Mesh_List(Teuchos::ParameterList* plist);
Teuchos::ParameterList get_Regions_List(Teuchos::ParameterList* plist);
Teuchos::Array<std::string> get_Variable_Macro(const std::string& macro_name, Teuchos::ParameterList* plist);

void init_global_info(Teuchos::ParameterList* plist);

void output_boundary_conditions(Teuchos::ParameterList* plist);
void check_AmanziInputVersion(Teuchos::ParameterList* plist);

static std::string phase_name;
static std::string phase_comp_name;
static Teuchos::Array<std::string> comp_names;
static std::map<std::string, int> comp_names_map;
static Teuchos::Array<std::string> mineral_names_;
static Teuchos::Array<std::string> sorption_site_names_;

static double constant_density;
static int spatial_dimension_;  
static bool need_dispersion_;

static std::string verbosity_level("low");
static int numproc_;

}
}
