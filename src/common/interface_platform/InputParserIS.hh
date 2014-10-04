#ifndef AMANZI_INPUT_PARSER_IS_HH_
#define AMANZI_INPUT_PARSER_IS_HH_

#include "boost/lambda/lambda.hpp"
#include "boost/bind.hpp"
#include "boost/lexical_cast.hpp"

#define  BOOST_FILESYTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/format.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace Amanzi {
namespace AmanziInput {

#define AMANZI_OLD_INPUT_VERSION_MAJOR 1
#define AMANZI_OLD_INPUT_VERSION_MINOR 2
#define AMANZI_OLD_INPUT_VERSION_MICRO 2

void InputParserIS_OutputBCs(Teuchos::ParameterList* plist);

class InputParserIS {
 public:
  InputParserIS() {
    flow_single_phase = false;
    verbosity_level = "low";
  };
  ~InputParserIS() {};

  // main members
  Teuchos::ParameterList Translate(Teuchos::ParameterList* plist, int numproc);
  void CheckCollectedRegions();

 private:
  // flow 
  Teuchos::ParameterList CreateFlowList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateFlowSrcList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateSS_FlowBC_List_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateFlowOperatorList_(const std::string& disc_method);
  Teuchos::ParameterList CreateWRM_List_(Teuchos::ParameterList* plist);

  // transport
  Teuchos::ParameterList CreateTransportList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateTransportSrcList_(Teuchos::ParameterList* plist);

  // chemistry
  Teuchos::ParameterList CreateChemistryList_(Teuchos::ParameterList* plist);

  // solvers and preconditioners
  Teuchos::ParameterList CreateSolversList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateNonlinearSolversList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreatePreconditionersList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateHypreAMG_List_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateDPC_List_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateBILU_List_(Teuchos::ParameterList* plist);

  // output 
  Teuchos::ParameterList CreateVisualizationDataList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateCheckpointDataList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateWalkaboutDataList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateObservationDataList_(Teuchos::ParameterList* plist);

  Teuchos::ParameterList CreateTimeMacro_(const std::string& macro, Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateCycleMacro_(const std::string& macro, Teuchos::ParameterList* plist);
  Teuchos::Array<std::string> CreateVariableMacro_(Teuchos::Array<std::string>& macros, Teuchos::ParameterList* plist);

  // state
  Teuchos::ParameterList CreateStateList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreatePartitionList_(Teuchos::ParameterList* plist);
  
  // MPC and state
  Teuchos::ParameterList CreateMPC_List_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateTimePeriodControlList_(Teuchos::ParameterList* plist);

  // mesh and geometry
  Teuchos::ParameterList CreateMeshList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CopyDomainList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CopyMeshList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CopyRegionsList_(Teuchos::ParameterList* plist);
  
  // other main members
  void CheckAmanziInputVersion_(Teuchos::ParameterList* plist);
  void InitGlobalInfo_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CreateVerbosityList_(const std::string& vlevel);
  Teuchos::ParameterList CreateAnalysisList_();
  Teuchos::Array<std::string> TranslateForms_(Teuchos::Array<std::string>& forms);

 private:
  std::string phase_name;
  std::string phase_comp_name;
  Teuchos::Array<std::string> comp_names;
  std::map<std::string, int> comp_names_map;
  Teuchos::Array<std::string> mineral_names_;
  Teuchos::Array<std::string> sorption_site_names_;

  double constant_density;
  int spatial_dimension_;  
  bool flow_single_phase;
  bool need_dispersion_;

  std::string verbosity_level;
  int numproc_;

 private:
  std::vector<std::string> vv_bc_regions;  // XML verification
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
