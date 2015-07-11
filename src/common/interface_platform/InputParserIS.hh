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

#include "VerboseObject.hh"
#include "InputParserIS_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

void InputParserIS_OutputBCs(Teuchos::ParameterList* plist);

class Phase {
 public:
  Phase() {};
  ~Phase() {};

 public:
  std::string name;
  std::string solute_name;  // We assume only one component
  std::vector<std::string> solute_comp_names;
};


class InputParserIS {
 public:
  InputParserIS() : vo_(NULL) {
    flow_single_phase = false;
    verbosity_level = "low";
    use_picard_ = USE_PICARD;
  };
  ~InputParserIS() {
    if (vo_ != NULL) delete vo_;
  };

  // main members
  Teuchos::ParameterList Translate(Teuchos::ParameterList* plist, int numproc);
  void CheckCollectedRegions();

 private:
  // flow 
  Teuchos::ParameterList CreateFlowList_(Teuchos::RCP<Teuchos::ParameterList>& plist, int time_regime = BOTH_REGIMES);
  Teuchos::ParameterList CreateFlowSrcList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateSS_FlowBC_List_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateFlowOperatorList_(
     const std::string& disc_method, const std::string& prec_method,
     const std::string& nonlinear_solver, const std::string& rel_perm);
  Teuchos::ParameterList CreateWRM_List_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreatePOM_List_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateInitializationList_(Teuchos::ParameterList& ini_list,
                                                   bool have_picard, Teuchos::ParameterList& pic_list);

  // transport
  Teuchos::ParameterList CreateTransportList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateTransportSrcList_(Teuchos::RCP<Teuchos::ParameterList>& plist);

  // chemistry
  Teuchos::ParameterList CreateChemistryList_(Teuchos::RCP<Teuchos::ParameterList>& plist);

  // solvers and preconditioners
  Teuchos::ParameterList CreateSolversList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateNonlinearSolversList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreatePreconditionersList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateHypreAMG_List_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateDPC_List_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateBILU_List_(Teuchos::RCP<Teuchos::ParameterList>& plist);

  // output 
  Teuchos::ParameterList CreateVisualizationDataList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateCheckpointDataList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateWalkaboutDataList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateObservationDataList_(Teuchos::RCP<Teuchos::ParameterList>& plist);

  Teuchos::ParameterList CreateTimeMacro_(const std::string& macro, Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateCycleMacro_(const std::string& macro, Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::Array<std::string> CreateVariableMacro_(Teuchos::Array<std::string>& macros, Teuchos::ParameterList* plist);

  // state
  Teuchos::ParameterList CreateStateList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreatePartitionList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  
  // MPC and state
  Teuchos::ParameterList CreateTimePeriodControlList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateCycleDriverList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  void CreatePKslist_(Teuchos::ParameterList& cycle_driver_list, Teuchos::ParameterList& pks_list);
  void RegisterPKlist_(Teuchos::ParameterList& pk_tree, Teuchos::ParameterList& pks_list);
  void FillPKslist_(Teuchos::RCP<Teuchos::ParameterList>& plist, Teuchos::ParameterList& pks_list);

  // mesh and geometry
  Teuchos::ParameterList CreateMeshList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CopyDomainList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CopyMeshList_(Teuchos::ParameterList* plist);
  Teuchos::ParameterList CopyRegionsList_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  
  // other main members
  void CheckAmanziInputVersion_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  void InitGlobalInfo_(Teuchos::RCP<Teuchos::ParameterList>& plist);
  Teuchos::ParameterList CreateVerbosityList_(const std::string& vlevel);
  Teuchos::ParameterList CreateAnalysisList_();
  Teuchos::Array<std::string> TranslateForms_(Teuchos::Array<std::string>& forms);
  void PrintUnused_(const Teuchos::ParameterList& p, VerboseObject* vo) const;

 private:
  std::vector<Phase> phases_;

  std::vector<std::string> comp_names_;  // aqueous components 
  std::vector<std::string> comp_names_all_;  // all components from aqueous and gaseous phases

  Teuchos::Array<std::string> mineral_names_;
  Teuchos::Array<std::string> sorption_site_names_;
  std::string chemistry_model_;

  double constant_density;
  int spatial_dimension_;  
  bool flow_single_phase, use_picard_, compressibility_;
  bool need_dispersion_;
  std::vector<std::string> transport_diagnostics_;

  std::string verbosity_level;
  int numproc_;

 private:
  std::vector<std::string> vv_bc_regions;  // XML verification
  std::vector<std::string> vv_src_regions;
  std::vector<std::string> vv_obs_regions;

 protected:
  VerboseObject* vo_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
