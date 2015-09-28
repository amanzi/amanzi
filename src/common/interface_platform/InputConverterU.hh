/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_INPUT_CONVERTER_UNSTRUCTURED_HH_
#define AMANZI_INPUT_CONVERTER_UNSTRUCTURED_HH_

// TPLs
#include "xercesc/dom/DOM.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

// Amanzi's
#include "VerboseObject.hh"

#include "InputConverter.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

typedef std::map<std::string, Teuchos::RCP<Teuchos::ParameterList> > PK;
typedef std::map<std::string, std::vector<std::string> > Tree;

class InputConverterU : public InputConverter {
 public:
  explicit InputConverterU(const std::string& input) :
      InputConverter(input), 
      vo_(NULL),
      flow_single_phase_(false),
      compressibility_(false),
      mesh_rectangular_(false),
      transport_permeability_(false),
      restart_(false),
      init_filename_(input) {};

  explicit InputConverterU(xercesc::DOMDocument* input) :
      InputConverter(input), 
      vo_(NULL),
      flow_single_phase_(false),
      compressibility_(false),
      mesh_rectangular_(false),
      transport_permeability_(false),
      restart_(false),
      init_filename_("") {};

  ~InputConverterU() { if (vo_ != NULL) delete vo_; }

  // main members
  Teuchos::ParameterList Translate(int rank, int num_proc);
  void SaveXMLFile(Teuchos::ParameterList& plist, std::string& filename);

 private:
  void VerifyXMLStructure_();
  void ParseSolutes_();
  void ParseModelDescription_();

  Teuchos::ParameterList TranslateVerbosity_();
  Teuchos::ParameterList TranslateMisc_();

  Teuchos::ParameterList TranslateMesh_();
  Teuchos::ParameterList TranslateRegions_();
  Teuchos::ParameterList TranslateOutput_();
  Teuchos::ParameterList TranslatePreconditioners_();
  Teuchos::ParameterList TranslateTrilinosML_();
  Teuchos::ParameterList TranslateHypreAMG_();
  Teuchos::ParameterList TranslateBILU_();
  Teuchos::ParameterList TranslateSolvers_();
  Teuchos::ParameterList TranslateState_();
  Teuchos::ParameterList TranslateMaterialsPartition_();
  Teuchos::ParameterList TranslateCycleDriver_();
  Teuchos::ParameterList TranslateCycleDriverNew_();
  Teuchos::ParameterList TranslateTimePeriodControls_();
  Teuchos::ParameterList TranslatePKs_(const Teuchos::ParameterList& cd_list);
  Teuchos::ParameterList TranslateDiffusionOperator_(
      const std::string& disc_methods, const std::string& pc_method,
      const std::string& nonlinear_solver, const std::string& extensions, bool gravity);
  Teuchos::ParameterList TranslateTimeIntegrator_(
      const std::string& err_options, const std::string& nonlinear_solver,
      bool modify_correction, const std::string& unstr_controls,
      double dt_cut_default, double dt_inc_default);
  Teuchos::ParameterList TranslateInitialization_(
      const std::string& unstr_controls);

  // -- flow
  Teuchos::ParameterList TranslateFlow_(const std::string& mode);
  Teuchos::ParameterList TranslateWRM_();
  Teuchos::ParameterList TranslatePOM_();
  Teuchos::ParameterList TranslateFlowMSM_();
  Teuchos::ParameterList TranslateFlowBCs_();
  Teuchos::ParameterList TranslateFlowSources_();

  // -- transport
  Teuchos::ParameterList TranslateTransport_();
  Teuchos::ParameterList TranslateTransportMSM_();
  Teuchos::ParameterList TranslateTransportBCs_();
  void TranslateTransportBCsGroup_(
      std::string& bcname, std::vector<std::string>& regions,
      xercesc::DOMNodeList* solutes, Teuchos::ParameterList& out_list);
  Teuchos::ParameterList TranslateTransportSources_();
  void TranslateTransportSourcesGroup_(
      std::string& srcname, std::vector<std::string>& regions,
      xercesc::DOMNodeList* solutes, xercesc::DOMNode* phase_l, Teuchos::ParameterList& out_list);

  // -- chemistry and energy
  Teuchos::ParameterList TranslateChemistry_();
  Teuchos::ParameterList TranslateEnergy_();
  Teuchos::ParameterList TranslateEnergyBCs_();

  void ProcessMacros_(const std::string& prefix, char* text_content,
                      Teuchos::ParameterList& mPL, Teuchos::ParameterList& outPL);

  void RegisterPKsList_(Teuchos::ParameterList& pk_tree, Teuchos::ParameterList& pks_list);

  Teuchos::ParameterList CreateAnalysis_();
  Teuchos::ParameterList CreateRegionAll_();
  std::string CreateBGDFile(std::string& filename);

  void FilterEmptySublists_(Teuchos::ParameterList& plist);

 private:
  int dim_;
  int rank_, num_proc_;
  std::vector<std::string> coords_;

  Tree tree_;
  Tree phases_;

  // global data
  std::map<std::string, std::string> pk_model_;
  std::map<std::string, bool> pk_master_;
  std::map<std::string, double> dt_cut_, dt_inc_;

  // global flow constants
  std::string flow_model_;  // global value
  bool flow_single_phase_;  // runtime value
  bool compressibility_;
  double rho_;

  // global mesh data
  bool mesh_rectangular_;

  // global transport and chemistry constants
  bool transport_permeability_;
  std::vector<std::string> comp_names_all_;

  // global state parameters
  // -- initialization filename, different from restart
  bool restart_;
  std::string init_filename_;

  // for analysis
  std::vector<std::string> transport_diagnostics_;

  std::vector<std::string> vv_bc_regions_;
  std::vector<std::string> vv_src_regions_;
  std::vector<std::string> vv_obs_regions_;

  Teuchos::ParameterList verb_list_;
  VerboseObject* vo_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
