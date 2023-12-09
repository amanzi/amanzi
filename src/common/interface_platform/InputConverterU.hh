/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Input Converter

*/

#ifndef AMANZI_INPUT_CONVERTER_UNSTRUCTURED_HH_
#define AMANZI_INPUT_CONVERTER_UNSTRUCTURED_HH_

// TPLs
#include "xercesc/dom/DOM.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

// Amanzi's
#include "Key.hh"
#include "VerboseObject.hh"

#include "InputConverter.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

typedef std::map<std::string, Teuchos::RCP<Teuchos::ParameterList>> PK;
typedef std::map<std::string, std::vector<std::string>> Tree;

struct Phase {
 public:
  Phase() : active(false){};

  bool active;
  std::string primary;
  std::string model;
  std::vector<std::string> dissolved;
};
typedef std::map<std::string, Phase> PhaseTree;


struct BCs {
 public:
  BCs() : mol_mass(-1.0){};
  BCs(double value) : mol_mass(value){};

  std::string type;
  std::vector<double> times, values, fluxes;
  std::vector<std::string> forms, formulas;

  std::string filename, xheader, yheader, variable;
  double mol_mass;
};


class InputConverterU : public InputConverter {
 public:
  explicit InputConverterU(const std::string& input_filename)
    : InputConverter(input_filename),
      multiphase_(false),
      isothermal_(true),
      gravity_on_(true),
      const_gravity_(GRAVITY_MAGNITUDE),
      const_atm_pressure_(ATMOSPHERIC_PRESSURE),
      fracture_network_(false),
      flow_single_phase_(false),
      compressibility_(false),
      beta_(0.0),
      mesh_rectangular_(false),
      use_transport_porosity_(false),
      use_transport_dispersion_(true),
      transport_permeability_(false),
      transport_implicit_(false),
      restart_(false),
      ic_time_flow_(0.0),
      ic_time_(0.0),
      output_prefix_(""),
      io_walkabout_(false),
      io_mesh_info_(false),
      vo_(NULL){};

  explicit InputConverterU(const std::string& input_filename,
                           xercesc::DOMDocument* input_doc,
                           const std::string& output_prefix)
    : InputConverter(input_filename, input_doc),
      multiphase_(false),
      isothermal_(true),
      gravity_on_(true),
      const_gravity_(GRAVITY_MAGNITUDE),
      const_atm_pressure_(ATMOSPHERIC_PRESSURE),
      fracture_network_(false),
      flow_single_phase_(false),
      compressibility_(false),
      beta_(0.0),
      mesh_rectangular_(false),
      use_transport_porosity_(false),
      use_transport_dispersion_(true),
      transport_permeability_(false),
      transport_implicit_(false),
      restart_(false),
      ic_time_flow_(0.0),
      ic_time_(0.0),
      output_prefix_(output_prefix),
      io_walkabout_(false),
      io_mesh_info_(false),
      vo_(NULL){};

  ~InputConverterU()
  {
    if (vo_ != NULL) delete vo_;
  }

  // main members
  Teuchos::ParameterList Translate(int rank, int num_proc);
  void SaveXMLFile(Teuchos::ParameterList& plist, std::string& filename);

 private:
  void VerifyXMLStructure_();
  void ParseSolutes_();
  void ParseModelDescription_();
  void ParseFractureNetwork_();
  void ModifyDefaultPhysicalConstants_();
  void ParseGlobalNumericalControls_();

  BCs ParseCondList_(DOMNode* node,
                     double vmin,
                     double vmax,
                     const std::string& unit,
                     bool is_bc = true);
  BCs ParseCondList_(DOMNode* node,
                     const std::string& bctype,
                     double vmin,
                     double vmax,
                     const std::string& unit,
                     bool is_bc = true);
  BCs ParseCondList_(std::vector<DOMNode*>& same_list,
                     const std::string& bctype,
                     double vmin,
                     double vmax,
                     const std::string& unit,
                     bool is_bc = true,
                     const std::string& filter_name = "");

  Teuchos::ParameterList TranslateVerbosity_();
  Teuchos::ParameterList TranslateUnits_();

  Teuchos::ParameterList TranslateMesh_();
  Teuchos::ParameterList TranslateRegions_();
  Teuchos::ParameterList TranslateOutput_();
  Teuchos::ParameterList TranslatePreconditioners_();
  Teuchos::ParameterList TranslateTrilinosML_();
  Teuchos::ParameterList TranslateHypreAMG_();
  Teuchos::ParameterList TranslateBILU_();
  Teuchos::ParameterList TranslateILU_();
  Teuchos::ParameterList
  TranslateLinearSolvers_(std::string tags, std::string method_default, std::string method_enforce);
  Teuchos::ParameterList TranslateSolvers_();
  Teuchos::ParameterList TranslateState_();
  Teuchos::ParameterList TranslateMaterialsPartition_();
  Teuchos::ParameterList TranslateCycleDriver_();
  Teuchos::ParameterList TranslateCycleDriverNew_();
  Teuchos::ParameterList TranslateTimePeriodControls_();
  Teuchos::ParameterList TranslatePKs_(Teuchos::ParameterList& glist);
  Teuchos::ParameterList TranslateDiffusionOperator_(const std::string& disc_methods,
                                                     const std::string& pc_method,
                                                     const std::string& nonlinear_solver,
                                                     const std::string& nonlinear_coef,
                                                     const std::string& extensions,
                                                     const std::string& domain,
                                                     bool gravity,
                                                     const std::string& pk = "flow");
  Teuchos::ParameterList TranslateTimeIntegrator_(const std::string& err_options,
                                                  const std::string& nonlinear_solver,
                                                  bool modify_correction,
                                                  const std::string& controls,
                                                  const std::string& linsolver,
                                                  double dt_cut_default,
                                                  double dt_inc_default);
  Teuchos::ParameterList TranslateInitialization_(const std::string& unstr_controls);

  // -- general
  Teuchos::ParameterList TranslateSources_(const std::string& domain, const std::string& pkname);

  // -- state
  void TranslateCommonContinuumFields_(const std::string& domain,
                                       Teuchos::ParameterList& out_ic,
                                       Teuchos::ParameterList& out_ev);

  void TranslateFieldEvaluator_(DOMNode* node,
                                const std::string& field,
                                const std::string& unit,
                                const std::string& reg_str,
                                const std::vector<std::string>& regions,
                                Teuchos::ParameterList& out_ic,
                                Teuchos::ParameterList& out_ev,
                                std::string data_key = "value",
                                std::string domain = "domain");
  void TranslateFieldIC_(DOMNode* node,
                         std::string field,
                         std::string unit,
                         const std::string& reg_str,
                         const std::vector<std::string>& regions,
                         Teuchos::ParameterList& out_ic,
                         std::string data_key = "value",
                         const std::vector<std::string>& components = { "cell" });

  void AddIndependentFieldEvaluator_(Teuchos::ParameterList& out_ev,
                                     const std::string& field,
                                     const std::string& region,
                                     const std::string& comp,
                                     double val);

  void AddSecondaryFieldEvaluator_(Teuchos::ParameterList& out_ev,
                                   const Key& field,
                                   const Key& key,
                                   const std::string& type,
                                   const std::string& eos_table_name);

  void AddConstantFieldInitialization_(Teuchos::ParameterList& out_ev,
                                       const std::string& field,
                                       const std::string& region,
                                       double val);

  // -- flow
  Teuchos::ParameterList
  TranslateFlow_(const std::string& mode, const std::string& domain, const std::string& pk_model);
  Teuchos::ParameterList TranslateWRM_(const std::string& pk_name);
  Teuchos::ParameterList TranslatePOM_(const std::string& domain);
  Teuchos::ParameterList TranslatePPM_(const std::string& domain);
  Teuchos::ParameterList TranslateFAM_(const std::string& domain);
  Teuchos::ParameterList TranslateFlowMSM_();
  Teuchos::ParameterList TranslateFlowBCs_(const std::string& domain);
  Teuchos::ParameterList TranslateFlowFractures_(const std::string& domain);

  // -- transport
  Teuchos::ParameterList TranslateTransport_(const std::string& domain);
  Teuchos::ParameterList TranslateMolecularDiffusion_();
  Teuchos::ParameterList TranslateTransportMSM_();
  Teuchos::ParameterList TranslateTransportBCs_(const std::string& domain);
  void TranslateTransportBCsGroup_(std::string& bcname,
                                   std::vector<std::string>& regions,
                                   xercesc::DOMNodeList* solutes,
                                   Teuchos::ParameterList& out_list);
  Teuchos::ParameterList TranslateTransportSources_();
  void TranslateTransportSourcesGroup_(std::string& srcname,
                                       std::vector<std::string>& regions,
                                       xercesc::DOMNodeList* solutes,
                                       xercesc::DOMNode* phase_l,
                                       Teuchos::ParameterList& out_list);
  void TranslateTransportGeochemistry_(DOMNode* node,
                                       std::string& bcname,
                                       std::vector<std::string>& regions,
                                       Teuchos::ParameterList& out_list);

  // -- backward compatibility
  void TranslateTransportBCsAmanziGeochemistry_(Teuchos::ParameterList& out_list);
  void TranslateStateICsAmanziGeochemistry_(Teuchos::ParameterList& out_list,
                                            std::string& constraint,
                                            std::vector<std::string>& regions,
                                            const std::string& domain);

  // -- chemistry and energy
  Teuchos::ParameterList TranslateChemistry_(const std::string& domain);
  Teuchos::ParameterList TranslateEnergy_(const std::string& domain, const std::string& pk_model);
  Teuchos::ParameterList TranslateEnergyBCs_(const std::string& domain);

  // -- multiphase
  bool multiphase_, isothermal_;
  Teuchos::ParameterList
  TranslateMultiphase_(const std::string& domain, Teuchos::ParameterList& state_list);
  Teuchos::ParameterList TranslateMultiphaseBCs_();

  // -- shallow water
  Teuchos::ParameterList TranslateShallowWater_(const std::string& domain);
  Teuchos::ParameterList TranslateShallowWaterBCs_();

  // -- mpc pks
  bool coupled_flow_, coupled_transport_, coupled_energy_, coupled_multiphase_;
  std::vector<std::string> fracture_regions_, surface_regions_;

  void ProcessMacros_(const std::string& prefix,
                      char* text_content,
                      Teuchos::ParameterList& mPL,
                      Teuchos::ParameterList& outPL);

  void PopulatePKTree_(Teuchos::ParameterList& pk_tree, const std::string pk_name);
  void RegisterPKsList_(Teuchos::ParameterList& pk_tree, Teuchos::ParameterList& pks_list);

  void FinalizeMPC_PKs_(Teuchos::ParameterList& glist);

  Teuchos::ParameterList CreateAnalysis_();
  Teuchos::ParameterList CreateRegionAll_();

  // -- complex functions
  void TranslateFunctionGaussian_(const std::vector<double>& data, Teuchos::ParameterList& bcfn);
  void TranslateFunctionGradient_(double refv,
                                  std::vector<double>& grad,
                                  std::vector<double>& refc,
                                  Teuchos::ParameterList& bcfn);

  void FilterEmptySublists_(Teuchos::ParameterList& plist);
  void MergeInitialConditionsLists_(Teuchos::ParameterList& plist, const std::string& chemistry);

  bool TranslateGenericMath_(const BCs& bcs, Teuchos::ParameterList& bcfn);

  // -- sort functions
  template <class Iterator>
  Iterator SelectUniqueEntries(Iterator first, Iterator last);

  // -- miscalleneous
  DOMNode* GetPKChemistryPointer_(bool& flag);
  bool FindNameInVector_(const std::string& name, const std::vector<std::string>& list);
  std::string CreateNameFromVector_(const std::vector<std::string>& list);
  bool WeightVolumeSubmodel_(const std::vector<std::string>& regions);
  std::string CreateUniqueName_(const Teuchos::Array<std::string>& list);
  void PrintStatistics_();
  bool HasSubmodel_(const std::string& model, const std::string& submodel);

 private:
  int dim_;
  int rank_, num_proc_;
  std::vector<std::string> coords_;

  Tree tree_;
  PhaseTree phases_;

  // global list
  Teuchos::RCP<Teuchos::ParameterList> glist_;

  // global data
  std::map<std::string, std::string> pk_domain_, pk_region_;
  std::map<std::string, std::vector<std::string>> pk_model_;
  std::map<std::string, bool> pk_master_;
  std::map<std::string, double> dt_cut_, dt_inc_;

  std::string eos_lookup_table_, eos_model_;

  // global physical constants prefixed with "const"
  bool gravity_on_;
  double const_gravity_;
  double const_atm_pressure_;

  // global physics
  bool fracture_network_;

  // global flow constants
  std::string flow_model_; // global value
  bool flow_single_phase_; // runtime value
  bool compressibility_;
  double rho_, beta_;

  // global mesh data
  bool mesh_rectangular_;
  std::map<std::string, int> region_type_; // flag for vofs

  // global transport and chemistry constants
  bool use_transport_porosity_, use_transport_dispersion_;
  bool transport_permeability_, transport_implicit_;
  std::vector<std::string> comp_names_all_;
  std::map<std::string, double> solute_molar_mass_;

  // global state parameters
  // -- initialization filename, different from restart
  bool restart_;
  std::string init_filename_;
  double ic_time_flow_, ic_time_;

  // global solvers parameters
  std::vector<std::pair<std::string, double>> gmres_solvers_;

  // global output parameters
  std::string output_prefix_;
  bool io_walkabout_, io_mesh_info_;

  // global names for visualization
  std::vector<std::string> material_regions_;
  std::vector<std::string> material_names_;
  std::vector<int> material_ids_;

  // for analysis
  std::vector<std::string> transport_diagnostics_;

  std::vector<std::string> vv_bc_regions_;
  std::vector<std::string> vv_src_regions_;
  std::vector<std::string> vv_obs_regions_;

  // for statistics
  Teuchos::ParameterList verb_list_;
  VerboseObject* vo_;
};


/* ******************************************************************
* Short functions
****************************************************************** */
inline bool
InputConverterU::HasSubmodel_(const std::string& model, const std::string& submodel)
{
  int n = pk_model_[model].size();
  for (int i = 0; i < n; ++i) {
    if (pk_model_[model][i] == submodel) return true;
  }
  return false;
}


/* ******************************************************************
* Selects unique entries and places them in [first, last)
****************************************************************** */
template <class Iterator>
Iterator
InputConverterU::SelectUniqueEntries(Iterator first, Iterator last)
{
  while (first != last) {
    Iterator next(first);
    last = std::remove(++next, last, *first);
    first = next;
  }
  return last;
}

} // namespace AmanziInput
} // namespace Amanzi

#endif
