#ifndef AMANZI_CHEMISTRY_STATE_NEW_HH_
#define AMANZI_CHEMISTRY_STATE_NEW_HH_

#include "State.hh"
#include "beaker.hh"

#ifdef ALQUIMIA_ENABLED
#include "Teuchos_RCP.hpp"
#include "ChemistryEngine.hh"
#endif

namespace Amanzi {
namespace AmanziChemistry {

class Chemistry_State {
 public:
  Chemistry_State(Teuchos::ParameterList& plist,
                  const std::vector<std::string>& component_names,
                  const Teuchos::RCP<State>& S);

  virtual ~Chemistry_State() {};

  void Setup();
  void AllocateAdditionalChemistryStorage(const Beaker::BeakerComponents&);
  void SetAuxDataNames(const std::vector<std::string>& aux_data_names);

  void Initialize();

  // access methods for state variables
  int number_of_aqueous_components(void) const { return number_of_aqueous_components_; }
  int number_of_ion_exchange_sites() const { return number_of_ion_exchange_sites_; }
  int number_of_sorption_sites() const { return number_of_sorption_sites_; }
  bool using_sorption() const { return using_sorption_; }
  bool using_sorption_isotherms() const { return using_sorption_isotherms_; }

#ifdef ALQUIMIA_ENABLED
  // Copies the chemistry state in the given cell to the given Alquimia containers.
   void CopyToAlquimia(const int cell_id,
                       AlquimiaMaterialProperties& mat_props,
                       AlquimiaState& state,
                       AlquimiaAuxiliaryData& aux_data);
  
  // Copies the chemistry state in the given cell to the given Alquimia containers, 
  // taking the aqueous components from the given multivector.
  void CopyToAlquimia(const int cell_id,
                      Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                      AlquimiaMaterialProperties& mat_props,
                      AlquimiaState& state,
                      AlquimiaAuxiliaryData& aux_data);

  // Copies the data in the given Alquimia containers to the given cell within the 
  // chemistry state. The aqueous component concentrations are placed into 
  // the aqueous_components multivector.
  void CopyFromAlquimia(const int cell_id,
                        const AlquimiaMaterialProperties& mat_props,
                        const AlquimiaState& state,
                        const AlquimiaAuxiliaryData& aux_data,
                        const AlquimiaAuxiliaryOutputData& aux_output,
                        Teuchos::RCP<const Epetra_MultiVector> aqueous_components);
#endif

 protected:
  void InitializeField_(Teuchos::ParameterList& ic_plist, std::string fieldname,
                        bool sane_default, double default_val);

  void SetupMineralNames_();
  void SetupSorptionSiteNames_();
  void ParseMeshBlocks_();
  void VerifyMineralogy_(const std::string& region_name,
                         const Teuchos::ParameterList& minerals_list);
  void VerifySorptionIsotherms_(const std::string& region_name,
                                const Teuchos::ParameterList& isotherms_list);
  void VerifySorptionSites_(const std::string& region_name,
                            const Teuchos::ParameterList& sorption_site_list);
  void RequireData_();
  void RequireAuxData_();
  
 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<State> S_;
  bool ghosted_;
  std::string name_;

  Teuchos::ParameterList plist_;

  int number_of_aqueous_components_;
  int number_of_minerals_;
  int number_of_ion_exchange_sites_;
  int number_of_sorption_sites_;
  bool using_sorption_;
  bool using_sorption_isotherms_;

  std::vector<std::string> compnames_;
  std::map<std::string,int> comp_name_id_map_;
  std::vector<std::string> mineral_names_;
  std::map<std::string, int> mineral_name_id_map_;
  std::vector<std::string> sorption_site_names_;
  std::map<std::string, int> sorption_site_name_id_map_;

  // Auxiliary data, maintained by Amanzi and updated
  int num_aux_data_;
  Teuchos::RCP<Epetra_MultiVector> aux_data_;
};

}
}

#endif
