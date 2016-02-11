/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Base class for chemical process kernels. It should be never
  instantiated.
*/
 
#ifndef AMANZI_CHEMISTRY_PK_HH_
#define AMANZI_CHEMISTRY_PK_HH_

#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "VerboseObject.hh"

// Amanzi
#ifdef ALQUIMIA_ENABLED
#include "ChemistryEngine.hh"
#endif
#include "Mesh.hh"
#include "PK.hh"
#include "State.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Chemistry_PK : public PK {
 public:
  Chemistry_PK();
  virtual ~Chemistry_PK() {};

  // required members for PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // Required members for chemistry interface
  // -- output of auxillary cellwise data from chemistry
  virtual Teuchos::RCP<Epetra_MultiVector> extra_chemistry_output_data() = 0;

  // Basic capabilities
  // -- get/set auxiliary tcc vector that now contains only aqueous components.
  Teuchos::RCP<Epetra_MultiVector> aqueous_components() { return aqueous_components_; } 
  void set_aqueous_components(Teuchos::RCP<Epetra_MultiVector> tcc) { aqueous_components_ = tcc; }

  // -- process various objects before/during setup phase
  void InitializeMinerals(Teuchos::RCP<Teuchos::ParameterList> plist);
  void InitializeSorptionSites(Teuchos::RCP<Teuchos::ParameterList> plist,
                               Teuchos::RCP<Teuchos::ParameterList> state_list);

  // -- access
#ifdef ALQUIMIA_ENABLED
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine() { return chem_engine_; }
#endif

  // -- output of error messages.
  void ErrorAnalysis(int ierr, std::string& internal_msg);

 private:
  void InitializeField_(std::string fieldname, double default_val);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<State> S_;
  std::string passwd_;

  int number_aqueous_components_;
  std::vector<std::string> comp_names_;
  Teuchos::RCP<Epetra_MultiVector> aqueous_components_;

  int number_minerals_;
  std::vector<std::string> mineral_names_;

  int number_sorption_sites_, number_total_sorbed_;
  std::vector<std::string> sorption_site_names_;
  bool using_sorption_, using_sorption_isotherms_;

  int number_free_ion_, number_ion_exchange_sites_;

#ifdef ALQUIMIA_ENABLED
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
#endif

  // verbosity object thatis not shared with common chemistry
  Teuchos::RCP<VerboseObject> vo_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
