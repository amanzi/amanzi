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
#include "Mesh.hh"
#include "PK.hh"
#include "State.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Chemistry_PK : public PK {
 public:
  Chemistry_PK();
  virtual ~Chemistry_PK() {};

  // required members
  virtual void Setup();
  virtual void Initialize();

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) = 0;
  virtual void CommitStep(double t_old, double t_new) = 0;

  // -- returns the (maximum) time step allowed for this chemistry PK.
  virtual double get_dt() = 0;

  // -- temporary members
  virtual void set_dt(double dt) {};
  virtual void CalculateDiagnostics() {};
  virtual std::string name() { return "chemistry"; }

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

  Teuchos::RCP<VerboseObject> vo_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
