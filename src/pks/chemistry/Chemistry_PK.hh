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
#include "State.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Chemistry_PK {
 public:
  Chemistry_PK();
  virtual ~Chemistry_PK() {};

  // required members
  virtual void Setup();
  virtual void Initialize();

  virtual void Advance(const double& dt,
                       Teuchos::RCP<Epetra_MultiVector> total_component_concentration) = 0;
  virtual void CommitState(const double& dt) = 0;

  // -- returns the (maximum) time step allowed for this chemistry PK.
  virtual double time_step() const = 0;

  // -- output of auxillary cellwise data from chemistry
  virtual Teuchos::RCP<Epetra_MultiVector> get_extra_chemistry_output_data() = 0;

  // -- process various objects
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
