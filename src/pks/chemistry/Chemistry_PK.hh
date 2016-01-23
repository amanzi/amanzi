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

#include "Chemistry_State.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Chemistry_PK {
 public:
  Chemistry_PK();
  virtual ~Chemistry_PK() {};

  // required members
  // -- initialization
  virtual void InitializeChemistry() = 0;

  virtual void Advance(const double& dt,
                       Teuchos::RCP<Epetra_MultiVector> total_component_concentration) = 0;
  virtual void CommitState(Teuchos::RCP<Chemistry_State> chem_state, const double& dt) = 0;

  // -- returns the (maximum) time step allowed for this chemistry PK.
  virtual double time_step() const = 0;

  // -- output of auxillary cellwise data from chemistry
  virtual Teuchos::RCP<Epetra_MultiVector> get_extra_chemistry_output_data() = 0;

  // Shared capabilities
  // -- register fields with the state
  virtual void Setup();
  
  // -- process various objects
  void InitializeMinerals(Teuchos::RCP<Teuchos::ParameterList> plist);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<State> S_;
  std::string passwd_;

  int number_minerals_;
  std::vector<std::string> mineral_names_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
