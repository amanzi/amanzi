/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Reduced multiphase model: H2O is only in liquid phase.
  Solution vectors: pressure, tcc, saturation.
*/

#ifndef AMANZI_MULTIPHASE_REDUCED_PK_HH_
#define AMANZI_MULTIPHASE_REDUCED_PK_HH_

// Amanzi

// Multiphase
#include "Multiphase_PK.hh"

namespace Amanzi {
namespace Multiphase {

class MultiphaseReduced_PK: public Multiphase_PK {
 public:
  MultiphaseReduced_PK(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln);

  ~MultiphaseReduced_PK() {};

  // interface for PKs
  virtual std::string name() override { return "multiphase reduced"; }

  // interface multiphase models
  virtual void InitializeSolution() override;
  virtual void PopulateBCs() override;

 private:
  int missed_bc_faces_;
};

}  // namespace Multiphase
}  // namespace Amanzi
#endif
