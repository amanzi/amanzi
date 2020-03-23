/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples chemistry PKs in matrix and fracture.
*/

#ifndef AMANZI_CHEMISTRY_MATRIX_FRACTURE_PK_HH_
#define AMANZI_CHEMISTRY_MATRIX_FRACTURE_PK_HH_

#include "Teuchos_RCP.hpp"

#include "PK_MPCWeak.hh"
#include "PK_Factory.hh"

namespace Amanzi {

class ChemistryMatrixFracture_PK : public PK_MPCWeak {
 public:
  ChemistryMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- setup
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- miscaleneous methods
  virtual std::string name() { return "coupled chemistry"; } 

 private:
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_, mesh_fracture_;

  // factory registration
  static RegisteredPKFactory<ChemistryMatrixFracture_PK> reg_;
};

}  // namespace Amanzi
#endif
