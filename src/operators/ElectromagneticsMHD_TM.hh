/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  TM mode of MHD: B field is on faces, E field is at nodes.
*/

#ifndef AMANZI_OPERATOR_ELECTROMAGNETICS_MHD_TM_HH_
#define AMANZI_OPERATOR_ELECTROMAGNETICS_MHD_TM_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi::Operators
#include "ElectromagneticsMHD.hh"

namespace Amanzi {
namespace Operators {

class ElectromagneticsMHD_TM : public ElectromagneticsMHD {
 public:
  ElectromagneticsMHD_TM(const Teuchos::RCP<Operator>& global_op)
    : ElectromagneticsMHD(global_op)
  {};

  ElectromagneticsMHD_TM(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : ElectromagneticsMHD(plist, mesh)
  {};

  // main virtual members
  // -- before solving the problem
  virtual void ModifyMatrices(CompositeVector& E, CompositeVector& B, double dt);
  virtual void ApplyBCs(bool primary, bool eliminate);

  // -- after solving the problem
  virtual void ModifyFields(CompositeVector& E, CompositeVector& B, double dt);

 private:
  void ApplyBCs_Node_(const Teuchos::Ptr<BCs>& bc_f,
                      const Teuchos::Ptr<BCs>& bc_v, bool primary, bool eliminate);
};

}  // namespace Operators
}  // namespace Amanzi

#endif


