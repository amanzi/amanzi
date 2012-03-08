/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Base interface for a general-purpose advection operator.
   ------------------------------------------------------------------------- */

#ifndef OPERATOR_ADVECTION_ADVECTION_HH_
#define OPERATOR_ADVECTION_ADVECTION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"
#include "composite_vector.hh"

namespace Amanzi {
namespace Operators {

class Advection {

public:

  Teuchos::RCP<const CompositeVector> flux() const { return flux_; }
  virtual void set_flux(Teuchos::RCP<const CompositeVector>& flux);

  int num_dofs() const { return num_dofs_; }
  virtual void set_num_dofs(int num_dofs);

  Teuchos::RCP<CompositeVector> field() { return field_; }

  virtual void Advect() = 0;

protected:
  int num_dofs_;
  Teuchos::RCP<const CompositeVector> flux_;
  Teuchos::RCP<CompositeVector> field_;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList advect_plist_;
}

} // namespace Operators
} // namespace Amanzi
