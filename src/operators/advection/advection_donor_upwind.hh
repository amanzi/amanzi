/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Donor upwind advection.
   ------------------------------------------------------------------------- */

#ifndef OPERATOR_ADVECTION_ADVECTION_DONOR_UPWIND_HH_
#define OPERATOR_ADVECTION_ADVECTION_DONOR_UPWIND_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_IntVector.h"

#include "Mesh.hh"
#include "composite_vector.hh"

#include "advection.hh"

namespace Amanzi {
namespace Operators {

class AdvectionDonorUpwind : public Advection {

public:

  AdvectionDonorUpwind(Teuchos::ParameterList& advect_plist,
                       Teuchos::RCP<AmanziMesh::Mesh> mesh);

  virtual void set_flux(Teuchos::RCP<const CompositeVector>& flux);
  virtual void Advect();

private:
  void IdentifyUpwindCells_();

  int f_min_;
  int f_count_;
}

} // namespace Operators
} // namespace Amanzi
