/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Simple factory of advection schemes.
*/

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_FACTORY_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "PDE_AdvectionUpwind.hh"
#include "PDE_AdvectionUpwindDFN.hh"
#include "PDE_AdvectionUpwindFracturedMatrix.hh"

namespace Amanzi {
namespace Operators {

class PDE_AdvectionUpwindFactory {
 public:
  PDE_AdvectionUpwindFactory() {}
  PDE_AdvectionUpwindFactory(Teuchos::ParameterList& oplist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : oplist_(oplist), mesh_(mesh)
  {}

  Teuchos::RCP<PDE_AdvectionUpwind> Create(const Teuchos::RCP<Operator>& global_op = Teuchos::null);

  // backward compatibility
  Teuchos::RCP<PDE_AdvectionUpwind>
  Create(Teuchos::ParameterList& oplist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  Teuchos::RCP<PDE_AdvectionUpwind>
  Create(Teuchos::ParameterList& oplist, const Teuchos::RCP<Operator>& global_op);

 private:
  Teuchos::ParameterList oplist_;
  const Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

} // namespace Operators
} // namespace Amanzi

#endif
