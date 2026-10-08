/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

// PDE_ElasticityFactory constructs objects which implement the interface for a PDE_Elasticity.

#ifndef AMANZI_OPERATOR_PDE_ELASTICITY_FACTORY_HH_
#define AMANZI_OPERATOR_PDE_ELASTICITY_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "PDE_Elasticity.hh"
#include "PDE_ElasticityFracturedMatrix.hh"

namespace Amanzi {
namespace Operators {

class PDE_ElasticityFactory {
 public:
  PDE_ElasticityFactory() {};
  ~PDE_ElasticityFactory() {};

  Teuchos::RCP<PDE_Elasticity> Create(Teuchos::ParameterList& oplist,
                                      const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  {
    bool fractured_matrix = oplist.sublist("schema").isParameter("fracture");

    if (fractured_matrix) {
      auto op = Teuchos::rcp(new PDE_ElasticityFracturedMatrix(oplist, mesh));
      op->Init(oplist);
      return op;

    } else {
      auto op = Teuchos::rcp(new PDE_Elasticity(oplist, mesh));
      op->Init(oplist);
      return op;
    }
  }
};

} // namespace Operators
} // namespace Amanzi

#endif
