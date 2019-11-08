/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_WITH_BND_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_WITH_BND_WITH_GRAVITY_HH_

#include <strings.h>

// TPLs
#include "Ifpack.h" 
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

// Operators
#include "PDE_DiffusionNLFVwithBndFaces.hh"
#include "PDE_DiffusionWithGravity.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionNLFVwithBndFacesGravity : public PDE_DiffusionNLFVwithBndFaces {
 public:
  PDE_DiffusionNLFVwithBndFacesGravity(Teuchos::ParameterList& plist,
                                       const Teuchos::RCP<Operator>& global_op) :
      PDE_DiffusionNLFVwithBndFaces(plist, global_op)
  {}  

  PDE_DiffusionNLFVwithBndFacesGravity(Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_DiffusionNLFVwithBndFaces(plist, mesh)
  {}

  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;

  // -- after solving the problem: postrocessing
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

};

}  // namespace Operators
}  // namespace Amanzi

#endif
