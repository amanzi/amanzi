/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_WITH_GRAVITY_HH_

#include <strings.h>

// TPLs
#include "Ifpack.h" 
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

// Operators
#include "PDE_DiffusionNLFV.hh"
#include "PDE_DiffusionWithGravity.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionNLFVwithGravity : public PDE_DiffusionNLFV {
 public:
  PDE_DiffusionNLFVwithGravity(Teuchos::ParameterList& plist,
          const Teuchos::RCP<Operator>& global_op) :
      PDE_DiffusionNLFV(plist, global_op)
  {}

  PDE_DiffusionNLFVwithGravity(Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh):
      PDE_DiffusionNLFV(plist, mesh)
  {}


  virtual void SetDensity(const Teuchos::RCP<const CompositeVector>& rho) override {
    Exceptions::amanzi_throw("PDE_DiffusionNLFVwithGravity does not support vector density.");
  }
  
  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;

  // -- after solving the problem: postrocessing
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeGravityFlux(int f) const override {
    Exceptions::amanzi_throw("PDE_DiffusionNLFVwithGravity::ComputeGravityFlux not implemented.");
    return 0.;
  };

  // virtual members from the base NLFV class
  // -- solution can be modified on boundary faces. This reflects specifics
  //    of nonlinear FV schemes.
  virtual double MapBoundaryValue_(int f, double u) override;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
