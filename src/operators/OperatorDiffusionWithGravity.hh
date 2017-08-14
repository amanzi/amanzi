/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_DIFFUSION_WITH_GRAVITY_HH_

#include "OperatorDiffusion.hh"

namespace Amanzi {
namespace Operators {

class OperatorDiffusionWithGravity : public virtual OperatorDiffusion {
 public:
  OperatorDiffusionWithGravity(const Teuchos::RCP<Operator>& global_op) :
      OperatorDiffusion(global_op),
      is_scalar_(false) {};

  OperatorDiffusionWithGravity(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      is_scalar_(false) {};

  OperatorDiffusionWithGravity(const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      is_scalar_(false) {};

  virtual ~OperatorDiffusionWithGravity() = default;
  
  virtual void SetGravity(const AmanziGeometry::Point& g) { g_ = g; }
  
  virtual void SetDensity(double rho) {
    is_scalar_ = true;
    rho_ = rho;
  }

  virtual void SetDensity(const Teuchos::RCP<const CompositeVector>& rho) {
    is_scalar_ = false;
    rho_cv_ = rho;
  }

 protected:
  bool is_scalar_;
  double rho_;
  Teuchos::RCP<const CompositeVector> rho_cv_;
  AmanziGeometry::Point g_;
};

}  // namespace Operators
}  // namespace Amanzi
  
#endif 
