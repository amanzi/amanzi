/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_WITH_GRAVITY_HH_

#include "PDE_Diffusion.hh"

namespace Amanzi {
namespace Operators {

class PDE_DiffusionWithGravity : public virtual PDE_Diffusion {
 public:
  PDE_DiffusionWithGravity(const Teuchos::RCP<Operator>& global_op)
    : PDE_Diffusion(global_op), is_scalar_(false){};

  PDE_DiffusionWithGravity(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(mesh), is_scalar_(false){};

  PDE_DiffusionWithGravity(const Teuchos::RCP<AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(mesh), is_scalar_(false){};

  virtual ~PDE_DiffusionWithGravity() = default;

  virtual void SetGravity(const AmanziGeometry::Point& g) { g_ = g; }

  virtual void SetDensity(double rho)
  {
    is_scalar_ = true;
    rho_ = rho;
  }

  virtual void SetDensity(const Teuchos::RCP<const CompositeVector>& rho)
  {
    is_scalar_ = false;
    if (rho->HasComponent("cell")) { rho_cv_ = rho; }
  }

  double GetDensity(int c)
  {
    if (is_scalar_) return rho_;
    return (*rho_cv_->ViewComponent("cell", true))[0][c];
  }

 protected:
  bool is_scalar_;
  double rho_;
  Teuchos::RCP<const CompositeVector> rho_cv_;
  AmanziGeometry::Point g_;
};

} // namespace Operators
} // namespace Amanzi

#endif
