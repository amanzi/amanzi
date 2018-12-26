/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_FRACTURED_MATRIX_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_FRACTURED_MATRIX_HH_

#include "Teuchos_RCP.hpp"

#include "PDE_DiffusionMFD.hh"

namespace Amanzi {
namespace Operators {

class PDE_DiffusionFracturedMatrix : public PDE_DiffusionMFD {
 public:
  PDE_DiffusionFracturedMatrix(Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_Diffusion(mesh),
      PDE_DiffusionMFD(plist, mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_FRACTURED_MATRIX;
  }

  // main interface members
  virtual void Init(Teuchos::ParameterList& plist) override;
};

}  // namespace Operators
}  // namespace Amanzi


#endif
