// PDE_DiffusionMFD: elliptic operators using the MFD family of discretizations.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_MFD_XMOF_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_MFD_XMOF_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "exceptions.hh"
#include "Tensor.hh"
#include "Point.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"

#include "BCs.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"
#include "xmof2D.h"

namespace Amanzi {
namespace Operators {

class PDE_DiffusionMFD_XMOF : public virtual PDE_DiffusionMFD {
 public:
 PDE_DiffusionMFD_XMOF(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<Operator>& global_op) :
   PDE_Diffusion(global_op),
   PDE_DiffusionMFD(plist,global_op)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_XMOF;
    InitDiffusion_(plist);
  }

  PDE_DiffusionMFD_XMOF(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    PDE_Diffusion(mesh),
    PDE_DiffusionMFD(plist, mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_XMOF;
    InitDiffusion_(plist);
  }

  PDE_DiffusionMFD_XMOF(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
    PDE_Diffusion(mesh),
    PDE_DiffusionMFD(plist, mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_XMOF;
    InitDiffusion_(plist);
  }

  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;
 

  void ConstructMiniMesh(const Teuchos::Ptr<const CompositeVector>& vol_fraction,
                         const Teuchos::Ptr<const CompositeVector>& centroids); 

 protected:
  void InitDiffusion_(Teuchos::ParameterList& plist);

  void ConstructBaseMesh_();
  void CreateMassMatrices_();


 protected:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<std::vector<XMOF2D::MeshConfig> > mesh_cfg_;
  Teuchos::RCP<std::vector<XMOF2D::CellsMatData> > mat_data_;
  Teuchos::RCP<std::vector<XMOF2D::XMOF_Reconstructor> > xmof_ir_;
  XMOF2D::IRTolerances ir_tolerances_;

  Teuchos::RCP<const CompositeVector>  vol_frac_, centroids_, sc_mat_coef_;
  std::vector<std::vector<WhetStone::DenseMatrix> >  Wff_cells_mm_;

};


}  // namespace Operators
}  // namespace Amanzi

#endif
