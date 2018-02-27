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
#include "PDE_Abstract.hh"
#include "PDE_DiffusionMFD.hh"
#include "xmof2D.h"

namespace Amanzi {
namespace Operators {

class PDE_DiffusionMFD_XMOF : public virtual PDE_Abstract {
 public:
 PDE_DiffusionMFD_XMOF(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<Operator>& global_op) :
   // PDE_Diffusion(global_op),
   // PDE_DiffusionMFD(plist,global_op)
   PDE_Abstract(plist, global_op),  factor_(1.0)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_XMOF;
    InitDiffusion_(plist);
  }

  PDE_DiffusionMFD_XMOF(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    // PDE_Diffusion(mesh),
    // PDE_DiffusionMFD(plist, mesh)
    PDE_Abstract(plist, mesh), factor_(1.0)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_XMOF;
    InitDiffusion_(plist);
  }

  // PDE_DiffusionMFD_XMOF(Teuchos::ParameterList& plist,
  //                  const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
  //   // PDE_Diffusion(mesh),
  //   // PDE_DiffusionMFD(plist, mesh)
  //   PDE_Abstract(plist, mesh)
  // {
  //   operator_type_ = OPERATOR_DIFFUSION_MFD_XMOF;
  //   InitDiffusion_(plist);
  // }

  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u, 
                              double dt);

  virtual void ApplyBCs(bool primary, bool eliminate);

  virtual void UpdateFlux(const Teuchos::Ptr<CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u, 
                          double dt);

  void ConstructMiniMesh(const Teuchos::Ptr<const CompositeVector>& vol_fraction,
                         const Teuchos::Ptr<const CompositeVector>& centroids); 

  void SetupMultiMatK(const Teuchos::RCP<std::vector<std::vector<WhetStone::Tensor> > > KMultiMat){KMultiMat_ = KMultiMat;}

  void WriteSpecialGMV(std::string filename, const Epetra_MultiVector& vol_frac_vec, const Epetra_MultiVector& solution);

 protected:
  void InitDiffusion_(Teuchos::ParameterList& plist);
  
  void ConstructBaseMesh_();
  void CreateMassMatrices_();
  void ApplyBCs_Mixed_(const Teuchos::Ptr<const BCs>& bc_trial,
                       const Teuchos::Ptr<const BCs>& bc_test,
                       bool primary, bool eliminate);


 protected:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<std::vector<XMOF2D::MeshConfig> > mesh_cfg_;
  Teuchos::RCP<std::vector<XMOF2D::CellsMatData> > mat_data_;
  Teuchos::RCP<std::vector<XMOF2D::XMOF_Reconstructor> > xmof_ir_;
  XMOF2D::IRTolerances ir_tolerances_;

  Teuchos::RCP<std::vector<std::vector<WhetStone::Tensor> > > KMultiMat_; 
  //Teuchos::RCP<CompositeVector>  vol_frac_, centroids_;
  Teuchos::RCP<Epetra_Vector> rhs_cells_;
  
  std::vector<std::vector<WhetStone::DenseMatrix> >  Wff_cells_mm_;
  std::vector<WhetStone::DenseMatrix> Abb_, Pb_, Pp_, T3_;
  
  double factor_;

};


}  // namespace Operators
}  // namespace Amanzi

#endif
