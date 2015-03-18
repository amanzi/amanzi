/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_MFD_HH_
#define AMANZI_OPERATOR_DIFFUSION_MFD_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "tensor.hh"
#include "Point.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"

#include "BCs.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusion.hh"

namespace Amanzi {
namespace Operators {

class OperatorDiffusionMFD : public OperatorDiffusion {
 public:
  OperatorDiffusionMFD(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<Operator>& global_op) :
      OperatorDiffusion(global_op),
      factor_(1.0)
  {
    InitDiffusion_(plist);
  }

  OperatorDiffusionMFD(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      factor_(1.0)
  {
    InitDiffusion_(plist);
  }

  OperatorDiffusionMFD(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      factor_(1.0)
  {
    InitDiffusion_(plist);
  }
  
  // main virtual members
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                     double rho, double mu);
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                     const Teuchos::RCP<const CompositeVector>& rho,
                     const Teuchos::RCP<const CompositeVector>& mu);
  virtual void Setup(const Teuchos::RCP<const CompositeVector>& k,
                     const Teuchos::RCP<const CompositeVector>& dkdp);
  using OperatorDiffusion::Setup;

  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
          const Teuchos::Ptr<const CompositeVector>& u);
  virtual void UpdateMatricesNewtonCorrection(
          const Teuchos::Ptr<const CompositeVector>& flux,
          const Teuchos::Ptr<const CompositeVector>& u);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);
  virtual void ApplyBCs(bool primary = true);
  virtual void ModifyMatrices(const CompositeVector& u);

  int nfailed_primary() { return nfailed_primary_; }
  void set_factor(double factor) { factor_ = factor; }
  
 protected:
  void InitDiffusion_(Teuchos::ParameterList& plist);
  void CreateMassMatrices_();

  void UpdateMatricesNodal_();
  void UpdateMatricesTPFA_();
  void UpdateMatricesMixed_(const Teuchos::Ptr<const CompositeVector>& flux);
  void UpdateMatricesMixedWithGrad_(const Teuchos::Ptr<const CompositeVector>& flux);

  void AddNewtonCorrectionCell_(const Teuchos::Ptr<const CompositeVector>& flux,
                                const Teuchos::Ptr<const CompositeVector>& u);

 protected:
  std::vector<WhetStone::DenseMatrix> Wff_cells_;
  int newton_correction_;
  bool exclude_primary_terms_;
  double factor_;
  bool scaled_constraint_;

  int mfd_primary_, mfd_secondary_, mfd_pc_primary_, mfd_pc_secondary_;
  int nfailed_primary_;

};

}  // namespace Operators
}  // namespace Amanzi


#endif


