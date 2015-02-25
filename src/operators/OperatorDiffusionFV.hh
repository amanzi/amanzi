/*
  This is the opeartors component of the Amanzi code.  

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_FV_HH_
#define AMANZI_OPERATOR_DIFFUSION_FV_HH_

#include <strings.h>

#include "Ifpack.h" 

#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"
#include "OperatorDiffusion.hh"


namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionFV : public OperatorDiffusion {
 public:

  OperatorDiffusionFV(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<Operator>& global_op) :
      OperatorDiffusion(global_op),
      gravity_(false)      
  {
    InitDiffusion_(plist);
  }

  OperatorDiffusionFV(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      gravity_(false)
  {
    InitDiffusion_(plist);
  }

  OperatorDiffusionFV(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      gravity_(false)
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
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);
  virtual void ApplyBCs(bool primary=true);
  virtual void ModifyMatrices(const CompositeVector& u) {}

  virtual void SetGravity(const AmanziGeometry::Point& g) {
    g_ = g;
    gravity_ = true;
  }
  virtual void SetVectorDensity(const Teuchos::RCP<const CompositeVector>& rho) { rho_cv_ = rho; }
  
  template <class Model> 
  double DeriveBoundaryFaceValue(int f, const CompositeVector& u, const Model& model);

  //access function
  const Epetra_Vector transmissibilities() { return *transmissibility_; }
  const Epetra_Vector gravity_terms() { return *gravity_term_; }

 protected:
  void ComputeTransmissibilities_();

  void AnalyticJacobian_(const CompositeVector& solution);

  void ComputeJacobianLocal_(
      int mcells, int f, int face_dir, int Krel_method,
      int bc_model, double bc_value,
      double *pres, double *dkdp_cell,
      WhetStone::DenseMatrix& Jpp);

  void InitDiffusion_(Teuchos::ParameterList& plist);

  
 protected:
  AmanziGeometry::Point g_;
  Teuchos::RCP<Epetra_Vector> gravity_term_;
  bool gravity_;

  Teuchos::RCP<Epetra_Vector> transmissibility_;

  int newton_correction_;
  bool exclude_primary_terms_;
  
};

}  // namespace Operators
}  // namespace Amanzi

// Description of templated function DeriveBoundaryFaceValue(f, u, model)
#include "FluxTPFABCfunc.hh"

#endif
