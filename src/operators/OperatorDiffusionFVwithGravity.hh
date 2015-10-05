/*
  This is the opeartors component of the Amanzi code.  

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_FV_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_DIFFUSION_FV_WITH_GRAVITY_HH_

#include <strings.h>

#include "Ifpack.h" 

#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"
#include "OperatorDiffusionFV.hh"


namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionFVwithGravity : public OperatorDiffusionFV {
 public:
  OperatorDiffusionFVwithGravity(Teuchos::ParameterList& plist,
                                 const Teuchos::RCP<Operator>& global_op,
                                 double rho, const AmanziGeometry::Point& g) :
      OperatorDiffusionFV(plist, global_op),
      rho_(rho),
      g_(g)
  {
    operator_type_ = OPERATOR_DIFFUSION_FV_GRAVITY;
    InitDiffusion_(plist);
  }

  OperatorDiffusionFVwithGravity(Teuchos::ParameterList& plist,
                                 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                 double rho, const AmanziGeometry::Point& g) :
      OperatorDiffusionFV(plist, mesh),
      rho_(rho),
      g_(g)
  {
    operator_type_ = OPERATOR_DIFFUSION_FV_GRAVITY;
    InitDiffusion_(plist);
  }

  OperatorDiffusionFVwithGravity(Teuchos::ParameterList& plist,
                                 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                 Teuchos::RCP<const CompositeVector> rho_cv,
                                 const AmanziGeometry::Point& g) :
      OperatorDiffusionFV(plist, mesh),
      rho_cv_(rho_cv),
      g_(g)
  {
    rho_ = 1.0;  // FIXME: vector is not supported
    operator_type_ = OPERATOR_DIFFUSION_FV_GRAVITY;
    InitDiffusion_(plist);
  }

  // main virtual members
  // -- setup
  using OperatorDiffusionFV::Setup;
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);
  virtual void Setup(const Teuchos::RCP<const CompositeVector>& k,
                     const Teuchos::RCP<const CompositeVector>& dkdp);

  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             double rho, const AmanziGeometry::Point& g) {
    rho_ = rho;
    g_ = g;
    Setup(K);
    Setup(k, dkdp);
  } 

  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             const Teuchos::RCP<const CompositeVector>& rho,
             const AmanziGeometry::Point& g) {
    rho_cv_ = rho;
    g_ = g;
    Setup(K);
    Setup(k, dkdp);
  }

  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
          const Teuchos::Ptr<const CompositeVector>& u);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);

  // -- modify an operator
  virtual void ApplyBCs(bool primary, bool eliminate);
  virtual void ModifyMatrices(const CompositeVector& u) {};
  virtual void ScaleMassMatrices(double s) {};

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeGravityFlux(int f) const;

  // access
  const CompositeVector& gravity_terms() { return *gravity_term_; }

 protected:
  virtual void ComputeJacobianLocal_(
      int mcells, int f, int face_dir, int Krel_method,
      int bc_model, double bc_value,
      double *pres, double *dkdp_cell,
      WhetStone::DenseMatrix& Jpp);

  void ComputeTransmissibility_(AmanziGeometry::Point* g, Teuchos::RCP<CompositeVector> g_cv);

  virtual void InitDiffusion_(Teuchos::ParameterList& plist);

 protected:
  double rho_;
  Teuchos::RCP<const CompositeVector> rho_cv_;

  AmanziGeometry::Point g_;
  Teuchos::RCP<CompositeVector> gravity_term_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif
