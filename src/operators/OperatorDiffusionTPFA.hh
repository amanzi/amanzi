/*
  This is the opeartors component of the Amanzi code.  

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_TPFA_HH_
#define AMANZI_OPERATOR_DIFFUSION_TPFA_HH_

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

class OperatorDiffusionTPFA : public OperatorDiffusion {
 public:
  OperatorDiffusionTPFA() {};
  OperatorDiffusionTPFA(Teuchos::RCP<const CompositeVectorSpace> cvs, 
                        Teuchos::ParameterList& plist, Teuchos::RCP<BCs> bc) 
      : OperatorDiffusion(cvs, plist, bc) {};
  OperatorDiffusionTPFA(const Operator& op, Teuchos::ParameterList& plist, Teuchos::RCP<BCs> bc) 
      : OperatorDiffusion(op, plist, bc) {};
  ~OperatorDiffusionTPFA() {};

  // re-implementation of basic operator virtual members
  void InitOperator(std::vector<WhetStone::Tensor>& K,
                    Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
                    double rho, double mu);

  void UpdateMatrices(Teuchos::RCP<const CompositeVector> flux, Teuchos::RCP<const CompositeVector> u);
  void ApplyBCs(); 
  void UpdateFlux(const CompositeVector& u, CompositeVector& flux);

  int Apply(const CompositeVector& X, CompositeVector& Y) const;
  int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const;
  void ComputeNegativeResidual(const CompositeVector& v, CompositeVector& r);

  void SetGravity(const AmanziGeometry::Point& g) { g_ = g; }
  void SetUpwind(int upwind_method) { upwind_ = upwind_method; }

  template <class Model> 
  double DeriveBoundaryFaceValue(int f, const CompositeVector& u, const Model& model);

  //access function
  const Epetra_Vector transmissibilities() { return *transmissibility_; }
  const Epetra_Vector gravity_terms() { return *gravity_term_; }

 private:
  void ComputeTransmissibilities_();

  void AnalyticJacobian_(const CompositeVector& solution);

  void ComputeJacobianLocal_(
      int mcells, int f, int face_dir, int Krel_method,
      int bc_model, double bc_value,
      double *pres, double *dkdp_cell,
      WhetStone::DenseMatrix& Jpp);
         
 private:
  AmanziGeometry::Point g_;
  Teuchos::RCP<Epetra_Vector> transmissibility_;
  Teuchos::RCP<Epetra_Vector> gravity_term_;
};

}  // namespace Operators
}  // namespace Amanzi

// Description of templated function DeriveBoundaryFaceValue(f, u, model)
#include "FluxTPFABCfunc.hh"

#endif
