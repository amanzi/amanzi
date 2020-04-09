/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Reaction operator. 
*/

#ifndef AMANZI_OPERATOR_PDE_REACTION_HH_
#define AMANZI_OPERATOR_PDE_REACTION_HH_

#include <string>

#include "Epetra_MultiVector.h"

// Amanzi
#include "BilinearForm.hh"
#include "VectorObjects.hh"

// Amanzi::Operators
#include "PDE_HelperDiscretization.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class PDE_Reaction : public PDE_HelperDiscretization {
 public:
  PDE_Reaction(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op) :
      K_(Teuchos::null),
      PDE_HelperDiscretization(global_op),
      coef_type_(CoefType::CONSTANT),
      static_matrices_initialized_(false) {
    InitReaction_(plist);
  }

  PDE_Reaction(Teuchos::ParameterList& plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh) : 
      K_(Teuchos::null),
      PDE_HelperDiscretization(mesh),
      coef_type_(CoefType::CONSTANT),
      static_matrices_initialized_(false) {
    InitReaction_(plist);
  }

  // required members 
  // -- setup
  void Setup(const Teuchos::RCP<Epetra_MultiVector>& K) { K_ = K; Kpoly_ = Teuchos::null; Kpoly_st_ = Teuchos::null; }

  template<typename T>
  void Setup(const Teuchos::RCP<std::vector<T> >& K, bool reset);

  // -- generate a linearized operator 
  virtual void UpdateMatrices();

  // -- new interface for pre-computed data  
  void UpdateMatrices(double t);

  // boundary conditions
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

 private:
  void InitReaction_(Teuchos::ParameterList& plist);
  void CreateStaticMatrices_();

 protected:
  Teuchos::RCP<const Epetra_MultiVector> K_;
  Teuchos::RCP<std::vector<WhetStone::Polynomial> > Kpoly_;
  Teuchos::RCP<std::vector<WhetStone::SpaceTimePolynomial> > Kpoly_st_;

  Teuchos::RCP<WhetStone::BilinearForm> mfd_;

 private:
  Schema global_schema_col_, global_schema_row_;
  Schema local_schema_col_, local_schema_row_;

  CoefType coef_type_;
  bool static_matrices_initialized_;
  std::vector<std::vector<WhetStone::DenseMatrix> > static_matrices_;
};


/* ******************************************************************
* Specialization of Setup
****************************************************************** */
template<>
inline
void PDE_Reaction::Setup<WhetStone::Polynomial>(
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& K, bool reset) {
  K_ = Teuchos::null;
  Kpoly_st_ = Teuchos::null;
  Kpoly_ = K;
  coef_type_ = CoefType::POLYNOMIAL;
}

template<>
inline
void PDE_Reaction::Setup<WhetStone::SpaceTimePolynomial>(
    const Teuchos::RCP<std::vector<WhetStone::SpaceTimePolynomial> >& K, bool reset) {
  K_ = Teuchos::null;
  Kpoly_st_ = K;
  coef_type_ = CoefType::VECTOR_SPACETIME_POLYNOMIAL;
  if (!static_matrices_initialized_ || reset) CreateStaticMatrices_();
}

}  // namespace Operators
}  // namespace Amanzi

#endif

