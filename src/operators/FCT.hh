/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Kontantin Lipnikov (lipnikov@lanl.gov)

  Flux Corrected Transport.
*/

#ifndef FCT_HH_
#define FCT_HH_

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

#include "CompositeVector.hh"
#include "LimiterCell.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace Operators {

class FCT {
 public:
  FCT(Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
      Teuchos::RCP<const AmanziMesh::Mesh> mesh1,
      Teuchos::RCP<const LimiterCell> limiter,
      Teuchos::RCP<const Epetra_MultiVector> field)
      : mesh0_(mesh0),
        mesh1_(mesh1),
        limiter_(limiter),
        field_(field),
        component_(0),
        weight0_(Teuchos::null),
        weight1_(Teuchos::null) {};
  ~FCT() {};

  // re-initialization for another pair field + component 
  void Init(Teuchos::RCP<const Epetra_MultiVector> field, int component,
            Teuchos::RCP<const Epetra_MultiVector> weight0 = Teuchos::null,
            Teuchos::RCP<const Epetra_MultiVector> weight1 = Teuchos::null) {
    field_ = field;
    component_ = component;
    weight0_ = weight0;
    weight1_ = weight1;
  }

  // flux is from 1st to 2nd cell in the list face->cells 
  void Compute(const CompositeVector& flux_lo,
               const CompositeVector& flux_ho, 
               const BCs& bc,
               CompositeVector& flux);

  // access
  double get_alpha_mean() { return alpha_mean_; }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_, mesh1_;
  Teuchos::RCP<const LimiterCell> limiter_;

  Teuchos::RCP<const Epetra_MultiVector> field_;
  int component_;

  Teuchos::RCP<const Epetra_MultiVector> weight0_, weight1_;

  double alpha_mean_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


