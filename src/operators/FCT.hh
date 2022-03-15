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
      Teuchos::RCP<const LimiterCell> lifting,
      Teuchos::RCP<const Epetra_MultiVector> field)
      : mesh0_(mesh0),
        mesh1_(mesh1),
        lifting_(lifting),
        field_(field) {};
  ~FCT() {};

  // flux is from 1st to 2nd cell in the list face->cells 
  void Compute(const CompositeVector& flux_lo,
               const CompositeVector& flux_ho, 
               const BCs& bc,
               CompositeVector& flux);

  // access
  double get_alpha_mean() { return alpha_mean_; }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_, mesh1_;
  Teuchos::RCP<const LimiterCell> lifting_;
  Teuchos::RCP<const Epetra_MultiVector> field_;

  double alpha_mean_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


