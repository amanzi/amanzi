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
  FCT(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
      Teuchos::RCP<const LimiterCell> lifting,
      Teuchos::RCP<const Epetra_MultiVector> field)
      : mesh_(mesh),
        lifting_(lifting),
        field_(field) {};
  ~FCT() {};

  void Compute(const CompositeVector& flux_lo,
               const CompositeVector& flux_ho, CompositeVector& flux);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const LimiterCell> lifting_;
  Teuchos::RCP<const Epetra_MultiVector> field_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


