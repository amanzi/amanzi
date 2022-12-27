/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_RECONSTRUCTION_HH_
#define AMANZI_RECONSTRUCTION_HH_

#include <vector>

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "CompositeVector.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "Point.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace Operators {

class Reconstruction {
 public:
  Reconstruction(){};
  Reconstruction(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : mesh_(mesh), field_(Teuchos::null), component_(0), weight_(WeightType::WT_CONSTANT){};
  virtual ~Reconstruction() = default;

  // main members
  virtual void Init(Teuchos::ParameterList& plist) = 0;

  virtual void Compute(const Teuchos::RCP<const Epetra_MultiVector>& field,
                       int component = 0,
                       const Teuchos::RCP<const BCs>& bc = Teuchos::null)
  {
    field_ = field;
    component_ = component;
  }

  virtual double getValue(int c, const AmanziGeometry::Point& p) = 0;
  virtual double getValueSlope(int c, const AmanziGeometry::Point& p) = 0;
  virtual WhetStone::Polynomial getPolynomial(int c) const = 0;

  // access function returns slope (gradient and higher-order derivatives)
  virtual Teuchos::RCP<CompositeVector> data() = 0;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const Epetra_MultiVector> field_;
  int component_;
  WeightType weight_;
};

} // namespace Operators
} // namespace Amanzi

#endif
