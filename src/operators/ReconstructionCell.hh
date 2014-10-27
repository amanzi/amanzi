/*
  This is the operator component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_RECONSTRUCTION_CELL_HH_
#define AMANZI_RECONSTRUCTION_CELL_HH_

#include <vector>

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "Reconstruction.hh"


namespace Amanzi {
namespace Operators {

class ReconstructionCell : public Reconstruction {  
 public:
  ReconstructionCell() {};
  ReconstructionCell(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : Reconstruction(mesh) {};
  ~ReconstructionCell() {};

  // main members for base class
  void Init(Teuchos::RCP<const Epetra_MultiVector> field);
  void Compute();
  void ApplyLimiter(Teuchos::RCP<Epetra_MultiVector> limiter);
  double getValue(int cell, const AmanziGeometry::Point& p);

  // main members
  double getValue(AmanziGeometry::Point& gradient, int cell, const AmanziGeometry::Point& p);

  // access
  Teuchos::RCP<CompositeVector> gradient() { return gradient_; }
 
 private:
  void PopulateLeastSquareSystem(AmanziGeometry::Point& centroid,
                                 double field_value,
                                 WhetStone::DenseMatrix& matrix,
                                 WhetStone::DenseVector& rhs);

  void PrintLeastSquareSystem(WhetStone::DenseMatrix& matrix,
                              WhetStone::DenseVector& rhs);

 private:
  int ncells_owned, nfaces_owned;
  int dim;
  Teuchos::RCP<CompositeVector> gradient_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
