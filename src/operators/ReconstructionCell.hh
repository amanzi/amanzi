/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  At the moment, we require the input field to have valid values
  in ghost cells. 
*/

#ifndef AMANZI_RECONSTRUCTION_CELL_HH_
#define AMANZI_RECONSTRUCTION_CELL_HH_

#include <vector>

#include "Epetra_IntVector.h"
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

  // save pointer to the already distributed field.
  virtual void Init(Teuchos::ParameterList& plist) override;

  // unlimited gradient
  // -- compute gradient and keep it internally
  virtual void ComputeGradient(const Teuchos::RCP<const Epetra_MultiVector>& field,
                               int component = 0) override {
    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    AmanziMesh::Entity_ID_List ids(ncells_wghost);
    for (int c = 0; c < ncells_wghost; ++c) ids[c] = c;
    ComputeGradient(ids, field, component);
  }

  // -- compute gradient only in specified cells
  void ComputeGradient(const AmanziMesh::Entity_ID_List& ids,
                       const Teuchos::RCP<const Epetra_MultiVector>& field, int component);

  // calculate value of a linear function at point p
  double getValue(int c, const AmanziGeometry::Point& p);
  double getValue(const AmanziGeometry::Point& gradient, int c, const AmanziGeometry::Point& p);

  // access
  Teuchos::RCP<CompositeVector> gradient() { return gradient_; }

 private:
  void PopulateLeastSquareSystem_(AmanziGeometry::Point& centroid,
                                  double field_value,
                                  WhetStone::DenseMatrix& matrix,
                                  WhetStone::DenseVector& rhs);

  // On intersecting manifolds, we extract neighboors living in the same manifold
  // using a smoothness criterion.
  void CellFaceAdjCellsNonManifold_(AmanziMesh::Entity_ID c,
                                    AmanziMesh::Parallel_type ptype,
                                    std::vector<AmanziMesh::Entity_ID>& cells) const;
 private:
  int dim, poly_order_;
  Teuchos::RCP<CompositeVector> gradient_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
