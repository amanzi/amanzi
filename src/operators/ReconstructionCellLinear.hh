/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Linear conservative reconstrution on cell data.
  Due to conservation, only gradient is needed to be stored.

  NOTE: At the moment, we require the input field to have valid
  values in ghost cells.
*/

#ifndef AMANZI_RECONSTRUCTION_CELL_LINEAR_HH_
#define AMANZI_RECONSTRUCTION_CELL_LINEAR_HH_

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

class ReconstructionCellLinear : public Reconstruction {
 public:
  ReconstructionCellLinear(){};

  ReconstructionCellLinear(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : Reconstruction(mesh), dim(mesh->getSpaceDimension()){};

  ReconstructionCellLinear(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                           Teuchos::RCP<CompositeVector>& gradient)
    : Reconstruction(mesh),
      dim(mesh->getSpaceDimension()),
      gradient_(gradient),
      gradient_c_(gradient->ViewComponent("cell")){};

  ~ReconstructionCellLinear(){};

  // save pointer to the already distributed field.
  virtual void Init(Teuchos::ParameterList& plist) override;

  // base class interface interface
  // -- compute gradient and keep it internally
  virtual void Compute(const Teuchos::RCP<const Epetra_MultiVector>& field,
                       int component = 0,
                       const Teuchos::RCP<const BCs>& bc = Teuchos::null) override
  {
    int ncells_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL);
    AmanziMesh::Entity_ID_List ids("ids", ncells_wghost);
    for (int c = 0; c < ncells_wghost; ++c) ids[c] = c;
    Compute(ids, field, component, bc);
  }

  // -- calculate value, deviation from mean, and full polynomial
  virtual double getValue(int c, const AmanziGeometry::Point& p) override;
  virtual double getValueSlope(int c, const AmanziGeometry::Point& p) override;
  virtual WhetStone::Polynomial getPolynomial(int c) const override;

  // -- access returns gradient, i.e. no mean value
  virtual Teuchos::RCP<CompositeVector> data() override { return gradient_; }

  // compute gradient only in specified cells
  void Compute(const AmanziMesh::Entity_ID_List& ids,
               const Teuchos::RCP<const Epetra_MultiVector>& field,
               int component,
               const Teuchos::RCP<const BCs>& bc = Teuchos::null);

 private:
  void PopulateLeastSquareSystem_(AmanziGeometry::Point& centroid,
                                  double field_value,
                                  WhetStone::DenseMatrix& matrix,
                                  WhetStone::DenseVector& rhs);

  // On intersecting manifolds, we extract neighboors living in the same manifold
  // using a smoothness criterion.
  void CellFaceAdjCellsManifold_(AmanziMesh::Entity_ID c,
                                 AmanziMesh::Parallel_type ptype,
                                 std::vector<AmanziMesh::Entity_ID>& cells) const;

 private:
  int dim;
  Teuchos::RCP<CompositeVector> gradient_;
  Teuchos::RCP<Epetra_MultiVector> gradient_c_;
};

} // namespace Operators
} // namespace Amanzi

#endif
