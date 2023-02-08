/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Polynomials conservative reconstrution on cell data.
  Due to conservation, only slope data have to be stored.
*/

#ifndef AMANZI_RECONSTRUCTION_CELL_POLYNOMIAL_HH_
#define AMANZI_RECONSTRUCTION_CELL_POLYNOMIAL_HH_

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

class ReconstructionCellPolynomial : public Reconstruction {
 public:
  ReconstructionCellPolynomial(){};
  ReconstructionCellPolynomial(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : Reconstruction(mesh), poly_(Teuchos::null){};
  ~ReconstructionCellPolynomial(){};

  // save pointer to the already distributed field.
  virtual void Init(Teuchos::ParameterList& plist) override;

  // base class interface
  // -- compute gradient and keep it internally
  virtual void Compute(const Teuchos::RCP<const Epetra_MultiVector>& field,
                       int component = 0,
                       const Teuchos::RCP<const BCs>& bc = Teuchos::null) override
  {
    int ncells_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
    AmanziMesh::Entity_ID_View ids("ids", ncells_wghost);
    for (int c = 0; c < ncells_wghost; ++c) ids[c] = c;
    Compute(ids, field, component, bc);
  }

  // -- calculate value, deviation from mean, and full polynomial
  virtual double getValue(int c, const AmanziGeometry::Point& p) override;
  virtual double getValueSlope(int c, const AmanziGeometry::Point& p) override;
  virtual WhetStone::Polynomial getPolynomial(int c) const override;

  // -- access returns full polynomial, including the mean value
  virtual Teuchos::RCP<CompositeVector> data() override { return poly_; }

  // compute gradient only in specified cells
  void Compute(const AmanziMesh::Entity_ID_View& ids,
               const Teuchos::RCP<const Epetra_MultiVector>& field,
               int component,
               const Teuchos::RCP<const BCs>& bc = Teuchos::null);

 private:
  void PopulateLeastSquareSystem_(WhetStone::DenseVector& coef,
                                  double field_value,
                                  WhetStone::DenseMatrix& matrix,
                                  WhetStone::DenseVector& rhs);

  // On intersecting manifolds, we extract neighboors living in the same manifold
  // using a smoothness criterion.
  void CellAllAdjCells_(AmanziMesh::Entity_ID c,
                        AmanziMesh::Parallel_kind ptype,
                        std::set<AmanziMesh::Entity_ID>& cells) const;

  void CellAllAdjFaces_(AmanziMesh::Entity_ID c,
                        const std::set<AmanziMesh::Entity_ID>& cells,
                        std::set<AmanziMesh::Entity_ID>& faces) const;

 private:
  int d_, degree_;
  Teuchos::RCP<CompositeVector> poly_, ortho_;
  Teuchos::RCP<Epetra_MultiVector> poly_c_, ortho_c_;
};

} // namespace Operators
} // namespace Amanzi

#endif
