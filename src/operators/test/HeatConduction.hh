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

#ifndef AMANZI_OPERATOR_HEAT_CONDUCTION_HH_
#define AMANZI_OPERATOR_HEAT_CONDUCTION_HH_

#include <string>
#include <vector>

// Amanzi
#include "Mesh_Algorithms.hh"
#include "CompositeVector.hh"

// Operators
#include "OperatorDefs.hh"
#include "Analytic03.hh"

namespace Amanzi {

class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh), ana_(mesh)
  {
    int dim = mesh_->getSpaceDimension();
    cvs_.SetMesh(mesh_);
    cvs_.SetGhosted(true);
    cvs_.AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    cvs_.AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    cvs_.AddComponent("grad", AmanziMesh::Entity_kind::CELL, dim);

    values_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
    derivatives_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
  }
  ~HeatConduction(){};

  // main members
  void UpdateValues(const CompositeVector& u,
                    const std::vector<int>& bc_model,
                    const std::vector<double>& bc_value)
  {
    Epetra_MultiVector& vcell = *values_->ViewComponent("cell", true);
    int ncells =
      mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL);

    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
      const WhetStone::Tensor& Kc = ana_.TensorDiffusivity(xc, 0.0);
      vcell[0][c] = Kc(0, 0);
    }

    // add boundary face component
    Epetra_MultiVector& values_f = *values_->ViewComponent("face", true);
    for (int f = 0; f != bc_model.size(); ++f) {
      if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);
        values_f[0][f] = Conduction(c, bc_value[f]);
      }
    }

    // add gradient component
    int dim = mesh_->getSpaceDimension();
    Epetra_MultiVector& vgrad = *values_->ViewComponent("grad", true);

    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
      AmanziGeometry::Point grad = ana_.ScalarTensorGradient(xc, 0.0);
      for (int i = 0; i < dim; i++) vgrad[i][c] = grad[i];
    }

    derivatives_->PutScalar(1.0);
  }

  // adds twin-component and over-writes face-components on discontinuity
  void UpdateValuesPostUpwind()
  {
    int dim = mesh_->getSpaceDimension();

    if (!values_->HasComponent("twin")) {
      cvs_.AddComponent("twin", AmanziMesh::Entity_kind::FACE, 1);
      Teuchos::RCP<CompositeVector> tmp =
        Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));

      *tmp->ViewComponent("cell") = *values_->ViewComponent("cell");
      *tmp->ViewComponent("face") = *values_->ViewComponent("face");
      *tmp->ViewComponent("grad") = *values_->ViewComponent("grad");
      values_ = tmp;
    }

    AmanziMesh::Entity_ID_List cells;
    Epetra_MultiVector& vcell = *values_->ViewComponent("cell", true);
    Epetra_MultiVector& vface = *values_->ViewComponent("face", true);
    Epetra_MultiVector& vtwin = *values_->ViewComponent("twin", true);
    Epetra_MultiVector& vgrad = *values_->ViewComponent("grad", true);

    vtwin = vface;
    int nfaces =
      mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL);

    for (int f = 0; f < nfaces; f++) {
      cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_type::ALL);
      int ncells = cells.size();

      if (ncells == 2) {
        int c1 = cells[0], c2 = cells[1];
        double v1 = vcell[0][c1];
        double v2 = vcell[0][c2];
        if (fabs(v1 - v2) > 2 * std::min(fabs(v1), fabs(v2))) {
          vface[0][f] = v1;
          vtwin[0][f] = v2;

          const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
          const AmanziGeometry::Point& xc1 = mesh_->getCellCentroid(c1);
          const AmanziGeometry::Point& xc2 = mesh_->getCellCentroid(c2);

          for (int i = 0; i < dim; ++i) {
            vface[0][f] += vgrad[i][c1] * (xf[i] - xc1[i]);
            vtwin[0][f] += vgrad[i][c2] * (xf[i] - xc2[i]);
          }
        }
      }
    }
  }

  // adds twin-component and over-writes face-components
  void UpdateValuesFaceTwin()
  {
    int dim = mesh_->getSpaceDimension();

    if (!values_->HasComponent("twin")) {
      cvs_.AddComponent("twin", AmanziMesh::Entity_kind::FACE, 1);
      Teuchos::RCP<CompositeVector> tmp =
        Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));

      *tmp->ViewComponent("cell") = *values_->ViewComponent("cell");
      *tmp->ViewComponent("face") = *values_->ViewComponent("face");
      *tmp->ViewComponent("grad") = *values_->ViewComponent("grad");
      values_ = tmp;
    }

    AmanziMesh::Entity_ID_List cells;
    Epetra_MultiVector& vcell = *values_->ViewComponent("cell", true);
    Epetra_MultiVector& vface = *values_->ViewComponent("face", true);
    Epetra_MultiVector& vgrad = *values_->ViewComponent("grad", true);
    Epetra_MultiVector& vtwin = *values_->ViewComponent("twin", true);

    int nfaces =
      mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL);
    AmanziGeometry::Point grad(dim), xc(dim);

    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);

      cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_type::ALL);
      int ncells = cells.size();

      int c = cells[0];
      for (int i = 0; i < dim; ++i) grad[i] = vgrad[i][c];
      xc = mesh_->getCellCentroid(c);
      vface[0][f] = vcell[0][c] + grad * (xf - xc);

      if (ncells == 2) {
        c = cells[1];
        for (int i = 0; i < dim; ++i) grad[i] = vgrad[i][c];
        xc = mesh_->getCellCentroid(c);
        vtwin[0][f] = vcell[0][c] + grad * (xf - xc);
      }
    }
  }

  double Conduction(int c, double T) const
  {
    const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana_.TensorDiffusivity(xc, 0.0);
    return Kc(0, 0);
  }

  Teuchos::RCP<CompositeVector> values() { return values_; }
  Teuchos::RCP<CompositeVector> derivatives() { return derivatives_; }

 private:
  CompositeVectorSpace cvs_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<CompositeVector> values_, derivatives_;
  mutable Analytic03 ana_;
};

} // namespace Amanzi

#endif
