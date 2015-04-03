/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_HEAT_CONDUCTION_HH_
#define AMANZI_OPERATOR_HEAT_CONDUCTION_HH_

#include <string>
#include <vector>

// Amanzi
#include "CompositeVector.hh"

// Operators
#include "Analytic03.hh"

namespace Amanzi{

// This class wraps scalar diffusion coefficient Analytic03.
class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh), ana_(mesh) { 
    int dim = mesh_->space_dimension();
    cvs_.SetMesh(mesh_);
    cvs_.SetGhosted(true);
    cvs_.SetComponent("cell", AmanziMesh::CELL, 1);
    cvs_.SetOwned(false);
    cvs_.AddComponent("face", AmanziMesh::FACE, 1);
    cvs_.AddComponent("grad", AmanziMesh::CELL, dim);

    values_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
    derivatives_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
  }
  ~HeatConduction() {};

  // main members
  void UpdateValues(const CompositeVector& u) { 
    Epetra_MultiVector& vcell = *values_->ViewComponent("cell", true); 
    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      const WhetStone::Tensor& Kc = ana_.Tensor(xc, 0.0);
      vcell[0][c] = Kc(0, 0);
    }

    // add gradient component
    int dim = mesh_->space_dimension();
    Epetra_MultiVector& vgrad = *values_->ViewComponent("grad", true); 

    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      AmanziGeometry::Point grad = ana_.ScalarTensorGradient(xc, 0.0);
      for (int i = 0; i < dim; i++) vgrad[i][c] = grad[i];
    }

    derivatives_->PutScalar(1.0);
  }

  void UpdateValuesPostUpwind() { 
    if (!values_->HasComponent("twin")) {
      cvs_.AddComponent("twin", AmanziMesh::FACE, 1);
      Teuchos::RCP<CompositeVector> tmp = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));

      *tmp->ViewComponent("cell") = *values_->ViewComponent("cell"); 
      *tmp->ViewComponent("face") = *values_->ViewComponent("face"); 
      *tmp->ViewComponent("grad") = *values_->ViewComponent("grad"); 
      values_ = tmp;
    }

    AmanziMesh::Entity_ID_List cells;
    Epetra_MultiVector& vcell = *values_->ViewComponent("cell", true); 
    Epetra_MultiVector& vface = *values_->ViewComponent("face", true); 
    Epetra_MultiVector& vtwin = *values_->ViewComponent("twin", true); 

    vtwin = vface;
    int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

    for (int f = 0; f < nfaces; f++) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();
      
      if (ncells == 2) {
        double v1 = vcell[0][cells[0]];
        double v2 = vcell[0][cells[1]];
        if (fabs(v1 - v2) > 2 * std::min(fabs(v1), fabs(v2))) {
          vface[0][f] = v1;
          vtwin[0][f] = v2;
        }  
      } 
    }
  }

  double Conduction(int c, double T) const {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana_.Tensor(xc, 0.0);
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

typedef double(HeatConduction::*ModelUpwindFn)(int c, double T) const; 

}  // namespace Amanzi

#endif
