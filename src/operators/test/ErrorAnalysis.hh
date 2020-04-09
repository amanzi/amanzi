/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Miscaleneous tools for error analysis.
*/

#ifndef AMANZI_OPERATOR_ERROR_ANALYSIS_HH_
#define AMANZI_OPERATOR_ERROR_ANALYSIS_HH_

#include <cmath>

#include "Epetra_MultiVector.h"

#include "Mesh.hh"

inline
void ComputeGradError(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                      Epetra_MultiVector& grad, Epetra_MultiVector& grad_exact,
                      double& err_int, double& err_glb, double& gnorm)
{
  int dim = mesh->space_dimension();
  int ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh->num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED);
  std::vector<int> flag(ncells_owned, 0);

  Amanzi::AmanziMesh::Entity_ID_List cells;
  double err_bnd(0.0);

  err_bnd = 0.0;
  for (int f = 0; f < nfaces_owned; ++f) {
    mesh->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);
    int c = cells[0];
    if (cells.size() == 1 && flag[c] == 0) {
      for (int i = 0; i < dim; ++i) {
        double tmp = grad[i][c] - grad_exact[i][c];
        err_bnd += tmp * tmp * mesh->cell_volume(c);
      }
      flag[c] = 1;
    }
  }

  gnorm = 0.0;
  err_glb = 0.0;
  for (int c = 0; c < ncells_owned; ++c) {
    double volume = mesh->cell_volume(c);
    for (int i = 0; i < dim; ++i) {
      double tmp = grad[i][c] - grad_exact[i][c];
      err_glb += tmp * tmp * volume;
      gnorm += grad_exact[i][c] * grad_exact[i][c] * volume;
    }
  }
  err_int = err_glb - err_bnd;

#ifdef HAVE_MPI
    double tmp = err_int;
    mesh->get_comm()->SumAll(&tmp, &err_int, 1);
    tmp = err_glb;
    mesh->get_comm()->SumAll(&tmp, &err_glb, 1);
    tmp = gnorm;
    mesh->get_comm()->SumAll(&tmp, &gnorm, 1);
#endif

  err_int = std::pow(err_int / gnorm, 0.5);
  err_glb = std::pow(err_glb / gnorm, 0.5);
  gnorm = std::pow(gnorm, 0.5);
}

#endif

