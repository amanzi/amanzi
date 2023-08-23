/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Miscaleneous tools for error analysis.
*/

#ifndef AMANZI_OPERATOR_ERROR_ANALYSIS_HH_
#define AMANZI_OPERATOR_ERROR_ANALYSIS_HH_

#include <cmath>

#include "Epetra_MultiVector.h"

#include "Mesh.hh"

inline void
ComputePolyError(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                 Epetra_MultiVector& poly,
                 Epetra_MultiVector& poly_exact,
                 double& err_int,
                 double& err_glb,
                 double& gnorm)
{
  int npoly = poly.NumVectors();
  int ncells_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                          Amanzi::AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE,
                                          Amanzi::AmanziMesh::Parallel_type::OWNED);
  std::vector<int> flag(ncells_owned, 0);

  Amanzi::AmanziMesh::Entity_ID_List cells;
  double err_bnd(0.0);

  err_bnd = 0.0;
  for (int f = 0; f < nfaces_owned; ++f) {
    cells = mesh->getFaceCells(f, Amanzi::AmanziMesh::Parallel_type::ALL);
    int c = cells[0];
    if (cells.size() == 1 && flag[c] == 0) {
      for (int i = 0; i < npoly; ++i) {
        double tmp = poly[i][c] - poly_exact[i][c];
        err_bnd += tmp * tmp * mesh->getCellVolume(c);
      }
      flag[c] = 1;
    }
  }

  gnorm = 0.0;
  err_glb = 0.0;
  for (int c = 0; c < ncells_owned; ++c) {
    double volume = mesh->getCellVolume(c);
    for (int i = 0; i < npoly; ++i) {
      double tmp = poly[i][c] - poly_exact[i][c];
      err_glb += tmp * tmp * volume;
      gnorm += poly_exact[i][c] * poly_exact[i][c] * volume;
    }
  }
  err_int = std::abs(err_glb - err_bnd);

#ifdef HAVE_MPI
  double tmp = err_int;
  mesh->getComm()->SumAll(&tmp, &err_int, 1);
  tmp = err_glb;
  mesh->getComm()->SumAll(&tmp, &err_glb, 1);
  tmp = gnorm;
  mesh->getComm()->SumAll(&tmp, &gnorm, 1);
#endif

  err_int = std::pow(err_int / gnorm, 0.5);
  err_glb = std::pow(err_glb / gnorm, 0.5);
  gnorm = std::pow(gnorm, 0.5);
}

#endif
