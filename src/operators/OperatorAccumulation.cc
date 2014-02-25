/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "OperatorAccumulation.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Adds time derivative ss * (u - u0) / dT.
****************************************************************** */
void OperatorAccumulation::UpdateMatrices(
    const CompositeVector& u0, const CompositeVector& ss, double dT)
{
  const Epetra_MultiVector& u0c = *u0.ViewComponent("cell");
  const Epetra_MultiVector& ssc = *u0.ViewComponent("cell");

  Epetra_MultiVector& diag = *diagonal_block_->ViewComponent("cell");
  Epetra_MultiVector& rhs = *rhs_->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume * ssc[0][c] / dT;
    diag[0][c] += factor;
    rhs[0][c] += factor * u0c[0][c];
  }
}

}  // namespace Operators
}  // namespace Amanzi
