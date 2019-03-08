/*
  Flow PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Darcy_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Estimate dT increase factor by comparing the 1st and 2nd order
* time approximations.
****************************************************************** */
double Darcy_PK::ErrorEstimate_(double* dt_factor)
{
  Epetra_MultiVector& p_cell = *solution->ViewComponent("cell");

  double tol, atol(1.0), rtol(1e-5), error, error_max(0.0), p(101325.0);
  double dt_factor_cell;

  *dt_factor = 100.0;
  for (int c = 0; c < ncells_owned; c++) {
    error = fabs((*pdot_cells)[c] - (*pdot_cells_prev)[c]) * dt_ / 2;
    // tol = rtol * fabs(p_cell[0][c]) + atol;
    tol = rtol * p + atol;

    dt_factor_cell = sqrt(tol / std::max(error, FLOW_DT_ADAPTIVE_ERROR_TOLERANCE));
    *dt_factor = std::min(*dt_factor, dt_factor_cell);

    error_max = std::max(error_max, error - tol);
  }

  *dt_factor *= FLOW_DT_ADAPTIVE_SAFETY_FACTOR;
  *dt_factor = std::min(*dt_factor, FLOW_DT_ADAPTIVE_INCREASE);
  *dt_factor = std::max(*dt_factor, FLOW_DT_ADAPTIVE_REDUCTION);

#ifdef HAVE_MPI
    double dt_tmp = *dt_factor;
    solution->Comm()->MinAll(&dt_tmp, dt_factor, 1);  // find the global minimum
 
    double error_tmp = error_max;
    solution->Comm()->MaxAll(&error_tmp, &error_max, 1);  // find the global maximum
#endif

  return error_max;
}

}  // namespace Flow
}  // namespace Amanzi


