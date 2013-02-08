/*
This is the flow component of the Amanzi code.

 
Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "Darcy_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Estimate dT increase factor by comparing the 1st and 2nd order
* time approximations.
****************************************************************** */
double Darcy_PK::ErrorEstimate(double* dTfactor)
{
  double tol, error, error_max = 0.0;
  double dTfactor_cell;

  *dTfactor = 100.0;
  for (int c = 0; c < ncells_owned; c++) {
    error = fabs((*pdot_cells)[c] - (*pdot_cells_prev)[c]) * dT / 2;
    tol = ti_specs->rtol * fabs((*solution)[c]) + ti_specs->atol;

    dTfactor_cell = sqrt(tol / std::max<double>(error, FLOW_DT_ADAPTIVE_ERROR_TOLERANCE));
    *dTfactor = std::min<double>(*dTfactor, dTfactor_cell);

    error_max = std::max<double>(error_max, error - tol);
  }

  *dTfactor *= FLOW_DT_ADAPTIVE_SAFETY_FACTOR;
  *dTfactor = std::min<double>(*dTfactor, FLOW_DT_ADAPTIVE_INCREASE);
  *dTfactor = std::max<double>(*dTfactor, FLOW_DT_ADAPTIVE_REDUCTION);

#ifdef HAVE_MPI
    double dT_tmp = *dTfactor;
    solution->Comm().MinAll(&dT_tmp, dTfactor, 1);  // find the global minimum
 
    double error_tmp = error_max;
    solution->Comm().MaxAll(&error_tmp, &error_max, 1);  // find the global maximum
#endif
  return error_max;
}

}  // namespace AmanziFlow
}  // namespace Amanzi


