/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Clip pressure using pressure threshold.
****************************************************************** */
void Richards_PK::ClipHydrostaticPressure(double pmin, Epetra_MultiVector& p)
{
  for (int c = 0; c < ncells_owned; c++) p[0][c] = std::max(p[0][c], pmin);
}


/* ******************************************************************
* Clip pressure using constant saturation.
****************************************************************** */
void Richards_PK::ClipHydrostaticPressure(double pmin, double s0, Epetra_MultiVector& p)
{
  for (int c = 0; c < ncells_owned; c++) {
    if (p[0][c] < pmin) {
      double pc = wrm_->second[(*wrm_->first)[c]]->capillaryPressure(s0);
      p[0][c] = atm_pressure_ - pc;
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi



