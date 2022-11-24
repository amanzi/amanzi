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
void
Richards_PK::ClipHydrostaticPressure(double pmin, Epetra_MultiVector& p)
{
  for (int c = 0; c < ncells_owned; c++) p[0][c] = std::max(p[0][c], pmin);
}


/* ******************************************************************
* Clip pressure using constant saturation.
****************************************************************** */
void
Richards_PK::ClipHydrostaticPressure(double pmin, double s0, Epetra_MultiVector& p)
{
  for (int c = 0; c < ncells_owned; c++) {
    if (p[0][c] < pmin) {
      double pc = wrm_->second[(*wrm_->first)[c]]->capillaryPressure(s0);
      p[0][c] = atm_pressure_ - pc;
    }
  }
}


/* ****************************************************************
* Plot water retention curves.
**************************************************************** */
void
Richards_PK::PlotWRMcurves_()
{
  int MyPID = mesh_->cell_map(false).Comm().MyPID();
  if (MyPID != 0) return;

  int ndata(1000);
  for (int n = 0; n < wrm_->second.size(); ++n) {
    std::ofstream ofile;
    std::string filename("wrm_" + std::to_string(n) + ".txt");
    ofile.open(filename.c_str());

    double sr = wrm_->second[n]->residualSaturation();
    double ds = (1.0 - sr) / ndata;

    for (int i = 0; i < ndata; i++) {
      double sat = sr + ds * (i + 0.5);
      double pc = wrm_->second[n]->capillaryPressure(sat);
      double krel = wrm_->second[n]->k_relative(pc);
      ofile << sat << " " << krel << " " << pc << std::endl;
    }
    ofile << std::endl;
    ofile.close();
  }
}

} // namespace Flow
} // namespace Amanzi
