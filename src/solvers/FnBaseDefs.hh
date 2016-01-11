/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FNBASE_DEFS_HH_
#define AMANZI_FNBASE_DEFS_HH_

namespace Amanzi {
namespace AmanziSolvers {
// contained in its own namespace for easier using
namespace FnBaseDefs {

// enum for ModifyCorrection control
enum ModifyCorrectionResult {
  CORRECTION_NOT_MODIFIED = 0,
  CORRECTION_MODIFIED = 1,
  CORRECTION_MODIFIED_LAG_BACKTRACKING
};

}  // namespace FnBaseDefs
}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
