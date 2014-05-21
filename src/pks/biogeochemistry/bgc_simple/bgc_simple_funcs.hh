/*
Main functions for biogeochemistry on a column.

Author: Ethan Coon (ecoon@lanl.gov)
        Chonggang Xu (cxu@lanl.gov)

License: BSD
*/


#ifndef ATS_BGC_SIMPLE_FUNCS_HH_
#define ATS_BGC_SIMPLE_FUNCS_HH_

#include <vector>

#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"

#include "utils.hh"
#include "PFT.hh"
#include "SoilCarbon.hh"

namespace Amanzi {
namespace BGC {

void advance(double t, double dt, double gridarea,
             const MetData& met,
             const Epetra_SerialDenseVector& SoilTArr,
             const Epetra_SerialDenseVector& SoilArr,
             const Epetra_SerialDenseVector& SoilDArr,
             const Epetra_SerialDenseVector& SoilThicknessArr,
             std::vector<Teuchos::RCP<PFT> >& pftarr,
             std::vector<Teuchos::RCP<SoilCarbon> >& soilcarr,
             Epetra_SerialDenseVector& SoilCO2Arr);

} // namespace
} // namespace

#endif

