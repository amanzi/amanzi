/*

Functor for evaluating QSat

Author: Ethan Coon (ecoon@lanl.gov)
        Chonggang Xu (cxu@lanl.gov)

Licencse: BSD
*/

#include <cmath>
#include "utils.hh"

namespace Amanzi {
namespace BGC {

double PermafrostDepth(const Epetra_SerialDenseVector& SoilTArr,
                       const Epetra_SerialDenseVector& SoilDArr,
                       double freeze_temp) {
  int i = 0;
  int nSoilLayers = SoilTArr.Length();
  while (SoilDArr[i] > 0. && i < nSoilLayers) i++;
  return SoilDArr[i];
}

// This function calculate the effect of temperature on biological process.
double TEffectsQ10(double Q10, double T, double refT) {
  return std::pow(Q10, 0.1 * (T - refT));
}

} // namespace
} // namespace
