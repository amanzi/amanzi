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
namespace BGCSimple {

// This function calculate the effect of temperature on biological process.
double TEffectsQ10(double Q10, double T, double refT) {
  return std::pow(Q10, 0.1 * (T - refT));
}


} // namespace
} // namespace
} // namespace
