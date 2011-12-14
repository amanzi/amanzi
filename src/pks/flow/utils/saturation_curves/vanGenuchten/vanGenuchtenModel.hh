#ifndef _VANGENUCHTENMODEL_HH_
#define _VANGENUCHTENMODEL_HH_

#include "SaturationCurve.hh"

namespace Amanzi {
class vanGenuchtenModel : public SaturationCurve {

public:
  vanGenuchtenModel(int meshblock, double lambda, double P0);

  // overridden from WaterRetentionBaseModel
  double s_star(double pc);

private:
  const double lambda_;     // van Genuchten m
  double nu_;           // van Genuchten n
  const double P0; // van Genuchten alpha
};
}
#endif
