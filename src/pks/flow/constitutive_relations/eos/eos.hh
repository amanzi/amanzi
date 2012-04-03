#ifndef _EOS_HH_
#define _EOS_HH_

// Equation of State model
namespace Amanzi {
class EOS {

public:
  virtual double density(double p, double T) = 0;
  virtual double viscosity(double p, double T) = 0;
};
}
#endif
