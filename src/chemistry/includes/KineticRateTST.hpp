/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __KINETIC_RATE_TST_HPP__

#define __KINETIC_RATE_TST_HPP__

/*******************************************************************************
**
**  Description: abstract base class for all kinetic rates
**
*******************************************************************************/
#include <vector>

#include "KineticRate.hpp"
#include "Species.hpp"
#include "StringTokenizer.hpp"

class Block; 

class KineticRateTST : public KineticRate
{
 public:
  KineticRateTST(void);
  ~KineticRateTST(void);

  void Update(const std::vector<Species> primarySpecies);
  void AddContributionToResidual(const double por_den_sat_vol, 
                                 std::vector<double> *residual);
                                 
  void AddContributionToJacobian(const std::vector<Species> primarySpecies,
                                 const double por_den_sat_vol,
                                 Block *J);
  void Display(void) const;

  void ParseParameters(StringTokenizer rate);

 protected:

 private:

};

#endif     /* __KINETIC_RATE_TST_HPP__ */

