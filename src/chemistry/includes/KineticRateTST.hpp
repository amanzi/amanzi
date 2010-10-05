/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __KINETIC_RATE_TST_HPP__

#define __KINETIC_RATE_TST_HPP__

/*******************************************************************************
**
**  Description: implementation of the TST rate law for mineral kinetics
**
**  R = k * A * Prod (a_i^m_i) * ( 1 - Q/Keq)^n
**
**  where:
**    R : reaction rate, [moles/sec]
**    Keq : equilibrium constant, [-]
**    Q : ion activity product, [-]
**    k : reaction rate constant, [moles m^2 s^-1]
**    A : reactive surface area, [m^2]
**    n : saturation state exponent [-]
**    a_i : activity of modifying species (primary or secondary)
**    m_i : exponent of modifying species
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

  void area(double set_area) { this->area_ = set_area; };
  double area(void) const { return this->area_; };
  void pK(double set_pK) { this->pK_ = set_pK; };
  double pK(void) const { return this->pK_; };
  void rate_constant(double set_rate_constant) { this->rate_constant_ = set_rate_constant; };
  double rate_constant(void) const { return this->rate_constant_; };
  void sat_state_exponent(double set_sat_state_exponent) { this->sat_state_exponent_ = set_sat_state_exponent; };
  double sat_state_exponent(void) const { return this->sat_state_exponent_; };

 protected:

 private:
  double area_;  // surface area [m^2]
  double pK_;  // pK [-]
  double rate_constant_;  // k, rate constant, [moles/m^2/sec]
  double sat_state_exponent_;  // n, saturation state exponent, [-]
  std::vector<SpeciesName> modifying_species_names;
  std::vector<double> modifying_exponents;
  std::vector<int> modifying_species_ids;
};

#endif     /* __KINETIC_RATE_TST_HPP__ */

