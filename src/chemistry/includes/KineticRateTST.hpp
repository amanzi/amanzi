/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __KINETIC_RATE_TST_HPP__

#define __KINETIC_RATE_TST_HPP__

/*******************************************************************************
**
**  Description: implementation of the TST rate law for mineral kinetics
**
**  R = k * A * Prod (a_i^m_i) * ( 1 - Q/Keq)
**
**  where:
**    R : reaction rate, [moles/sec]
**    Keq : equilibrium constant, [-]
**    Q : ion activity product, [-]
**    k : reaction rate constant, [moles m^2 s^-1]
**    A : reactive surface area, [m^2]
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

  void Setup(std::string reaction, StringTokenizer reaction_data,
             const SpeciesArray primary_species);
  void Update(const SpeciesArray primary_species);
  void AddContributionToResidual(const double por_den_sat_vol, 
                                 std::vector<double> *residual);
                                 
  void AddContributionToJacobian(const SpeciesArray primary_species,
                                 const double por_den_sat_vol,
                                 Block *J);
  void Display(void) const;

  void ParseParameters(StringTokenizer rate);

  /*
  ** end of KineticRate inherited interface
  */


 protected:

  Species mineral(void) const { return this->mineral_; };

  void area(double set_area) { this->area_ = set_area; };
  double area(void) const { return this->area_; };
  void pK(double set_pK) { this->pK_ = set_pK; };
  double pK(void) const { return this->pK_; };
  void rate_constant(double set_rate_constant) { this->rate_constant_ = set_rate_constant; };
  double rate_constant(void) const { return this->rate_constant_; };
  void sat_state_exponent(double set_sat_state_exponent) { this->sat_state_exponent_ = set_sat_state_exponent; };
  double sat_state_exponent(void) const { return this->sat_state_exponent_; };


 private:
  double area_;  // surface area [m^2]
  double pK_;  // pK [-]
  double rate_constant_;  // k, rate constant, [moles/m^2/sec]
  double sat_state_exponent_;  // n, saturation state exponent, [-]
  Species mineral_;
  std::vector<SpeciesName> modifying_species_names;
  std::vector<double> modifying_exponents;
  std::vector<int> modifying_primary_ids;
  std::vector<double> modifying_primary_exponents;
  std::vector<int> modifying_secondary_ids;
  std::vector<double> modifying_secondary_exponents;
};

#endif     /* __KINETIC_RATE_TST_HPP__ */

