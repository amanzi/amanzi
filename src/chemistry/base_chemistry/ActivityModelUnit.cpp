#include <cmath>
#include <iostream>
#include "Activity.hpp"

Activity::Activity() : I_(0.),
                       debyeA(0.4939),
                       debyeB(0.3253),
                       debyeBdot(0.0374)
{
} // end Activity constructor


Activity::~Activity() 
{
  
} // end Activity destructor

void Activity::calculateIonicStrength(
                     std::vector<Species> primarySpecies,
                     std::vector<AqueousEquilibriumComplex> secondarySpecies)
{

  // I = 0.5 * sum_i(m_i*z_i^2)

  I_ = 0.;

  // primary species
  for (std::vector<Species>::iterator i=primarySpecies.begin();
       i != primarySpecies.end(); i++) {
    I_ += i->molality() * i->charge() * i->charge();
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::iterator i = secondarySpecies.begin();
       i != secondarySpecies.end(); i++) {
    I_ += i->molality() * i->charge() * i->charge();
  }

  I_ *= 0.5;

} // end calculateIonicStrength()

void Activity::calculateActivityCoefficients(
                     std::vector<Species> &primarySpecies,
                     std::vector<AqueousEquilibriumComplex> &secondarySpecies)
{

 
  // primary species
  for (std::vector<Species>::iterator i=primarySpecies.begin();
       i != primarySpecies.end(); i++) {
    double gamma = calculateActivityCoefficient(i->charge(),
                                                i->ion_size_parameter());
    i->act_coef(gamma);
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::iterator i = secondarySpecies.begin();
       i != secondarySpecies.end(); i++) {
    double gamma = calculateActivityCoefficient(i->charge(),
                                                i->ion_size_parameter());
    i->act_coef(gamma);
  }

} // end calculateActivityCoefficients()

double Activity::calculateActivityCoefficient(double charge, 
                                              double ion_size_parameter)
{
  // log(gamma_i) = - A * z_i^2 * sqrt(I) / (1 + a0 * B * sqrt(I)) + Bdot * I

  double sqrt_I = std::sqrt(I_);
  double log_gamma = -debyeA * charge * charge * sqrt_I / 
                     (1 + ion_size_parameter * debyeB * sqrt_I) +
                     debyeBdot * I_;
  return std::exp(log_to_ln(log_gamma));

} // end calculateActivityCoefficient()

void Activity::display(void) const
{
} // end display()
