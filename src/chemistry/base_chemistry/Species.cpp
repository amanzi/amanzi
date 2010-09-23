#include <cmath>
#include <iostream>
#include "Species.hpp"

Species::Species() 
  : molality_(0.0), 
    activity_(1.0),
    act_coef_(1.0),
    ln_molality_(0.0),
    ln_activity_(0.0),
    ln_act_coef_(0.0),
    identifier_(0),
    charge_(0),
    gram_molecular_weight_(0.0),
    ion_size_parameter_(0.0),
    name_("")
{
//  ActivityCoefficient* activityCoefficient;
  // end Species()
}


Species::Species(SpeciesId id, SpeciesName name, double charge, double mol_wt, 
                 double size)
  : molality_(0.0), 
    activity_(1.0),
    act_coef_(1.0),
    ln_molality_(0.0),
    ln_activity_(0.0),
    ln_act_coef_(0.0),
    identifier_(id),
    charge_(charge),
    gram_molecular_weight_(mol_wt),
    ion_size_parameter_(size),
    name_(name)
{
  //  ActivityCoefficient* activityCoefficient;
  // end Species()
}

Species::~Species() {
  // end ~Species()
}

void Species::update(const double molality) 
{
  molality_ = molality;
  act_coef_ = 1.; // to be replaced with call to activiti coef function
  activity_ = act_coef_ * molality_;
  ln_molality_ = std::log(molality_);
  ln_act_coef_ = std::log(act_coef_);
  ln_activity_ = ln_molality_ + ln_act_coef_;
  // end update()
}

void Species::update(void)
{
  act_coef_ = 1.; // to be replaced with call to activiti coef function
  activity_ = act_coef_ * molality_;
  ln_molality_ = std::log(molality_);
  ln_act_coef_ = std::log(act_coef_);
  ln_activity_ = ln_molality_ + ln_act_coef_;
  // end update()
}

void Species::display(void) const
{
  std::cout << "    " << get_name();
  std::cout << std::endl;
  std::cout << "        charge = " << get_charge() << std::endl;
  std::cout << "        mol wt = " << get_gram_molecular_weight() << std::endl;
  
  // end display()
}
