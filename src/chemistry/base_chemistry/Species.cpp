#include "Species.hpp"

Species::Species() {

  molality_ = 0.;
  activity_ = 1.;
  act_coef_ = 1.;
  ln_molality_ = 0.;
  ln_activity_ = 0.;
  ln_act_coef_ = 0.;
 
//  ActivityCoefficient* activityCoefficient;
  
  identifier_ = 0;
  charge_ = 0.;
  gram_molecular_weight_ = 0.;
  ion_size_parameter_ = 0.;
  name_ = "";

}

Species::Species(SpeciesId id, SpeciesName name, double charge, double mol_wt, 
                 double size) {

  molality_ = 0.;
  activity_ = 0.;
  act_coef_ = 1.;
  ln_molality_ = 0.;
  ln_activity_ = 0.;
  ln_act_coef_ = 0.;
 
//  ActivityCoefficient* activityCoefficient;
  
  identifier_ = 0;
  charge_ = charge;
  gram_molecular_weight_ = mol_wt;
  ion_size_parameter_ = size;
  name_ = name;

}

void Species::update(const double molality) {

  molality_ = molality;
  act_coef_ = 1.; // to be replaced with call to activiti coef function
  activity_ = act_coef_*molality_;
  ln_molality_ = log(molality_);
  ln_act_coef_ = log(act_coef_);
  ln_activity_ = ln_molality_+ln_act_coef_;

}

void Species::update(void) {

  act_coef_ = 1.; // to be replaced with call to activiti coef function
  activity_ = act_coef_*molality_;
  ln_molality_ = log(molality_);
  ln_act_coef_ = log(act_coef_);
  ln_activity_ = ln_molality_+ln_act_coef_;

}

Species::~Species() {
}
