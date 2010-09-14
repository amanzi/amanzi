#ifndef __Species_hpp__
#define __Species_hpp__

#include "Teuchos_RCP.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"


// Base class for species

class Species {

public:
  ~Species();

  typdef std::string SpeciesName;
  typedef int SpeciesId;

  // update(): calculate the new activity coefficient, set the molarity,
  // activity and associated log values. Need to look at different
  // ActivityCoefficient models and determine what the most generic
  // version of this interface will require.
  void update(const double s_molarity, const double s_ionic_strength); 

  // accessor methods
  double molarity(void) const;
  double activity(void) const;
  double act_coeff(void) const;

  // log10 values
  double log_molarity(void) const;
  double log_activity(void) const;
  double log_act_coeff(void) const;

  SpeciesId identifier(void) const;
  int charge(void) const;
  double gram_molecular_weight(void) const;
  SpeciesName name(void) const;

protected:
  Species(const double s_charge, const double s_GMW, const ActivityCoefficient* s_activityCoefficient, SpeciesName s_name);

  double molarity_;
  double activity_;
  double act_coeff_;
  double log_molarity_;
  double log_activity_;
  double log_actCoeff_;
  
  ActivityCoefficient* activityCoefficient;


private:
  SpeciesId identifier_;
  int charge_;
  double gram_molecular_weight_;
  SpeciesName name_;

};

inline double Species::molarity(void) const { return this->molarity_; }
inline double Species::activity(void) const { return this->activity_; }
inline double Species::act_coeff(void) const { return this->act_coeff_; }

inline double Species::log_molarity(void) const { return this->molarity_; }
inline double Species::log_activity(void) const { return this->activity_; }
inline double Species::log_act_coeff(void) const { return this->act_coeff_; }

inline Species::SpeciesId Species::identifier(void) const { return this->identifier_; }
inline int Species::change(void) const { return this->charge_; }
inline double Species::gram_molecular_weight(void) const { return this->gram_molecular_weight_; }
inline Species::SpeciesName Species::name(void) const { return this->name_; }

#endif

