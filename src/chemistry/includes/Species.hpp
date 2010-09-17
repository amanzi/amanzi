#ifndef __Species_hpp__
#define __Species_hpp__

// Base class for species

#include <string>
#include <math.h>

typedef std::string SpeciesName;
typedef int SpeciesId;

class Species {

public:
  Species();
  ~Species();


  Species(SpeciesId id, SpeciesName name, double charge, double mol_wt, 
          double size);

  // update(): calculate the new activity coefficient, set the molarity,
  // activity and associated log values. Need to look at different
  // ActivityCoefficient models and determine what the most generic
  // version of this interface will require.
  void update(const double s_molarity, const double s_ionic_strength); 
  void update(void);
  void update(const double molality);

  // accessor methods
  double get_molality(void) const { return this->molality_; }
  double get_activity(void) const { return this->activity_; }
  double get_act_coef(void) const { return this->act_coef_; }

  // natural log versions
  double get_ln_molality(void) const { return this->ln_molality_; }
  double get_ln_activity(void) const { return this->ln_activity_; }
  double get_ln_act_coef(void) const { return this->ln_act_coef_; }

  SpeciesId get_identifier(void) const { return this->identifier_; }
  double get_charge(void) const { return this->charge_; }
  double get_gram_molecular_weight(void) const { return this->gram_molecular_weight_; }
  double get_ion_size_parameter(void) const { return this->ion_size_parameter_; }
  SpeciesName get_name(void) const { return this->name_; }

  void set_molality(double d) { this->molality_ = d; }
  void set_activity(double d) { this->activity_ = d; }
  void set_act_coef(double d) { this->act_coef_ = d; }

  void set_ln_molality(double d) { this->ln_molality_ = d; }
  void set_ln_activity(double d) { this->ln_activity_ = d; }
  void set_ln_act_coef(double d) { this->ln_act_coef_ = d; }

  void set_identifier(SpeciesId i) { this->identifier_ = i; }
  void set_charge(double d) { this->charge_ = d; }
  void set_gram_molecular_weight(double d) { this->gram_molecular_weight_ = d; }
  void set_ion_size_parameter(double d) { this->ion_size_parameter_ = d; }
  void set_name(SpeciesName name) { this->name_ = name; }

protected:
//  Species(const double s_charge, const double s_GMW, const ActivityCoefficient* s_activityCoefficient, SpeciesName s_name);

  double molality_;
  double activity_;
  double act_coef_;
  double ln_molality_;
  double ln_activity_;
  double ln_act_coef_;
  
//  ActivityCoefficient* activityCoefficient;


private:
  SpeciesId identifier_;
  double charge_;
  double gram_molecular_weight_;
  double ion_size_parameter_;
  SpeciesName name_;

};

#endif

