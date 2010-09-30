#ifndef __Species_hpp__
#define __Species_hpp__

// Base class for species

#include <string>

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
  void update(const double molarity, const double ionic_strength); 
  void update(void);
  void update(const double molality);

  // accessor methods
  double molality(void) const { return this->molality_; }
  double activity(void) const { return this->activity_; }
  double act_coef(void) const { return this->act_coef_; }

  // natural log versions
  double ln_molality(void) const { return this->ln_molality_; }
  double ln_activity(void) const { return this->ln_activity_; }
  double ln_act_coef(void) const { return this->ln_act_coef_; }

  SpeciesId identifier(void) const { return this->identifier_; }
  double charge(void) const { return this->charge_; }
  double gram_molecular_weight(void) const { return this->gram_molecular_weight_; }
  double ion_size_parameter(void) const { return this->ion_size_parameter_; }
  SpeciesName name(void) const { return this->name_; }

  void molality(double d) { this->molality_ = d; }
  void activity(double d) { this->activity_ = d; }
  void act_coef(double d) { this->act_coef_ = d; }

  void ln_molality(double d) { this->ln_molality_ = d; }
  void ln_activity(double d) { this->ln_activity_ = d; }
  void ln_act_coef(double d) { this->ln_act_coef_ = d; }

  void identifier(SpeciesId i) { this->identifier_ = i; }
  void charge(double d) { this->charge_ = d; }
  void gram_molecular_weight(double d) { this->gram_molecular_weight_ = d; }
  void ion_size_parameter(double d) { this->ion_size_parameter_ = d; }
  void name(SpeciesName name) { this->name_ = name; }

  void display(void) const;

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
  double charge_; // why is this a double rather than int...?
  double gram_molecular_weight_;
  double ion_size_parameter_;
  SpeciesName name_;

};

#endif // __Species_hpp__

