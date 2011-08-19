/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_AQUEOUS_SPECIES_HH_
#define AMANZI_CHEMISTRY_AQUEOUS_SPECIES_HH_

//
// Class for describing aqueous species. 
//
// Notes: will this class need to be subclassed (make the private mutators protected)?
// 

#include <string>
#include <vector>
#include <ostream>

namespace amanzi {
namespace chemistry {

class AqueousSpecies {
 public:
  typedef std::vector<AqueousSpecies> AqueousSpeciesArray;
  typedef std::string AqueousSpeciesName;
  typedef int AqueousSpeciesId;  // unsigned int?

  AqueousSpecies();  // this is only present for stl containers, don't use it
  AqueousSpecies(const AqueousSpeciesId& id, const AqueousSpeciesName& name, 
                 const double charge, const double mol_wt,
                 const double size);
  virtual ~AqueousSpecies();

  // update(): calculate the new activity coefficient, set the molarity,
  // activity and associated log values. 
  void Update(const double molarity, const double ionic_strength);
  virtual void Update(void);
  virtual void Update(const double molality);
  void Display(std::ostream& output) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

  //
  // accessor methods
  //

  // concentrations
  double molality(void) const {
    return this->molality_;
  }
  double activity(void) const {
    return this->activity_;
  }
  double activity_coeff(void) const {
    return this->activity_coeff_;
  }

  // natural log versions
  double ln_molality(void) const {
    return this->ln_molality_;
  }
  double ln_activity(void) const {
    return this->ln_activity_;
  }
  double ln_activity_coeff(void) const {
    return this->ln_activity_coeff_;
  }

  // invariant data
  AqueousSpeciesId identifier(void) const {
    return this->identifier_;
  }
  double charge(void) const {
    return this->charge_;
  }
  double gram_molecular_weight(void) const {
    return this->gram_molecular_weight_;
  }
  double ion_size_parameter(void) const {
    return this->ion_size_parameter_;
  }
  AqueousSpeciesName name(void) const {
    return this->name_;
  }


  //
  // mutator methods
  //

  // these should only be used by the activity coefficient model
  void set_activity_coeff(double d) {
    this->activity_coeff_ = d;
  }
  void ln_activity_coeff(double d) {
    this->ln_activity_coeff_ = d;
  }

 protected:

 private:
  //
  // mutator methods
  //

  // these are dangerous, should only be used internally. Use the
  // public update() to ensure that all related data gets updated!
  void set_molality(double d) {
    this->molality_ = d;
  }
  void set_activity(double d) {
    this->activity_ = d;
  }
  void set_ln_molality(double d) {
    this->ln_molality_ = d;
  }
  void set_ln_activity(double d) {
    this->ln_activity_ = d;
  }

  AqueousSpeciesId identifier_;
  double charge_;  // why is this a double rather than int...? to avoid casting?
  double gram_molecular_weight_;  // [grams/mole]
  double ion_size_parameter_;  // angstroms
  AqueousSpeciesName name_;

  double molality_;  // [moles/kg]
  double activity_;  // [moles/kg]
  double activity_coeff_;  // [-]
  double ln_molality_;
  double ln_activity_;
  double ln_activity_coeff_;


};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_AQUEOUS_SPECIES_HH_
