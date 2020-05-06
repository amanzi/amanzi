/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  Base class for species
*/

#ifndef AMANZI_CHEMISTRY_SPECIES_HH_
#define AMANZI_CHEMISTRY_SPECIES_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziChemistry {

typedef std::string SpeciesName;
typedef int SpeciesId;  // unsigned int?
// typedef std::vector<Species> SpeciesArray;
/*  SpeciesArray is actually defined at the end of the file because we
**  can't use it in the typdef until it is declared.... Put all these
**  in the class to get around this...?  */

class Species {
 public:
  Species();  // this is only present for stl containers, don't use it
  Species(SpeciesId id, SpeciesName name, double charge, double mol_wt,
          double size);
  virtual ~Species() {};

  // update(): calculate the new activity coefficient, set the molarity,
  // activity and associated log values. Need to look at different
  // ActivityCoefficient models and determine what the most generic
  // version of this interface will require.
  void update(const double molarity, const double ionic_strength);
  virtual void update(void);
  virtual void update(const double molality);

  // accessor methods for calculated values
  double molality(void) const {
    return this->molality_;
  }
  double activity(void) const {
    return this->activity_;
  }
  double act_coef(void) const {
    return this->act_coef_;
  }

  // natural log versions
  double ln_molality(void) const {
    return this->ln_molality_;
  }
  double ln_activity(void) const {
    return this->ln_activity_;
  }
  double ln_act_coef(void) const {
    return this->ln_act_coef_;
  }

  // access invariant data
  SpeciesId identifier(void) const {
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
  SpeciesName name(void) const {
    return this->name_;
  }

  // these should only be used by the activity coefficient model
  void act_coef(double d) {
    this->act_coef_ = d;
  }
  void ln_act_coef(double d) {
    this->ln_act_coef_ = d;
  }

  void display(void) const;
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

 protected:
  // Species(const double s_charge, const double s_GMW,
  //         const ActivityCoefficient* s_activityCoefficient, SpeciesName s_name);

  // these are dangerous, should only be used internally. Use the
  // public update() to ensure that all related data gets updated!
  void molality(double d) {
    this->molality_ = d;
  }
  void activity(double d) {
    this->activity_ = d;
  }
  void ln_molality(double d) {
    this->ln_molality_ = d;
  }
  void ln_activity(double d) {
    this->ln_activity_ = d;
  }

  double molality_;
  double activity_;
  double act_coef_;
  double ln_molality_;
  double ln_activity_;
  double ln_act_coef_;

  // ActivityCoefficient* activityCoefficient;

 private:
  SpeciesId identifier_;
  double charge_;  // why is this a double rather than int...?
  double gram_molecular_weight_;
  double ion_size_parameter_;  // angstroms
  SpeciesName name_;

  // these are data and should not be changed during a
  // simulation. Remove this interface and require correct data be
  // provided during initialization! Probably means we need to ditch
  // the empty constructor as well.
  // void identifier(SpeciesId i) { this->identifier_ = i; }
  // void charge(double d) { this->charge_ = d; }
  // void gram_molecular_weight(double d) { this->gram_molecular_weight_ = d; }
  // void ion_size_parameter(double d) { this->ion_size_parameter_ = d; }
  // void name(SpeciesName name) { this->name_ = name; }
};

typedef std::vector<Species> SpeciesArray;

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
