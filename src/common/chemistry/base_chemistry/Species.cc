/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ben Andre
*/

/*
  Chemistry

  Base class for species
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "VerboseObject.hh"

#include "exceptions.hh"
#include "errors.hh"

#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

Species::Species()
  : molality_(1.e-9),
    activity_(1.0),
    act_coef_(1.0),
    ln_molality_(0.0),
    ln_activity_(0.0),
    ln_act_coef_(0.0),
    identifier_(0),
    charge_(0),
    gram_molecular_weight_(0.0),
    ion_size_parameter_(0.0),
    name_(""){};


Species::Species(int id, const std::string& name, const Teuchos::ParameterList& plist)
  : molality_(1.e-9),
    activity_(1.0),
    act_coef_(1.0),
    ln_molality_(0.0),
    ln_activity_(0.0),
    ln_act_coef_(0.0),
    identifier_(id),
    name_(name)
{
  ion_size_parameter_ = 0.0;
  if (plist.isParameter("ion size parameter")) {
    ion_size_parameter_ = plist.get<double>("ion size parameter");
  }

  charge_ = 0.0;
  if (plist.isParameter("charge")) { charge_ = plist.get<int>("charge"); }

  gram_molecular_weight_ = plist.get<double>("gram molecular weight");

  if (identifier() < 0 || gram_molecular_weight() < 0.0 || ion_size_parameter() < 0.0) {
    std::ostringstream oss;
    oss << "Invalid species data, id = " << identifier() << std::endl
        << "   gram molecular weight = " << gram_molecular_weight() << std::endl
        << "      ion size parameter = " << ion_size_parameter() << std::endl;
    Exceptions::amanzi_throw(Errors::Message(oss.str()));
  }
}


void
Species::update(double molality)
{
  molality_ = molality;
  // note that activity coefficient not updated
  activity_ = act_coef_ * molality_;
  ln_molality_ = std::log(molality_);
  ln_act_coef_ = std::log(act_coef_);
  ln_activity_ = ln_molality_ + ln_act_coef_;
}


void
Species::update()
{
  activity_ = act_coef_ * molality_;
  ln_molality_ = std::log(molality_);
  ln_act_coef_ = std::log(act_coef_);
  ln_activity_ = ln_molality_ + ln_act_coef_;
}


void
Species::Display(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << name() << std::fixed << std::setprecision(2) << std::setw(10)
          << charge() << std::setprecision(5) << std::setw(10) << gram_molecular_weight()
          << std::setprecision(2) << std::setw(10) << ion_size_parameter() << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void
Species::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << name() << std::scientific << std::setprecision(5) << std::setw(15)
          << molality() << std::setw(15) << act_coef() << std::setw(15) << activity() << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

} // namespace AmanziChemistry
} // namespace Amanzi
