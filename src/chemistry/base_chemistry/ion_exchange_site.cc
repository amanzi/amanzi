/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "ion_exchange_site.hh"

#include <cmath>

#include <iostream>
#include <iomanip>

namespace amanzi {
namespace chemistry {

IonExchangeSite::IonExchangeSite()
    : Species(),
      cation_exchange_capacity_(0.0) {
}  // end IonExchangeSite constructor


IonExchangeSite::IonExchangeSite(const SpeciesName exchanger_name,
                                 const SpeciesId exchanger_id,
                                 const double exchanger_charge,
                                 const std::string exchanger_location,
                                 const double mol_wt, const double size)
    : Species(exchanger_id, exchanger_name, exchanger_charge, mol_wt, size),
      cation_exchange_capacity_(0.0),
      location_(exchanger_location) {
}  // end IonExchangeSite constructor

IonExchangeSite::~IonExchangeSite() {
}  // end IonExchangeSite destructor


void IonExchangeSite::update(void) {
  // dummy update function, activity/conc/act coeff always equal to one?
  molality_ = 1.0;
  act_coef_ = 1.0;
  activity_ = act_coef() * molality();
  ln_molality_ = std::log(molality());
  ln_act_coef_ = std::log(act_coef());
  ln_activity_ = ln_molality() + ln_act_coef();
}  // end Update()

void IonExchangeSite::update(const double in_molality) {
  // dummy update function, activity/conc/act coeff always equal to one?
  molality_ = 1.0;
  act_coef_ = 1.0;
  activity_ = act_coef() * molality();
  ln_molality_ = std::log(molality());
  ln_act_coef_ = std::log(act_coef());
  ln_activity_ = ln_molality() + ln_act_coef();
}  // end Update()


void IonExchangeSite::Display(void) const {
  std::cout << std::setw(15) << name()
            << std::setw(15) << location()
            << std::setw(10) << charge()
            << std::setw(10) << cation_exchange_capacity()
            << std::endl;
}  // end Display()

void IonExchangeSite::DisplayResultsHeader(void) const {
  std::cout << std::setw(15) << "Name"
            << std::setw(15) << "CEC"
            << std::endl;
}  // end DisplayResultsHeader()

void IonExchangeSite::DisplayResults(void) const {
  std::cout << std::setw(15) << name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << cation_exchange_capacity()
            << std::endl;
}  // end DisplayResults()

}  // namespace chemistry
}  // namespace amanzi
