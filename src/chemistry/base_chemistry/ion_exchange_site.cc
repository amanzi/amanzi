/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "ion_exchange_site.hh"

#include <cmath>

#include <iostream>
#include <iomanip>

namespace amanzi {
namespace chemistry {

IonExchangeSite::IonExchangeSite()
    : name_(),
      cation_exchange_capacity_(0.0),
      mineral_name_("bulk"),
      charge_(0.0) {
}  // end IonExchangeSite constructor

IonExchangeSite::IonExchangeSite(const IonxSiteName in_name)
    : name_(in_name),
      cation_exchange_capacity_(0.0),
      mineral_name_("bulk"),
      charge_(0.0) {
}  // end IonExchangeSite constructor

IonExchangeSite::IonExchangeSite(const IonxSiteName name,
                                 const double charge,
                                 const std::string location)
    : name_(name),
      cation_exchange_capacity_(0.0),
      mineral_name_(location),
      charge_(charge) {
}  // end IonExchangeSite constructor

IonExchangeSite::~IonExchangeSite() {
}  // end IonExchangeSite destructor

void IonExchangeSite::Display(void) const {
  std::cout << std::setw(15) << name()
            << std::setw(20) << mineral_name()
            << std::setw(10) << std::fixed << charge()
            << std::setw(10) << std::scientific << cation_exchange_capacity()
            << std::fixed << std::endl;
}  // end Display()

void IonExchangeSite::DisplayResultsHeader(void) const {
  std::cout << std::setw(15) << "Name"
            << std::setw(15) << "CEC"
            << std::endl;
}  // end DisplayResultsHeader()

void IonExchangeSite::DisplayResults(void) const {
  std::cout << std::setw(15) << name()
            << std::setw(15) << cation_exchange_capacity()
            << std::endl;
}  // end DisplayResults()

}  // namespace chemistry
}  // namespace amanzi
