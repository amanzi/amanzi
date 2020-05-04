/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  Base class for ion exchange sites (e.g. X- in standard geochemistry notation)
*/

#include "ion_exchange_site.hh"

#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziChemistry {

IonExchangeSite::IonExchangeSite()
    : name_(),
      cation_exchange_capacity_(0.0),
      mineral_name_("bulk"),
      charge_(0.0) {
}

IonExchangeSite::IonExchangeSite(const IonxSiteName in_name)
    : name_(in_name),
      cation_exchange_capacity_(0.0),
      mineral_name_("bulk"),
      charge_(0.0) {
}

IonExchangeSite::IonExchangeSite(const IonxSiteName name,
                                 const double charge,
                                 const std::string location)
    : name_(name),
      cation_exchange_capacity_(0.0),
      mineral_name_(location),
      charge_(charge) {
}


void IonExchangeSite::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::setw(20) << mineral_name()
          << std::setw(10) << std::fixed << charge()
          << std::setw(10) << std::scientific << cation_exchange_capacity()
          << std::fixed << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}


void IonExchangeSite::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << "Name"
          << std::setw(15) << "CEC"
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}


void IonExchangeSite::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::setw(15) << cation_exchange_capacity()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
