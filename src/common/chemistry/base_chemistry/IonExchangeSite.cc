/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  Base class for ion exchange sites (e.g. X- in standard geochemistry notation)
*/

#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "Teuchos_ParameterList.hpp"
#include "VerboseObject.hh"

#include "IonExchangeSite.hh"

namespace Amanzi {
namespace AmanziChemistry {

IonExchangeSite::IonExchangeSite()
  : name_(), mineral_name_("bulk"), cation_exchange_capacity_(0.0), charge_(0.0)
{}


IonExchangeSite::IonExchangeSite(const std::string& name, const Teuchos::ParameterList& plist)
  : name_(name), cation_exchange_capacity_(0.0)
{
  charge_ = plist.get<int>("charge");
  mineral_name_ = plist.get<std::string>("location");
}


void
IonExchangeSite::Display(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << name_ << std::setw(20) << get_mineral_name() << std::setw(10)
          << std::fixed << charge_ << std::setw(10) << std::scientific << cation_exchange_capacity_
          << std::fixed << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void
IonExchangeSite::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << "Name" << std::setw(15) << "CEC" << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void
IonExchangeSite::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << name_ << std::setw(15) << cation_exchange_capacity_ << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

} // namespace AmanziChemistry
} // namespace Amanzi
