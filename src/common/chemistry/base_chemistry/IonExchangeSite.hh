/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  Base class for ion exchange sites (e.g. X- in standard notation)
*/

#ifndef AMANZI_CHEMISTRY_ION_EXCHANGE_SITE_HH_
#define AMANZI_CHEMISTRY_ION_EXCHANGE_SITE_HH_

#include <cmath>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziChemistry {

class IonExchangeSite {
 public:
  IonExchangeSite();
  IonExchangeSite(const std::string& name);
  IonExchangeSite(const std::string& name, const double charge, const std::string& location);
  virtual ~IonExchangeSite() {};

  virtual void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  // access and mutators
  void set_cation_exchange_capacity(double value) { cation_exchange_capacity_ = value; }
  double get_cation_exchange_capacity() const { return cation_exchange_capacity_; }

  void set_charge(const double value) { charge_ = value; }
  double get_charge() const { return charge_; }

  std::string get_name() const { return name_; }

  void set_mineral_name(const std::string& name) { mineral_name_ = name; }
  std::string get_mineral_name() const { return mineral_name_; }

 protected:
  std::string name_;
  double cation_exchange_capacity_;  // units...
  std::string mineral_name_;
  double charge_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
