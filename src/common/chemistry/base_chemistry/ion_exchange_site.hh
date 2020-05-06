/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  Base class for ion exchange sites (e.g. X- in standard geochemistry notation)
*/

#ifndef AMANZI_CHEMISTRY_IONEXCHANGESITE_HH_
#define AMANZI_CHEMISTRY_IONEXCHANGESITE_HH_

#include <cmath>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziChemistry {

typedef std::string IonxSiteName;
typedef int IonxSiteId; 

class IonExchangeSite {
 public:
  IonExchangeSite();
  IonExchangeSite(const IonxSiteName in_name);
  IonExchangeSite(const IonxSiteName in_name, const double charge, const std::string location);
  virtual ~IonExchangeSite() {};

  virtual void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  // mutators
  void set_cation_exchange_capacity(const double in_value) {
    cation_exchange_capacity_ = in_value;
  };

  void set_charge(const double value) {
    charge_ = value;
  };

  void set_name(const IonxSiteName in_name) {
    name_ = in_name;
  };

  void set_mineral_name(const std::string name) {
    mineral_name_ = name;
  };

  // accessors
  IonxSiteName name(void) const {
    return name_;
  };

  double cation_exchange_capacity(void) const {
    return cation_exchange_capacity_;
  };

  std::string mineral_name(void) const {
    return mineral_name_;
  };

  double charge(void) const {
    return charge_;
  };

 protected:
  IonxSiteName name_;
  double cation_exchange_capacity_;  // units...
  std::string mineral_name_;
  double charge_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
