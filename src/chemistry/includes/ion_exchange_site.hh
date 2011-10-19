/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_IONEXCHANGESITE_HH_
#define AMANZI_CHEMISTRY_IONEXCHANGESITE_HH_

/*
**  Base class for ion exchange sites (e.g. X- in standard geochemistry notation)
**
**
*/
#include <cmath>

#include <string>
#include <vector>

#include "block.hh"

namespace amanzi {
namespace chemistry {

typedef std::string IonxSiteName;
typedef int IonxSiteId; 

class IonExchangeSite {
 public:
  IonExchangeSite();
  IonExchangeSite(const IonxSiteName in_name);
  virtual ~IonExchangeSite();

  virtual void Display(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

  void set_cation_exchange_capacity(const double in_value) {
    cation_exchange_capacity_ = in_value;
  };
  void set_name(const IonxSiteName in_name) {
    name_ = in_name;
  };
  double cation_exchange_capacity(void) const {
    return cation_exchange_capacity_;
  };
  IonxSiteName name(void) const {
    return name_;
  };

 protected:
  IonxSiteName name_;
  double cation_exchange_capacity_;  // units...

 private:
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_IONEXCHANGESITE_HH_
