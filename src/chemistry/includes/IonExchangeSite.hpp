/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __IonExchangeSite_hpp__
#define __IonExchangeSite_hpp__

/*
**  Base class for ion exchange sites (e.g. X- in standard geochemistry notation)
**
**
*/

#include <string>
#include <vector>
#include <cmath>

#include "Species.hpp"
#include "Block.hpp"

class IonExchangeSite : public Species {

 public:
  IonExchangeSite();
  IonExchangeSite(const SpeciesName exchanger_name,
                  const SpeciesId exchanger_id,
                  const double exchanger_charge, 
                  const std::string exchanger_location,
                  const double mol_wt,
                  const double size);

  virtual ~IonExchangeSite();

  virtual void update(void);
  virtual void update(const double in_molality);
  virtual void Display(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

  void set_cation_exchange_capacity(const double in_value) { this->cation_exchange_capacity_ = in_value; };
  double cation_exchange_capacity(void) const { return this->cation_exchange_capacity_; };

  std::string location(void) const { return this->location_; };

 protected:
  double cation_exchange_capacity_;  // units...
  std::string location_;

 private:
};

#endif  // __IonExchangeSite_hpp__
