/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Reaction_hpp__
#define __Reaction_hpp__

/*
** Description: Base class for reactions (aqueous equilibrium
** complexes, minerals)
*/

#include <string>
#include <vector>
#include <cmath>

#include "Species.hpp"

// forward declarations
class Block;

class Reaction {
 public:
  virtual ~Reaction();

  // update molalities
  virtual void Update(const std::vector<Species> primary_species) = 0;
  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double>& total) = 0;
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<Species> primary_species,
                                       Block* dtotal) = 0;

  virtual void Display(void) const = 0;

  void set_name(const std::string& s) { this->name_ = s; };
  std::string name(void) const { return this->name_; };

 protected:
  Reaction();
  Reaction(const std::string in_name);

  double log_to_ln(double d) { return d*2.30258509299; }
  double ln_to_log(double d) { return d*0.434294481904; }

 private:
  std::string name_;
};

#endif  // __Reaction_hpp__
