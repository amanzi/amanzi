/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_REACTION_HH_
#define AMANZI_CHEMISTRY_REACTION_HH_

//
// Base class for describing reactions
//
// Questions:
//
// Are we considering this to be a reaction of the form:
//
// (1) aA + bB <---> cC + dD
//
// or
// 
// (2) aA = cC + dD - bB
//
// Or are we starting with (1) and using the basis swapping
// transformation to generate (2)?
// How are we going to deal with water in a reaction, is it an
// aqueous species? is it always present?
//
// How are we going to handle temperature dependance? Store and update
// temperature, or pass it in as a parameter to function calls?
//
//

#include <string>
#include <vector>
#include <ostream>

#include "aqueous_species.hh"

namespace amanzi {
namespace chemistry {

class Reaction {
 public:
  Reaction();  // this is only present for stl containers, don't use it
  Reaction(const AqueousSpeciesArray& aqueous_species,
           const std::vector<double>& stoichiometries,
           const std::vector<double>& log10_Keq_temperature);
  virtual ~Reaction();

  //
  // public interface
  //
  virtual void Update(const AqueousSpecies::AqueousSpeciesArray& primary_species, 
                      const AqueousSpecies& water_species) = 0;
  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double> *total) = 0;
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<AqueousSpecies>& primary_species,
                                       Block* dtotal) = 0;

  void Display(std::ostream& output) const;


 protected:

 private:

  AqueousSpecies::AqueousSpeciesArray aqueous_species_;
  std::vector<double> stoichiometries_;
  std::vector<double> log10_Keq_temperature_;  // array of temperature dependent log10 Keq

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_REACTION_HH_
