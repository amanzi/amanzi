/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_GENERAL_REACTION_HH_
#define AMANZI_CHEMISTRY_GENERAL_REACTION_HH_

//
// Class implementing general reactions
//
//

#include <string>
#include <vector>
#include <ostream>

#include "aqueous_species.hh"
#include "reaction.hh"

namespace amanzi {
namespace chemistry {

class GeneralReaction : public Reaction {
 public:
  GeneralReaction();
  GeneralReaction(const AqueousSpeciesArray& aqueous_species,
                  const std::vector<double>& forward_stoichiometries,
                  const std::vector<double>& backward_stoichiometries,
                  const std::vector<double>& k_forward_temperature,
                  const std::vector<double>& k_backward_temperature);
  virtual ~GeneralReaction();

  //
  // public interface
  //
  virtual void Update(const AqueousSpecies::AqueousSpeciesArray& primary_species, 
                      const AqueousSpecies& water_species);
  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<AqueousSpecies>& primary_species,
                                       Block* dtotal);

  void Display(std::ostream& output) const;

 protected:

 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_GENERAL_REACTION_HH_
