/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_MINERAL_REACTION_HH_
#define AMANZI_CHEMISTRY_MINERAL_REACTION_HH_

//
// Class implementing mineral reactions
//
//

#include <string>
#include <vector>
#include <ostream>

#include "aqueous_species.hh"
#include "mineral_species.hh"
#include "reaction.hh"

namespace amanzi {
namespace chemistry {

class MineralReaction : public Reaction {
 public:
  typedef std::vector<MineralReaction> MineralReactionArray;

  MineralReaction();  // this is only present for stl containers, don't use it
  MineralReaction(const MineralSpecies& mineral,
                  const AqueousSpeciesArray& aqueous_species,
                  const std::vector<double>& stoichiometries,
                  const std::vector<double>& log10_Keq_temperature,
                  const std::string& rate_name);
  virtual ~MineralReaction();

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
  ReactionRate* reaction_rate_;

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_MINERAL_REACTION_HH_
