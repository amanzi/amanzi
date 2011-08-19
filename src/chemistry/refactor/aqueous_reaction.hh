/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_AQUEOUS_REACTION_HH_
#define AMANZI_CHEMISTRY_AQUEOUS_REACTION_HH_

//
// Class describing aqueous equilibrium reactions.
//
//
//

#include <vector>
#include <ostream>

#include "aqueous_species.hh"
#include "reaction.hh"


namespace amanzi {
namespace chemistry {

// forward declarations from chemistry
class Block;

class AqueousReaction : public Reaction {
 public:
  AqueousReaction();
  AqueousReaction(const AqueousReactionName name,
                  const AqueousReactionId aqueous_reaction_id,
                  const std::vector<AqueousSpeciesName>& species,
                  const std::vector<double>& stoichiometries,
                  const std::vector<AqueousSpeciesId>& species_ids,
                  const double h2o_stoich,
                  const std::vector<double>& logKeq_temperature);

  virtual ~AqueousReaction();

  // update molalities
  virtual void Update(const std::vector<AqueousSpecies>& primary_species);
  // add stoichiometric contribution of complex to totalmembermember
  virtual void AddContributionToTotal(std::vector<double> *total) = 0;
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<AqueousSpecies>& primary_species,
                                       Block* dtotal) = 0;

  int num_components(void) const {
    return this->num_components_in_rxncomp_;
  };

  // do we really want this mutator (and do we want it to be public),
  // or do we want something like set_logK(temperature)?
  void set_logK(const double in_logK) {
    this->logK_ = in_logK;
  };

  std::vector<AqueousSpeciesName> species_names(void) const {
    return this->species_names_;
  };
  std::vector<AqueousSpeciesId> species_ids(void) const {
    return this->species_ids_;
  };
  std::vector<double> stoichiometry(void) const {
    return this->stoichiometry_;
  };

  virtual void Display(std::ostream& output) const;

 protected:
 private:
  double logK(void) const {
    return this->logK_;
  };
  double lnK(void) const {
    return this->lnK_;
  };
  double lnQK(void) const {
    return this->lnQK_;
  };

  double log_to_ln(double d) {
    return d * 2.30258509299;
  }
  double ln_to_log(double d) {
    return d * 0.434294481904;
  }

  int num_components_in_rxn_;  // # components in reaction
  double h2o_stoich_;                  // stoichiometry of water in equation
  double lnK_;                         // log value of equlibrium constant
  double lnQK_;                        // store lnQK for derivatives later
  double logK_;

  // should not be able to change num_components after it is set in the constructor...?
  void set_num_components_in_rxn(const int nun_components_in_rxn) {
    this->num_components_in_rxn_ = num_components_in_rxn;
  };
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_AQUEOUS_REACTION_HH_
