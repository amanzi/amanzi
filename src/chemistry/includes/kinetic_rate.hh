/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_KINETIC_RATE_HH_

#define AMANZI_CHEMISTRY_KINETIC_RATE_HH_

/*******************************************************************************
 **
 **  Description: abstract base class for all kinetic rates
 **
 *******************************************************************************/
#include <vector>
#include <string>

#include "species.hh"
#include "secondary_species.hh"
#include "mineral.hh"
#include "string_tokenizer.hh"
#include "verbosity.hh"

namespace amanzi {
namespace chemistry {

class Block;

class KineticRate {
 public:
  virtual ~KineticRate(void);

  virtual void Setup(const SecondarySpecies& reaction,
                     const StringTokenizer& reaction_data,
                     const SpeciesArray& primary_species) = 0;
  virtual void Update(const SpeciesArray& primary_species,
                      const std::vector<Mineral>& minerals) = 0;
  virtual void AddContributionToResidual(const std::vector<Mineral>& minerals,
                                         const double por_den_sat_vol,
                                         std::vector<double> *residual) = 0;
  virtual void AddContributionToJacobian(const SpeciesArray& primary_species,
                                         const std::vector<Mineral>& minerals,
                                         const double por_den_sat_vol,
                                         Block* J) = 0;
  virtual void Display(void) const = 0;

  virtual void ParseParameters(const StringTokenizer& rate_parameters) = 0;

  void SetSpeciesIds(const SpeciesArray& species,
                     const std::string& species_type,
                     const std::vector<SpeciesName>& in_names,
                     const std::vector<double>& in_stoichiometry,
                     std::vector<SpeciesId>* out_ids,
                     std::vector<double>* out_stoichiometry);

  void DisplayReaction(void) const;

  void set_verbosity(const Verbosity s_verbosity) {
    this->verbosity_ = s_verbosity;
  };
  Verbosity verbosity(void) const {
    return this->verbosity_;
  };
  std::string name(void) const {
    return this->name_;
  };
  SpeciesId identifier(void) const {
    return this->identifier_;
  };

 protected:
  KineticRate(void);

  void set_name(const std::string in_name) {
    this->name_ = in_name;
  };
  void set_identifier(const SpeciesId in_id) {
    this->identifier_ = in_id;
  };

  std::vector<SpeciesName> reactant_names;
  std::vector<double> reactant_stoichiometry;
  std::vector<SpeciesId> reactant_ids;

 private:
  Verbosity verbosity_;
  std::string name_;
  SpeciesId identifier_;
};

}  // namespace chemistry
}  // namespace amanzi
#endif     /* AMANZI_CHEMISTRY_KINETIC_RATE_HH_ */
