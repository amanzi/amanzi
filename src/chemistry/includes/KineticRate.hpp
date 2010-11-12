/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __KINETIC_RATE_HPP__

#define __KINETIC_RATE_HPP__

/*******************************************************************************
**
**  Description: abstract base class for all kinetic rates
**
*******************************************************************************/
#include <vector>

#include "Species.hpp"
#include "SecondarySpecies.hpp"
#include "StringTokenizer.hpp"
#include "Verbosity.hpp"

class Block; 

class KineticRate
{
 public:
  ~KineticRate(void);

  virtual void Setup(const SecondarySpecies& reaction,
                     const StringTokenizer reaction_data,
                     const SpeciesArray primary_species) = 0;
  virtual void Update(const SpeciesArray primary_species) = 0;
  virtual void AddContributionToResidual(const double por_den_sat_vol,
                                         std::vector<double> *residual) = 0;
  virtual void AddContributionToJacobian(const SpeciesArray primary_species,
                                         const double por_den_sat_vol,
                                         Block *J) = 0;
  virtual void Display(void) const = 0;

  virtual void ParseParameters(const StringTokenizer rate) = 0;

  void SetSpeciesIds(const SpeciesArray species,
                     const std::string species_type,
                     const std::vector<SpeciesName> in_names,
                     const std::vector<double> in_stoichiometry,
                     std::vector<SpeciesId>* out_ids,
                     std::vector<double>* out_stoichiometry);

  void DisplayReaction(void) const;

  void verbosity(const Verbosity s_verbosity) { this->verbosity_ = s_verbosity; };
  Verbosity verbosity(void) const { return this->verbosity_; };

  void name(const std::string in_name) { this->name_ = in_name; };
  std::string name(void) const { return this->name_; };

  void identifier(const SpeciesId in_id) { this->identifier_ = in_id; };
  SpeciesId identifier(void) const { return this->identifier_; };

 protected:
  KineticRate(void);

  std::vector<SpeciesName> reactant_names;
  std::vector<double> reactant_stoichiometry;
  std::vector<SpeciesId> reactant_ids;

 private:
  Verbosity verbosity_;
  std::string name_;
  SpeciesId identifier_;
};

#endif     /* __KINETIC_RATE_HPP__ */

