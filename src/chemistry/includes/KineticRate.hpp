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
#include "StringTokenizer.hpp"
#include "Verbosity.hpp"

class Block; 

class KineticRate
{
 public:
  ~KineticRate(void);

  virtual void Setup(const std::string reaction, 
                     const StringTokenizer reaction_data,
                     const std::vector<Species> primary_species) = 0;
  virtual void Update(const std::vector<Species> primary_species) = 0;
  virtual void AddContributionToResidual(const double por_den_sat_vol,
                                         std::vector<double> *residual) = 0;
  virtual void AddContributionToJacobian(const std::vector<Species> primary_species,
                                         const double por_den_sat_vol,
                                         Block *J) = 0;
  virtual void Display(void) const = 0;

  virtual void ParseParameters(const StringTokenizer rate) = 0;

  void ParseReaction(const std::string rxn_string);

  void DisplayReaction(void) const;

  void verbosity(const Verbosity s_verbosity) { this->verbosity_ = s_verbosity; };
  Verbosity verbosity(void) const { return this->verbosity_; };

 protected:
  KineticRate(void);

  std::vector<SpeciesName> reactant_names;
  std::vector<double> reactant_stoichiometery;
  std::vector<SpeciesId> reactant_ids;
  std::vector<SpeciesName> product_names;
  std::vector<double> product_stoichiometery;
  std::vector<SpeciesId> product_ids;

 private:
  Verbosity verbosity_;
};

#endif     /* __KINETIC_RATE_HPP__ */

