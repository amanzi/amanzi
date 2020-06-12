/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Implementation of the TST rate law for mineral kinetics
 
    R = k * A * Prod (a_i^m_i) * ( 1 - Q/Keq)
 
  where:
    R : reaction rate, [moles/sec]
    Keq : equilibrium constant, [-]
    Q : ion activity product, [-]
    k : reaction rate constant, [moles m^2 s^-1]
    A : reactive surface area, [m^2]
    a_i : activity of modifying species (primary or secondary)
    m_i : exponent of modifying species
*/

#ifndef AMANZI_CHEMISTRY_KINETIC_RATE_TST_HH_
#define AMANZI_CHEMISTRY_KINETIC_RATE_TST_HH_

#include <vector>

#include "VerboseObject.hh"

#include "kinetic_rate.hh"
#include "species.hh"
#include "secondary_species.hh"
#include "mineral.hh"
#include "string_tokenizer.hh"

namespace Amanzi {
namespace AmanziChemistry {

class MatrixBlock;

class KineticRateTST : public KineticRate {
 public:
  KineticRateTST(void);
  virtual ~KineticRateTST(void) {};

  void Setup(const SecondarySpecies& reaction,
             const StringTokenizer& reaction_data,
             const SpeciesArray& primary_species);
  void Update(const SpeciesArray& primary_species,
              const std::vector<Mineral>& minerals);
  void AddContributionToResidual(const std::vector<Mineral>& minerals,
                                 const double bulk_volume,
                                 std::vector<double> *residual);

  void AddContributionToJacobian(const SpeciesArray& primary_species,
                                 const std::vector<Mineral>& minerals,
                                 const double bulk_volume_vol,
                                 MatrixBlock* J);
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;

  void ParseParameters(const StringTokenizer& rate);

  /*
  ** end of KineticRate inherited interface
  */

 protected:
  void area(double set_area) {
    this->area_ = set_area;
  };
  double area(void) const {
    return this->area_;
  };
  void log_Keq(double set_log_Keq) {
    this->log_Keq_ = set_log_Keq;
  };
  double log_Keq(void) const {
    return this->log_Keq_;
  };

  void rate_constant(double rate_constant) {
    this->rate_constant_ = rate_constant;
  };
  double rate_constant(void) const {
    return this->rate_constant_;
  };

  void log10_rate_constant(double log10_k) {
    this->log10_rate_constant_ = log10_k;
  };
  double log10_rate_constant(void) const {
    return this->log10_rate_constant_;
  };

  void sat_state_exponent(double set_sat_state_exponent) {
    this->sat_state_exponent_ = set_sat_state_exponent;
  };
  double sat_state_exponent(void) const {
    return this->sat_state_exponent_;
  };

  void Q_over_Keq(const double QK) {
    this->Q_over_Keq_ = QK;
  };
  double Q_over_Keq(void) const {
    return this->Q_over_Keq_;
  };

  void modifying_term(const double mod) {
    this->modifying_term_ = mod;
  };
  double modifying_term(void) const {
    return this->modifying_term_;
  };

 private:
  double area_;  // surface area [m^2]
  double log_Keq_;  // log_Keq [-]
  double rate_constant_;  // k, rate constant, [moles/m^2/sec]
  double log10_rate_constant_;  // log10(k),
  double sat_state_exponent_;  // n, saturation state exponent, [-]

  double Q_over_Keq_;
  double modifying_term_;

  // length is the number of modifying species (primary and secondary
  // in the same list!)
  std::vector<SpeciesName> modifying_species_names;
  std::vector<double> modifying_exponents;
  std::vector<int> modifying_primary_ids;  // not needed?
  std::vector<int> modifying_secondary_ids;  // not needed?

  // these vectors have length primary_species.size() so that dot
  // products are easier. any unneeded values are set to zero.
  std::vector<double> primary_stoichiometry;
  std::vector<double> modifying_primary_exponents;
  // length secondary_species.size()
  std::vector<double> modifying_secondary_exponents;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
