/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Implementation of the TST rate law for mineral kinetics
 
    R = k * A * Prod (a_i^m_i) * (1 - Q/Keq)
 
    Q = Prod_p (a_p^{nu_p})
 
  where:
    R : reaction rate, [moles/sec]
    Keq : equilibrium constant, [-]
    Q : ion activity product, [-]
    a_p : activity of primary species
    nu_p : stoichiometric coefficient of primary species
    k : reaction rate constant, [moles m^2 s^-1]
    A : reactive surface area, [m^2]
    a_m : activity of modifying species
    mu_m : exponent of modifying species

  Residual:
    r_i = nu_i R
 
  Jacobian contributions:
    dR/dC_j = k * A * (da_j/dC_j) *
           ( (1-Q/Keq) * (mu_j * a_j^{mu_j - 1}) * Prod_{m!=j}(a_m^{mu_m}) -
           Prod_{m}(a_m^{mu_m})/Keq * (nu_j * a_j^{nu_j - 1}) * Prod_{p!=j}(a_p^{nu_p}) )

    J_ij = nu_i * dR/dCj
 
  Notes:
   - R is the dissolution rate for the mineral, so positive when
     dissolution, negative when precipitation. dCalcite/dt = -R, dCa/dt = R
 
   - assume that the "reactants" list as defined in KineticRate
     consists solely of the mineral and the coefficient is always
     one, so they are ignored in the ion activity product. The
     mineral is only accounted for in the accumulation / consumption
     of solid mass.
 
   - assume that the stoichiometry, nu, and exponents, mu, are zero
     for any species that does not take part in the reaction!
 
   - assume that the mineral reactions are always written so the
     products consist of only primary species.
 
   - assume that the modifying species are only primary species
 
   - ignoring the (1-(Q/Keq)^n) and (1-(Q/Keq)^n)^p forms for now.

   - TODO(bandre): need to calculate the area (DONE) or read from the file.

   - TODO(bandre): DONE: need to obtain log_Keq from the mineral object rather than read from a file.
 
   - TODO(bandre): DONE: units of A and k are not consistent with the input
     file. need to pick a set!

   - TODO(bandre): where should the mineral mass get updated at....?
*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

// TPLs
#include "boost/algorithm/string.hpp"

// Chemistry
#include "ChemistryUtilities.hh"
#include "MatrixBlock.hh"
#include "SecondarySpecies.hh"

#include "KineticRateTST.hh"

namespace Amanzi {
namespace AmanziChemistry {

KineticRateTST::KineticRateTST()
  : KineticRate(),
    area_(0.0),
    log_Keq_(0.0),
    rate_constant_(0.0),
    log10_rate_constant_(0.0),
    Q_over_Keq_(1.0),
    modifying_term_(1.0) 
{};


void KineticRateTST::Setup(const Mineral& reaction,
                           double rate,
                           const std::string& modifiers,
                           const SpeciesArray& primary_species)
{
  // break the reaction string into reactants and products
  set_name(reaction.name());
  set_identifier(reaction.identifier());

  // copy the reactant species, ids and stoichiometry from the reaction species
  reactant_names = reaction.species_names();
  reactant_stoichiometry = reaction.stoichiometry();
  reactant_ids = reaction.species_ids();

  std::string species_type("primary");
  SetSpeciesIds(primary_species, species_type,
                reactant_names, reactant_stoichiometry,
                &reactant_ids, &primary_stoichiometry);

  // extract the rate parameters from the data string, including the
  // modifying species data. adds data to the mineral species!
  log10_rate_constant_ = rate;  // SI units!
  rate_constant_ = std::pow(10.0, log10_rate_constant_);

  double coeff;
  std::string mod_name;
  std::istringstream iss(modifiers);
  while (iss >> mod_name || !iss.eof()) {
    iss >> coeff;
    modifying_species_names.push_back(mod_name);
    modifying_exponents.push_back(coeff);
  }

  // determine the species ids for the modifying species
  SetSpeciesIds(primary_species, species_type,
                modifying_species_names, modifying_exponents,
                &modifying_primary_ids, &modifying_primary_exponents);
}


/* ******************************************************************
* TBW
****************************************************************** */
void KineticRateTST::Update(const SpeciesArray& primary_species,
                            const std::vector<Mineral>& minerals)
{
  // update surace area from the minerals
  area_ = minerals.at(identifier()).specific_surface_area();
  log_Keq_ = minerals.at(identifier()).logK();

  // calculate the Q/K term
  double lnQ = 0.0;
  for (int p = 0; p < primary_species.size(); p++) {
    lnQ += primary_stoichiometry.at(p) * primary_species.at(p).ln_activity();
  }
  double Q = std::exp(lnQ);
  double Keq = std::pow(10.0, log_Keq_);
  Q_over_Keq_ = Q / Keq;

  // calculate the modifying primary species term:
  double ln_mod_term = 0.0;
  for (int p = 0; p < primary_species.size(); p++) {
    ln_mod_term += modifying_primary_exponents.at(p) * primary_species.at(p).ln_activity();
  }
  modifying_term_ = std::exp(ln_mod_term);
}


/* ******************************************************************
* NOTE: residual has units of moles/sec
****************************************************************** */
void KineticRateTST::AddContributionToResidual(const std::vector<Mineral>& minerals,
                                               double bulk_volume,
                                               std::vector<double> *residual)
{
  // Calculate saturation state term: 1-Q/K
  double sat_state = 1.0 - Q_over_Keq_;

  // Calculate overall rate [moles/s]
  // area = [m^2 mnrl / m^3 bulk]
  // volume = [m^3 bulk]
  // rate constant = [moles / m^2 mnrl / s]
  // modifing term = [-]
  // saturation state = [-]
  double rate = area_ * bulk_volume * rate_constant_ * modifying_term_ * sat_state;

  // only add contribution if precipitating or mineral exists
  if (minerals.at(identifier()).volume_fraction() > 0.0 || sat_state < 0.0) {
    for (int p = 0; p < reactant_stoichiometry.size(); p++) {
      int id = reactant_ids.at(p);
      residual->at(id) -= reactant_stoichiometry.at(p) * rate;
    }
  }
  // store the reaction rate so we can use it to update volume fractions
  // need [moles/sec/m^3 bulk], so we divide by volume!
  set_reaction_rate(rate / bulk_volume);
}


/* ******************************************************************
* TBW
****************************************************************** */
void KineticRateTST::AddContributionToJacobian(const SpeciesArray& primary_species,
                                               const std::vector<Mineral>& minerals,
                                               double bulk_volume,
                                               MatrixBlock* J)
{
  /*
  ** Evaluate the dR/dC terms for this rate and add to the appropriate
  ** location in the jacobian. See file description above for the jacobian.
  **
  ** NOTE: jacobian has units of kg_water/sec, dR/dC = (moles/s) / (moles/kg_water)
  */

  // double dadC = 1.0;  // da_i/dC_j = nu_ij * m_i/m_j ; da_j/dC_j = m_j/m_j = 1

  double Keq = std::pow(10.0, log_Keq_);
  double one_minus_QK = 1.0 - Q_over_Keq_;  // (1-Q/Keq)
  double area_rate_constant = area_ * bulk_volume * rate_constant_;  // k*A

  // only add contribution if precipitating or mineral exists
  if (minerals.at(identifier()).volume_fraction() > 0.0 || one_minus_QK < 0.0) {
    // the contribution to every row of the jacobian has the same dR/dci
    // term. Each colum scaled by the stoichiometric coefficient.
    std::vector<double> dRdC_row(primary_species.size());  // primary_species.size() == J.getSize()!

    for (int p = 0; p < primary_species.size(); p++) {
      // temp_modifying_term = Prod_{m!=j}(a_m^{mu_m})
      double temp_modifying_term = modifying_term_;
      double remove_modifying = modifying_primary_exponents.at(p) * primary_species.at(p).ln_activity();
      remove_modifying = std::exp(remove_modifying);
      temp_modifying_term /= remove_modifying;

      // modifying_deriv = (mu_j * a_j^{mu_j - 1})
      double modifying_deriv = (modifying_primary_exponents.at(p) - 1.0) *
          primary_species.at(p).ln_activity();
      modifying_deriv = std::exp(modifying_deriv) * modifying_primary_exponents.at(p);

      // temp_Q = Prod_{p!=j}(a_p^{nu_p}
      double temp_Q = primary_stoichiometry.at(p) * primary_species.at(p).ln_activity();
      temp_Q = std::exp(temp_Q);
      temp_Q = Q_over_Keq_ * Keq / temp_Q;

      // primary_deriv = (nu_j * a_j^{nu_j - 1})
      double primary_deriv = (primary_stoichiometry.at(p) - 1.0) *
          primary_species.at(p).ln_activity();
      primary_deriv = std::exp(primary_deriv) * primary_stoichiometry.at(p);

      // modifying_term_ / Keq = Prod_{m}(a_m^{mu_m})/Keq *

      double dRdC = area_rate_constant *
          (one_minus_QK * modifying_deriv * temp_modifying_term -
           (modifying_term_ / Keq) * primary_deriv * temp_Q);
      dRdC_row[p] = -dRdC;  // where does the neg sign come from...?
    }

    // J_ij = nu_i * dR/dCj
    for (int i = 0; i < J->size(); i++) {
      for (int j = 0; j < J->size(); j++) {
        J->AddValue(i, j, dRdC_row.at(j) * primary_stoichiometry.at(i));
      }
    }
  }
}


/* ******************************************************************
* TBW
****************************************************************** */
void KineticRateTST::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  vo->Write(Teuchos::VERB_HIGH, "    Rate law: TST\n");
  DisplayReaction(vo);

  std::stringstream message;
  message << "    Parameters:" << std::endl;
  message << "      mineral = " << name() << std::endl;
  message << "      mineral id = " << identifier() << std::endl;
  message << "      log10_rate constant = " << log10_rate_constant_ << std::endl;
  message << "      rate constant = " << std::scientific << rate_constant_<< std::endl;
  message << "      rate modifiers: " << std::endl;
  for (int mod = 0; mod < modifying_species_names.size(); mod++) {
    message << "        ";
    message << "{ " << modifying_species_names.at(mod) << " }";
    message << "^" << modifying_exponents.at(mod) << " " << std::endl;
  }
  message << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
