/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Implementation of the TST rate law for mineral kinetics
 
    R = k * A * Prod (a_i^m_i) * ( 1 - Q/Keq)
 
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
#include <sstream>
#include <iostream>
#include <sstream>

// TPLs
#include "boost/algorithm/string.hpp"

// Chemistry
#include "secondary_species.hh"
#include "matrix_block.hh"
#include "string_tokenizer.hh"
#include "chemistry_utilities.hh"
#include "chemistry_exception.hh"

#include "kinetic_rate_tst.hh"

namespace Amanzi {
namespace AmanziChemistry {

KineticRateTST::KineticRateTST(void)
    : KineticRate(),
      area_(0.0),
      log_Keq_(0.0),
      rate_constant_(0.0),
      log10_rate_constant_(0.0),
      sat_state_exponent_(0.0),
      Q_over_Keq_(1.0),
      modifying_term_(1.0) {
  this->modifying_species_names.clear();
  this->modifying_exponents.clear();
  this->modifying_primary_ids.clear();
  this->modifying_secondary_ids.clear();

  this->primary_stoichiometry.clear();
  this->modifying_primary_exponents.clear();

  this->modifying_secondary_exponents.clear();
}


void KineticRateTST::Setup(const SecondarySpecies& reaction,
                           const StringTokenizer& reaction_data,
                           const SpeciesArray& primary_species) {
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
  ParseParameters(reaction_data);

  // determine the species ids for the modifying species
  SetSpeciesIds(primary_species, species_type,
                modifying_species_names, modifying_exponents,
                &modifying_primary_ids, &modifying_primary_exponents);


  //   std::cout << std::endl;
  //   for (unsigned int i = 0; i < modifying_primary_ids.size(); i++) {
  //     std::cout << "  Modifier: " << std::endl;
  //     std::cout << "        id: " << modifying_primary_ids.at(i) << std::endl;
  //     std::cout << "      name: " << primary_species.at(modifying_primary_ids.at(i)).name() << std::endl;
  //     std::cout << "    coeff: " << modifying_primary_exponents.at(i) << std::endl;
  //   }
}


void KineticRateTST::Update(const SpeciesArray& primary_species,
                            const std::vector<Mineral>&  minerals) {
  // update surace area from the minerals
  area(minerals.at(identifier()).specific_surface_area());
  log_Keq(minerals.at(identifier()).logK());
  //   std::cout << "area : " << area() << "  ";
  //   std::cout << std::endl;
  // calculate the Q/K term
  double lnQ = 0.0;
  for (unsigned int p = 0; p < primary_species.size(); p++) {
    lnQ += primary_stoichiometry.at(p) * primary_species.at(p).ln_activity();
    /*
    if (debug()) {
      std::stringstream message;
      message << "  Update: p: " << p << "  coeff: " << primary_stoichiometry.at(p)
              << "  ln_a: " << primary_species.at(p).ln_activity() << std::endl;
      vo->Write(Teuchos::VERB_EXTREME, message);
    }
    */
  }
  double Q = std::exp(lnQ);
  double Keq = std::pow(10.0, log_Keq());
  Q_over_Keq(Q / Keq);

  // calculate the modifying primary species term:
  double ln_mod_term = 0.0;
  for (unsigned int p = 0; p < primary_species.size(); p++) {
    ln_mod_term += modifying_primary_exponents.at(p) * primary_species.at(p).ln_activity();
  }
  modifying_term(std::exp(ln_mod_term));

  // calculate the modifying secondary species term:
  /*
  if (debug()) {
    std::stringstream message;
    message << "  Update:\n"
            << "    Rate: " << name()
            << "\n    lnQ = " << lnQ 
            << "\n     Q = " << Q
            << "\n     Keq = " << Keq
            << "\n    Q/K = " << Q_over_Keq()
            << "\n    lnQK = " << std::log(Q_over_Keq()) << "\n"
            << std::scientific
            << "    area: " << area() << " [m^2]\n"
            << "    rate_constant: " << rate_constant() << " [moles/m^2/s]\n"
            << "    1-Q/K: " << 1.0 - Q_over_Keq() << " [--]\n"
            << std::fixed
            << std::endl;
    vo->Write(Teuchos::VERB_EXTREME, message);
  }
  */
}


void KineticRateTST::AddContributionToResidual(const std::vector<Mineral>&  minerals,
                                               const double bulk_volume,
                                               std::vector<double> *residual) {
  /*
  ** NOTE: residual has units of moles/sec
  */

  // Calculate saturation state term: 1-Q/K
  double sat_state = 1.0 - Q_over_Keq();

  // Calculate overall rate [moles/s]
  // area = [m^2 mnrl / m^3 bulk]
  // volume = [m^3 bulk]
  // rate constant = [moles / m^2 mnrl / s]
  // modifing term = [-]
  // saturation state = [-]
  double rate = area() * bulk_volume * rate_constant() * modifying_term() * sat_state;

  // only add contribution if precipitating or mineral exists
  if (minerals.at(identifier()).volume_fraction() > 0.0 || sat_state < 0.0) {
    // add or subtract from the residual....
    for (unsigned int p = 0; p < reactant_stoichiometry.size(); p++) {
      int reactant_id = reactant_ids.at(p);
      residual->at(reactant_id) -= reactant_stoichiometry.at(p) * rate;
      /*
      if (debug()) {
        std::stringstream message;
        message << "  AddToResidual p: " << p
                << "  coeff: " << reactant_stoichiometry.at(p)
                << "  rate: " << rate / bulk_volume << "  redsidual: " << residual->at(p) << std::endl;
        vo->Write(Teuchos::VERB_EXTREME, message);
      }
      */
    }
  }
  // store the reaction rate so we can use it to update volume fractions
  // need [moles/sec/m^3 bulk], so we divide by volume!
  set_reaction_rate(rate / bulk_volume);
}


void KineticRateTST::AddContributionToJacobian(const SpeciesArray& primary_species,
                                               const std::vector<Mineral>&  minerals,
                                               const double bulk_volume,
                                               MatrixBlock* J) {
  /*
  ** Evaluate the dR/dC terms for this rate and add to the appropriate
  ** location in the jacobian. See file description above for the jacobian.
  **
  ** NOTE: jacobian has units of kg_water/sec, dR/dC = (moles/s) / (moles/kg_water)
  */

  // double dadC = 1.0;  // da_i/dC_j = nu_ij * m_i/m_j ; da_j/dC_j = m_j/m_j = 1

  double Keq = std::pow(10.0, log_Keq());
  double one_minus_QK = 1.0 - Q_over_Keq();  // (1-Q/Keq)
  double area_rate_constant = area() * bulk_volume * rate_constant();  // k*A

  // only add contribution if precipitating or mineral exists
  if (minerals.at(identifier()).volume_fraction() > 0.0 || one_minus_QK < 0.0) {
    // the contribution to every row of the jacobian has the same dR/dci
    // term. Each colum scaled by the stoichiometric
    // coefficient.
    std::vector<double> dRdC_row(primary_species.size());  // primary_species.size() == J.getSize()!

    for (unsigned int p = 0; p < primary_species.size(); p++) {
      // temp_modifying_term = Prod_{m!=j}(a_m^{mu_m})
      double temp_modifying_term = modifying_term();
      double remove_modifying = modifying_primary_exponents.at(p) *
          primary_species.at(p).ln_activity();
      remove_modifying = std::exp(remove_modifying);
      temp_modifying_term /= remove_modifying;

      // modifying_deriv = (mu_j * a_j^{mu_j - 1})
      double modifying_deriv = (modifying_primary_exponents.at(p) - 1.0) *
          primary_species.at(p).ln_activity();
      modifying_deriv = std::exp(modifying_deriv) * modifying_primary_exponents.at(p);

      // temp_Q = Prod_{p!=j}(a_p^{nu_p}
      double temp_Q = primary_stoichiometry.at(p) * primary_species.at(p).ln_activity();
      temp_Q = std::exp(temp_Q);
      temp_Q = Q_over_Keq() * Keq / temp_Q;

      // primary_deriv = (nu_j * a_j^{nu_j - 1})
      double primary_deriv = (primary_stoichiometry.at(p) - 1.0) *
          primary_species.at(p).ln_activity();
      primary_deriv = std::exp(primary_deriv) * primary_stoichiometry.at(p);

      // modifying_term() / Keq = Prod_{m}(a_m^{mu_m})/Keq *

      double dRdC = area_rate_constant *
          (one_minus_QK * modifying_deriv * temp_modifying_term -
           (modifying_term() / Keq) * primary_deriv * temp_Q);
      dRdC_row[p] = -dRdC;  // where does the neg sign come from...?
      /*
      if (debug()) {
        std::stringstream message;
        message << "J_row_contrib: p: " << p
                << std::scientific << "\tA*k: " << area_rate_constant
                << "\t1-Q/K: " << one_minus_QK
                << "\tmod_deriv: " << modifying_deriv
                << "\ttemp_mod_term: " << temp_modifying_term
                << "\tmod_term: " << modifying_term()
                << "\tKeq: " << Keq
                << "\tmod_term/Keq: " << modifying_term() / Keq
                << "\tprimary_deriv: " << primary_deriv
                << "\ttemp_Q: " << temp_Q
                << "\trow: " << dRdC_row.at(p)
                << std::endl;
        vo->Write(Teuchos::VERB_EXTREME, message);
      }
      */
    }

    // J_ij = nu_i * dR/dCj
    for (int i = 0; i < J->size(); i++) {
      for (int j = 0; j < J->size(); j++) {
        J->AddValue(i, j, dRdC_row.at(j) * primary_stoichiometry.at(i));
        // std::cout << dRdC_row.at(j) * primary_stoichiometry.at(i) << " ";
      }
    }
  }
}


void KineticRateTST::ParseParameters(const StringTokenizer& reaction_data) {
  std::string space(" ");
  StringTokenizer st;

  std::vector<std::string>::const_iterator field = reaction_data.begin();

  for (; field != reaction_data.end(); field++) {
    st.tokenize(*field, space);

    if (st.at(0) == "log10_rate_constant") {
      double value = std::atof(st.at(1).c_str());
      // what units do we have [moles/cm^2/sec] or [moles/m^2/sec]? We
      // need to set [moles/m^2/sec]!
      std::string units = st.at(2);
      if (boost::iequals(units, "moles/cm^2/sec")) {
        // add 4 in log10 space to convert cm^-2 m^-2
        value += 4.0;
      } else if (boost::iequals(units, "moles/m^2/sec")) {
        // no change
      } else {
        std::ostringstream error_stream;
        error_stream << "KineticRateTST::ParseParameters(" << name()
                     << "): unknown units for parameter 'log10_rate_constant': '"
                     << st.at(2) 
                     << "'. Valid units are 'moles/m^2/sec' and 'moles/cm^2/sec'.\n";
        Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
      }
      log10_rate_constant(value);
      rate_constant(std::pow(10.0, log10_rate_constant()));
    } else if (st.at(0) == "rate_constant") {
      std::ostringstream error_stream;
      error_stream << "KineticRateTST::ParseParameters(" << name()
                   << "): use of parameter 'rate_constant' has not been "
                   << "implemented yet. Please use 'log10_rate_constant'.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));

    } else {
      // assume we are dealing with the list of rate modifying species
      for (unsigned int modifier = 0; modifier < st.size(); modifier++) {
        this->modifying_species_names.push_back(st.at(modifier));
        modifier++;  // increment to get the exponent of this modifier
        this->modifying_exponents.push_back(std::atof(st.at(modifier).c_str()));
      }
    }
    if (0) {
      std::cout << "  Field: " << *field << std::endl;
    }
  }
}


void KineticRateTST::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  vo->Write(Teuchos::VERB_HIGH, "    Rate law: TST\n");
  this->DisplayReaction(vo);
  std::stringstream message;
  message << "    Parameters:" << std::endl;
  message << "      mineral = " << name() << std::endl;
  message << "      mineral id = " << identifier() << std::endl;
  message << "      log10_rate constant = " << log10_rate_constant() << std::endl;
  message << "      rate constant = " << std::scientific << rate_constant() << std::endl;
  message << "      rate modifiers: " << std::endl;
  for (unsigned int mod = 0; mod < this->modifying_species_names.size(); mod++) {
    message << "        ";
    message << "{ " << this->modifying_species_names.at(mod) << " }";
    message << "^" << this->modifying_exponents.at(mod) << " " << std::endl;
  }
  message << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
