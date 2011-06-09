/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*******************************************************************************
**
**  Description: implementation of the TST rate law for mineral kinetics
**
**  R = k * A * Prod_m (a_m^{mu_m}) * (1 - Q/Keq)
**
**  Q = Prod_p (a_p^{nu_p})
**
**  where:
**    R : reaction rate, [moles/sec]
**    Keq : equilibrium constant, [-]
**    Q : ion activity product, [-]
**    a_p : activity of primary species
**    nu_p : stoichiometric coefficient of primary species
**    k : reaction rate constant, [moles m^2 s^-1]
**    A : reactive surface area, [m^2]
**    a_m : activity of modifying species
**    mu_m : exponent of modifying species
**
**  Residual:
**
**    r_i = nu_i R
**
**  Jacobian contributions:
**
**  dR/dC_j = k * A * (da_j/dC_j) * 
**           ( (1-Q/Keq) * (mu_j * a_j^{mu_j - 1}) * Prod_{m!=j}(a_m^{mu_m}) - 
**             Prod_{m}(a_m^{mu_m})/Keq * (nu_j * a_j^{nu_j - 1}) * Prod_{p!=j}(a_p^{nu_p}) )
**
**  J_ij = nu_i * dR/dCj
**
**
**  Notes:
**
**    - R is the dissolution rate for the mineral, so positive when
**    dissolution, negative when precipitation. dCalcite/dt = -R, dCa/dt = R
**
**    - assume that the "reactants" list as defined in KineticRate
**    consists solely of the mineral and the coefficient is always
**    one, so they are ignored in the ion activity product. The
**    mineral is only accounted for in the accumulation / consumption
**    of solid mass.
**
**    - assume that the stoichiometry, nu, and exponents, mu, are zero
**    for any species that does not take part in the reaction!
**
**    - assume that the mineral reactions are always written so the
**    products consist of only primary species.
**
**    - assume that the modifying species are only primary species
**
**    - ignoring the (1-(Q/Keq)^n) and (1-(Q/Keq)^n)^p forms for now.
**
**    - TODO: need to calculate the area (DONE) or read from the file.
**
**    - TODO: DONE: need to obtain log_Keq from the mineral object rather than read from a file. 
**
**    - TODO: DONE: units of A and k are not consistent with the input
**    file. need to pick a set!
**
**    - TODO: where should the mineral mass get updated at....?
**
*******************************************************************************/
#include <cmath>
#include <cstdlib>

#include <iostream>

#include "secondary_species.hh"
#include "kinetic_rate_tst.hh"
#include "block.hh"
#include "string_tokenizer.hh"
#include "verbosity.hh"

KineticRateTST::KineticRateTST(void)
    : KineticRate(),
      area_(0.0),
      log_Keq_(0.0),
      rate_constant_(0.0),
      log10_rate_constant_(0.0),
      sat_state_exponent_(0.0),
      Q_over_Keq_(1.0),
      modifying_term_(1.0)
{
  this->modifying_species_names.clear();
  this->modifying_exponents.clear();
  this->modifying_primary_ids.clear();
  this->modifying_secondary_ids.clear();

  this->primary_stoichiometry.clear();
  this->modifying_primary_exponents.clear();

  this->modifying_secondary_exponents.clear();
}  // end KineticRateTST constructor

KineticRateTST::~KineticRateTST(void)
{
}  // end KineticRateTST destructor


void KineticRateTST::Setup(const SecondarySpecies& reaction,
                           const StringTokenizer& reaction_data,
                           const SpeciesArray& primary_species)
{
  // break the reaction string into reactants and products
  set_name(reaction.name());
  set_identifier(reaction.identifier());

  // copy the reactant species, ids and stoichiometry from the reaction species
  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "  KineticRateTST::Setup(): Searching for reactant species ids..." 
              << std::endl;
  }
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
  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "  KineticRateTST::Setup(): Searching for modifying species ids..." 
              << std::endl;
  }
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
}  // end Setup()



void KineticRateTST::Update(const SpeciesArray& primary_species,
                            const std::vector<Mineral>&  minerals)
{
  // update surace area from the minerals
  area(minerals.at(identifier()).surface_area());
  log_Keq(minerals.at(identifier()).logK());
//   std::cout << "area : " << area() << "  ";
//   std::cout << std::endl;
  // calculate the Q/K term
  double lnQ = 0.0;
  for (unsigned int p = 0; p < primary_species.size(); p++) {
    lnQ += primary_stoichiometry.at(p) * primary_species.at(p).ln_activity();
    if (verbosity() == kDebugMineralKinetics) {
      std::cout << "  Update: p: " << p << "  coeff: " << primary_stoichiometry.at(p)
                << "  ln_a: " << primary_species.at(p).ln_activity() << std::endl;
    }
  }
  double Q = std::exp(lnQ);
  double Keq = std::pow(10.0, log_Keq());
  Q_over_Keq(Q/Keq);

  // calculate the modifying primary species term:
  double ln_mod_term = 0.0;
  for (unsigned int p = 0; p < primary_species.size(); p++) {
    ln_mod_term += modifying_primary_exponents.at(p) * primary_species.at(p).ln_activity();
  }
  modifying_term(std::exp(ln_mod_term));

  // calculate the modifying secondary species term:

  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "  Update:\n"
              << "  Rate: " << name()
              << "\n  lnQ = " << lnQ << "\n   Q = " << Q 
              << "\n   Keq = " << Keq 
              << "\n  Q/K = " << Q_over_Keq() 
              << "\n  lnQK = " << std::log(Q_over_Keq()) << "\n"
              << std::scientific
              << "  area: " << area() << " [m^2]\n"
              << "  rate_constant: " << rate_constant() << " [moles/m^2/s]\n"
              << "  1-Q/K: " << 1.0-Q_over_Keq() << " [--]\n"
              << std::fixed
              << std::endl;
  }
}  // end Update()

void KineticRateTST::AddContributionToResidual(const std::vector<Mineral>&  minerals,
                                               const double por_den_sat_vol, 
                                               std::vector<double> *residual)
{
  /*
  ** NOTE: residual has units of moles/sec
  */
  static_cast<void>(por_den_sat_vol);

  // Calculate saturation state term: 1-Q/K 
  double sat_state = 1.0 - Q_over_Keq();

  // Calculate overall rate 
  double rate = area() * rate_constant() * modifying_term() * sat_state;

  // only add contribution if precipitating or mineral exists
  if (minerals.at(identifier()).volume_fraction() > 0. || sat_state < 0.) {

    // add or subtract from the residual....
    for (unsigned int p = 0; p < reactant_stoichiometry.size(); p++) {
      int reactant_id = reactant_ids.at(p);
      (*residual)[reactant_id] -= reactant_stoichiometry.at(p) * rate;
      //if (true) {
      if (verbosity() == kDebugMineralKinetics) {
        std::cout << "  AddToResidual p: " << p
                << "  coeff: " << reactant_stoichiometry.at(p)
                << "  rate: " << rate << "  redsidual: " << residual->at(p) << std::endl;
      }
    }
  }

  // TODO: updating the mineral mass.....

}  // end AddContributionToResidual()

void KineticRateTST::AddContributionToJacobian(const SpeciesArray& primary_species,
                                               const std::vector<Mineral>&  minerals,
                                               const double por_den_sat_vol,
                                               Block *J)
{
  /*
  ** Evaluate the dR/dC terms for this rate and add to the appropriate
  ** location in the jacobian. See file description above for the jacobian.
  **
  ** NOTE: jacobian has units of kg_water/sec, dR/dC = (moles/s) / (moles/kg_water)
  */
  static_cast<void>(por_den_sat_vol);

  // double dadC = 1.0;  // da_i/dC_j = nu_ij * m_i/m_j ; da_j/dC_j = m_j/m_j = 1

  double Keq = std::pow(10.0, log_Keq());
  double one_minus_QK = 1.0 - Q_over_Keq();  // (1-Q/Keq)
  double area_rate_constant = area() * rate_constant();  // k*A

  // only add contribution if precipitating or mineral exists
  if (minerals.at(identifier()).volume_fraction() > 0. || one_minus_QK < 0.) {

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
      dRdC_row[p] = -dRdC; // where does the neg sign come from...?
      if (verbosity() == kDebugMineralKinetics) {
        std::cout << "J_row_contrib: p: " << p
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
        
      }
    }

    // J_ij = nu_i * dR/dCj
    for (int i = 0; i < J->getSize(); i++) {
      for (int j = 0; j < J->getSize(); j++) {
        J->addValue(i, j, dRdC_row.at(j) * primary_stoichiometry.at(i));
          //std::cout << dRdC_row.at(j) * primary_stoichiometry.at(i) << " ";
      }
      //std::cout << std::endl;
    }
  }

}  // end AddContributionToJacobian()

void KineticRateTST::ParseParameters(const StringTokenizer& reaction_data)
{
  std::string space(" ");
  StringTokenizer st;

  std::vector<std::string>::const_iterator field = reaction_data.begin();

  for (; field != reaction_data.end(); field++) {
    st.tokenize(*field, space);

    if (st.at(0) == "log10_rate_constant") {
      log10_rate_constant(std::atof(st.at(1).c_str()));
      rate_constant(std::pow(10.0, log10_rate_constant()));
    }
    else {
      // assume we are dealing with the list of rate modifying species
      for (unsigned int modifier = 0; modifier < st.size(); modifier++) {
        this->modifying_species_names.push_back(st.at(modifier));
        modifier++; // increment to get the exponent of this modifier
        this->modifying_exponents.push_back(std::atof(st.at(modifier).c_str()));
      }      
    }
    if (0) {
      std::cout << "  Field: " << *field << std::endl;
    }
  }
}  // end ParseParameters()

void KineticRateTST::Display(void) const
{
  std::cout << "    Rate law: TST" << std::endl;
  this->DisplayReaction();
  std::cout << "    Parameters:" << std::endl;
  std::cout << "      mineral = " << name() << std::endl;
  std::cout << "      mineral id = " << identifier() << std::endl;
  std::cout << "      log10_rate constant = " << log10_rate_constant() << std::endl;
  std::cout << "      rate constant = " << std::scientific << rate_constant() << std::endl;
  std::cout << "      rate modifiers: " << std::endl;
  std::cout << "        ";
  for (unsigned int mod = 0; mod < this->modifying_species_names.size(); mod++) {
    std::cout << "{ " << this->modifying_species_names.at(mod) << " }";
    std::cout << "^" << this->modifying_exponents.at(mod) << " " << std::endl;
  }
  std::cout << std::endl;
}  // end Display()
