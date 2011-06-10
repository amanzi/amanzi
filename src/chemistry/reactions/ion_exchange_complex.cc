/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
**  Description: Implementation of ion exchange reactions as
**  equilibrium half reactions. Using the activity of exchange complex
**  is the equivalent fraction of exchangable cations (exchange sites)
**
**  Assumptions:
**
**    - Assume that the reactions are always going to be written as
**    one complex equals one primary plus one exchange site:
**
**        1.0 complex = double primary + double exchange_site
**
**    - can we assume that the exchanger activity is always one...?
**
**  Examples:
**
**  Ca0.5X = 0.5Ca++ + X- ; Keq = double
**
**  CaX_2 = Ca++ + 2X- ; Keq = double
**
**  NaX = Na+ + X- ; Keq = double
**
**  Mass Action Law:
**
**    K_eq = a_p^nu_p * a_x^nu_x / a_c^nu_c
**
**    where:
**
**      subscript p = primary species
**
**      subscript x = exchange site
**
**      subscript c = ion exchange complex
**
**      a_ = activity of species
**
**      nu_ = stoichiometric coefficient of species
**
**      Keq = half reaction equilibrium constant
**
**    a_c = eq_c / CEC
**
**        = m_c * nu_x / CEC
**
**    where:
**
**      a_c = activity of the exchanger-cation = equivalent fraction of c
**
**      eq_c = equivalents of complex c, equivalents per MASS? VOLUME?
**             equivalents are equivalent number of exchange sites occupied
**             NaX, Ca_{0.5}X : eq_c = conc of complex, complex occupies one exchange site
**             CaX_2 : eq_c = 2*conc of complex, complex occupies two exchange sites.
**           = nu_x * m_c
**
**      CEC = cation exchange capacity, equivalents per MASS? VOLUME?
**
**  Residual:
**
**    R(a_p, a_x, a_c) = -ln_Keq +
**                       nu_p * ln(a_p) +
**                       nu_x * ln(a_x) -
**                       nu_c * (ln(m_c) + ln(nu_x) - ln(CEC))
**
**
*/

#include <iostream>
#include <iomanip>

#include "ion_exchange_complex.hh"
#include "block.hh"


IonExchangeComplex::IonExchangeComplex()
    : Species() {
}  // end IonExchangeComplex() constructor

IonExchangeComplex::IonExchangeComplex(
    const SpeciesName in_name,
    const SpeciesId in_id,
    const SpeciesName in_primary_name,
    const double in_primary_stoichiometry,
    const SpeciesId in_primary_id,
    const SpeciesName in_exchange_site_name,
    const double in_exchange_site_stoichiometry,
    const SpeciesId in_exchange_site_id,
    const double in_h2o_stoich,
    const double in_log_Keq)
    : Species(in_id, in_name,
              0.0, 0.0, 0.0),
      complex_stoich_coeff_(1.0),
      primary_name_(in_primary_name),
      primary_stoichiometry_(in_primary_stoichiometry),
      primary_id_(in_primary_id),
      exchange_site_name_(in_exchange_site_name),
      exchange_site_stoichiometry_(in_exchange_site_stoichiometry),
      exchange_site_id_(in_exchange_site_id),
      h2o_stoich_(in_h2o_stoich),
      log_Keq_(in_log_Keq),
      ln_Keq_(0.0),
      ln_QKeq_(0.0) {
  double Keq = std::pow(10.0, log_Keq());
  ln_Keq_ = std::log(Keq);
}  // end IonExchangeComplex costructor


IonExchangeComplex::~IonExchangeComplex() {
}  // end IonExchangeComplex() destructor


void IonExchangeComplex::update(void) {
  /* update function for species....
  **
  ** assumes that activity coefficient and molality have already been
  ** set externally so that the activity is correct. Ugly hack....
  **
  ** activity = equivalent fraction of complex
  **          = molality * stoich_coeff_exchanger / cation_exchange_capacity
  */

  activity_ = act_coef() / molality();
  ln_molality_ = std::log(molality());
  ln_act_coef_ = std::log(act_coef());
  ln_activity_ = std::log(activity());
}  // end Update()

void IonExchangeComplex::update(const double molality) {
  // dummy function
  static_cast<void>(molality);
}  // end update()


void IonExchangeComplex::Update(const std::vector<Species> primary_species,
                                const std::vector<IonExchangeSite> exchange_sites) {
  /* The Q/K calculated in this function is:
  **   cPX <==> pP + xX
  **   K = a_P^p * a_X^x / a_PX^c
  **   c = 1.0
  **   a_PX = Q/K = a_P^p * a_X^x / K
  **
  ** The activity coefficient of the exchange complex is
  **   a_PX = gamma_PX * m_PX
  **   a_PX = m_PX * x / CEC
  **   gamma_PX = x / CEC
  */

  double ln_QK = -ln_Keq();
  // + p * ln(a_P), standard aqueous activities for the primary
  ln_QK += primary_stoichiometry() * primary_species.at(primary_id()).ln_activity();

  // + x * ln(a_X)
  ln_QK += exchange_site_stoichiometry() * exchange_sites.at(exchange_site_id()).ln_activity();

  set_ln_QKeq(ln_QK);

  // activity coefficient
  // gamma_PB = x / CEC
  double CEC = exchange_sites.at(exchange_site_id()).cation_exchange_capacity();
  double activity_coefficient = exchange_site_stoichiometry() / CEC;
  act_coef(activity_coefficient);

  //  m_PX = QK / gamma_PX
  molality(std::exp(ln_QK) / act_coef());

  update();  // update the Species portion of this object!
}  // end Update()

void IonExchangeComplex::AddContributionToTotal(std::vector<double> &total) {
  // only the primary species are in total, so we only need to update a single entry!
  total[primary_id()] += primary_stoichiometry() * molality();
}  // end addContributionToTotal()

void IonExchangeComplex::AddContributionToDTotal(const std::vector<Species> primary_species,
                                                 Block* dtotal) {
  static_cast<void>(primary_species);
  static_cast<void>(dtotal);
}  // end addContributionToDTotal()


/*
**
**  Display functions
**
*/
void IonExchangeComplex::display(void) const {
  DisplayReaction();
  std::cout << "      log_Keq: " << log_Keq()
            << std::endl;
}  // end Display()

void IonExchangeComplex::Display(void) const {
  DisplayReaction();
  std::cout << std::setw(40) << " "
            << std::setw(10) << log_Keq()
            << std::endl;
}  // end Display()

void IonExchangeComplex::DisplayReaction(void) const {
  std::cout << "    " << name() << " = "
            << primary_stoichiometry() << " " << primary_name()
            << exchange_site_stoichiometry() << " "
            << exchange_site_name();
  std::cout << std::endl;
}  // end DisplayReaction()


void IonExchangeComplex::DisplayResultsHeader(void) const {
  std::cout << std::setw(15) << "Name"
            << std::setw(15) << "Molarity"
            << std::setw(15) << "Activity"
            << std::endl;
}  // end DisplayResultsHeader()

void IonExchangeComplex::DisplayResults(void) const {
  std::cout << std::setw(15) << name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << molality()
            << std::setw(15) << activity()
            << std::endl;
}  // end DisplayResults()



// void IonExchangeComplex::DisplayReaction_arrays(void) const
// {
//   std::cout << "    " << name() << " = ";
//   if (species_names_.size() > 0) {
//     for (unsigned int i = 0; i < species_names_.size(); i++) {
//       std::cout << stoichiometry_[i] << " " << species_names_[i];
//       if (i < species_names_.size() - 1) {
//         std::cout << " + ";
//       }
//     }
//   }

//   if (exchange_site_names_.size() > 0) {
//     std::cout << " + ";
//     for (unsigned int i = 0; i < exchange_site_names_.size(); i++) {
//       std::cout << exchange_site_stoichiometries_[i] << " "
//                 << exchange_site_names_[i];
//       if (i < exchange_site_names_.size() - 1) {
//         std::cout << " + ";
//       }
//     }
//   }
//   std::cout << std::endl;
// }  // end DisplayReaction_arrays()
