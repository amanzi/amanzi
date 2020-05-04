/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Description: Implementation of ion exchange reactions as
  equilibrium half reactions. Using the activity of exchange complex
  is the equivalent fraction of exchangable cations (exchange sites)

  Assumptions:
    - Assume that the reactions are always going to be written as
      one complex equals one primary plus one exchange site:

        1.0 complex = double primary + double exchange_site

    - can we assume that the exchanger activity is always one...?

  Examples:
    Ca0.5X = 0.5Ca++ + X- ; Keq = double

    CaX_2 = Ca++ + 2X- ; Keq = double

    NaX = Na+ + X- ; Keq = double

  Mass Action Law:
    K_eq = a_p^nu_p * a_x^nu_x / a_c^nu_c

  where:
    subscript p = primary species
    subscript x = exchange site
    subscript c = ion exchange complex

    a_ = activity of species
    nu_ = stoichiometric coefficient of species
    Keq = half reaction equilibrium constant

    a_c = eq_c / CEC
        = m_c * nu_x / CEC

  where:
    a_c = activity of the exchanger-cation = equivalent fraction of c

    eq_c = equivalents of complex c, equivalents per MASS? VOLUME?
           equivalents are equivalent number of exchange sites occupied
           NaX, Ca_{0.5}X : eq_c = conc of complex, complex occupies one exchange site
           CaX_2 : eq_c = 2*conc of complex, complex occupies two exchange sites.
         = nu_x * m_c

    CEC = cation exchange capacity, equivalents per MASS? VOLUME?

  Residual:
    R(a_p, a_x, a_c) = -ln_Keq +
                       nu_p * ln(a_p) +
                       nu_x * ln(a_x) -
                       nu_c * (ln(m_c) + ln(nu_x) - ln(CEC))
*/

#include <iostream>
#include <iomanip>
#include <sstream>

#include "VerboseObject.hh"

#include "ion_exchange_complex.hh"

namespace Amanzi {
namespace AmanziChemistry {

IonExchangeComplex::IonExchangeComplex(
    const IonxComplexName in_name,
    const IonxComplexId in_id,
    const SpeciesName in_primary_name,
    const SpeciesId in_primary_id,
    const double in_K)
    : name_(in_name),
      id_(in_id),
      primary_name_(in_primary_name),
      primary_id_(in_primary_id),
      K_(in_K),
      concentration_(0.),
      X_(0.) {
}


/*
**
**  Display functions
**
*/
void IonExchangeComplex::display(const Teuchos::Ptr<VerboseObject> vo) const {
  DisplayReaction(vo);
  std::cout << "      K: " << K() << std::endl;
}


void IonExchangeComplex::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  DisplayReaction(vo);
  std::stringstream message;
  message << std::setw(40) << " "
          << std::setw(10) << K()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}


void IonExchangeComplex::DisplayReaction(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << "    " << name() << " = "
          << primary_name()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}


void IonExchangeComplex::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << "Name"
          << std::setw(15) << "X"
          << std::setw(15) << "Concentration"
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}


void IonExchangeComplex::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::scientific << std::setprecision(5)
          << std::setw(15) << X()
          << std::setw(15) << concentration()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
