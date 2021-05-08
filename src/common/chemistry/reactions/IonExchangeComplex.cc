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

#include "IonExchangeComplex.hh"

namespace Amanzi {
namespace AmanziChemistry {

IonExchangeComplex::IonExchangeComplex(
    const std::string& name, int id,
    const Teuchos::ParameterList& plist,
    const std::vector<Species>& primary_species)
  : name_(name),
    id_(id),
    concentration_(0.0),
    X_(0.0)
{
  double coeff;

  K_ = plist.get<double>("equilibrium constant");
  std::istringstream iss(plist.get<std::string>("reaction"));
  iss >> coeff;
  iss >> primary_name_;
  iss >> coeff;
  iss >> site_name_;

  primary_id_ = -999;
  for (int i = 0; i < primary_species.size(); ++i) {
    if (primary_name_ == primary_species[i].name()) {
      primary_id_ = i;
      break;
    }
  }
}


void IonExchangeComplex::Display(const Teuchos::Ptr<VerboseObject> vo) const
{
  DisplayReaction(vo);
  std::stringstream message;
  message << std::setw(40) << " "
          << std::setw(10) << K()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void IonExchangeComplex::DisplayReaction(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << "    " << name() << " = "
          << primary_name()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void IonExchangeComplex::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << "Name"
          << std::setw(15) << "X"
          << std::setw(15) << "Concentration"
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void IonExchangeComplex::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const
{
  std::stringstream message;
  message << std::setw(15) << name()
          << std::scientific << std::setprecision(5)
          << std::setw(15) << X()
          << std::setw(15) << concentration()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
