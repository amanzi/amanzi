/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>

#include <iostream>

#include "KineticRateTST.hpp"
#include "StringTokenizer.hpp"
#include "Verbosity.hpp"

KineticRateTST::KineticRateTST(void)
    : KineticRate()
{
}  // end KineticRateTST constructor

KineticRateTST::~KineticRateTST(void)
{
}  // end KineticRateTST destructor

void KineticRateTST::Update(const std::vector<Species> primarySpecies)
{

}  // end Update()

void KineticRateTST::AddContributionToResidual(const double por_den_sat_vol, 
                                               std::vector<double> *residual)
{

}  // end AddContributionToResidual()
void KineticRateTST::AddContributionToJacobian(const std::vector<Species> primarySpecies,
                                               const double por_den_sat_vol,
                                               Block *J)
{

}  // end AddContributionToJacobian()

void KineticRateTST::ParseParameters(StringTokenizer rate)
{
  double area = 0.0;
  double pK = 0.0;
  double k = 0.0;
  double sat_state_exponent = 0.0;
  std::vector<SpeciesName> modifying_species_names;
  std::vector<double> modifying_exponents;
  std::vector<int> modifying_species_ids;

  std::string space(" ");
  StringTokenizer st;

  std::vector<std::string>::iterator field = rate.begin();
  field++; // field[0] is the reaction string, handle it elsewhere
  field++; // field[1] is the rate name, handled elsewhere

  for (; field != rate.end(); field++) {
    st.tokenize(*field, space);
    if (!(st.at(0).compare("area"))) {
      area = std::atof(st.at(1).c_str());
    } else if (!(st.at(0).compare("pK"))) {
      pK = std::atof(st.at(1).c_str());
    } else if (!(st.at(0).compare("k"))) {
      k = std::atof(st.at(1).c_str());
    } else if (!(st.at(0).compare("n"))) {
      sat_state_exponent = std::atof(st.at(1).c_str());
    } else {
      // assume we are dealing with the list of rate modifying speces
      for (unsigned int modifier = 0; modifier < st.size(); modifier++) {
        modifying_species_names.push_back(st.at(modifier));
        modifier++; // increment to get the exponent of this modifier
        modifying_exponents.push_back(std::atof(st.at(modifier).c_str()));
      }      
    }
    if (verbosity() == kDebugMineralKinetics) {
      std::cout << "  Field: " << *field << std::endl;
    }
  }
  if (verbosity() == kDebugMineralKinetics) {
    std::cout << "area = " << area << std::endl;
    std::cout << "pK = " << pK << std::endl;
    std::cout << "k = " << k << std::endl;
    std::cout << "n = " << sat_state_exponent << std::endl;
    std::cout << "rate modifiers: " << std::endl;
    for (unsigned int mod = 0; mod < modifying_species_names.size(); mod++) {
      std::cout << "  { " << modifying_species_names.at(mod) << " }";
      std::cout << "^" << modifying_exponents.at(mod) << std::endl;
    }
  }
}  // end ParseParameters()

void KineticRateTST::Display(void) const
{
  std::cout << "TST rate: " << std::endl;
}  // end Display()
