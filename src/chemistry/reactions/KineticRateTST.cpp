/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>

#include <iostream>

#include "KineticRateTST.hpp"
#include "StringTokenizer.hpp"
#include "Verbosity.hpp"

KineticRateTST::KineticRateTST(void)
    : KineticRate(),
      area_(0.0),
      pK_(0.0),
      rate_constant_(0.0),
      sat_state_exponent_(0.0)
{
  this->modifying_species_names.clear();
  this->modifying_exponents.clear();
  this->modifying_species_ids.clear();
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
  // assume that the "reactants" as defined in KineticRate consists
  // solely of the mineral and the coefficient is always one, so they
  // are ignored in the ion activity product. The mineral is only
  // accounted for in the accumulation / consumption of solid mass.

  // assume that the mineral reactions are always written so the
  // products consist of only primary species.

  // Calculate K/Q 

  // Calculate modifying term based on primary species

  // Calculate modifying term based on secondary species

  // Calculate overall rate 
  
}  // end AddContributionToResidual()

void KineticRateTST::AddContributionToJacobian(const std::vector<Species> primarySpecies,
                                               const double por_den_sat_vol,
                                               Block *J)
{

}  // end AddContributionToJacobian()

void KineticRateTST::ParseParameters(StringTokenizer rate)
{
  std::string space(" ");
  StringTokenizer st;

  std::vector<std::string>::iterator field = rate.begin();
  field++; // field[0] is the reaction string, handle it elsewhere
  field++; // field[1] is the rate name, handled elsewhere

  for (; field != rate.end(); field++) {
    st.tokenize(*field, space);
    if (!(st.at(0).compare("area"))) {
      area(std::atof(st.at(1).c_str()));
    } else if (!(st.at(0).compare("pK"))) {
      pK(std::atof(st.at(1).c_str()));
    } else if (!(st.at(0).compare("k"))) {
      rate_constant(std::atof(st.at(1).c_str()));
    } else if (!(st.at(0).compare("n"))) {
      sat_state_exponent(std::atof(st.at(1).c_str()));
    } else {
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
  // TODO: now need to determine if the modifying species are primary
  // or secondary and identify the indicies in the proper species
  // list....
}  // end ParseParameters()

void KineticRateTST::Display(void) const
{
  std::cout << "  Rate law: TST" << std::endl;
  this->DisplayReaction();
  std::cout << "  Parameters:" << std::endl;
  std::cout << "    area = " << area() << " [m^2]" << std::endl;
  std::cout << "    pK = " << pK() << " [-]" << std::endl;
  std::cout << "    rate constant = " << rate_constant() << " [moles/m^2/sec]" << std::endl;
  std::cout << "    n = " << sat_state_exponent() << " [-]" << std::endl;
  std::cout << "    rate modifiers: " << std::endl;
  for (unsigned int mod = 0; mod < this->modifying_species_names.size(); mod++) {
    std::cout << "      { " << this->modifying_species_names.at(mod) << " }";
    std::cout << "^" << this->modifying_exponents.at(mod) << std::endl;
  }
  std::cout << std::endl;
}  // end Display()
