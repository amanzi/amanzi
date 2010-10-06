/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*******************************************************************************
**
**  Description: implementation of the TST rate law for mineral kinetics
**
**  R = k * A * Prod (a_i^m_i) * (1 - Q/Keq)
**
**  where:
**    R : reaction rate, [moles/sec]
**    Keq : equilibrium constant, [-]
**    Q : ion activity product, [-]
**    k : reaction rate constant, [moles m^2 s^-1]
**    A : reactive surface area, [m^2]
**    a_i : activity of modifying species
**    m_i : exponent of modifying species
**
**
**  Notes:
**
**    - assume that the "reactants" list as defined in KineticRate
**    consists solely of the mineral and the coefficient is always
**    one, so they are ignored in the ion activity product. The
**    mineral is only accounted for in the accumulation / consumption
**    of solid mass.
**
**    - assume that the mineral reactions are always written so the
**    products consist of only primary species.
**
**    - assume that the modifying species are only primary species
**
**    - ignoring the (1-(Q/Keq)^n) and (1-(Q/Keq)^n)^p forms for now.
**
*******************************************************************************/

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
      sat_state_exponent_(0.0),
      mineral_()
{
  this->modifying_species_names.clear();
  this->modifying_exponents.clear();
  this->modifying_species_ids.clear();
}  // end KineticRateTST constructor

KineticRateTST::~KineticRateTST(void)
{
}  // end KineticRateTST destructor


void KineticRateTST::Setup(const Verbosity s_verbosity, std::string reaction, 
                           StringTokenizer reaction_data)
{
  verbosity(s_verbosity);
  ParseReaction(reaction);
  // first element in KineticRate.reactants is the mineral
  SpeciesName mineral_name(reactants_names.at(0));
  mineral_.name(mineral_name);
  ParseParameters(reaction_data);
  
}  // end Setup()


void KineticRateTST::Update(const std::vector<Species> primarySpecies)
{

}  // end Update()

void KineticRateTST::AddContributionToResidual(const double por_den_sat_vol, 
                                               std::vector<double> *residual)
{

  // Calculate K/Q 

  // Calculate modifying term based on primary species

  // Calculate overall rate 
  
}  // end AddContributionToResidual()

void KineticRateTST::AddContributionToJacobian(const std::vector<Species> primarySpecies,
                                               const double por_den_sat_vol,
                                               Block *J)
{

}  // end AddContributionToJacobian()

void KineticRateTST::ParseParameters(StringTokenizer reaction_data)
{
  std::string space(" ");
  StringTokenizer st;

  std::vector<std::string>::iterator field = reaction_data.begin();

  for (; field != reaction_data.end(); field++) {
    st.tokenize(*field, space);

    if (!(st.at(0).compare("gmw"))) {
      mineral_.gram_molecular_weight(std::atof(st.at(1).c_str()));
    }
    else if (!(st.at(0).compare("molar_density"))) {
      mineral_.molar_density(std::atof(st.at(1).c_str()));
    }
    else if (!(st.at(0).compare("rate_constant"))) {
      rate_constant(std::atof(st.at(1).c_str()));
    }
    else if (!(st.at(0).compare("area"))) {
      area(std::atof(st.at(1).c_str()));
    } 
    else if (!(st.at(0).compare("pK"))) {
      pK(std::atof(st.at(1).c_str()));
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
  // TODO: now need to identify the indicies in the primary species
  // list....
}  // end ParseParameters()

void KineticRateTST::Display(void) const
{
  std::cout << "  Rate law: TST" << std::endl;
  this->DisplayReaction();
  std::cout << "  Parameters:" << std::endl;
  std::cout << "    mineral = " << mineral_.name() << std::endl;
  std::cout << "    molar density = " << mineral_.molar_density() 
            << " [cm^3/mole]" << std::endl;
  std::cout << "    gram molecular weight = " 
            << mineral_.gram_molecular_weight() << " [grams/mole]" 
            << std::endl;
  std::cout << "    rate constant = " << rate_constant() 
            << " [moles/m^2/sec]" << std::endl;
  std::cout << "    area = " << area() << " [m^2]" << std::endl;
  std::cout << "    pK = " << pK() << " [-]" << std::endl;
  std::cout << "    rate modifiers: " << std::endl;
  for (unsigned int mod = 0; mod < this->modifying_species_names.size(); mod++) {
    std::cout << "      { " << this->modifying_species_names.at(mod) << " }";
    std::cout << "^" << this->modifying_exponents.at(mod) << std::endl;
  }
  std::cout << std::endl;
}  // end Display()
