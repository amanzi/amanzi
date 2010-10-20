/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <iomanip>

#include "IonExchangeComplex.hpp"
#include "SecondarySpecies.hpp"

IonExchangeComplex::IonExchangeComplex() 
    : SecondarySpecies()
{
} // end IonExchangeComplex() constructor

IonExchangeComplex::IonExchangeComplex(
    const SpeciesName complex_name, 
    SpeciesId complex_id,
    std::vector<SpeciesName> complex_primary_species,
    std::vector<double> complex_primary_stoichiometries,
    std::vector<SpeciesId> complex_primary_ids,
    std::vector<SpeciesName> complex_exchange_site_names,
    std::vector<double> complex_exchange_site_stoichiometries,
    std::vector<SpeciesId> complex_exchange_site_ids,
    const double complex_h2o_stoich,
    const double complex_log_Keq)
    : SecondarySpecies(complex_name, complex_id, 
                       complex_primary_species, complex_primary_stoichiometries, 
                       complex_primary_ids,
                       complex_h2o_stoich, 0.0, 0.0, 0.0, 
                       complex_log_Keq)
{
  exchange_site_names_ = complex_exchange_site_names;
  exchange_site_stoichiometries_ = complex_exchange_site_stoichiometries;
  exchange_site_ids_ = complex_exchange_site_ids;
}  // end IonExchangeComplex costructor


IonExchangeComplex::~IonExchangeComplex() 
{
} // end IonExchangeComplex() destructor


void IonExchangeComplex::Update(const std::vector<Species> primary_species,
                                const std::vector<IonExchangeSite> exchange_sites) 
{
  static_cast<void>(primary_species);
  static_cast<void>(exchange_sites);
} // end update()

void IonExchangeComplex::AddContributionToTotal(std::vector<double> &total) 
{
  static_cast<void>(total);
} // end addContributionToTotal()

void IonExchangeComplex::AddContributionToDTotal(const std::vector<Species> primary_species,
                                                 Block *dtotal) 
{
  static_cast<void>(primary_species);
  static_cast<void>(dtotal);
} // end addContributionToDTotal()


/*
**
**  Display functions
**
*/
void IonExchangeComplex::display(void) const
{
  DisplayReaction();
  std::cout << "      log_K: " << logK()
            << std::endl;
} // end Display()

void IonExchangeComplex::Display(void) const
{
  DisplayReaction();
  std::cout << std::setw(40) << " " 
            << std::setw(10) << logK()
            << std::endl;
} // end Display()

void IonExchangeComplex::DisplayReaction(void) const
{
  std::cout << "    " << name() << " = ";
  if (species_names_.size() > 0) {
    for (unsigned int i = 0; i < species_names_.size(); i++) {
      std::cout << stoichiometry_[i] << " " << species_names_[i];
      if (i < species_names_.size() - 1) {
        std::cout << " + ";
      }
    }
  }

  if (exchange_site_names_.size() > 0) {
    std::cout << " + ";
    for (unsigned int i = 0; i < exchange_site_names_.size(); i++) {
      std::cout << exchange_site_stoichiometries_[i] << " " 
                << exchange_site_names_[i];
      if (i < exchange_site_names_.size() - 1) {
        std::cout << " + ";
      }
    }
  }
  std::cout << std::endl;
} // end DisplayReaction()
