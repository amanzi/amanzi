/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <iomanip>

#include "Mineral.hpp"
#include "SecondarySpecies.hpp"

Mineral::Mineral() 
    : SecondarySpecies(),
      molar_density_(0.0)
{
} // end Mineral() constructor

Mineral::Mineral(const SpeciesName in_name, 
                 const SpeciesId in_id,
                 std::vector<SpeciesName> in_species,
                 std::vector<double> in_stoichiometries,
                 std::vector<int> in_species_ids,
                 const double in_h2o_stoich, 
                 const double in_mol_wt, 
                 const double in_logK, 
                 const double in_molar_density)
    : SecondarySpecies(in_name, in_id, 
                       in_species, in_stoichiometries, in_species_ids,
                       in_h2o_stoich, 0., in_mol_wt, 0., in_logK),
      molar_density_(in_molar_density)
{
}  // end Mineral costructor


Mineral::~Mineral() 
{
} // end Mineral() destructor

/*
**
**  these functions are only needed if mineral equilibrium is added.
**
*/
void Mineral::Update(const std::vector<Species> primary_species) 
{
  double lnQK = -lnK_;
  for (int i = 0; i < ncomp_; i++) {
    lnQK += stoichiometry_[i] * primary_species[species_ids_[i]].ln_activity();
  }
  lnQK_ = lnQK;
  static_cast<void>(primary_species);
} // end update()

void Mineral::AddContributionToTotal(std::vector<double> &total) 
{
  static_cast<void>(total);
} // end addContributionToTotal()

void Mineral::AddContributionToDTotal(const std::vector<Species> primary_species,
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
void Mineral::Display(void) const
{
  std::cout << "    " << name() << " = ";
  for (int i = 0; i < (int)species_names_.size(); i++) {
    std::cout << stoichiometry_[i] << " " << species_names_[i];
    if (i < (int)species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
  std::cout << std::setw(40) << " " 
            << std::setw(10) << logK_
            << std::setw(15) << molar_density()
            << std::setw(10) << gram_molecular_weight()
            << std::endl;
} // end Display()

void Mineral::DisplayResultsHeader(void) const
{
  std::cout << std::setw(15) << "Name" 
            << std::setw(15) << "Q/K"
            << std::setw(15) << "SI"
            << std::endl;
} // end DisplayResultsHeader()

void Mineral::DisplayResults(void) const
{
  std::cout << std::setw(15) << name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << Q_over_K()
            << std::fixed << std::setprecision(3)
            << std::setw(15) << saturation_index()
            << std::endl;
} // end DisplayResults()
