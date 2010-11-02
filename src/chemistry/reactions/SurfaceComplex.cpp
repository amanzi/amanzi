#include <iomanip>

#include "SurfaceComplex.hpp"

SurfaceComplex::SurfaceComplex()
{
  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();
  logK_array_.clear();

} // end SurfaceComplex() constructor

SurfaceComplex::SurfaceComplex(const SpeciesName name, 
                            const SpeciesId id,
                            std::vector<SpeciesName>species,
                            std::vector<double>stoichiometries,
                            std::vector<int>species_ids,
                            const double h2o_stoich, 
                            const double free_site_stoich,
                            const double charge, 
                            const double logK) 
                            : name_(name),
                              identifier_(id),
                              charge_(charge),
                              surface_concentration_(0.),
                              free_site_stoichiometry_(free_site_stoich),
                              h2o_stoichiometry_(h2o_stoich),
                              lnK_(log_to_ln(logK)),
                              lnQK_(0.),
                              logK_(logK)
{

  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();
  logK_array_.clear();

  set_ncomp(static_cast<int>(stoichiometries.size()));

  // species names
  for (std::vector<SpeciesName>::const_iterator i = species.begin(); 
       i != species.end(); i++) {
    species_names_.push_back(*i);
  } 
  // species stoichiometries
  for (std::vector<double>::const_iterator i = stoichiometries.begin();
       i != stoichiometries.end(); i++) {
    stoichiometry_.push_back(*i);
  }
  // species ids
  for (std::vector<int>::const_iterator i = species_ids.begin();
       i != species_ids.end(); i++) {
    species_ids_.push_back(*i);
  }
} // end SurfaceComplex() constructor

SurfaceComplex::~SurfaceComplex() 
{
} // end SurfaceComplex() destructor

void SurfaceComplex::Update(const std::vector<Species> primarySpecies,
                            const SurfaceSite& surface_site) 
{

  double lnQK_temp = -lnK_;

  // Need to consider activity of water

  for (int i = 0; i < ncomp_; i++) {
    lnQK_temp += stoichiometry_[i] * 
                 primarySpecies[species_ids_[i]].ln_activity();
  }
  lnQK_temp += free_site_stoichiometry() *
               surface_site.ln_free_site_concentration();
  set_lnQK(lnQK_temp);
  set_surface_concentration(std::exp(lnQK()));

} // end Update()

void SurfaceComplex::AddContributionToTotal(std::vector<double> *total) 
{
  for (int i = 0; i < ncomp_; i++) {
    (*total)[species_ids_[i]] += stoichiometry_[i] * surface_concentration(); 
  }
} // end AddContributionToTotal()

void SurfaceComplex::AddContributionToDTotal(
                                   const std::vector<Species> primarySpecies,
                                   Block *dtotal) 
{
} // end AddContributionToDTotal()

void SurfaceComplex::Display(void) const
{
/*
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
            << std::setw(10) << charge()
            << std::setw(10) << gram_molecular_weight()
            << std::endl;
*/
} // end Display()
