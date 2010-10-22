#include <iomanip>

#include "SurfaceSite.hpp"

SurfaceSite::SurfaceSite() 
{
  minerals_.clear();
} 

SurfaceSite::~SurfaceSite() 
{
}

// Add a pointer to mineral list
void SurfaceSite::AddMineral(Mineral *mineral) {
  minerals_.push_back(mineral);
}

// Sum the total site concentration based on minerals 
double SurfaceSite::SiteDensity(void) const
{
  double sum = 0.;
  for (std::vector<Mineral*>::const_iterator mineral_ptr=minerals_.begin();
       mineral_ptr != minerals_.end(); mineral_ptr++)
    // surface_area [m^2 mineral / m^3 mineral]
    // volume_fraction [m^3 mineral / m^3 bulk]
    sum += (*mineral_ptr)->surface_area()*(*mineral_ptr)->volume_fraction();
  // sum [m^2 mineral / m^3 bulk]
  // molar_surface_density [moles sites / m^2 mineral]
  // return [moles sites / m^3 bulk]
  return sum*molar_surface_density();
} // end SiteDensity()

void SurfaceSite::display(void) const
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
  std::cout << "        logK = " << logK_ << std::endl;
  std::cout << "        charge = " << charge() << std::endl;
  std::cout << "        mol wt = " << gram_molecular_weight() << std::endl;
*/
} // end display()

void SurfaceSite::Display(void) const
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
