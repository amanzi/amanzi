/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  Class for storing surface site data for surface complexation
*/

#include <sstream>
#include <iostream>
#include <iomanip>

#include "VerboseObject.hh"

#include "surface_site.hh"

namespace Amanzi {
namespace AmanziChemistry {

SurfaceSite::SurfaceSite()
    : name_(""),
      identifier_(0),
      charge_(0.),
      molar_density_(0.),
      molar_surface_density_(0.),
      free_site_concentration_(0.),
      ln_free_site_concentration_(0.) {
  // minerals_.clear();
}

SurfaceSite::SurfaceSite(const SpeciesName name,
                         const SpeciesId id,
                         const double molar_density)
    : name_(name),
      identifier_(id),
      charge_(0.),
      molar_density_(molar_density),
      molar_surface_density_(0.),
      // initialize to 10% of molar_density
      free_site_concentration_(0.1 * molar_density),
      ln_free_site_concentration_(std::log(0.1 * molar_density)) {
  // minerals_.clear();
}


// Add a pointer to mineral list
void SurfaceSite::AddMineral(Mineral* mineral) {
  // minerals_.push_back(mineral);
}


void SurfaceSite::UpdateSiteDensity(const double site_density) {
  // needs to change for minerals....
  set_molar_density(site_density);
} 


// Sum the total site concentration based on minerals
double SurfaceSite::SiteDensity(void) const {
  /* For now, we are skipping the use of minerals - geh
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
  */
  return molar_density();
}


void SurfaceSite::display(void) const {
  std::cout << "    " << name() << std::endl;
  std::cout << "        site density = " << molar_density() << std::endl;
}


void SurfaceSite::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::setw(15) << std::scientific << molar_density()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}


void SurfaceSite::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << "Site Name"
          << std::setw(15) << "Free Conc."
          << std::endl;
  message << std::setw(15) << " "
          << std::setw(15) << "[mol/m^3 bulk]"
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}


void SurfaceSite::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::scientific << std::setprecision(5)
          << std::setw(15) << free_site_concentration()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message);
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
