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

#include "Teuchos_ParameterList.hpp"
#include "VerboseObject.hh"

#include "SurfaceSite.hh"

namespace Amanzi {
namespace AmanziChemistry {

SurfaceSite::SurfaceSite()
  : name_(""),
    identifier_(0),
    charge_(0.0),
    molar_density_(0.0),
    molar_surface_density_(0.0),
    free_site_concentration_(0.0),
    ln_free_site_concentration_(0.0) {
  // minerals_.clear();
}

SurfaceSite::SurfaceSite(const std::string& name, int id,
                         const Teuchos::ParameterList& plist)
  : name_(name),
    identifier_(id),
    charge_(0.0),
    molar_surface_density_(0.0)
{
  molar_density_ = plist.get<double>("density");

  // initialize to 10% of molar_density
  free_site_concentration_ = 0.1 * molar_density_;
  ln_free_site_concentration_ = std::log(0.1 * molar_density_);
}


// Add a pointer to mineral list
void SurfaceSite::AddMineral(Mineral* mineral) {
  // minerals_.push_back(mineral);
}


void SurfaceSite::UpdateSiteDensity(double site_density) {
  // needs to change for minerals....
  molar_density_  = site_density;
} 


// Sum the total site concentration based on minerals
double SurfaceSite::SiteDensity() const {
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


void SurfaceSite::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::setw(15) << std::scientific << molar_density()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void SurfaceSite::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << "Site Name"
          << std::setw(15) << "Free Conc."
          << std::endl;
  message << std::setw(15) << " "
          << std::setw(15) << "[mol/m^3 bulk]"
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void SurfaceSite::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(15) << name()
          << std::scientific << std::setprecision(5)
          << std::setw(15) << free_site_concentration()
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
