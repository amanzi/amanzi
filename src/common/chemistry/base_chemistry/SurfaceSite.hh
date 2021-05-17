/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  Class for storing surface site data for surface complexation
*/

#ifndef AMANZI_CHEMISTRY_SURFACE_SITE_HH_
#define AMANZI_CHEMISTRY_SURFACE_SITE_HH_

#include <vector>

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations
class Mineral;

class SurfaceSite {
 public:
  SurfaceSite();
  SurfaceSite(const std::string& name, int id,
              const Teuchos::ParameterList& plist);
  ~SurfaceSite() {};

  // adds a pointer to mineral list
  void AddMineral(Mineral* mineral);
  double SiteDensity() const;
  void UpdateSiteDensity(const double site_density);

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  // These are changing two variables.
  void set_free_site_concentration(double d) {
    free_site_concentration_ = d;
    ln_free_site_concentration_ = std::log(d);
  }
  void set_ln_free_site_concentration(double d) {
    free_site_concentration_ = std::exp(d);
    ln_free_site_concentration_ = d;
  }
  void set_molar_density(double d) { molar_density_ = d; }
  void set_molar_surface_density(double d) { molar_surface_density_ = d; }

  std::string name() const { return name_; }
  int identifier() const { return identifier_; }
  double charge() const { return charge_; }

  double free_site_concentration() const { return free_site_concentration_; }
  double ln_free_site_concentration() const { return ln_free_site_concentration_; }

  double molar_density() const { return molar_density_; }
  double molar_surface_density() const { return molar_surface_density_; }

 private:
  std::string name_;
  int identifier_;
  double charge_;

  double molar_density_;  // [moles sites / m^3 mineral]
  double molar_surface_density_;  // [moles sites / m^2 mineral]
  double free_site_concentration_;  // [moles sites / m^3 bulk]
  double ln_free_site_concentration_;  // [-]

  // for future use
  // a list of pointers to minerals assocaited with the site
  // std::vector<Mineral*> minerals_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
