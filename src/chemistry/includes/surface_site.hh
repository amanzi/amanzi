#ifndef AMANZI_CHEMISTRY_SURFACESITE_HH_
#define AMANZI_CHEMISTRY_SURFACESITE_HH_

// Class for storing surface site data for surface complexation

#include <vector>

#include "mineral.hh"

// forward declarations
class Block;

class SurfaceSite {

 public:
  SurfaceSite();
  SurfaceSite(const SpeciesName name, 
              const SpeciesId id,
              const double molar_density);
  ~SurfaceSite();

  // adds a pointer to mineral list
  void AddMineral(Mineral *mineral);
  double SiteDensity(void) const;

  void display(void) const;
  void Display(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

  // TODO: should all these really be public? should not be changing the
  // name, id, etc after creation!
  void set_name(SpeciesName name) { this->name_ = name; }
  void set_identifier(SpeciesId i) { this->identifier_ = i; }
  void set_charge(double d) { this->charge_ = d; }

  // TODO: I'd like to keep set_xyz just as a mutator for a single
  // variable. Since these are changing two variables, can we call
  // them update_free_site_concentration, update_ln_free_site....?
  void set_free_site_concentration(const double d) 
    { this->free_site_concentration_ = d; 
      this->ln_free_site_concentration_ = std::log(d); }
  void set_ln_free_site_concentration(const double d) 
    { this->free_site_concentration_ = std::exp(d); 
      this->ln_free_site_concentration_ = d; }
  void set_molar_density(const double d) { this->molar_density_ = d; }
  void set_molar_surface_density(const double d) 
    { this->molar_surface_density_ = d; }

  SpeciesName name(void) const { return this->name_; }
  SpeciesId identifier(void) const { return this->identifier_; }
  double charge(void) const { return this->charge_; }

  double free_site_concentration(void) const 
    { return this->free_site_concentration_; }
  double ln_free_site_concentration(void) const 
    { return this->ln_free_site_concentration_; }
  double molar_density(void) const { return this->molar_density_; }
  double molar_surface_density(void) const 
    { return this->molar_surface_density_; }

 protected:

 private:

  SpeciesName name_;
  SpeciesId identifier_;
  double charge_;

  double molar_density_; // [moles sites / m^3 mineral]
  double molar_surface_density_; // [moles sites / m^2 mineral]
  double free_site_concentration_; // [moles sites / m^3 bulk]
  double ln_free_site_concentration_; // [-]

  // for future use
  // a list of pointers to minerals assocaited with the site
 // std::vector<Mineral*> minerals_;

};

#endif // AMANZI_CHEMISTRY_SURFACESITE_HH_
