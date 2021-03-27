/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for aqueous equilibrium complexation reaction
*/

#ifndef AMANZI_CHEMISTRY_SURFACECOMPLEX_HH_
#define AMANZI_CHEMISTRY_SURFACECOMPLEX_HH_

#include <vector>

#include "species.hh"
#include "surface_site.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class SurfaceComplex {
 public:
  SurfaceComplex();
  SurfaceComplex(const SpeciesName name,
                 const SpeciesId id,
                 const std::vector<SpeciesName>& species,
                 const std::vector<double>& stoichiometries,
                 const std::vector<int>& species_ids,
                 const double h2o_stoich,
                 const double free_site_stoichiometry,
                 const double charge,
                 const double logK);
  SurfaceComplex(const SpeciesName name,
                 const SpeciesId id,
                 const std::vector<SpeciesName>& species,
                 const std::vector<double>& stoichiometries,
                 const std::vector<int>& species_ids,
                 const double h2o_stoich,
                 SpeciesName free_site_name,
                 const double free_site_stoichiometry,
                 SpeciesId free_site_id,
                 const double charge,
                 const double logK);
  ~SurfaceComplex() {};

  // update molalities
  void Update(const std::vector<Species>& primarySpecies,
              const SurfaceSite& surface_site);
  // add stoichiometric contribution of complex to total
  void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to dtotal
  void AddContributionToDTotal(const std::vector<Species>& primarySpecies,
                               MatrixBlock* dtotal);

  void display(const Teuchos::Ptr<VerboseObject> vo) const;
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  void set_name(SpeciesName name) { this->name_ = name; }
  void set_identifier(SpeciesId i) { this->identifier_ = i; }
  void set_charge(double d) { this->charge_ = d; }

  void set_free_site_stoichiometry(const double d) { this->free_site_stoichiometry_ = d; }
  void set_lnQK(const double d) { this->lnQK_ = d; };
  void set_logK(const double d) { this->logK_ = d; };
  void set_ncomp(const int i) { this->ncomp_ = i; };
  void set_surface_concentration(const double d) { this->surface_concentration_ = d; };

  SpeciesName name(void) const { return this->name_; }
  SpeciesId identifier(void) const { return this->identifier_; }
  double charge(void) const { return this->charge_; }

  double free_site_stoichiometry(void) const { return this->free_site_stoichiometry_; }
  double stoichiometry(const int i) const { return this->stoichiometry_[i]; }
  double lnQK(void) const { return this->lnQK_; };
  double logK(void) const { return this->logK_; };
  int ncomp(void) const { return this->ncomp_; };
  int species_id(const int i) const { return this->species_ids_[i]; };
  double surface_concentration(void) const { return this->surface_concentration_; };

 private:
  double log_to_ln(double d) {
    return d * 2.302585092994046;
  }
  double ln_to_log(double d) {
    return d * 0.434294481903252;
  }

  SpeciesName name_;
  SpeciesId identifier_;
  double charge_;

  double surface_concentration_;  // units? ?[mol/m^3 bulk]?

  int ncomp_;  // # components in reaction
  std::vector<SpeciesName> species_names_;
  std::vector<SpeciesId> species_ids_;       // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn
  SpeciesName free_site_name_;
  double free_site_stoichiometry_;     // stoichiometry of free site in rxn
  SpeciesId free_site_id_;
  std::vector<double> logK_array_;     // for temperature dep. logK
  double h2o_stoichiometry_;           // stoichiometry of water in equation
  double lnK_;                         // log value of equlibrium constant
  double lnQK_;                        // store lnQK for derivatives later
  double logK_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
