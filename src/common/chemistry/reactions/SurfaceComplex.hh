/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for surface equilibrium complexation reaction
*/

#ifndef AMANZI_CHEMISTRY_SURFACECOMPLEX_HH_
#define AMANZI_CHEMISTRY_SURFACECOMPLEX_HH_

#include <string>
#include <vector>

#include "Species.hh"
#include "SurfaceSite.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class SurfaceComplex {
 public:
  SurfaceComplex();
  SurfaceComplex(const std::string& name,
                 const int id,
                 const std::vector<std::string>& species,
                 const std::vector<double>& stoichiometries,
                 const std::vector<int>& species_ids,
                 const double h2o_stoich,
                 const double free_site_stoichiometry,
                 const double charge,
                 const double logK);

  SurfaceComplex(const std::string& name,
                 const int id,
                 const std::vector<std::string>& species,
                 const std::vector<double>& stoichiometries,
                 const std::vector<int>& species_ids,
                 const double h2o_stoich,
                 const std::string& free_site_name,
                 const double free_site_stoichiometry,
                 int free_site_id,
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

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  void set_identifier(int id) { identifier_ = id; }
  void set_charge(double d) { charge_ = d; }

  void set_free_site_stoichiometry(double d) { free_site_stoichiometry_ = d; }
  void set_ncomp(const int i) { ncomp_ = i; };
  void set_surface_concentration(double d) { surface_concentration_ = d; };

  std::string name() const { return name_; }
  int identifier() const { return identifier_; }
  double charge() const { return charge_; }

  double free_site_stoichiometry() const { return free_site_stoichiometry_; }
  double stoichiometry(int i) const { return stoichiometry_[i]; }
  double lnQK() const { return lnQK_; };
  double logK() const { return logK_; };
  int ncomp() const { return ncomp_; };
  int species_id(int i) const { return species_ids_[i]; };
  double surface_concentration() const { return surface_concentration_; };

 private:
  std::string name_;
  int identifier_;
  double charge_;

  double surface_concentration_;  // units? ?[mol/m^3 bulk]?

  int ncomp_;  // # components in reaction
  std::vector<std::string> species_names_;
  std::vector<int> species_ids_;  // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn
  std::string free_site_name_;
  double free_site_stoichiometry_;  // stoichiometry of free site in rxn
  int free_site_id_;
  std::vector<double> logK_array_;  // for temperature dep. logK
  double h2o_stoichiometry_;  // stoichiometry of water in equation
  double lnK_;  // log value of equlibrium constant
  double lnQK_;  // store lnQK for derivatives later
  double logK_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
