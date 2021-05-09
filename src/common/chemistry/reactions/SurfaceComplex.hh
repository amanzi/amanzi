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

#include "Teuchos_ParameterList.hpp"

#include "FunctionTabular.hh"

#include "Species.hh"
#include "SurfaceSite.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class SurfaceComplex {
 public:
  SurfaceComplex() {};
  SurfaceComplex(const std::string& name, int id,
                 const std::vector<Species>& primary_species,
                 const std::vector<SurfaceSite>& surface_sites,
                 const Teuchos::ParameterList& plist);
  ~SurfaceComplex() {};

  // update molalities
  void Update(const std::vector<Species>& primary_species,
              const SurfaceSite& surface_site);

  // update temperature dependent quantities
  void UpdateTemperatureDependentCoefs(double T);

  // add stoichiometric contribution of complex to total
  void AddContributionToTotal(std::vector<double> *total);

  // add derivative of total with respect to free-ion to dtotal
  void AddContributionToDTotal(const std::vector<Species>& primarySpecies, MatrixBlock* dtotal);

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  std::string name() const { return name_; }

  int free_site_id() const { return free_site_id_; }
  double free_site_stoichiometry() const { return free_site_stoichiometry_; }
  double stoichiometry(int i) const { return stoichiometry_[i]; }

  int ncomp() const { return ncomp_; };
  int species_id(int i) const { return species_ids_[i]; };
  double surface_concentration() const { return surface_concentration_; };

 private:
  std::string name_;
  int identifier_;
  double charge_;

  double surface_concentration_;  // units? ?[mol/m^3 bulk]?

  int ncomp_;  // numebr components in reaction
  std::vector<std::string> species_names_;
  std::vector<int> species_ids_;  // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn

  std::string free_site_name_;
  double free_site_stoichiometry_;  // stoichiometry of free site in rxn
  int free_site_id_;

  double h2o_stoichiometry_;

  double lnK_;  // log value of equlibrium constant
  double lnQK_;  // store lnQK for derivatives later
  double logK_;
  Teuchos::RCP<FunctionTabular> func_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
