/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for surface complexation reaction

  Notes:
  - Each instance of this class should contain a single unique
    surface site (e.g. >FeOH) and ALL surface complexes associated
    with that site!
*/

#ifndef AMANZI_CHEMISTRY_SURFACE_COMPLEXATION_RXN_HH_
#define AMANZI_CHEMISTRY_SURFACE_COMPLEXATION_RXN_HH_

#include <vector>

#include "SurfaceComplex.hh"
#include "SurfaceSite.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class SurfaceComplexationRxn {
 public:
  SurfaceComplexationRxn();
  SurfaceComplexationRxn(SurfaceSite* surface_sites,
                         const std::vector<SurfaceComplex>& surface_complexes);
  explicit SurfaceComplexationRxn(SurfaceSite surface_sites);
  ~SurfaceComplexationRxn() {};

  // add complexes to the reaction
  void AddSurfaceComplex(SurfaceComplex surface_complex);
  void UpdateSiteDensity(double d);
  double GetSiteDensity() const { return surface_site_.at(0).molar_density(); }
  int SiteId() const { return surface_site_.at(0).identifier(); }

  double free_site_concentration() const {
    return surface_site_.at(0).free_site_concentration();
  }

  void set_free_site_concentration(const double value) {
    surface_site_.at(0).set_free_site_concentration(value);
  }

  // update sorbed concentrations
  void Update(const std::vector<Species>& primarySpecies);

  // add stoichiometric contribution of complex to sorbed total
  void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to sorbed dtotal
  void AddContributionToDTotal(const std::vector<Species>& primarySpecies,
                               MatrixBlock* dtotal);
  // If the free site stoichiometry in any of the surface complexes
  // is not equal to 1., we must use Newton's method to solve for
  // the free site concentration.  This function determines if this
  // is the case.
  void SetNewtonSolveFlag();

  void display(const Teuchos::Ptr<VerboseObject> vo) const;
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplaySite(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayComplexes(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

 private:
  std::vector<SurfaceComplex> surface_complexes_;
  std::vector<SurfaceSite> surface_site_;
  bool use_newton_solve_;

  // std::vector<double> dSx_dmi_;  // temporary storage for derivative calculations
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
