/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for mineral reaction, should be written with the mineral
  as the reactant:

    Calcite = 1.0 Ca++ + 1.0 HCO3- -1.0 H+
*/
 
#ifndef AMANZI_CHEMISTRY_MINERAL_HH_
#define AMANZI_CHEMISTRY_MINERAL_HH_

#include <cmath>
#include <vector>

#include "species.hh"
#include "secondary_species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class Mineral : public SecondarySpecies {
 public:
  Mineral();
  Mineral(const SpeciesName name,
          SpeciesId mineral_id,
          const std::vector<SpeciesName>& species,
          const std::vector<double>& stoichiometries,
          const std::vector<SpeciesId>& species_ids,
          const double h2o_stoich, const double mol_wt,
          const double logK, const double molar_volume,
          const double specific_surface_area);
  ~Mineral() {};

  // update molalities
  using SecondarySpecies::Update;
  virtual void Update(const std::vector<Species>& primary_species, const Species& water_species);
  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<Species>& primary_species,
                                       MatrixBlock* dtotal);

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  double Q_over_K(void) const {
    return std::exp(this->lnQK_);
  };
  double saturation_index(void) const {
    return std::log10(Q_over_K());
  };  // SI = log10(Q/Keq)

  double molar_volume(void) const {
    return this->molar_volume_;
  }

  double specific_surface_area(void) const { return this->specific_surface_area_; }
  void set_specific_surface_area(const double d) { this->specific_surface_area_ = d; }

  void UpdateSpecificSurfaceArea(void);

  double volume_fraction(void) const { return this->volume_fraction_; }
  void set_volume_fraction(const double d) { this->volume_fraction_ = d; }

  void UpdateVolumeFraction(const double rate, const double delta_time);

 private:
  double molar_volume_;     // [m^3 / moles]
  double specific_surface_area_;  // [m^2 mineral / m^3 bulk]
  double volume_fraction_;   // [m^3 mineral / m^3 bulk]
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
