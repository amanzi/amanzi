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
#include <string>
#include <vector>

#include "SecondarySpecies.hh"
#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class Mineral : public SecondarySpecies {
 public:
  Mineral();
  Mineral(int id,
          const std::string& name,
          const Teuchos::ParameterList& plist,
          const std::vector<Species>& primary_species);
  ~Mineral(){};

  // update molalities
  void Update(const std::vector<Species>& primary_species, const Species& water_species);

  // add stoichiometric contribution of complex to total
  void AddContributionToTotal(std::vector<double>* total);

  // add derivative of total with respect to free-ion to dtotal
  void AddContributionToDTotal(const std::vector<Species>& primary_species, MatrixBlock* dtotal);

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  double Q_over_K() const { return std::exp(lnQK_); };
  double saturation_index() const { return std::log10(Q_over_K()); }; // SI = log10(Q/Keq)

  double specific_surface_area() const { return specific_surface_area_; }
  void set_specific_surface_area(double d) { specific_surface_area_ = d; }

  double molar_volume() const { return molar_volume_; }

  // not supported yet
  void UpdateSpecificSurfaceArea(){};

  double volume_fraction() const { return volume_fraction_; }
  void set_volume_fraction(double d) { volume_fraction_ = d; }

  void UpdateVolumeFraction(double rate, double dt);

 private:
  double molar_volume_;          // [m^3 / moles]
  double specific_surface_area_; // [m^2 mineral / m^3 bulk]
  double volume_fraction_;       // [m^3 mineral / m^3 bulk]
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
