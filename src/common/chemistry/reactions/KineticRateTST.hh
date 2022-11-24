/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Implementation of the TST rate law for mineral kinetics
 
    R = k * A * Prod (a_i^m_i) * ( 1 - Q/Keq)
 
  where:
    R : reaction rate, [moles/sec]
    Keq : equilibrium constant, [-]
    Q : ion activity product, [-]
    k : reaction rate constant, [moles m^2 s^-1]
    A : reactive surface area, [m^2]
    a_i : activity of modifying species (primary or secondary)
    m_i : exponent of modifying species
*/

#ifndef AMANZI_CHEMISTRY_KINETIC_RATE_TST_HH_
#define AMANZI_CHEMISTRY_KINETIC_RATE_TST_HH_

#include <string>
#include <vector>

#include "VerboseObject.hh"

#include "KineticRate.hh"
#include "Mineral.hh"
#include "SecondarySpecies.hh"
#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

class MatrixBlock;

class KineticRateTST : public KineticRate {
 public:
  KineticRateTST();
  virtual ~KineticRateTST(){};

  void Setup(const Mineral& reaction,
             double rate,
             const std::string& modifiers,
             const SpeciesArray& primary_species);

  virtual void
  Update(const SpeciesArray& primary_species, const std::vector<Mineral>& minerals) override;

  virtual void AddContributionToResidual(const std::vector<Mineral>& minerals,
                                         double bulk_volume,
                                         std::vector<double>* residual) override;

  virtual void AddContributionToJacobian(const SpeciesArray& primary_species,
                                         const std::vector<Mineral>& minerals,
                                         double bulk_volume_vol,
                                         MatrixBlock* J) override;

  virtual void Display(const Teuchos::Ptr<VerboseObject> vo) const override;

 private:
  double area_;                // surface area [m^2]
  double log_Keq_;             // log_Keq [-]
  double rate_constant_;       // k, rate constant, [moles/m^2/sec]
  double log10_rate_constant_; // log10(k),

  double Q_over_Keq_;
  double modifying_term_;

  // length is the number of modifying species (primary and secondary
  // in the same list!)
  std::vector<std::string> modifying_species_names;
  std::vector<double> modifying_exponents;
  std::vector<int> modifying_primary_ids;   // not needed?
  std::vector<int> modifying_secondary_ids; // not needed?

  // these vectors have length primary_species.size() so that dot
  // products are easier. any unneeded values are set to zero.
  std::vector<double> primary_stoichiometry;
  std::vector<double> modifying_primary_exponents;
  // length secondary_species.size()
  std::vector<double> modifying_secondary_exponents;
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
