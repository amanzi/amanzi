/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Base class for secondary species (aqueous equilibrium complexes,
  minerals.
*/
 
#ifndef AMANZI_CHEMISTRY_SECONDARY_SPECIES_HH_
#define AMANZI_CHEMISTRY_SECONDARY_SPECIES_HH_

#include <vector>

#include "species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class SecondarySpecies : public Species {
 public:
  SecondarySpecies();
  SecondarySpecies(const SpeciesName name,
                   const SpeciesId secondary_id,
                   const std::vector<SpeciesName>& species,
                   const std::vector<double>& stoichiometries,
                   const std::vector<SpeciesId>& species_ids,
                   const double h2o_stoich,
                   const double charge,
                   const double mol_wt,
                   const double size,
                   const double logK);

  virtual ~SecondarySpecies() {};

  // update molalities
  virtual void Update(const std::vector<Species>& primary_species);
  // add stoichiometric contribution of complex to totalmembermember
  virtual void AddContributionToTotal(std::vector<double> *total) = 0;
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<Species>& primary_species,
                                       MatrixBlock* dtotal) = 0;

  int ncomp() const { return ncomp_; };

  void set_logK(const double in_logK) { logK_ = in_logK; };
  double logK() const { return logK_; };
  double lnK() const { return lnK_; };
  double lnQK() const { return lnQK_; };

  std::vector<SpeciesName> species_names() const { return species_names_; };
  std::vector<SpeciesId> species_ids() const { return species_ids_; };
  std::vector<double> stoichiometry() const { return stoichiometry_; };

  virtual void Display() const;

 protected:
  int ncomp_;  // # components in reaction
  std::vector<SpeciesName> species_names_;
  std::vector<SpeciesId> species_ids_;  // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn
  std::vector<double> logK_array_;  // for temperature dep. logK
  double h2o_stoich_;  // stoichiometry of water in equation
  double lnK_;  // log value of equlibrium constant
  double lnQK_;  // store lnQK for derivatives later
  double logK_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif