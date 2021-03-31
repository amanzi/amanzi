/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for aqueous equilibrium complexation reaction
*/

#ifndef AMANZI_CHEMISTRY_AQUEOUS_EQUILIBRIUM_COMPLEX_HH_
#define AMANZI_CHEMISTRY_AQUEOUS_EQUILIBRIUM_COMPLEX_HH_

#include <vector>

#include "secondary_species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class AqueousEquilibriumComplex : public SecondarySpecies {
 public:
  AqueousEquilibriumComplex() : SecondarySpecies() {};
  AqueousEquilibriumComplex(const std::string& name,
                            const int id,
                            const std::vector<std::string>& species,
                            const std::vector<double>& stoichiometry,
                            const std::vector<int>& species_ids,
                            const double h2o_stoich,
                            const double charge,
                            const double mol_wt,
                            const double size,
                            const double logK);
  ~AqueousEquilibriumComplex() {};

  // update molalities
  using SecondarySpecies::Update;
  virtual void Update(const std::vector<Species>& primary_species, const Species& water_species);

  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double>* total);

  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(
          const std::vector<Species>& primary_species, MatrixBlock* dtotal);

  using SecondarySpecies::Display;
  void display(const Teuchos::Ptr<VerboseObject> vo) const;
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
