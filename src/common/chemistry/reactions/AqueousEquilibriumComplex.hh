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

#include "SecondarySpecies.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class AqueousEquilibriumComplex : public SecondarySpecies {
 public:
  AqueousEquilibriumComplex() : SecondarySpecies() {};
  AqueousEquilibriumComplex(int id, const std::string& name,
                            const Teuchos::ParameterList& plist,
                            const std::vector<Species>& primary_species);
  ~AqueousEquilibriumComplex() {};

  // update molalities
  virtual void Update(const std::vector<Species>& primary_species, const Species& water_species);

  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double>* total);

  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(
          const std::vector<Species>& primary_species, MatrixBlock* dtotal);

  using SecondarySpecies::Display;
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
