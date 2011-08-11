/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_AQUEOUS_EQUILIBRIUM_COMPLEX_HH_
#define AMANZI_CHEMISTRY_AQUEOUS_EQUILIBRIUM_COMPLEX_HH_

// Class for aqueous equilibrium complexation reaction

#include <vector>

#include "secondary_species.hh"

// forward declarations from outside chemistry

namespace amanzi {
namespace chemistry {

// forward declarations from chemistry
class Block;

class AqueousEquilibriumComplex : public SecondarySpecies {
 public:
  AqueousEquilibriumComplex();
  AqueousEquilibriumComplex(const SpeciesName name,
                            const SpeciesId id,
                            const std::vector<SpeciesName>& species,
                            const std::vector<double>& stoichiometry,
                            const std::vector<SpeciesId>& species_ids,
                            const double h2o_stoich,
                            const double charge,
                            const double mol_wt,
                            const double size,
                            const double logK);
  ~AqueousEquilibriumComplex();

  // update molalities
  virtual void Update(const std::vector<Species>& primary_species, const Species& water_species);
  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double>* total);
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<Species>& primary_species,
                                       Block* dtotal);

  void display(void) const;
  void Display(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

 protected:

 private:
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_AQUEOUS_EQUILIBRIUM_COMPLEX_HH_
