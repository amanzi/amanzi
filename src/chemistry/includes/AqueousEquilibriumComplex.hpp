#ifndef __AqueousEquilibriumComplex_hpp__
#define __AqueousEquilibriumComplex_hpp__

#include <string>
#include <vector>
#include <cmath>

#include "Species.hpp"
#include "SecondarySpecies.hpp"
#include "Block.hpp"

// Class for aqueous equilibrium complexation reaction

class AqueousEquilibriumComplex : public SecondarySpecies {

 public:
  AqueousEquilibriumComplex();
  AqueousEquilibriumComplex(std::string s);
  AqueousEquilibriumComplex(const SpeciesName name, 
                            const SpeciesId id,
                            std::vector<SpeciesName> species,
                            std::vector<double> stoichiometries,
                            std::vector<int> species_ids,
                            const double h2o_stoich, 
                            const double charge, 
                            const double mol_wt,
                            const double size, 
                            const double logK);
  ~AqueousEquilibriumComplex();

  // update molalities
  void Update_kludge(const std::vector<Species>primarySpecies);
  // add stoichiometric contribution of complex to total
  void AddContributionToTotal(std::vector<double> &total);
  // add derivative of total with respect to free-ion to dtotal
  void AddContributionToDTotal(const std::vector<Species> primarySpecies,
                               Block *dtotal);

  void display(void) const;
  void Display(void) const;

 protected:

 private:

};

#endif // __AqueousEquilibriumComplex_hpp__
