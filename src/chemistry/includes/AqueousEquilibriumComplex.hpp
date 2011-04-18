#ifndef __AqueousEquilibriumComplex_hpp__
#define __AqueousEquilibriumComplex_hpp__

// Class for aqueous equilibrium complexation reaction

#include <vector>

#include "SecondarySpecies.hpp"

class Block;

class AqueousEquilibriumComplex : public SecondarySpecies {

 public:
  AqueousEquilibriumComplex();
  AqueousEquilibriumComplex(const SpeciesName name, 
                            const SpeciesId id,
                            std::vector<SpeciesName> species,
                            std::vector<double> stoichiometry,
                            std::vector<SpeciesId> species_ids,
                            const double h2o_stoich, 
                            const double charge, 
                            const double mol_wt,
                            const double size, 
                            const double logK);
  ~AqueousEquilibriumComplex();

  // update molalities
  virtual void Update(const std::vector<Species>& primary_species);
  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<Species> primary_species,
                               Block* dtotal);
  
  void display(void) const;
  void Display(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

 protected:

 private:


};

#endif // __AqueousEquilibriumComplex_hpp__
