#ifndef __AqueousEquilibriumComplex_hpp__
#define __AqueousEquilibriumComplex_hpp__

#include <string>
#include <vector>
#include <math.h>

#include "Species.hpp"
#include "Block.hpp"

// Class for aqueous equilibrium complexation reaction

class AqueousEquilibriumComplex : public Species {

public:
  AqueousEquilibriumComplex();
  AqueousEquilibriumComplex(std::string s);
  AqueousEquilibriumComplex(SpeciesName name, 
                            std::vector<SpeciesName> species,
                            std::vector<double> stoichiometries,
                            std::vector<int> species_ids,
                            double h2o_stoich, double charge, double mol_wt,
                            double size, double logK);
  ~AqueousEquilibriumComplex();

  // update molalities
  void update(const std::vector<Species>primarySpecies);
  // add stoichiometric contribution of complex to total
  void addContributionToTotal(std::vector<double> &total);
  // add derivative of total with respect to free-ion to dtotal
  void addContributionToDTotal(const std::vector<Species> primarySpecies,
                               Block *dtotal);

  void display(void) const;

protected:

private:

  int ncomp_; // # components in reaction
  std::vector<SpeciesName> species_names_;
  std::vector<int> species_ids_;       // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn
  std::vector<double> logK_array_;     // for temperature dep. logK
  double h2o_stoich_;
  double lnK_;
  double lnQK_;                        // store lnQK for derivatives later
  double logK_;
  
  double log_to_ln(double d) { return d*2.30258509299; }
  double ln_to_log(double d) { return d*0.434294481904; }

};

#endif
