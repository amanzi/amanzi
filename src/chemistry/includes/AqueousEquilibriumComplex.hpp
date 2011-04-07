#ifndef __AqueousEquilibriumComplex_hpp__
#define __AqueousEquilibriumComplex_hpp__

// Class for aqueous equilibrium complexation reaction

#include <vector>

#include "Species.hpp"

class Block;

class AqueousEquilibriumComplex : public Species {

 public:
  AqueousEquilibriumComplex();
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
  void Update(const std::vector<Species> primary_species);
  // add stoichiometric contribution of complex to total
  void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to dtotal
  void AddContributionToDTotal(const std::vector<Species> primary_species,
                               Block* dtotal);
  
  void set_ncomp(const int ncomp) { this->ncomp_ = ncomp; };
  int ncomp(void) const { return this->ncomp_; };

  void set_logK(const double logK) { this->logK_ = logK; };
  double logK(void) const { return this->logK_; };

  void display(void) const;
  void Display(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

 protected:

 private:

  double log_to_ln(double d) { return d*2.30258509299; }
  double ln_to_log(double d) { return d*0.434294481904; }

  int ncomp_; // # components in reaction
  std::vector<SpeciesName> species_names_;
  std::vector<SpeciesId> species_ids_;       // ids of primary species in rxn
  std::vector<double> stoichiometry_;  // stoich of primary species in rxn
  std::vector<double> logK_array_;     // for temperature dep. logK
  double h2o_stoichiometry_;           // stoichiometry of water in equation
  double lnK_;                         // log value of equlibrium constant
  double lnQK_;                        // store lnQK for derivatives later
  double logK_;

};

#endif // __AqueousEquilibriumComplex_hpp__
