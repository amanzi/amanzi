#ifndef __Geochemistry_hpp__
#define __Geochemistry_hpp__

#include "Species.hpp"
#include "AqueousEquilibriumComplex.hpp"
//#include "MineralReaction.hpp"
//#include "GasExchange.hpp"
//#include "IonExchange.hpp"
//#include "SurfaceComplexation.hpp"
#include "Block.hpp"
#include "LU.hpp"

#include <vector>

// Driver class for evalating geochemical related processes at a
// single computational node

class Geochemistry {

public:
  Geochemistry();
  ~Geochemistry();

  // inheriting classes setup the species, etc
  virtual void setup(std::vector<double> *total);

  // speciate for free-ion concentrations
  int speciate(std::vector<double> target_totals);

  void initializeMolalities(double initial_molality);
  // update activities, equilibrium complex concentrations, etc.
  void updateChemistry(void);
  void calculateTotal(std::vector<double> &total);
  // calculate block of Jacobian corresponding to derivatives of total with
  // respect to free-ion
  void calculateDTotal(Block *dtotal);
  void scaleRHSAndJacobian(double *rhs, Block *J);

  void set_ncomp(int i) { this->ncomp_ = i; }
  void addPrimarySpecies(Species s);
  void addAqueousEquilibriumComplex(AqueousEquilibriumComplex c);

  void display(void) const;
  void print_results(void) const;
  
  int ncomp(void) const { return this->ncomp_; };
  void ncomp(const int s_ncomp) { this->ncomp_ = s_ncomp; };

  void verbose(const int s_verbose) { this->verbose_ = s_verbose; };
  int verbose(void) const { return this->verbose_; };
private:
  int verbose_;
  int ncomp_;                // # basis species
  std::vector<double> totals_;
  
  std::vector<Species> primarySpecies_; // list of primary species
  Species water_;

//  std::vector<ActivityCoefficient*> activityCoefficients_;

  std::vector<AqueousEquilibriumComplex> aqComplexRxns_; // list of aqueous equilibrium complexation reactions
//  vector<MineralReaction*> mineralRxns_;
//  vector<GasExchange*> gasRxns_;
//  vector<IonExchange*> ionExchangeRxns_;
//  vector<SurfaceComplexation*> surfaceComplexationRxns_;

};

#endif
