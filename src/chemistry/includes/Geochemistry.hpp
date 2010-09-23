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

  // speciate for free-ion concentrations
  int speciate(std::vector<double> totals_);

  // solve a chemistry step
  int react(std::vector<double> total, double volume, double porosity, 
            double saturation, double dt);

  void initializeMolalities(double initial_molality);
  // update activities, equilibrium complex concentrations, etc.
  void updateChemistry(void);
  void calculateTotal(std::vector<double> &total);
  // calculate block of Jacobian corresponding to derivatives of total with
  // respect to free-ion
  void calculateDTotal(Block *dtotal);
  // accumulation terms
  void calculateAccumulation(std::vector<double> total, double *residual,
                             double volume, double porosity,
                             double saturation, double dt);
  void calculateAccumulationDerivative(Block *dtotal,Block *J,
                                       double volume, double porosity,
                                       double saturation, double dt);

  void scaleRHSAndJacobian(double *rhs, Block *J);

  void addPrimarySpecies(Species s);
  void addAqueousEquilibriumComplex(AqueousEquilibriumComplex c);

  void set_ncomp(int i) { this->ncomp_ = i; }
  
  int get_ncomp(void) const { return this->ncomp_; }

private:

  int ncomp_;                // # basis species

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
