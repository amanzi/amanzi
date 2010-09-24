#ifndef __Beaker_hpp__
#define __Beaker_hpp__

#include "Species.hpp"
#include "AqueousEquilibriumComplex.hpp"
//#include "MineralReaction.hpp"
//#include "GasExchange.hpp"
//#include "IonExchange.hpp"
//#include "SurfaceComplexation.hpp"
#include "Block.hpp"
#include "LU.hpp"

#include <vector>

// Driver class for evaluating geochemical related processes at a
// single computational node

class Beaker {

 public:
  Beaker();
  ~Beaker();

  // resizes matrix and vectors for nonlinear system
  void resize();
  void resize(int ncomp);

  // inheriting classes setup the species, etc
  virtual void setup(std::vector<double> &total);

  // speciate for free-ion concentrations
  int speciate(std::vector<double> target_total);

  // solve a chemistry step
  int react(std::vector<double> &total, double volume, double porosity, 
            double saturation, double dt);

  void initializeMolalities(double initial_molality);
  void initializeMolalities(std::vector<double> initial_molalities);
  // update activities, equilibrium complex concentrations, etc.
  void updateEquilibriumChemistry(void);
  void calculateTotal(void);
  void calculateTotal(std::vector<double> &total);
  // calculate block of Jacobian corresponding to derivatives of total with
  // respect to free-ion
  void calculateDTotal(void);
  void calculateDTotal(Block *dtotal);
  // accumulation terms
  void calculateAccumulation(std::vector<double> &residual);
  void calculateAccumulation(std::vector<double> total, 
                             std::vector<double> &residual);
  void calculateAccumulationDerivative(Block *J);
  void calculateAccumulationDerivative(Block *dtotal, Block *J);

  void updateMolalitiesWithTruncation(std::vector<double> &update, 
                                      std::vector<double> &prev_solution,
                                      double max_change);
  void calculateMaxRelChangeInMolality(std::vector<double> prev_molal, 
                                       double &max_rel_change);

  void scaleRHSAndJacobian(double *rhs, Block *J);
  void scaleRHSAndJacobian(std::vector<double> &rhs, Block *J);
  void solveLinearSystem(Block *A, std::vector<double> &b);

  void addPrimarySpecies(Species s);
  void addAqueousEquilibriumComplex(AqueousEquilibriumComplex c);


  void display(void) const;
  void print_results(void) const;
  void print_linear_system(string s, Block *A, std::vector<double> vector);

  void ncomp(int i) { this->ncomp_ = i; }
  int ncomp(void) const { return this->ncomp_; }

  void updateParameters(double por, double sat, double vol, double dt);
 
  void porosity(double d) { this->porosity_ = d; }
  void saturation(double d) { this->saturation_ = d; }
  void volume(double d) { this->volume_ = d; }
  void dt(double d) { this->dt_ = d; }

  double porosity(void) const { return this->porosity_; }
  double saturation(void) const { return this->saturation_; }
  double volume(void) const { return this->volume_; }
  double dt(void) const { return this->dt_; }

  void accumulation_coef(double por, double sat, double vol, double dtt) 
                         { this->accumulation_coef_ = por*sat*vol*1000./dtt; }
  double accumulation_coef(void) const { return this->accumulation_coef_; }

  void verbose(const int s_verbose) { this->verbose_ = s_verbose; };
  int verbose(void) const { return this->verbose_; };

private:
  int verbose_;
  int ncomp_;                   // # basis species
  std::vector<double> total_;  // total component concentrations of basis species
  Block *dtotal_;      // matrix that holds derivative of total concentration w/respec to free-ion

  // common parameters among reactions
  double porosity_;          // [m^3 pore / m^3 bulk]
  double saturation_;        // [m^3 water / m^3 pore] 
  double volume_;            // cell volume [m^3 bulk]
  double dt_;                // time step size [seconds]
  // accumulation_coef_ = porosity*saturation*volume*1000./dt [L water/sec]
  // units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = (L water/sec)
  double accumulation_coef_; 
  
  std::vector<Species> primarySpecies_; // list of primary species
  Species water_;

//  std::vector<ActivityCoefficient*> activityCoefficients_;

  std::vector<AqueousEquilibriumComplex> aqComplexRxns_; // list of aqueous equilibrium complexation reactions
//  vector<MineralReaction*> mineralRxns_;
//  vector<GasExchange*> gasRxns_;
//  vector<IonExchange*> ionExchangeRxns_;
//  vector<SurfaceComplexation*> surfaceComplexationRxns_;

  // solver data structures

  std::vector<double> fixed_residual; // fixed (time t) portion of residual
  std::vector<double> residual;       // entire residual
  std::vector<double> prev_molal;     // previous molality of primary species

  std::vector<double> rhs;            // right-hand-side of system
  std::vector<int> indices;           // array for pivoting in LU
  Block *J;                           // Jacobian

};

#endif // __Beaker_hpp__
