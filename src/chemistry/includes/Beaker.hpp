#ifndef __Beaker_hpp__
#define __Beaker_hpp__

#include "Species.hpp"
#include "AqueousEquilibriumComplex.hpp"
#include "GeneralRxn.hpp"
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

  void addPrimarySpecies(Species s);
  void addAqueousEquilibriumComplex(AqueousEquilibriumComplex c);
  void addGeneralRxn(GeneralRxn r);

  // speciate for free-ion concentrations
  int speciate(std::vector<double> target_total, double water_density);

  // solve a chemistry step
  int react(std::vector<double> &total, double porosity, double saturation, 
            double density, double volume, double dt);

  void initializeMolalities(double initial_molality);
  void initializeMolalities(std::vector<double> initial_molalities);

  // equilibrium chemistry
  // update activities, equilibrium complex concentrations, etc.
  void updateEquilibriumChemistry(void);
  void calculateTotal(void);
  void calculateTotal(std::vector<double> &total);
  // calculate block of Jacobian corresponding to derivatives of total with
  // respect to free-ion
  void calculateDTotal(void);
  void calculateDTotal(Block *dtotal);
  // kinetic chemistry
  void updateKineticChemistry(void);
  void addKineticChemistryToResidual(std::vector<double> &residual);
  void addKineticChemistryToJacobian(Block *J);
  // accumulation terms
  void addAccumulation(std::vector<double> &residual);
  void addAccumulation(std::vector<double> total, 
                       std::vector<double> &residual);
  void addAccumulationDerivative(Block *J);
  void addAccumulationDerivative(Block *J, Block *dtotal);
  void calculateFixedAccumulation(std::vector<double> total,
                                  std::vector<double> &fixed_accumulation);
  // residual and Jacobian
  void calculateResidual(std::vector<double> &residual, 
                         std::vector<double> fixed_residual);
  void calculateJacobian(Block *J);

  // utilities for updating solution, convergence checks
  void updateMolalitiesWithTruncation(std::vector<double> &update, 
                                      std::vector<double> &prev_solution,
                                      double max_change);
  void calculateMaxRelChangeInMolality(std::vector<double> prev_molal, 
                                       double &max_rel_change);
  // solvers
  void scaleRHSAndJacobian(double *rhs, Block *J);
  void scaleRHSAndJacobian(std::vector<double> &rhs, Block *J);
  void solveLinearSystem(Block *A, std::vector<double> &b);

  void display(void) const;
  void print_results(void) const;
  void print_results(double time) const;
  void print_linear_system(string s, Block *A, std::vector<double> vector);

  void ncomp(int i) { this->ncomp_ = i; }
  int ncomp(void) const { return this->ncomp_; }

  // update discretization and flow parameters
  // water_density [kg/m^3]
  // volume [m^3]
  void updateParameters(double por, double sat, double den, double vol, double dt);
 
  void porosity(double d) { this->porosity_ = d; }
  void saturation(double d) { this->saturation_ = d; }
  // updates both water density variables
  void water_density_kg_m3(double d) { this->water_density_kg_m3_ = d; 
                                       this->water_density_kg_L_ = d/1000.; }
  void water_density_kg_L(double d) { this->water_density_kg_m3_ = d*1000.; 
                                      this->water_density_kg_L_ = d; }
  void volume(double d) { this->volume_ = d; }
  void dt(double d) { this->dt_ = d; }
  // calculates the coefficient in aqueous portion of accumulation term
  void accumulation_coef(double por, double sat, double vol, double dtt) 
                         { this->accumulation_coef_ = por*sat*vol*1000./dtt; }
  // calculates product of porosity,saturation,water_density[kg/m^3],volume
  void por_sat_den_vol(double por, double sat, double den, double vol) 
                       { this->por_sat_den_vol_ = por*sat*den*vol; }

  double porosity(void) const { return this->porosity_; }
  double saturation(void) const { return this->saturation_; }
  double water_density_kg_m3(void) const { return this->water_density_kg_m3_; }
  double water_density_kg_L(void) const { return this->water_density_kg_L_; }
  double volume(void) const { return this->volume_; }
  double dt(void) const { return this->dt_; }

  double accumulation_coef(void) const { return this->accumulation_coef_; }
  double por_sat_den_vol(void) const { return this->por_sat_den_vol_; }

  void verbose(const int s_verbose) { this->verbose_ = s_verbose; };
  int verbose(void) const { return this->verbose_; };

private:
  int verbose_;
  int ncomp_;                   // # basis species
  std::vector<double> total_;  // total component concentrations of basis species
  Block *dtotal_;      // matrix that holds derivative of total concentration w/respec to free-ion

  // common parameters among reactions
  double porosity_;            // [m^3 pore / m^3 bulk]
  double saturation_;          // [m^3 water / m^3 pore] 
  double water_density_kg_m3_; // [kg water / m^3 water]
  double water_density_kg_L_;  // [kg water / L water]
  double volume_;              // cell volume [m^3 bulk]
  double dt_;                  // time step size [seconds]
  // accumulation_coef_ = porosity*saturation*volume*1000./dt [L water/sec]
  // units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = (L water/sec)
  double accumulation_coef_;  
  // por_sat_den_vol_ = porosity * saturation * water_density * volume [kg water]
  double por_sat_den_vol_; 
  
  std::vector<Species> primarySpecies_; // list of primary species
  Species water_;

//  std::vector<ActivityCoefficient*> activityCoefficients_;

  std::vector<AqueousEquilibriumComplex> aqComplexRxns_; // list of aqueous equilibrium complexation reactions
  std::vector<GeneralRxn> generalKineticRxns_; //list of general kinetic reactions
//  vector<MineralReaction*> mineralRxns_;
//  vector<GasExchange*> gasRxns_;
//  vector<IonExchange*> ionExchangeRxns_;
//  vector<SurfaceComplexation*> surfaceComplexationRxns_;

  // solver data structures
  std::vector<double> fixed_accumulation; // fixed (time t) portion of accumulation term
  std::vector<double> residual;       // entire residual
  std::vector<double> prev_molal;     // previous molality of primary species

  std::vector<double> rhs;            // right-hand-side of system
  std::vector<int> indices;           // array for pivoting in LU
  Block *J;                           // Jacobian

};

#endif // __Beaker_hpp__
