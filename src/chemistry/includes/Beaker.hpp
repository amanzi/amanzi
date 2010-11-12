#ifndef __Beaker_hpp__
#define __Beaker_hpp__

// Driver class for evaluating geochemical related processes at a
// single computational node

#include <string>

#include "ActivityModel.hpp"
#include "ActivityModelFactory.hpp"
#include "AqueousEquilibriumComplex.hpp"
#include "Block.hpp"
#include "GeneralRxn.hpp"
#include "IonExchangeComplex.hpp"
#include "IonExchangeSite.hpp"
#include "KineticRate.hpp"
#include "LU.hpp"
#include "Mineral.hpp"
#include "MineralKineticsFactory.hpp"
#include "Species.hpp"
#include "SurfaceComplexationRxn.hpp"
#include "Verbosity.hpp"

#include <vector>

class KineticRate;

class Beaker {

 public:
  Beaker();
  virtual ~Beaker();

  struct BeakerComponents {
    std::vector<double> free_ion;
    std::vector<double> minerals;
    std::vector<double> ion_exchange_sites;
    std::vector<double> total;
    std::vector<double> total_sorbed;
  };

  struct BeakerParameters {
    // input files
    std::string thermo_database_file;
    // solver parameters
    double tolerance;
    unsigned int max_iterations;
    // models
    std::string activity_model_name;
    // physical parameters
    double porosity; // [-]
    double saturation; // [-]
    double water_density; // [kg/m^3]
    double volume; // [m^3]
  };

  // resizes matrix and vectors for nonlinear system
  void resize();
  void resize(int ncomp);

  // inheriting classes setup the species, etc
  virtual void Setup(const Beaker::BeakerComponents& components,
                     const Beaker::BeakerParameters& parameters);
  void SetupActivityModel(std::string model);
  void VerifyComponentSizes(const Beaker::BeakerComponents& components);
  void SetComponents(const Beaker::BeakerComponents& components);

  void addPrimarySpecies(Species s);
  void AddIonExchangeSite(IonExchangeSite exchanger);
  void AddIonExchangeComplex(IonExchangeComplex exchange_complex);
  void addAqueousEquilibriumComplex(AqueousEquilibriumComplex c);
  void addMineral(Mineral m);
  void AddMineralKineticRate(KineticRate* rate);
  void addGeneralRxn(GeneralRxn r);
  void addSurfaceComplexationRxn(SurfaceComplexationRxn r);

  bool HaveKinetics(void) const;

  // speciate for free-ion concentrations
  int Speciate(const BeakerComponents& components, const BeakerParameters& parameters);

  // solve a chemistry step
  int ReactionStep(BeakerComponents* components, const BeakerParameters& parameters, 
		   double dt);

  void updateActivityCoefficients();
  void initializeMolalities(double initial_molality);
  void initializeMolalities(std::vector<double> initial_molalities);

  // equilibrium chemistry
  // update activities, equilibrium complex concentrations, etc.
  void updateEquilibriumChemistry(void);
  void calculateTotal(void);
  void calculateTotal(std::vector<double> *total,
                      std::vector<double> *total_sorbed);
  // calculate block of Jacobian corresponding to derivatives of total with
  // respect to free-ion
  void calculateDTotal(void);
  void calculateDTotal(Block *dtotal, Block *dtotal_sorbed);
  // kinetic chemistry
  void updateKineticChemistry(void);
  void addKineticChemistryToResidual(std::vector<double> *residual);
  void addKineticChemistryToJacobian(Block *J);
  // accumulation terms
  void addAccumulation(std::vector<double> *residual);
  void addAccumulation(std::vector<double> total, 
                       std::vector<double> total_sorbed, 
                       std::vector<double> *residual);
  void addAccumulationDerivative(Block *J);
  void addAccumulationDerivative(Block *J, Block *dtotal, Block *dtotal_sorbed);
  void calculateFixedAccumulation(std::vector<double> total,
                                  std::vector<double> total_sorbed,
                                  std::vector<double> *fixed_accumulation);
  // residual and Jacobian
  void calculateResidual(std::vector<double> *residual, 
                         std::vector<double> fixed_residual);
  void calculateJacobian(Block *J);

  // utilities for updating solution, convergence checks
  void updateMolalitiesWithTruncation(std::vector<double> &update, 
                                      std::vector<double> &prev_solution,
                                      double max_change);
  double calculateMaxRelChangeInMolality(std::vector<double> prev_molal);
  // solvers
  void scaleRHSAndJacobian(double *rhs, Block *J);
  void scaleRHSAndJacobian(std::vector<double> &rhs, Block *J);
  void solveLinearSystem(Block *A, std::vector<double> &b);

  virtual void display(void) const;
  void print_results(void) const;
  void print_results(double time) const;
  void print_linear_system(string s, Block *A, std::vector<double> vector);

  void Display(void) const;
  void DisplayParameters(void) const;
  void DisplayPrimary(void) const;
  void DisplayAqueousEquilibriumComplexes(void) const;
  void DisplayMinerals(void) const;
  void DisplayMineralKinetics(void) const;
  void DisplayIonExchangeSites(void) const;
  void DisplayIonExchangeComplexes(void) const;

  void DisplayTotalColumnHeaders(void) const;
  void DisplayTotalColumns(const double time, const std::vector<double>& total) const;
  void DisplayResults(void) const;

  void ncomp(int i) { this->ncomp_ = i; }
  int ncomp(void) const { return this->ncomp_; }

  double tolerance(void) const { return this->tolerance_; }
  unsigned int max_iterations(void) const { return this->max_iterations_; }
  double porosity(void) const { return this->porosity_; }
  double saturation(void) const { return this->saturation_; }
  double water_density_kg_m3(void) const { return this->water_density_kg_m3_; }
  double water_density_kg_L(void) const { return this->water_density_kg_L_; }
  double volume(void) const { return this->volume_; }
  double dt(void) const { return this->dt_; }

  double aqueous_accumulation_coef(void) const 
    { return this->aqueous_accumulation_coef_; }
  double sorbed_accumulation_coef(void) const 
    { return this->sorbed_accumulation_coef_; }
  double por_sat_den_vol(void) const { return this->por_sat_den_vol_; }

  virtual void verbosity(const Verbosity s_verbosity) { this->verbosity_ = s_verbosity; };
  virtual Verbosity verbosity(void) const { return this->verbosity_; };

  BeakerParameters GetDefaultParameters(void) const;
  BeakerParameters GetCurrentParameters(void) const;
  void SetParameters(const BeakerParameters& parameters);

protected:
  // update discretization and flow parameters
  // water_density [kg/m^3]
  // volume [m^3]
  void updateParameters(const BeakerParameters& parameters, double dt);

  void tolerance(double value) { this->tolerance_ = value; }
  void max_iterations(unsigned int value) { this->max_iterations_ = value; }
  void porosity(double d) { this->porosity_ = d; }
  void saturation(double d) { this->saturation_ = d; }
  // updates both water density variables
  void water_density_kg_m3(double d) { this->water_density_kg_m3_ = d; 
                                       this->water_density_kg_L_ = d/1000.; }
  void water_density_kg_L(double d) { this->water_density_kg_m3_ = d*1000.; 
                                      this->water_density_kg_L_ = d; }
  void volume(double d) { this->volume_ = d; }
  void dt(double d) { this->dt_ = d; }
  void aqueous_accumulation_coef(double d) 
    { this->aqueous_accumulation_coef_ = d; }
  void sorbed_accumulation_coef(double d) 
    { this->sorbed_accumulation_coef_ = d; }
  void por_sat_den_vol(double d) { this->por_sat_den_vol_ = d; }
  // calculates the coefficient in aqueous portion of accumulation term
  void update_accumulation_coefficients(void); 
  // calculates product of porosity,saturation,water_density[kg/m^3],volume
  void update_por_sat_den_vol(void);


  const std::vector<Species>& primary_species(void) const { return this->primarySpecies_; };
  const std::vector<Mineral>& minerals(void) const { return this->minerals_; };
  const std::vector<IonExchangeSite>& ion_exchange_sites(void) const { return this->ion_exchange_sites_; };
  const std::vector<double>& total(void) const { return this->total_; };
  const std::vector<double>& total_sorbed(void) const { return this->total_sorbed_; };

private:
  Verbosity verbosity_;
  double tolerance_;
  unsigned int max_iterations_;
  int ncomp_;                   // # basis species

  // Aqueous phase total component concentrations for basis species
  std::vector<double> total_; // [mol/L]
  // Matrix block containing derivative of total concentration w/respec to 
  // free-ion
  Block *dtotal_; // [kg water/sec]

  // Sorbed phase total component concentrations for basis species
  std::vector<double> total_sorbed_;  // [mol/m^3 bulk]
  // Matrix block containing derivative of total sorbed concentration 
  // w/respec to free-ion
  Block *dtotal_sorbed_; // [kg water/sec]

  // common parameters among reactions
  double porosity_;            // [m^3 pore / m^3 bulk]
  double saturation_;          // [m^3 water / m^3 pore] 
  double water_density_kg_m3_; // [kg water / m^3 water]
  double water_density_kg_L_;  // [kg water / L water]
  double volume_;              // cell volume [m^3 bulk]
  double dt_;                  // time step size [seconds]
  // aqueous_accumulation_coef_ = porosity*saturation*volume*1000./dt [L water/sec]
  // units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = (L water/sec)
  double aqueous_accumulation_coef_;  
  // sorbed_accumulation_coef_ = volume/dt [m^3 bulk/sec]
  double sorbed_accumulation_coef_;  
  // por_sat_den_vol_ = porosity * saturation * water_density * volume [kg water]
  double por_sat_den_vol_; 
  
  Species water_;
  std::vector<Species> primarySpecies_; // list of primary species
  std::vector<Mineral> minerals_; // list of mineral species
  std::vector<IonExchangeSite> ion_exchange_sites_;

  ActivityModel *activity_model_;

  std::vector<AqueousEquilibriumComplex> aqComplexRxns_; // list of aqueous equilibrium complexation reactions
  std::vector<GeneralRxn> generalKineticRxns_; //list of general kinetic reactions
  std::vector<KineticRate*> mineral_rates_;
//  vector<GasExchange*> gasRxns_;
  std::vector<IonExchangeComplex> ion_exchange_rxns_;
  std::vector<SurfaceComplexationRxn> surfaceComplexationRxns_;

  // solver data structures
  std::vector<double> fixed_accumulation; // fixed (time t) portion of accumulation term
  std::vector<double> residual;       // entire residual [mol/sec]
  std::vector<double> prev_molal;     // previous molality of primary species

  std::vector<double> rhs;            // right-hand-side of system
  std::vector<int> indices;           // array for pivoting in LU
  Block *J;                           // Jacobian [kg water/sec]

  
  static const double tolerance_default;
  static const unsigned int max_iterations_default;
  static const double porosity_default;
  static const double saturation_default;
  static const double water_density_kg_m3_default;
  static const double volume_default;

};

#endif // __Beaker_hpp__
