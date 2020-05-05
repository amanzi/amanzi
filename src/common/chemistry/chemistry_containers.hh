/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_CHEMISTRY_CONTAINERS_HH_
#define AMANZI_CHEMISTRY_CHEMISTRY_CONTAINERS_HH_

/*******************************************************************************
 **
 ** All the containers used within the chemistry code to pass data
 ** between different modules should be defined here.
 **
 ** NOTEs: 
 **
 ** - Containers can contain other containers, i.e. the options
 **   container for the chemistry coordinator will have an options
 **   container for the speciation evalautor. To keep the compiler
 **   happy, the low level objects must come first.
 **
 ** - To avoid circular dependancies, do not include any of the main
 **   chemistry classes in the containers in this file. All structs
 **   should contain only int, double, bool, strings and vectors of
 **   those.
 **
 ** TL;DR: The order of the containers in this file is important. 
 ** Containers should only include POD, vectors, strings and other
 ** containers.
 **
 ******************************************************************************/

#include <vector>
#include <string>


namespace Amanzi {
namespace AmanziChemistry {

/*******************************************************************************
 **
 ** Solver Options
 **
 ******************************************************************************/
struct NonlinearSolverOptions {
  NonlinearSolverOptions()
      : linear_solver_name(strings::kLU),
        system_size(0),
        under_relax(false),
        under_relax_coeff(1.0),
        max_iterations(50),
        tolerance_function(1.0e-10),
        tolerance_solution(1.0e-10) {
  }  
  std::string linear_solver_name;
  int system_size;
  bool under_relax;
  double under_relax_coeff;
  int max_iterations;
  double tolerance_function;
  double tolerance_solution;
};

struct ODESolverOptions {
  ODESolverOptions()
      : linear_solver_name(strings::kLU),
        system_size(0),
        weighting_factor(0.5),
        max_substeps(1) {  
  }
  std::string linear_solver_name;
  int system_size;
  double weighting_factor;  // 1=fully implicit, 0.5=crank-nicholson, 0=fully explicit
  int max_substeps;
};

/* TODO(bandre): do we want seperate status objects for nonlinear and ode solvers? */
struct SolverStatus {
  SolverStatus() {
    Reset();
  }
  void Reset(void) {
    num_rhs_evaluations = 0;
    num_jacobian_evaluations = 0;
    num_iterations = 0;
    converged = false;
  }
  unsigned int num_rhs_evaluations;
  unsigned int num_jacobian_evaluations;
  unsigned int num_iterations;  // or time substeps for ode solvers
  bool converged;
};

struct OutputOptions {
  OutputOptions()
      : use_stdout(true),
        file_name("") {
    verbosity_levels.clear();
  }
  bool use_stdout;
  std::string file_name;
  std::vector<std::string> verbosity_levels;
};

/*******************************************************************************
 **
 ** Database : format independent containers for info from
 ** thermodynamic databases, kinetics, etc
 **
 ******************************************************************************/
struct AqueousSpeciesData {
  AqueousSpeciesData(void) 
      : name(""),
        charge(0.0),
        gmw(0.0),
        ion_size(0.0) {}
  void Clear(void) {
    this->name = "";
    this->charge = 0.0;
    this->gmw = 0.0;
    this->ion_size = 0.0;
  };
  std::string name;
  double charge;  // [-]
  double gmw;  // [grams/mole]
  double ion_size;  // [angstroms]
};

struct AqueousReactionData {
  AqueousReactionData(void)
      : aqueous_species() {
    stoich_coeff.clear();
    species_names.clear();
    log_K_temperature.clear();
  }
  void Clear(void) {
    this->aqueous_species.Clear();
    this->stoich_coeff.clear();
    this->species_names.clear();
    this->log_K_temperature.clear();
  }
  AqueousSpeciesData aqueous_species;
  std::vector<double> stoich_coeff;
  std::vector<std::string> species_names;
  std::vector<double> log_K_temperature;
};

struct MineralSpeciesData {
  MineralSpeciesData(void) 
      : name(""),
        molar_volume(0.0),
        gmw(0.0),
        specific_surface_area(0.0),
        bulk_surface_area(0.0) {}
  void Clear(void) {
    this->name = "";
    this->molar_volume = 0.0;
    this->gmw = 0.0;
    this->specific_surface_area = 0.0;
    this->bulk_surface_area = 0.0;
  }
  std::string name;
  double molar_volume;  // [cm^3/mole]
  double gmw;  // [grams/mole]
  double specific_surface_area;  // [m^2/gram]
  double bulk_surface_area;  // [m^2/m^3]
};

struct MineralReactionData {
  MineralReactionData(void)
      : mineral_species() {
    stoich_coeff.clear();
    species_names.clear();
    log_K_temperature.clear();
  }
  void Clear(void) {
    this->mineral_species.Clear();
    this->stoich_coeff.clear();
    this->species_names.clear();
    this->log_K_temperature.clear();
  }
  MineralSpeciesData mineral_species;
  std::vector<double> stoich_coeff;
  std::vector<std::string> species_names;
  std::vector<double> log_K_temperature;
};

struct GasSpeciesData {
  GasSpeciesData(void) 
      : name(""),
        molar_volume(0.0),
        gmw(0.0) {}
  void Clear(void) {
    this->name = "";
    this->molar_volume = 0.0;
    this->gmw = 0.0;
  }
  std::string name;
  double molar_volume;  // [cm^3/mole]
  double gmw;  // [grams/mole]
};

struct GasReactionData {
  GasReactionData(void)
      : gas_species() {
    stoich_coeff.clear();
    species_names.clear();
    log_K_temperature.clear();
  }
  void Clear(void) {
    gas_species.Clear();
  this->stoich_coeff.clear();
  this->species_names.clear();
  this->log_K_temperature.clear();
} // end InitializeGasReaction
  GasSpeciesData gas_species;
  std::vector<double> stoich_coeff;
  std::vector<std::string> species_names;
  std::vector<double> log_K_temperature;
};

struct AqueousKineticsTSTData {
  AqueousKineticsTSTData(void)
      : rate(0.0) {
    prefactor_species.clear();
    prefactor_exponents.clear();
  }
  double rate;  // [mol/kgw/year]
  std::vector<std::string> prefactor_species;
  std::vector<double> prefactor_exponents;
};

struct AqueousKineticsIrreversibleData {
  AqueousKineticsIrreversibleData(void)
      : rate(0.0) {
    prefactor_species.clear();
    prefactor_exponents.clear();
  }
  double rate;  // [mol/kgw/year]
  std::vector<std::string> prefactor_species;
  std::vector<double> prefactor_exponents;
};

struct AqueousKineticsData {
  AqueousKineticsData(void)
      : name(""),
        rate_type(""),
        log_10_K(0.0) {
    stoich_coeff.clear();
    species_names.clear();
    tst_rates.clear();
  }
  std::string name;
  std::string rate_type;
  double log_10_K;  
  std::vector<double> stoich_coeff;
  std::vector<std::string> species_names;
  std::vector<AqueousKineticsTSTData> tst_rates;
  std::vector<AqueousKineticsIrreversibleData> irreversible_rates;
};

struct MineralKineticsTSTData {
  // see process models doc, Eq 5.129
  MineralKineticsTSTData(void)
      : label(""),
        mineral_name(""),
        rate(0.0),
        activation_energy(0.0),  // default, exp(0.0) = 1.0
        affinity_exponent_M(1.0),  // default x^1.0 = x
        affinity_exponent_n(1.0) {  // default y^1.0 = y
    prefactor_species.clear();
    prefactor_exponents.clear();
  }
  std::string label;
  std::string mineral_name;
  double rate;  // 25C [mol/kgw/year]
  double activation_energy;  // [kcal/mol]
  double affinity_exponent_M;
  double affinity_exponent_n;  // (1-(Q/K)^M)^n
  std::vector<std::string> prefactor_species;
  std::vector<double> prefactor_exponents;
};

struct MineralKineticsData {
  MineralKineticsData(void) {
    tst_rates.clear();
    //irreversible_rates.clear();
  }
  std::vector<MineralKineticsTSTData> tst_rates;
  //std::vector<MineralKineticsIrreversibleData> irreversible_rates;
};


struct DatabaseOptions {
  DatabaseOptions() 
      : file_name(""),
        type("") {
    aqueous_species_names.clear();
    aqueous_complex_names.clear();
    mineral_species_names.clear();
    gas_species_names.clear();
    aqueous_kinetics_names.clear();
    mineral_kinetics_names.clear();
  }
  std::string file_name;
  std::string type;
  std::vector<std::string> aqueous_species_names;
  std::vector<std::string> aqueous_complex_names;
  std::vector<std::string> mineral_species_names;
  std::vector<std::string> gas_species_names;
  std::vector<std::string> aqueous_kinetics_names;
  std::vector<std::string> mineral_kinetics_names;
};

/*******************************************************************************
 **
 ** Chemistry Module Options
 **
 ******************************************************************************/
struct ActivityModelParameters {
  ActivityModelParameters() 
      : activity_model_name(strings::kDebyeHuckelBdot),
        database_file_name(""),
        interpolator_name("") {
    temperature_points.clear();
    debye_huckel_a_points.clear();
    debye_huckel_b_points.clear();
    debye_huckel_bdot_points.clear();
  }
  // models
  std::string activity_model_name;
  // database file with model specific parameters
  std::string database_file_name;
  // parameters that come from the general database file
  std::string interpolator_name;
  std::vector<double> temperature_points;
  std::vector<double> debye_huckel_a_points;
  std::vector<double> debye_huckel_b_points;
  std::vector<double> debye_huckel_bdot_points;
};

struct ReactionRateParameters {
  ReactionRateParameters()
      : reaction_rate_type(""),
        rate_constant_25C(0.0),
        activation_energy(0.0) {
    prefactor_names.clear();
  }
  std::string reaction_rate_type;
  double rate_constant_25C;
  double activation_energy;  // [kcal/mole]
  std::vector<std::string> prefactor_names;
  
};

/*******************************************************************************
 **
 ** Chemical Constraint / Condition
 **
 ******************************************************************************/
struct ConcentrationConstraint {
  ConcentrationConstraint()
      : name(""),
        type(""),
        value(0.0) {}
  std::string name;
  std::string type;
  double value;
};

struct MineralConstraint {
  MineralConstraint()
      : name(""),
        volume_fraction(0.0),
        surface_area_type(""),
        surface_area_value(0.0) {}
  std::string name;
  double volume_fraction;
  std::string surface_area_type;
  double surface_area_value;
};

struct ChemicalCondition {
  ChemicalCondition()
      : name(""),
        temperature(25.0) {
    conc_constraints.clear();
    mineral_constraints.clear();
  }
  std::string name;
  double temperature;  // [celsius]
  std::vector<ConcentrationConstraint> conc_constraints;
  std::vector<MineralConstraint> mineral_constraints;
};

/*******************************************************************************
 **
 ** Evaluator Options
 **
 ******************************************************************************/
struct SpeciationOptions {
  SpeciationOptions()
      : system_size(0),
        nonlinear_solver_name(strings::kNewtonRaphson),
        nonlinear_solver_options() {}
  int system_size;
  std::string nonlinear_solver_name;
  NonlinearSolverOptions nonlinear_solver_options;
};

struct OperatorSplittingOptions {
  OperatorSplittingOptions() 
      : system_size(0),
        ode_solver_name("Weighted Euler"),
        ode_solver_options(),
        nonlinear_solver_name(""),
        nonlinear_solver_options() {}
  int system_size;
  std::string ode_solver_name;
  ODESolverOptions ode_solver_options;
  std::string nonlinear_solver_name;
  NonlinearSolverOptions nonlinear_solver_options;
};

struct GlobalImplicitOptions {
  GlobalImplicitOptions() 
      : name("") {}
  std::string name;
};


/*******************************************************************************
 **
 ** Chemistry Coordinator Options and Parameters
 **
 ******************************************************************************/
struct CoordinatorOptions {
  CoordinatorOptions()
      : activity_model_parameters(),
        system_size(0),
        speciation_options(),
        evaluator_type_name(strings::kOperatorSplittingNR),
        os_options(),
        gi_options() {

  }
  ActivityModelParameters activity_model_parameters;
  int system_size;
  SpeciationOptions speciation_options;
  std::string evaluator_type_name;
  OperatorSplittingOptions os_options;
  GlobalImplicitOptions gi_options;
};


struct PrimaryComponentsNames {
  std::vector<std::string> aqueous_species;  // primaries
  std::vector<std::string> aqueous_equilibrium_reactions;  // secondaries
  std::vector<std::string> minerals;
  std::vector<std::string> ion_exchange_sites;
};

struct PrimaryComponents {
  std::vector<double> free_ion;  // [molality]
  std::vector<double> mineral_volume_frac;  // [volume fractions]
  std::vector<double> ion_exchange_sites;  // [CEC]
  std::vector<double> total;  // [molarity]
  std::vector<double> total_sorbed;
};


struct CellPhysicalParameters {
  CellPhysicalParameters() 
      : porosity(1.0), 
        saturation(1.0), 
        water_density_kg_m3(1000.0),
        volume(1.0),
        temperature(25.0) {
  }
  double porosity;  // [-]
  double saturation;  // [-]
  double water_density_kg_m3;  // [kg/m^3]
  double volume;  // [m^3]
  double temperature;  // [celsius]
};


}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif  // AMANZI_CHEMISTRY_CHEMISTRY_CONTAINERS_HH_
