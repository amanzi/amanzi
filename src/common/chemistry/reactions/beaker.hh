/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Driver class for evaluating geochemical related processes at a
  single computational node
*/

#ifndef AMANZI_CHEMISTRY_BEAKER_HH_
#define AMANZI_CHEMISTRY_BEAKER_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "VerboseObject.hh"

#include "activity_model.hh"
#include "aqueous_equilibrium_complex.hh"
#include "general_rxn.hh"
#include "radioactive_decay.hh"
#include "ion_exchange_rxn.hh"
#include "kinetic_rate.hh"
#include "mineral.hh"
#include "species.hh"
#include "sorption_isotherm_rxn.hh"
#include "surface_complexation_rxn.hh"
#include "lu_solver.hh"
#include "matrix_block.hh"
#include "chemistry_utilities.hh"

namespace Amanzi {
namespace AmanziChemistry {

class KineticRate;

class Beaker {
 public:
  Beaker(const Teuchos::Ptr<VerboseObject> vo);
  virtual ~Beaker();

  struct BeakerState {
    // TODO(bandre): move all "state" variables (porosity, density, 
    // volume, mineral ssa, isotherms, etc) into a single struct.
    std::vector<double> total;  // molarity
    std::vector<double> total_sorbed;
    std::vector<double> free_ion;  // molality
    std::vector<double> primary_activity_coeff;
    std::vector<double> secondary_activity_coeff;
    std::vector<double> mineral_volume_fraction;  // volume fractions
    std::vector<double> mineral_specific_surface_area;  // [m^2 mineral/ m^3 bulk]
    std::vector<double> ion_exchange_sites;  // CEC
    std::vector<double> ion_exchange_ref_cation_conc;  // [?]
    std::vector<double> surface_site_density;
    std::vector<double> surface_complex_free_site_conc;  // [moles sites / m^3 bulk]
    std::vector<double> isotherm_kd;
    std::vector<double> isotherm_freundlich_n;
    std::vector<double> isotherm_langmuir_b;
  };

  struct BeakerParameters {
    std::string thermo_database_file;

    // solver parameters
    double tolerance;
    unsigned int max_iterations;

    // models
    std::string activity_model_name;
    // Name of the Pitzer virial coefficients database
    std::string pitzer_database;
    // Name of the approach for J's functions for the Pitzer model
    std::string jfunction_pitzer;

    // physical parameters
    double porosity;  // [-]
    double saturation;  // [-]
    double water_density;  // [kg/m^3]
    double volume;  // [m^3]
  };

  struct SolverStatus {
    unsigned int num_rhs_evaluations;
    unsigned int num_jacobian_evaluations;
    unsigned int num_newton_iterations;
    bool converged;
  };

  //
  // public routines that may need to be accessed by the driver
  //

  // inheriting classes setup the species, etc
  virtual void Setup(const BeakerState& state,
                     const BeakerParameters& parameters);
  void CopyBeakerToState(BeakerState* state);

  BeakerParameters GetDefaultParameters() const;
  // take parameters object and map the data into chemistry object
  void SetParameters(const BeakerParameters& parameters);
  void CopyState(const BeakerState& from, BeakerState* to) { *to = from; }

  void GetPrimaryNames(std::vector<std::string>* names) const;
  int GetPrimaryIndex(const std::string& name) const;

  bool HaveKinetics() const;

  // speciate for free-ion concentrations
  int Speciate(BeakerState* state, const BeakerParameters& parameters);

  // solve a chemistry step
  int ReactionStep(BeakerState* state, const BeakerParameters& parameters, double dt);

  // enforce constraint
  int EnforceConstraint(BeakerState* state, const BeakerParameters& parameters,
                        const std::vector<std::string>& names,
                        const std::vector<double>& values);

  // i/o
  void Display() const;
  void DisplayComponents(const BeakerState& state) const;
  void DisplayTotalColumnHeaders(const bool display_free) const;
  void DisplayTotalColumns(const double time, 
                           const BeakerState& total,
                           const bool display_free) const;
  void DisplayResults() const;

  void print_results() const;
  void print_results(double time) const;

  // access
  // int ncomp() const { return ncomp_; }

  double water_density_kg_m3() const { return water_density_kg_m3_; }
  double water_density_kg_L() const { return water_density_kg_L_; }

  double sorbed_accumulation_coef() const { return sorbed_accumulation_coef_; }
  double por_sat_den_vol() const { return por_sat_den_vol_; }

  SolverStatus status() const { return status_; }

  const std::vector<Mineral>& minerals() const { return minerals_; }
  const std::vector<Species>& primary_species() const { return primary_species_; }
  const std::vector<IonExchangeRxn>& ion_exchange_rxns() const { return ion_exchange_rxns_; }

  const std::vector<double>& total() const { return total_; }
  const std::vector<double>& total_sorbed() const { return total_sorbed_; }

 protected:
  // resizes matrix and vectors for nonlinear system
  void ResizeInternalMemory(const int size);

  void SetupActivityModel(std::string model, std::string pitzer_database, std::string jfunction_pitzer);
  void VerifyState(const BeakerState& state) const;
  void CopyStateToBeaker(const BeakerState& state);

  void AddPrimarySpecies(const Species& s);
  void AddIonExchangeRxn(const IonExchangeRxn& ionx_rxn);
  void AddIonExchangeComplex(const int irxn, const IonExchangeComplex& ionx_complex);
  void AddAqueousEquilibriumComplex(const AqueousEquilibriumComplex& c);
  void AddMineral(const Mineral& m);
  void AddMineralKineticRate(KineticRate* rate);
  void AddGeneralRxn(const GeneralRxn& r);
  void AddRadioactiveDecayRxn(const RadioactiveDecay& r);
  void AddSurfaceComplexationRxn(const SurfaceComplexationRxn& r);
  void AddSorptionIsothermRxn(const SorptionIsothermRxn& r);

  void set_ncomp(int i) { ncomp_ = i; }
  void set_tolerance(double value) { tolerance_ = value; }
  void set_max_iterations(int value) { max_iterations_ = value; }

  void set_porosity(double value) { porosity_ = value; }
  void set_saturation(double value) { saturation_ = value; }

  // updates both water density variables
  void water_density_kg_m3(double d) {
    water_density_kg_m3_ = d;
    water_density_kg_L_ = d / 1000.;
  }
  void water_density_kg_L(double d) {
    water_density_kg_m3_ = d * 1000.;
    water_density_kg_L_ = d;
  }

  void set_volume(double vol) { volume_ = vol; }
  void set_dt(double dt) { dt_ = dt; }
  void set_aqueous_accumulation_coef(double coef) { aqueous_accumulation_coef_ = coef; }
  void set_sorbed_accumulation_coef(double coef) { sorbed_accumulation_coef_ = coef; }
  void por_sat_den_vol(double coef) { por_sat_den_vol_ = coef; }

  void set_use_log_formulation(bool value) { use_log_formulation_ = value; }

 private:
  void CheckChargeBalance(const std::vector<double>& aqueous_totals) const;

  // update discretization and flow parameters
  // water_density [kg/m^3]
  // volume [m^3]
  void UpdateParameters(const BeakerParameters& parameters, double dt);

  void UpdateActivityCoefficients();
  void UpdateKineticMinerals();
  void InitializeMolalities(double initial_molality);
  void InitializeMolalities(const std::vector<double>& initial_molalities);

  // equilibrium chemistry
  // update activities, equilibrium complex concentrations, etc.
  void UpdateEquilibriumChemistry();
  void CalculateTotal();

  // calculate block of Jacobian corresponding to derivatives of total with
  // respect to free-ion
  void CalculateDTotal();

  // kinetic chemistry
  void UpdateKineticChemistry();
  void AddKineticChemistryToResidual();
  void AddKineticChemistryToJacobian();

  // accumulation terms
  void AddAccumulation(const std::vector<double>& total,
                       const std::vector<double>& total_sorbed,
                       std::vector<double> *residual);
  void AddAccumulationDerivative(MatrixBlock* J, MatrixBlock* dtotal, MatrixBlock* dtotal_sorbed);
  void CalculateFixedAccumulation(const std::vector<double>& total,
                                  const std::vector<double>& total_sorbed,
                                  std::vector<double> *fixed_accumulation);
  // residual and Jacobian
  void CalculateResidual();
  void CalculateJacobian();

  // utilities for updating solution, convergence checks
  void UpdateMolalitiesWithTruncation(const double max_change);
  void CalculateMaxRelChangeInMolality(double* max_rel_change, int* max_rel_index);
  void ValidateSolution();

  // solvers
  void ScaleRHSAndJacobian();

  // output
  void DisplayParameters() const;
  void DisplayPrimary() const;
  void DisplayAqueousEquilibriumComplexes() const;
  void DisplayGeneralKinetics() const;
  void DisplayRadioactiveDecayRxns() const;
  void DisplayMinerals() const;
  void DisplayMineralKinetics() const;
  void DisplayIonExchangeRxns() const;
  void DisplayIonExchangeSites() const;
  void DisplayIonExchangeComplexes() const;
  void DisplaySurfaceSites() const;
  void DisplaySurfaceComplexes() const;
  void DisplaySorptionIsotherms() const;

 protected:
  Teuchos::Ptr<VerboseObject> vo_;

 private:
  double tolerance_;
  unsigned int max_iterations_;
  int ncomp_;                   // # basis species

  // Aqueous phase total component concentrations for basis species
  std::vector<double> total_;  // [mol/L]
  // Matrix block containing derivative of total concentration w/respec to
  // free-ion
  MatrixBlock dtotal_;  // [kg water/sec]

  // Sorbed phase total component concentrations for basis species
  std::vector<double> total_sorbed_;  // [mol/m^3 bulk]
  // Matrix block containing derivative of total sorbed concentration
  // w/respec to free-ion
  MatrixBlock dtotal_sorbed_;  // [kg water/sec]

  // common parameters among reactions
  double porosity_;  // [m^3 pore / m^3 bulk]
  double saturation_;  // [m^3 water / m^3 pore]
  double water_density_kg_m3_;  // [kg water / m^3 water]
  double water_density_kg_L_;  // [kg water / L water]
  double volume_;  // cell volume [m^3 bulk]
  double dt_; 
  // aqueous_accumulation_coef_ = porosity * saturation * volume * 1000 / dt [L water/sec]
  // units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = (L water/sec)
  double aqueous_accumulation_coef_;
  // sorbed_accumulation_coef_ = volume/dt [m^3 bulk/sec]
  double sorbed_accumulation_coef_;
  // por_sat_den_vol_ = porosity * saturation * water_density * volume [kg water]
  double por_sat_den_vol_;

  Species water_;
  std::vector<Species> primary_species_;  // list of primary species
  std::vector<Mineral> minerals_;  // list of mineral species

  ActivityModel* activity_model_;

  std::vector<AqueousEquilibriumComplex> aqComplexRxns_;  // list of aqueous equilibrium complexation reactions
  std::vector<GeneralRxn> generalKineticRxns_;  // list of general kinetic reactions
  std::vector<RadioactiveDecay> radioactive_decay_rxns_;  // list of radioactive decay rxns
  std::vector<KineticRate*> mineral_rates_;
  // vector<GasExchange*> gasRxns_;
  std::vector<IonExchangeRxn> ion_exchange_rxns_;
  std::vector<SurfaceComplexationRxn> surfaceComplexationRxns_;
  std::vector<SorptionIsothermRxn> sorption_isotherm_rxns_;

  // solver data structures
  std::vector<double> fixed_accumulation_;  // fixed (time t) portion of accumulation term
  std::vector<double> residual_;  // entire residual [mol/sec]
  std::vector<double> prev_molal_;  // previous molality of primary species

  std::vector<double> rhs_;
  MatrixBlock jacobian_;  // Jacobian [kg water/sec]

  SolverStatus status_;
  void ResetStatus();
  static const double tolerance_default_;
  static const unsigned int max_iterations_default_;
  static const double porosity_default_;
  static const double saturation_default_;
  static const double water_density_kg_m3_default_;
  static const double volume_default_;

  LUSolver lu_solver_;

  bool use_log_formulation_;

  std::vector<double> sorption_isotherm_params_;
};


// non-member functions
inline
void Display(const Beaker::BeakerState& state,
             const std::string& message,
             const Teuchos::RCP<VerboseObject>& vo) {
  vo->Write(Teuchos::VERB_HIGH, message);
  utilities::PrintVector("Totals", state.total, vo, 16, true);
  utilities::PrintVector("Total sorbed", state.total_sorbed, vo, 16, true);
  utilities::PrintVector("Free Ion", state.free_ion, vo, 16, true);
  utilities::PrintVector("Primary activity coeff", state.primary_activity_coeff, vo, 16, true);
  utilities::PrintVector("Secondary activity coeff", state.secondary_activity_coeff, vo, 16, true);
  utilities::PrintVector("Mineral VF", state.mineral_volume_fraction, vo, 16, true);
  utilities::PrintVector("Mineral SSA", state.mineral_specific_surface_area, vo, 16, true);
  utilities::PrintVector("Ion Exchange Sites", state.ion_exchange_sites, vo, 16, true);
  utilities::PrintVector("Ion Exchange ref cation conc", state.ion_exchange_ref_cation_conc, vo, 16, true);
  utilities::PrintVector("Surface Site Density", state.surface_site_density, vo, 16, true);
  utilities::PrintVector("Surface Complex free site conc", state.surface_complex_free_site_conc, vo, 16, true);
  utilities::PrintVector("Kd", state.isotherm_kd, vo, 16, true);
  utilities::PrintVector("Freundlich n", state.isotherm_freundlich_n, vo, 16, true);
  utilities::PrintVector("Langmuir b", state.isotherm_langmuir_b, vo, 16, true);
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
