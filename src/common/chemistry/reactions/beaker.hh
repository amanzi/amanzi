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

  struct BeakerComponents {
    // TODO(bandre): rename to BeakerState and move all "state"
    // variables (porosity, density, volume, mineral ssa, isotherms,
    // etc) into a single struct.
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

    void Display(const std::string& message, const Teuchos::RCP<VerboseObject>& vo) const {
      vo->Write(Teuchos::VERB_HIGH, message);
      utilities::PrintVector("Totals", total, vo, 16, true);
      utilities::PrintVector("Total sorbed", total_sorbed, vo, 16, true);
      utilities::PrintVector("Free Ion", free_ion, vo, 16, true);
      utilities::PrintVector("Primary activity coeff", primary_activity_coeff, vo, 16, true);
      utilities::PrintVector("Secondary activity coeff", secondary_activity_coeff, vo, 16, true);
      utilities::PrintVector("Mineral VF", mineral_volume_fraction, vo, 16, true);
      utilities::PrintVector("Mineral SSA", mineral_specific_surface_area, vo, 16, true);
      utilities::PrintVector("Ion Exchange Sites", ion_exchange_sites, vo, 16, true);
      utilities::PrintVector("Ion Exchange ref cation conc", ion_exchange_ref_cation_conc, vo, 16, true);
      utilities::PrintVector("Surface Site Density", surface_site_density, vo, 16, true);
      utilities::PrintVector("Surface Complex free site conc", surface_complex_free_site_conc, vo, 16, true);
      utilities::PrintVector("Kd", isotherm_kd, vo, 16, true);
      utilities::PrintVector("Freundlich n", isotherm_freundlich_n, vo, 16, true);
      utilities::PrintVector("Langmuir b", isotherm_langmuir_b, vo, 16, true);
    }
    void DumpCfg(const std::vector<std::string>& names, const Teuchos::RCP<VerboseObject>& vo) const {
      std::stringstream message;
      message << std::scientific << std::setprecision(12);
      message << "\n[total]\n";
      for (int i = 0; i < total.size(); ++i) {
        message << names.at(i) << " = " << total.at(i) << "\n";
      }
      message << "\n[total_sorbed]\n";
      for (int i = 0; i < total_sorbed.size(); ++i) {
        message << names.at(i) << " = " << total_sorbed.at(i) << "\n";
      }
      message << "\n[free_ion]\n";
      for (int i = 0; i < free_ion.size(); ++i) {
        message << names.at(i) << " = " << free_ion.at(i) << "\n";
      }
      vo->Write(Teuchos::VERB_HIGH, message.str());
    }
  };

  struct BeakerParameters {
    // input files
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
  virtual void Setup(const Beaker::BeakerComponents& components,
                     const Beaker::BeakerParameters& parameters);
  void CopyBeakerToComponents(Beaker::BeakerComponents* components);

  BeakerParameters GetDefaultParameters(void) const;
  BeakerParameters GetCurrentParameters(void) const;
  void SetParameters(const BeakerParameters& parameters);
  void CopyComponents(const Beaker::BeakerComponents& from,
                      Beaker::BeakerComponents* to);

  void GetPrimaryNames(std::vector<std::string>* names) const;
  int GetPrimaryIndex(const std::string& name) const;

  void Display(void) const;
  void DisplayComponents(const BeakerComponents& components) const;
  void DisplayTotalColumnHeaders(const bool display_free) const;
  void DisplayTotalColumns(const double time, 
                           const BeakerComponents& total,
                           const bool display_free) const;
  void DisplayResults(void) const;


  bool HaveKinetics(void) const;

  // speciate for free-ion concentrations
  int Speciate(BeakerComponents* components,
               const BeakerParameters& parameters);

  // solve a chemistry step
  int ReactionStep(BeakerComponents* components, const BeakerParameters& parameters,
                   double dt);

  void print_results(void) const;
  void print_results(double time) const;

  int ncomp(void) const { return this->ncomp_; }
  double tolerance(void) const { return this->tolerance_; }
  unsigned int max_iterations(void) const { return this->max_iterations_; }

  double porosity(void) const { return this->porosity_; }
  double saturation(void) const { return this->saturation_; }

  double water_density_kg_m3(void) const { return this->water_density_kg_m3_; }
  double water_density_kg_L(void) const { return this->water_density_kg_L_; }
  double volume(void) const { return this->volume_; }
  double dt(void) const { return this->dt_; }

  double aqueous_accumulation_coef(void) const { return this->aqueous_accumulation_coef_; }
  double sorbed_accumulation_coef(void) const { return this->sorbed_accumulation_coef_; }
  double por_sat_den_vol(void) const { return this->por_sat_den_vol_; }

  SolverStatus status(void) const { return this->status_; }

  const std::vector<Mineral>& minerals(void) const { return this->minerals_; }
  const std::vector<Species>& primary_species(void) const { return this->primary_species_; }
  const std::vector<IonExchangeRxn>& ion_exchange_rxns(void) const { return this->ion_exchange_rxns_; }

  const std::vector<double>& total(void) const { return this->total_; }
  const std::vector<double>& total_sorbed(void) const { return this->total_sorbed_; }

  const std::vector<double>& fixed_accumulation(void) const { return this->fixed_accumulation_; }
  const std::vector<double>& prev_molal(void) const { return this->prev_molal_; }

 protected:
  // resizes matrix and vectors for nonlinear system
  void ResizeInternalMemory(const int size);

  void SetupActivityModel(std::string model, std::string pitzer_database, std::string jfunction_pitzer);
  void VerifyComponentSizes(const Beaker::BeakerComponents& components) const;
  void CopyComponentsToBeaker(const Beaker::BeakerComponents& components);

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

  void ncomp(int i) { this->ncomp_ = i; }
  void tolerance(double value) { this->tolerance_ = value; }
  void max_iterations(unsigned int value) { this->max_iterations_ = value; }

  void porosity(double d) { this->porosity_ = d; }
  void saturation(double d) { this->saturation_ = d; }

  // updates both water density variables
  void water_density_kg_m3(double d) {
    this->water_density_kg_m3_ = d;
    this->water_density_kg_L_ = d / 1000.;
  }
  void water_density_kg_L(double d) {
    this->water_density_kg_m3_ = d * 1000.;
    this->water_density_kg_L_ = d;
  }

  void volume(double d) { this->volume_ = d; }
  void dt(double d) { this->dt_ = d; }
  void aqueous_accumulation_coef(double d) { this->aqueous_accumulation_coef_ = d; }
  void sorbed_accumulation_coef(double d) { this->sorbed_accumulation_coef_ = d; }
  void por_sat_den_vol(double d) { this->por_sat_den_vol_ = d; }

  // calculates the coefficient in aqueous portion of accumulation term
  void update_accumulation_coefficients(void);
  // calculates product of porosity,saturation,water_density[kg/m^3],volume
  void update_por_sat_den_vol(void);

  void set_use_log_formulation(const bool value) { use_log_formulation_ = value; }
  bool use_log_formulation(void) const { return use_log_formulation_; }

 protected:
  Teuchos::Ptr<VerboseObject> vo_;

 private:
  void CheckChargeBalance(const std::vector<double>& aqueous_totals) const;

  // update discretization and flow parameters
  // water_density [kg/m^3]
  // volume [m^3]
  void UpdateParameters(const BeakerParameters& parameters, double dt);

  void UpdateActivityCoefficients(void);
  void UpdateKineticMinerals(void);
  void InitializeMolalities(double initial_molality);
  void InitializeMolalities(const std::vector<double>& initial_molalities);

  // equilibrium chemistry
  // update activities, equilibrium complex concentrations, etc.
  void UpdateEquilibriumChemistry(void);
  void CalculateTotal(void);

  // calculate block of Jacobian corresponding to derivatives of total with
  // respect to free-ion
  void CalculateDTotal(void);

  // kinetic chemistry
  void UpdateKineticChemistry(void);
  void AddKineticChemistryToResidual(void);
  void AddKineticChemistryToJacobian(void);

  // accumulation terms
  void AddAccumulation(const std::vector<double>& total,
                       const std::vector<double>& total_sorbed,
                       std::vector<double> *residual);
  void AddAccumulationDerivative(MatrixBlock* J, MatrixBlock* dtotal, MatrixBlock* dtotal_sorbed);
  void CalculateFixedAccumulation(const std::vector<double>& total,
                                  const std::vector<double>& total_sorbed,
                                  std::vector<double> *fixed_accumulation);
  // residual and Jacobian
  void CalculateResidual(void);
  void CalculateJacobian(void);

  // utilities for updating solution, convergence checks
  void UpdateMolalitiesWithTruncation(const double max_change);
  void CalculateMaxRelChangeInMolality(double* max_rel_change, int* max_rel_index);
  void ValidateSolution(void);

  // solvers
  void ScaleRHSAndJacobian(void);

  // output
  void DisplayParameters(void) const;
  void DisplayPrimary(void) const;
  void DisplayAqueousEquilibriumComplexes(void) const;
  void DisplayGeneralKinetics(void) const;
  void DisplayRadioactiveDecayRxns(void) const;
  void DisplayMinerals(void) const;
  void DisplayMineralKinetics(void) const;
  void DisplayIonExchangeRxns(void) const;
  void DisplayIonExchangeSites(void) const;
  void DisplayIonExchangeComplexes(void) const;
  void DisplaySurfaceSites(void) const;
  void DisplaySurfaceComplexes(void) const;
  void DisplaySorptionIsotherms(void) const;

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
  double porosity_;            // [m^3 pore / m^3 bulk]
  double saturation_;          // [m^3 water / m^3 pore]
  double water_density_kg_m3_;  // [kg water / m^3 water]
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
  std::vector<Species> primary_species_;  // list of primary species
  std::vector<Mineral> minerals_;  // list of mineral species

  ActivityModel* activity_model_;

  std::vector<AqueousEquilibriumComplex> aqComplexRxns_;  // list of aqueous equilibrium complexation reactions
  std::vector<GeneralRxn> generalKineticRxns_;  // list of general kinetic reactions
  std::vector<RadioactiveDecay> radioactive_decay_rxns_;  // list of radioactive decay rxns
  std::vector<KineticRate*> mineral_rates_;
  //  vector<GasExchange*> gasRxns_;
  std::vector<IonExchangeRxn> ion_exchange_rxns_;
  std::vector<SurfaceComplexationRxn> surfaceComplexationRxns_;
  std::vector<SorptionIsothermRxn> sorption_isotherm_rxns_;

  // solver data structures
  std::vector<double> fixed_accumulation_;  // fixed (time t) portion of accumulation term
  std::vector<double> residual_;       // entire residual [mol/sec]
  std::vector<double> prev_molal_;     // previous molality of primary species

  std::vector<double> rhs_;            // right-hand-side of system
  MatrixBlock jacobian_;              // Jacobian [kg water/sec]

  SolverStatus status_;
  void ResetStatus(void);
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

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
