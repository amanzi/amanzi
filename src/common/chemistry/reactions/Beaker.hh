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

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "VerboseObject.hh"

#include "ActivityModel.hh"
#include "AqueousEquilibriumComplex.hh"
#include "BeakerState.hh"
#include "BeakerParameters.hh"
#include "ChemistryUtilities.hh"
#include "GeneralRxn.hh"
#include "IonExchangeRxn.hh"
#include "KineticRate.hh"
#include "LUSolver.hh"
#include "MatrixBlock.hh"
#include "Mineral.hh"
#include "RadioactiveDecay.hh"
#include "Species.hh"
#include "SorptionIsothermRxn.hh"
#include "SurfaceComplexationRxn.hh"

namespace Amanzi {
namespace AmanziChemistry {

class KineticRate;

class Beaker {
 public:
  Beaker(const Teuchos::Ptr<VerboseObject> vo);
  virtual ~Beaker();

  struct SolverStatus {
    int num_rhs_evaluations;
    int num_jacobian_evaluations;
    int num_newton_iterations;
    bool converged;
  };

  // inheriting classes setup the species, etc
  virtual void Initialize(const BeakerState& state,
                          const BeakerParameters& parameters);

  // we only copy data allocate by Amanzi state
  void CopyStateToBeaker(const BeakerState& state);
  // we copy all and allocate memory as needed
  void CopyBeakerToState(BeakerState* state);
  void CopyState(const BeakerState& from, BeakerState* to) { *to = from; }

  int GetPrimaryIndex(const std::string& name) const;

  // speciate for free-ion concentrations
  int Speciate(BeakerState* state);

  // solve a chemistry step
  int ReactionStep(BeakerState* state,
                   const BeakerParameters& parameters, double dt);

  // enforce constraint
  int EnforceConstraint(BeakerState* state, const BeakerParameters& parameters,
                        const std::vector<std::string>& names,
                        const std::vector<double>& values);

  // i/o
  void Display() const;
  void DisplayComponents(const BeakerState& state) const;
  void DisplayTotalColumns(const double time, 
                           const BeakerState& total,
                           const bool display_free) const;
  void DisplayResults() const;

  // access
  SolverStatus status() const { return status_; }

  const std::vector<Mineral>& minerals() const { return minerals_; }
  const std::vector<Species>& primary_species() const { return primary_species_; }
  const std::vector<AqueousEquilibriumComplex>& secondary_species() { return aqComplexRxns_; }
  const std::vector<IonExchangeRxn>& ion_exchange_rxns() const { return ion_exchange_rxns_; }

  const std::vector<double>& total() const { return total_; }
  const std::vector<double>& total_sorbed() const { return total_sorbed_; }
  std::vector<SorptionIsothermRxn> sorption_isotherm_rxns() const { return sorption_isotherm_rxns_; }

 protected:
  // resizes matrix and vectors for nonlinear system
  void ResizeInternalMemory();

  void SetupActivityModel(std::string model, std::string pitzer_database, std::string jfunction_pitzer);
  void VerifyState(const BeakerState& state) const;

  void AddIonExchangeRxn(const IonExchangeRxn& ionx_rxn);
  void AddIonExchangeComplex(int irxn, const IonExchangeComplex& ionx_complex);
  void AddAqueousEquilibriumComplex(const AqueousEquilibriumComplex& c);
  void AddMineral(const Mineral& m);
  void AddMineralKineticRate(KineticRate* rate);
  void AddGeneralRxn(const GeneralRxn& r);
  void AddRadioactiveDecayRxn(const RadioactiveDecay& r);
  void AddSurfaceComplexationRxn(const SurfaceComplexationRxn& r);
  void AddSorptionIsothermRxn(const SorptionIsothermRxn& r);

  void set_dt(double dt) { dt_ = dt; }

 private:
  void CheckChargeBalance_(const std::vector<double>& aqueous_totals) const;

  void UpdateActivityCoefficients();
  void UpdateKineticMinerals();
  void InitializeMolalities_(double initial_molality);
  void InitializeMolalities_(const std::vector<double>& initial_molalities);

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

  void ResetStatus();

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

  Species water_;
  std::vector<Species> primary_species_; 
  std::vector<Mineral> minerals_;

 private:
  double tolerance_;
  unsigned int max_iterations_;
  int ncomp_;  // numabe of primaryspecies

  std::vector<double> total_;  // aqueous tcc for primaries [mol/L]
  MatrixBlock dtotal_;  // derivaties wrt free-ion [kg water/sec]

  // Sorbed phase total component concentrations for basis species
  std::vector<double> total_sorbed_;  // [mol/m^3 bulk]
  MatrixBlock dtotal_sorbed_;  // derivaties wrt free-ion [kg water/sec]

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
  double sorbed_accumulation_coef_;  // [m^3 bulk/sec]
  double por_sat_den_vol_;

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
  LUSolver lu_solver_;

  bool use_log_formulation_;

  std::vector<double> sorption_isotherm_params_;
};


// non-member functions
inline
void Display(const BeakerState& state,
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
