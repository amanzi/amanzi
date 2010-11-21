/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Chemistry_PK_hpp__
#define __Chemistry_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "State.hpp"
#include "Chemistry_State.hpp"
#include "Beaker.hpp"
#include "ChemistryException.hpp"
#include "Verbosity.hpp"

// Chemistry Process Kernel Interface

class Chemistry_PK {
 public:


  Chemistry_PK (Teuchos::ParameterList &param_list,
                Teuchos::RCP<Chemistry_State> chem_state);

  ~Chemistry_PK ();

  void advance(const double& delta_time,
               Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star);
  void commit_state ( Teuchos::RCP<Chemistry_State> chem_state, const double& delta_time);
  ChemistryException::Status status(void) const { return this->status_; };
  Teuchos::RCP<Epetra_MultiVector> get_total_component_concentration(void) const;


  Verbosity verbosity(void) const { return this->verbosity_; };
  void set_verbosity(const Verbosity verbosity) { this->verbosity_ = verbosity; };

  void set_max_time_step(const double mts) { this->max_time_step_ = mts; };
  double max_time_step(void) const { return this->max_time_step_; };

  int number_aqueous_components(void) const { return this->number_aqueous_components_; };
  void set_number_aqueous_components(const int nac) { this->number_aqueous_components_ = nac; };

  int number_minerals(void) const { return this->number_minerals_; };
  void set_number_minerals(const int nm) { this->number_minerals_ = nm; };

  int number_ion_exchange_sites(void) const { return this->number_ion_exchange_sites_; };
  void set_number_ion_exchange_sites(const int ies) { this->number_ion_exchange_sites_ = ies; };

  int number_sorption_sites(void) const { return this->number_sorption_sites_; };
  void set_number_sorption_sites(const int scs) { this->number_sorption_sites_ = scs; };

 protected:
  void set_status(ChemistryException::Status status) { this->status_ = status; };
 private:
  ChemistryException::Status status_;
  Verbosity verbosity_;
  double max_time_step_;
  // auxilary state for process kernel
  Teuchos::RCP<Chemistry_State> chemistry_state_;

  // parameter list
  Teuchos::ParameterList parameter_list_;

  Beaker* chem_;
  Beaker::BeakerParameters beaker_parameters_;
  Beaker::BeakerComponents beaker_components_;

  double current_time_;
  double saved_time_;
  int number_aqueous_components_;
  int number_minerals_;
  int number_ion_exchange_sites_;
  int number_sorption_sites_;

  Teuchos::RCP<const Epetra_Vector> current_porosity_;
  Teuchos::RCP<const Epetra_Vector> current_water_saturation_;
  Teuchos::RCP<const Epetra_Vector> current_water_density_;
  Teuchos::RCP<const Epetra_SerialDenseVector> current_volume_;

  Teuchos::RCP<Epetra_MultiVector> aqueous_components_;
  Teuchos::RCP<Epetra_MultiVector> minerals_;
  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites_;
  Teuchos::RCP<Epetra_MultiVector> sorption_sites_;

  Teuchos::RCP<Epetra_Vector> saved_porosity_;
  Teuchos::RCP<Epetra_Vector> saved_water_saturation_;
  Teuchos::RCP<Epetra_Vector> saved_water_density_;

  void XMLParameters(void);
  void LocalPhysicalState(void);
  void LocalInitialConditions(void);
  void ExtractInitialCondition(const std::string& type,
                               const std::string& keyword,
                               const int number_to_find,
                               const int block,
                               const Teuchos::ParameterList& mesh_block_list,
                               const int mesh_block_ID,
                               Teuchos::RCP<Epetra_MultiVector> data);
  void set_const_values_for_block(const std::vector<double>& values,
                                  const int num_values,
                                  Teuchos::RCP<Epetra_MultiVector>& vec,
                                  const int mesh_block_id);
  void set_cell_value_in_mesh_block(const double value,
                                    Epetra_Vector& vec,
                                    const int mesh_block_id);
  void SizeBeakerComponents(void);
  void CopyCellToBeakerComponents(int cell_id,
                                  Teuchos::RCP<const Epetra_MultiVector> aqueous_components);
  void CopyBeakerComponentsToCell(int cell_id);
  void CopyStateToBeakerParameters(int cell_id);

};

#endif
