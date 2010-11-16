/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Chemistry_PK_hpp__
#define __Chemistry_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"
#include "Beaker.hpp"
#include "Verbosity.hpp"

// Chemistry Process Kernel Interface

class Chemistry_PK {
 public:

  enum ChemistryStatus {
    kChemistryOK, kChemistryError
  };

  Chemistry_PK (Teuchos::ParameterList &param_list, 
                Teuchos::RCP<Chemistry_State> chem_state);

  ~Chemistry_PK ();

  void advance(const double& delta_time,
               Teuchos::RCP<Epetra_MultiVector> total_component_concentration_star);
  void commit_state ( Teuchos::RCP<Chemistry_State> chem_state, const double& delta_time);
  ChemistryStatus status(void) const { return this->status_; };

  Verbosity verbosity(void) const { return this->verbosity_; };
  void set_verbosity(const Verbosity verbosity) { this->verbosity_ = verbosity; };

  unsigned int num_aqueous_components(void) const { return this->num_aqueous_components_; };
  void num_aqueous_components(const unsigned int nac) { this->num_aqueous_components_ = nac; };

 protected:
  void set_status(const ChemistryStatus& status) { this->status_ = status; };
 private:
  ChemistryStatus status_;
  Verbosity verbosity_;
  // auxilary state for process kernel
  Teuchos::RCP<Chemistry_State> chemistry_state_;

  // parameter list
  Teuchos::ParameterList parameter_list_;

  Beaker* chem_;
  Beaker::BeakerParameters beaker_parameters_;
  Beaker::BeakerComponents beaker_components_;

  Teuchos::RCP<const Epetra_Vector> current_porosity_;
  Teuchos::RCP<const Epetra_Vector> current_water_saturation_;
  Teuchos::RCP<const Epetra_Vector> current_water_density_;

  Teuchos::RCP<Epetra_Vector> saved_porosity_;
  Teuchos::RCP<Epetra_Vector> saved_water_saturation_;
  Teuchos::RCP<Epetra_Vector> saved_water_density_;

  double current_time_;
  double saved_time_;
  unsigned int num_aqueous_components_;

  void XMLParameters(void);
  void LocalPhysicalState(void);
  void LocalInitialConditions(void);

};

#endif
