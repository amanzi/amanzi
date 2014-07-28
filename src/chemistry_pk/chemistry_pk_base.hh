/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_PK_BASE_HH_
#define AMANZI_CHEMISTRY_PK_BASE_HH_

#include "Teuchos_RCP.hpp"
#include "Chemistry_State.hh"

// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;

namespace Amanzi {
namespace AmanziChemistry {

// Abstract base class for Chemistry PKs.
class Chemistry_PK_Base {
 public:

  Chemistry_PK_Base() {}
  virtual ~Chemistry_PK_Base() {}

  virtual void InitializeChemistry(void) = 0;

  virtual void Advance(const double& delta_time,
                       Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star) = 0;
  virtual void CommitState(Teuchos::RCP<Chemistry_State> chem_state, const double& delta_time) = 0;
  virtual Teuchos::RCP<Epetra_MultiVector> get_total_component_concentration(void) const = 0;

  // Returns the maximum time step allowed for this chemistry PK.
  virtual double max_time_step(void) const = 0;

  // Ben: the following two routines provide the interface for
  // output of auxillary cellwise data from chemistry
  virtual Teuchos::RCP<Epetra_MultiVector> get_extra_chemistry_output_data() = 0;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif // AMANZI_CHEMISTRY_PK_BASE_HH_
