/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/
 
#ifndef AMANZI_CHEMISTRY_PK_HH_
#define AMANZI_CHEMISTRY_PK_HH_

#include "Teuchos_RCP.hpp"
#include "Chemistry_State.hh"

// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;

namespace Amanzi {
namespace AmanziChemistry {

// Abstract base class for Chemistry PKs.
class Chemistry_PK {
 public:
  Chemistry_PK() {};
  virtual ~Chemistry_PK() {};

  virtual void InitializeChemistry() = 0;

  virtual void Advance(const double& delta_time,
                       Teuchos::RCP<Epetra_MultiVector> total_component_concentration_star) = 0;
  virtual void CommitState(Teuchos::RCP<Chemistry_State> chem_state, const double& delta_time) = 0;

  // Returns the (maximum) time step allowed for this chemistry PK.
  virtual double time_step() const = 0;

  // Ben: the following two routines provide the interface for
  // output of auxillary cellwise data from chemistry
  virtual Teuchos::RCP<Epetra_MultiVector> get_extra_chemistry_output_data() = 0;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
