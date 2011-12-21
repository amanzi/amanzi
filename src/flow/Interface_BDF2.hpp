/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __INTERFACE_BDF2_HPP__
#define __INTERFACE_BDF2_HPP__

#include "Epetra_Vector.h"

#include "BDF2_fnBase.hpp"

namespace Amanzi {
namespace AmanziFlow {

class Interface_BDF2 : public BDF2::fnBase {
 public:
  Interface_BDF2() {}; 
  ~Interface_BDF2() {}; 

  // required methods from the base class
  void initialize(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params);

  void fun(const double t, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& f);
  void precon(const Epetra_Vector& u, Epetra_Vector& Pu);
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du);
  void update_precon(const double t, const Epetra_Vector& up, const double h, int& errc);

  bool is_admissible(const Epetra_Vector& up) { return true; }
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
