/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __INTERFACE_NOX_HPP__
#define __INTERFACE_NOX_HPP__

#include "exceptions.hh"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "BDF2_Dae.hpp"


namespace Amanzi {
namespace AmanziFlow {

class Interface_NOX : public NOX::Epetra::Interface::Required,
                      public NOX::Epetra::Interface::Jacobian,
                      public NOX::Epetra::Interface::Preconditioner {
 public:
  Interface_NOX(BDF2::fnBase* FPK, const Epetra_Vector& uprev, double time_, double dt) : 
		FPK_(FPK), u0(uprev), lag_prec_(1), lag_count_(0) {
			time = time_;
			deltaT = dt;
			fun_eval = 0;
			fun_eval_time = 0;
		};
  ~Interface_NOX() {};

  // required interface members
  bool computeF(const Epetra_Vector& x, Epetra_Vector& f, FillType flag);
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& J) { assert(false); }
  bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* params);
  void printTime();

  inline void setPrecLag(int lag_prec) { lag_prec_ = lag_prec;}
  inline void resetPrecLagCounter() { lag_count_ = 0; }
  inline int getPrecLag() const { return lag_prec_; }
  inline int getPrecLagCounter() const { return lag_count_; }

 private:
  BDF2::fnBase* FPK_;
  const Epetra_Vector& u0;	// value at the previous time step

  double deltaT, time;		// time step
  int lag_prec_;  // the preconditioner is lagged this many times before it is recomputed
  int lag_count_; // this counts how many times the preconditioner has been lagged
  int fun_eval;
  double  fun_eval_time;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
