#ifndef _BDF2_FNBASE_HPP_
#define _BDF2_FNBASE_HPP_

namespace BDF2 {

// this is the interface definition for the BDF2 class
// the nonlinear functional, preconditioner, and error 
// functions must be derived from this class to be 
// usable with BDF2::Dae

class fnBase {
 public:
  // computes the non-linear functional f = f(t,u,udot) 
  virtual void fun(const double t, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& f, const double dT=0.0) = 0;

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(const Epetra_Vector& u, Epetra_Vector& Pu) = 0;

  // computes a norm on u-du and returns the result 
  virtual double enorm(const Epetra_Vector& u, const Epetra_Vector& du) = 0;

  // updates the preconditioner
  virtual void update_precon(const double t, const Epetra_Vector& up, const double h, int& errc) = 0;
   
  // check the admissibility of a solution
  // override with the actual admissibility check
  virtual bool is_admissible(const Epetra_Vector& up) { return true; }
};

}  // namespace BDF2

#endif  // _BDF2_FNBASE_HPP_
