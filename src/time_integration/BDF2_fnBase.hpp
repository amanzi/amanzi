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
    virtual void fun(double t, Epetra_Vector& u, Epetra_Vector& udot, Epetra_Vector& f) = 0;

    // applies preconditioner to u and returns the result in Pu
    virtual void precon(Epetra_Vector& u, Epetra_Vector& Pu) = 0;

    // computes a norm on u-du and returns the result 
    virtual double enorm(Epetra_Vector& u, Epetra_Vector& du) = 0;

    // updates the preconditioner
    virtual void update_precon(double t, Epetra_Vector& up, double h, int& errc) = 0;
   
    // check the admissibility of a solution
    // override with the actual admissibility check
    bool is_admissible(Epetra_Vector& up) { return true; }

  };

}

#endif  // _BDF2_FNBASE_HPP_
