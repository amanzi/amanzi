#ifndef _BDF2_FNBASE_HPP_
#define _BDF2_FNBASE_HPP_

namespace BDF2 {

  class fnBase {
    
  public:
        
    virtual void fun(double t, Epetra_Vector& u, Epetra_Vector& udot, Epetra_Vector& f) = 0;
    virtual void precon(Epetra_Vector& u, Epetra_Vector& Pu) = 0;
    virtual double enorm(Epetra_Vector& u, Epetra_Vector& du) = 0;
    virtual void update_precon(double t, Epetra_Vector& up, double h, int& errc) = 0;
   
  };

}

#endif  // _BDF2_FNBASE_HPP_
