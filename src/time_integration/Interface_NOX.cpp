/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Daniil Svyatskiy (version 3) (dasvyat@lanl.gov)
*/

// #include "dbc.hh"

// #include "Flow_PK.hpp"
#include "Interface_NOX.hpp"
#include "Timer.hh"


namespace Amanzi {
namespace AmanziFlow {


bool Interface_NOX::computeF(const Epetra_Vector& x,
                                   Epetra_Vector& f,
                                   FillType flag) {
//   for (int i=0; i<x.MyLength(); i++) f[i] = x[i]*x[i] - 9;

    Epetra_Vector u_tmp(x);

//     u_tmp = x;
    u_tmp.Update(-1.0/deltaT, u0, 1.0/deltaT);

    // evaluate nonlinear functional
    Timer t1;
    t1.start();
    FPK->fun(time, x, u_tmp, f, deltaT);
    t1.stop();
//     std::cout << "evaluate nonlinear functional: " << t1 << std::endl;
    fun_eval++;
    fun_eval_time += t1.getTime();

//     printTime();

    return true;
}


/* ******************************************************************
* NOT USED
****************************************************************** */
bool Interface_NOX::computePreconditioner(
  const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* params) {



  lag_count_++;
  lag_count_ %= lag_prec_;
  int errc;

//   std::cout << lag_count_ << " " << lag_prec_<<"\n";
  if (lag_count_==0) {
    std::cout << "computePreconditioner\n";
    FPK->update_precon(time, x, deltaT, errc);
  }


  return true;
}

void Interface_NOX::printTime() {

    std::cout << "Function evalution was called " << fun_eval << "times" << std::endl;
    std::cout << "Function evalution took " << fun_eval_time << "seconds" << std::endl;

}



}  // namespace AmanziFlow
}  // namespace Amanzi
