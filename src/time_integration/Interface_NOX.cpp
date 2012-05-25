/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

// #include "dbc.hh"

// #include "Flow_PK.hpp"
#include "Interface_NOX.hpp"
#include "Timer.hh"

namespace Amanzi {
namespace AmanziFlow {


bool Interface_NOX::computeF(const Epetra_Vector& x, Epetra_Vector& f, FillType flag)
{
//   for (int i=0; i<x.MyLength(); i++) f[i] = x[i]*x[i] - 9;
  
    
  
    Epetra_Vector u_tmp(x);
     
//     u_tmp = x;
    u_tmp.Update(-1.0/deltaT, u0, 1.0/deltaT);

    // evaluate nonlinear functional
    Timer t1;
    t1.start();
    FPK_->fun(time, x, u_tmp, f, deltaT);
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
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* params)
{
 
  lag_count_++;
  lag_count_ %= lag_prec_;
  
  FPK_->compute_precon(time, deltaT, x, M, params);

  return true;
}

void Interface_NOX::printTime(){
	
	std::cout<<"Function evalution was called "<<fun_eval<<"times"<<std::endl;
	std::cout<<"Function evalution took "<<fun_eval_time<<"seconds"<<std::endl;
	
}

}  // namespace AmanziFlow
}  // namespace Amanzi


 