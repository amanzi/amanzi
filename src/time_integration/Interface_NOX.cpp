/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Daniil Svyatskiy (version 3) (dasvyat@lanl.gov)
*/

// #include "dbc.hh"

// #include "Flow_PK.hpp"
#include "NOX_Abstract_Vector.H"
#include "NOX_Common.H"
#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_Vector.H"
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
* Update preconditioner for JFNK method
****************************************************************** */
bool Interface_NOX::computePreconditioner(
  const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* params) {
  lag_count_++;
  lag_count_ %= lag_prec_;
  int errc;

  if (lag_count_ == 0) {
    FPK->update_precon(time, x, deltaT, errc);
  }


  return true;
}

void Interface_NOX::printTime() {
    std::cout << "Function evalution was called " << fun_eval << "times" << std::endl;
    std::cout << "Function evalution took " << fun_eval_time << "seconds" << std::endl;
}

NOX::StatusTest::StatusType PK_enorm::
checkStatus(const NOX::Solver::Generic& problem,
            NOX::StatusTest::CheckType checkType) {
        
        if (checkType == NOX::StatusTest::None) {
                status = NOX::StatusTest::Unevaluated;
                normF = -1.0;
                return status;
        }

        // On the first iteration, the old and current solution are the same so
        // we should return the test as unconverged until there is a valid 
        // old solution (i.e. the number of iterations is greater than zero).
        int niters = problem.getNumIterations();
        if (niters == 0) {
                status = NOX::StatusTest::Unconverged;
                normF = -1.0;
                return status;
        } 

        // Check that F exists!
        if (!problem.getSolutionGroup().isF()) {
                status = NOX::StatusTest::Unconverged;
                normF = -1.0;
                return status;
        } 

        const NOX::Abstract::Vector& oldSoln = problem.getPreviousSolutionGroup().getX();
        const NOX::Abstract::Vector& curSoln = problem.getSolutionGroup().getX();
        
        if (Teuchos::is_null(updateVectorPtr)) 
                updateVectorPtr = curSoln.clone();

        updateVectorPtr->update(1.0, curSoln, -1.0, oldSoln, 0.0); 
        
        const Epetra_Vector& u =(dynamic_cast<const NOX::Epetra::Vector&>(curSoln)).getEpetraVector();
        
        const Epetra_Vector& du = (dynamic_cast<const NOX::Epetra::Vector&>(*updateVectorPtr)).getEpetraVector();
        
        normF = FPK->enorm(u, du);
  
        status = ((normF != -1) && (normF < specifiedTolerance)) ? NOX::StatusTest::Converged : NOX::StatusTest::Unconverged;

        return status;
}

NOX::StatusTest::StatusType PK_enorm::getStatus() const {
  return status;
}

ostream& PK_enorm::print(ostream& stream, int indent) const {
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "PK_enorm = " << NOX::Utils::sciformat(normF,3);
  stream << " < " << NOX::Utils::sciformat(specifiedTolerance, 3);
  stream << "\n";
  stream << endl;

  return stream;
}



}  // namespace AmanziFlow
}  // namespace Amanzi
