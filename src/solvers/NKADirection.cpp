#include "NKADirection.H"
#include "Teuchos_ENull.hpp"

NKADirection::NKADirection(const Teuchos::RCP<NOX::GlobalData>& gd, 
			   Teuchos::ParameterList& param, 
			   const NOX::Abstract::Vector& initvec)
{
  // get the parameters the parameter list
  
  int mvec    = param.get<int>("maxv", 10);
  double vtol = param.get<double>("vtol", 1e-4);
  
  state = new nka (mvec, vtol, initvec);
  
  reset(gd, param);
  
};
  

NKADirection::~NKADirection() 
{
  delete state;
};


bool NKADirection::reset(const Teuchos::RCP<NOX::GlobalData> &gd, 
			 Teuchos::ParameterList &param)
{
  globalDataPtr = gd;
  utils = gd->getUtils();
  paramPtr = &param;
  
  state->nka_restart();
  
  return true;
  
};

  
bool NKADirection::compute (NOX::Abstract::Vector &dir, 
			    NOX::Abstract::Group &soln, 
			    const NOX::Solver::Generic  &solver)
{
  NOX::Abstract::Group::ReturnType status;
  
  // evaluate the nonlinear functional
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) {
    throwError("compute", "Unable to compute F");
  }

  Teuchos::RCP<NOX::Abstract::Vector> fptr 
    = soln.getF().clone(NOX::DeepCopy);
  
  Teuchos::RCP<NOX::Abstract::Vector> precond_fptr 
    = (*fptr).clone(NOX::ShapeCopy);

  status 
    = soln.applyRightPreconditioning(false,
				     paramPtr->sublist("Linear Solver"), 
				     *fptr, *precond_fptr);
  
  if (status != NOX::Abstract::Group::Ok) {
    throwError("compute", "Unable to apply preconditioner"); 
  } 

  state->nka_correction(dir, precond_fptr);
  dir.scale(-1.0);

  return true;

};
  
bool NKADirection::compute (NOX::Abstract::Vector  &dir, 
			    NOX::Abstract::Group  &soln, 
			    const NOX::Solver::LineSearchBased  &solver)
{
  return NOX::Direction::Generic::compute( dir, soln, solver );
};

  
void NKADirection::throwError(const string& functionName, 
			      const string& errorMsg)
{
  if (utils->isPrintType(NOX::Utils::Error))
    utils->err() << "NKADirection::" << functionName 
		 << " - " << errorMsg << endl;
  throw "NOX Error";
}
