#include "RichardsModelEvaluator.hpp"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#  include "mpi.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI

#include "Epetra_Vector.h"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


RichardsModelEvaluator::RichardsModelEvaluator(RichardsProblem *problem, 
					       Teuchos::ParameterList &plist, 
					       const Epetra_Map &map,
					       Teuchos::RCP<const Flow_State> FS) 
  : problem_(problem), D(problem->Matrix()),  map_(map), plist_(plist),
    FS_(FS)
{
  this->setLinePrefix("RichardsModelEvaluator");
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist_,this);
  
}

void RichardsModelEvaluator::initialize(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params)
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab  
  
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
    {
      *out << "initialize o.k." << std::endl;
    }
  

}

// Overridden from BDF2::fnBase

void RichardsModelEvaluator::fun(const double t, const Epetra_Vector& u, 
				 const Epetra_Vector& udot, Epetra_Vector& f) 
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab  
  
 
  // compute F(u)
  problem_->ComputeF(u, f);
 

  Epetra_Vector *uc     = problem_->CreateCellView(u);  
  Epetra_Vector *udotc  = problem_->CreateCellView(udot);
  Epetra_Vector *fc     = problem_->CreateCellView(f);
  
  // compute S'(p)
  Epetra_Vector dS (problem_->CellMap());
  problem_->dSofP(*uc, dS);
  const Epetra_Vector& phi = FS_->porosity();
  double rho;
  problem_->GetFluidDensity(rho);

  // assume that porosity is piecewise constant
  
  dS.Multiply(0.0,dS,phi,rho);
  
  dS.PutScalar(1.0);

  dS.Multiply(1.0,dS,*(problem_->cell_vols()),0.0);

  // on the cell unknowns compute f=f+dS*udotc*rho*phi
  fc->Multiply(1.0,dS,*udotc,1.0);
  
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
    {
      *out << "fun o.k." << std::endl;
    }

}

void RichardsModelEvaluator::precon(const Epetra_Vector& X, Epetra_Vector& Y)
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab  


  (problem_->Precon()).ApplyInverse(X, Y);

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
    {
      *out << "precon o.k." << std::endl;
    }


}

void RichardsModelEvaluator::update_precon(const double t, const Epetra_Vector& up, const double h, int& errc)
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab  


  problem_->ComputePrecon(up,h);

  errc = 0;

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
    {
      *out << "update_precon done" << std::endl;
    }


}



double RichardsModelEvaluator::enorm(const Epetra_Vector& u, const Epetra_Vector& du)
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab  


  // simply use 2-norm of the difference for now
  
  Epetra_Vector *u_cell  = problem_->CreateCellView(u);
  Epetra_Vector *du_cell = problem_->CreateCellView(du);
  
  double atol = 0.00001;
  double rtol = 0.0;

  double en = 0.0;
  for (int j=0; j<u.MyLength(); j++)
    {
      double tmp = abs(du[j])/(atol+rtol*abs(u[j]));
      en = std::max<double>(en, tmp);
    }
  
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
    {
      *out << "enorm done" << std::endl;
    }

  return  en;

}


bool RichardsModelEvaluator::is_admissible(const Epetra_Vector& up)
{
  return true;
}
