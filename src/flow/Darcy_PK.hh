/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef AMANZI_DARCY_PK_HH_
#define AMANZI_DARCY_PK_HH_

#include "Teuchos_RCP.hpp"

#include "Epetra_Vector.h"
#include "AztecOO.h"

#include "Mesh.hh"
#include "Point.hh"
#include "tensor.hh"

#include "Flow_PK.hh"
#include "Matrix_MFD.hh"
#include "TI_Specs.hh"


namespace Amanzi {
namespace AmanziFlow {

class Darcy_PK : public Flow_PK {
 public:
  Darcy_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S);
  ~Darcy_PK();

  // main methods
  void InitPK();
  void InitSteadyState(double T0, double dT0);
  void InitTransient(double T0, double dT0);
  void InitPicard(double T0) {};  // not used yet.
  void InitNextTI(double T0, double dT0, TI_Specs ti_specs);

  double CalculateDt() { return dT_desirable_; }
  int Advance(double dT); 
  int AdvanceToSteadyState(double T0, double dT0);
  void InitializeAuxiliaryData();
  void InitializeSteadySaturated();

  void CommitState(Teuchos::RCP<State> S);

  // methods required for time integration
  void fun(const double T, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& rhs, double dT = 0.0) {};
  void precon(const Epetra_Vector& u, Epetra_Vector& Hu) {};
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du) { return 0.0; }
  void update_precon(const double T, const Epetra_Vector& up, const double h, int& errc) {};
  void update_norm(double rtol, double atol) {};

  // other main methods
  void AddTimeDerivativeSpecificStorage(Epetra_MultiVector& p_cells, double dTp, Matrix_MFD* matrix_operator);
  void AddTimeDerivativeSpecificYield(Epetra_MultiVector& p_cells, double dTp, Matrix_MFD* matrix_operator);
  void UpdateSpecificYield();

  double ErrorEstimate(double* dTfactor);

  // linear solvers
  void SolveFullySaturatedProblem(double T, CompositeVector& u);
  void SolveFullySaturatedProblem(double T, const Epetra_Vector& rhs, Epetra_Vector& u);
  int ApllyPrecInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) { Y = X; return 1; }
  void AssembleMatrixMFD();

  // control methods (for unit tests)
  void ResetParameterList(const Teuchos::ParameterList& dp_list_new) { dp_list_ = dp_list_new; }

 private:
  Teuchos::ParameterList dp_list_;

  Teuchos::RCP<Matrix_MFD> matrix_;

  int error_control_;
  double dT_desirable_;

  Teuchos::RCP<CompositeVector> solution;  // next pressure state
  Teuchos::RCP<Epetra_Vector> pdot_cells_prev;  // time derivative of pressure
  Teuchos::RCP<Epetra_Vector> pdot_cells;
 
  Teuchos::RCP<Epetra_IntVector> upwind_cell, downwind_cell;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

