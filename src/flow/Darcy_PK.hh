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
  void fun(const double Told, double Tnew, 
           Teuchos::RCP<CompositeVector> u, Teuchos::RCP<CompositeVector> udot, 
           Teuchos::RCP<CompositeVector> rhs) {};
  void precon(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<CompositeVector> Hu) {};
  double enorm(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<const CompositeVector> du) { 
    return 0.0; 
  }
  void update_precon(double T, Teuchos::RCP<const CompositeVector> up, double h) {};
  void update_norm(double rtol, double atol) {};
  bool is_admissible(Teuchos::RCP<const CompositeVector> up) { 
   return false; 
  }
  bool modify_predictor(double h, Teuchos::RCP<CompositeVector> up) {
    return false;
  }
  bool modify_correction(double h, Teuchos::RCP<const CompositeVector> res,
                         Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<CompositeVector> du) {
    return false;
  }
  void changed_solution() {};

  // other main methods
  void AddTimeDerivativeSpecificStorage(Epetra_MultiVector& p, double dTp, Matrix_MFD* matrix_operator);
  void AddTimeDerivativeSpecificYield(Epetra_MultiVector& p, double dTp, Matrix_MFD* matrix_operator);
  void UpdateSpecificYield();

  double ErrorEstimate(double* dTfactor);

  // linear solvers
  void SolveFullySaturatedProblem(double T, CompositeVector& u);
  void SolveFullySaturatedProblem(double T, const CompositeVector& rhs, CompositeVector& u);
  int ApllyPrecInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) { Y = X; return 1; }
  void AssembleMatrixMFD();

  // control methods (for unit tests)
  void ResetParameterList(const Teuchos::ParameterList& dp_list_new) { dp_list_ = dp_list_new; }

  // access
  Teuchos::RCP<FlowMatrix> matrix() { return matrix_; }
  Teuchos::RCP<CompositeVector> get_solution() { return solution; }

 private:
  Teuchos::ParameterList dp_list_;

  Teuchos::RCP<FlowMatrix> matrix_;

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

