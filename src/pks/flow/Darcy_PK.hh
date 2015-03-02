/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef AMANZI_DARCY_PK_HH_
#define AMANZI_DARCY_PK_HH_

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "FnBaseDefs.hh"
#include "Operator.hh"
#include "OperatorDiffusion.hh"
#include "OperatorAccumulation.hh"

#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

class Darcy_PK : public Flow_PK {
 public:
  Darcy_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
           const std::string& pk_list_name,
           Teuchos::RCP<State> S);
  ~Darcy_PK();

  // main PK methods
  void Setup();
  void Initialize();
  bool Advance(double dT, double &dT_actual); 

  void set_dt(double dt){dT = dt; dT_desirable_ = dT;}
  double get_dt() { return dT_desirable_; }

  void CommitState(double dt, const Teuchos::Ptr<State>& S);
  void CalculateDiagnostics(const Teuchos::Ptr<State>& S);

  // main flow methods
  void InitTimeInterval();
  void InitializeAuxiliaryData();

  // methods required for time integration
  void Functional(const double Told, double Tnew, 
                  Teuchos::RCP<CompositeVector> u, Teuchos::RCP<CompositeVector> udot, 
                  Teuchos::RCP<CompositeVector> rhs) {};
  void ApplyPreconditioner(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<CompositeVector> Hu) {};
  double ErrorNorm(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<const CompositeVector> du) { 
    return 0.0; 
  }
  void UpdatePreconditioner(double T, Teuchos::RCP<const CompositeVector> up, double h) {};
  void update_norm(double rtol, double atol) {};
  bool IsAdmissible(Teuchos::RCP<const CompositeVector> up) { 
   return false; 
  }
  bool ModifyPredictor(double dT, Teuchos::RCP<const CompositeVector> u0, Teuchos::RCP<CompositeVector> u) {
    return false;
  }
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dT, Teuchos::RCP<const CompositeVector> res,
                       Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<CompositeVector> du) {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }
  void ChangedSolution() {};

  // linear solvers
  void SolveFullySaturatedProblem(double T, CompositeVector& u);
  int ApllyPrecInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) { Y = X; return 1; }

 private:
  void UpdateSpecificYield_();
  double ErrorEstimate_(double* dTfactor);

 private:
  Teuchos::RCP<Teuchos::ParameterList> dp_list_;
  Teuchos::RCP<Operators::Operator> op_;
  Teuchos::RCP<Operators::OperatorDiffusion> op_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> op_acc_;
  Teuchos::RCP<Operators::BCs> op_bc_;

  int error_control_;
  double dT_desirable_, dTmax_, dTfactor_;
  std::vector<std::pair<double, double> > dT_history_;  // statistics

  std::string preconditioner_name_, solver_name_;
  bool initialize_with_darcy_;
  int num_itrs_;

  Teuchos::RCP<CompositeVector> solution;  // next pressure state
  Teuchos::RCP<Epetra_Vector> pdot_cells_prev;  // time derivative of pressure
  Teuchos::RCP<Epetra_Vector> pdot_cells;

  Teuchos::RCP<CompositeVector> specific_yield_copy_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif

