/*
  Navier Stokes PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_NAVIER_STOKES_PK_HH_
#define AMANZI_NAVIER_STOKES_PK_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BDF1_TI.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "PK_PhysicalBDF.hh"
#include "primary_variable_field_evaluator.hh"
#include "State.hh"
#include "TreeVector.hh"
#include "Units.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace NavierStokes {

class NavierStokes_PK {
 public:
  NavierStokes_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  NavierStokes_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const std::string& pk_list_name,
                  Teuchos::RCP<State> S,
                  const Teuchos::RCP<TreeVector>& soln);

  ~NavierStokes_PK() {};

  // methods required for PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual double get_dt() { return dt_; }
  virtual void set_dt(double dt) { dt_ = dt; dt_desirable_ = dt_; }

  // methods required for time integration interface
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  void Functional(const double t_old, double t_new, 
                  Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
                  Teuchos::RCP<TreeVector> f);
  double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);
  
  // -- management of the preconditioner
  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu);
  void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt);

  // other methods
  // -- access
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae() { return bdf1_dae_; }

 public:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> ns_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<const Teuchos::ParameterList> linear_solver_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_owned, nnodes_wghost;

  double dt_, dt_next_, dt_desirable_;

 protected:
  Teuchos::RCP<TreeVector> soln_;

  Teuchos::RCP<PrimaryVariableFieldEvaluator> pressure_eval_, fluid_velocity_eval_;
 
 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<State> S_;
  std::string passwd_;
  int dim;

  // time integrators
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae_;

  // io
  Utils::Units units_;
  Teuchos::RCP<VerboseObject> vo_;

  // factory registration
  static RegisteredPKFactory<NavierStokes_PK> reg_;
};

}  // namespace NavierStokes
}  // namespace Amanzi

#endif
