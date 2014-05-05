/*
  Authors: Ethan Coon (ecoon@lanl.gov)

  Decorator for using a Solver with JFNK as the preconditioner.
*/

#ifndef AMANZI_JFNK_SOLVER_
#define AMANZI_JFNK_SOLVER_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "SolverFnBase.hh"
#include "SolverDefs.hh"
#include "Solver.hh"
#include "SolverFnBaseJF.hh"
#include "SolverNKA.hh"
#include "SolverNKA_BT.hh"
#include "SolverNKA_BT_ATS.hh"
#include "SolverNewton.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class SolverJFNK : public Solver<Vector,VectorSpace> {

 public:
  SolverJFNK(Teuchos::ParameterList& plist) :
      plist_(plist) {}

  SolverJFNK(Teuchos::ParameterList& plist,
	     const Teuchos::RCP<SolverFnBase<Vector> >& fn,
	     const VectorSpace& map) :
      plist_(plist) {
    Init(fn, map);
  }

  void Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
            const VectorSpace& map);

  int Solve(const Teuchos::RCP<Vector>& u) {
    return solver_->Solve(u);
  }

  double residual() {
    return solver_->residual();
  }

  int num_itrs() {
    return solver_->num_itrs();
  }

  void set_pc_lag(double pc_lag) {
    return solver_->set_pc_lag(pc_lag);
  }

  int returned_code() {
    return solver_->returned_code();
  }
  
 protected:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<SolverFnBaseJF<Vector,VectorSpace> > jf_fnbase_;
  Teuchos::RCP<SolverFnBase<Vector> > wrapped_fnbase_;
  Teuchos::RCP<Solver<Vector,VectorSpace> > solver_;
};


/* ******************************************************************
* Public Init method.
****************************************************************** */
template<class Vector, class VectorSpace>
void
SolverJFNK<Vector,VectorSpace>::Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
        const VectorSpace& map)
{
  // create the nonlinear solver
  Teuchos::ParameterList& slist = plist_.sublist("nonlinear solver");

  // NOTE that we recreate much of SolverFactory here to avoid recursion.
  if (slist.isParameter("solver type")) {
    std::string type = slist.get<std::string>("solver type");

    if (type == "nka") {
      if (!slist.isSublist("nka parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"nka parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList nka_list = slist.sublist("nka parameters");
      solver_ = Teuchos::rcp(new SolverNKA<Vector,VectorSpace>(nka_list));

    } else if (type == "Newton") {
      if (!slist.isSublist("Newton parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"Newton parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList newton_list = slist.sublist("Newton parameters");
      solver_ = Teuchos::rcp(new SolverNewton<Vector,VectorSpace>(newton_list));

    } else if (type == "nka_bt") {
      if (!slist.isSublist("nka_bt parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"nka_bt parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList nka_list = slist.sublist("nka_bt parameters");
      solver_ = Teuchos::rcp(new SolverNKA_BT<Vector,VectorSpace>(nka_list));

    } else if (type == "nka_bt_ats") {
      if (!slist.isSublist("nka_bt_ats parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"nka_bt_ats parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList nka_list = slist.sublist("nka_bt_ats parameters");
      solver_ = Teuchos::rcp(new SolverNKA_BT_ATS<Vector,VectorSpace>(nka_list));

    } else {
      Errors::Message msg("SolverFactory: wrong value of parameter `\"solver type`\"");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg("SolverFactory: parameter `\"solver type`\" is missing");
    Exceptions::amanzi_throw(msg);
  }

  // save the base fn
  wrapped_fnbase_ = fn;

  // wrap the base fn in a JF fn
  jf_fnbase_ = Teuchos::rcp(new SolverFnBaseJF<Vector,VectorSpace>(plist_,wrapped_fnbase_, map));

  // Init the nonlinear method with this wrapped SolverFnBase
  solver_->Init(jf_fnbase_, map);
  
}

} // namespace
} // namespace

#endif
