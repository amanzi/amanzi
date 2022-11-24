#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "exceptions.hh"
#include "FnBaseDefs.hh"
#include "SolverJFNK.hh"
#include "SolverFnBase.hh"

using namespace Amanzi;

SUITE(SOLVERS)
{
  class LinearFn : public AmanziSolvers::SolverFnBase<Epetra_Vector> {
   public:
    LinearFn() {}
    LinearFn(const Teuchos::RCP<Epetra_Map>& map) : map_(map)
    {
      rhs_ = Teuchos::rcp(new Epetra_Vector(*map_));
      soln_ = Teuchos::rcp(new Epetra_Vector(*map_));
    };
    ~LinearFn(){};
    LinearFn(const LinearFn& other) : map_(other.map_) {}

    Teuchos::RCP<LinearFn> Clone() const { return Teuchos::rcp(new LinearFn(*this)); }

    // computes the non-linear functional r = F(u)
    virtual void
    Residual(const Teuchos::RCP<Epetra_Vector>& u, const Teuchos::RCP<Epetra_Vector>& r)
    {
      Apply(*u, *r);
      r->Update(-1., *rhs_, 1.);
    }

    // preconditioner toolkit
    virtual int ApplyPreconditioner(const Teuchos::RCP<const Epetra_Vector>& r,
                                    const Teuchos::RCP<Epetra_Vector>& Pr)
    {
      return ApplyInverse(*r, *Pr);
    }
    virtual void UpdatePreconditioner(const Teuchos::RCP<const Epetra_Vector>& u) {}

    // error analysis
    virtual double ErrorNorm(const Teuchos::RCP<const Epetra_Vector>& u,
                             const Teuchos::RCP<const Epetra_Vector>& du)
    {
      double res;
      du->NormInf(&res);
      return res;
    }

    virtual int Apply(const Epetra_Vector& v, Epetra_Vector& mv) const
    {
      for (int i = 0; i < 5; i++) mv[i] = 2 * v[i];
      for (int i = 1; i < 5; i++) mv[i] -= v[i - 1];
      for (int i = 0; i < 4; i++) mv[i] -= v[i + 1];
      return 0;
    }
    virtual int ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
    {
      for (int i = 0; i < 5; i++) hv[i] = v[i];
      return 0;
    }

    void SetSolution(const Epetra_Vector& soln)
    {
      *soln_ = soln;
      Apply(soln, *rhs_);
    }

    bool CheckSolution(const Epetra_Vector& soln)
    {
      for (int i = 0; i < 5; i++) { CHECK_CLOSE(soln[i], (*soln_)[i], 1e-6); }
      return true;
    }

    virtual const Epetra_Map& DomainMap() const { return *map_; }
    virtual const Epetra_Map& RangeMap() const { return *map_; }

   private:
    Teuchos::RCP<Epetra_Map> map_;
    Teuchos::RCP<Epetra_Vector> rhs_;
    Teuchos::RCP<Epetra_Vector> soln_;
  };

  TEST(JFNK_SOLVER_LINEAR)
  {
    std::cout << "Checking JFNK solver on linear problem..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(5, 0, *comm));

    // create the SolverFnBase object
    Teuchos::RCP<LinearFn> fn = Teuchos::rcp(new LinearFn(map));

    Epetra_Vector soln(*map);
    soln[0] = -0.1666666666666;
    soln[1] = 0.6666666666666;
    soln[2] = 0.5;
    soln[3] = 0.33333333333333;
    soln[4] = 0.16666666666666;
    fn->SetSolution(soln);

    // create the JFNK object
    Teuchos::ParameterList jfnk_list;
    jfnk_list.sublist("nonlinear solver").set("solver type", "Newton");
    jfnk_list.sublist("nonlinear solver")
      .sublist("Newton parameters")
      .sublist("verbose object")
      .set("verbosity level", "extreme");
    jfnk_list.sublist("nonlinear solver")
      .sublist("Newton parameters")
      .set("monitor", "monitor residual");
    jfnk_list.sublist("JF matrix parameters").set("finite difference epsilon", 1e-7);
    jfnk_list.sublist("linear operator").set("iterative method", "pcg");
    jfnk_list.sublist("linear operator")
      .sublist("pcg parameters")
      .sublist("verbose object")
      .set("verbosity level", "extreme");

    Teuchos::RCP<AmanziSolvers::SolverJFNK<Epetra_Vector, Epetra_Map>> jfnk =
      Teuchos::rcp(new AmanziSolvers::SolverJFNK<Epetra_Vector, Epetra_Map>(jfnk_list));
    jfnk->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*map));
    jfnk->Solve(u);

    fn->CheckSolution(*u);
    CHECK_EQUAL(1, jfnk->num_itrs());

    delete comm;
  };
}
