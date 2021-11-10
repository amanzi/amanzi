/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_BDF1_SOLVER_FNBASE_
#define AMANZI_BDF1_SOLVER_FNBASE_

#include "Teuchos_RCP.hpp"

#include "SolverFnBase.hh"
#include "BDFFnBase.hh"

#include "AmanziDebug.hh"


namespace Amanzi {

template <class Vector>
class BDF1_SolverFnBase : public AmanziSolvers::SolverFnBase<Vector> {
 public:
  BDF1_SolverFnBase(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                    const Teuchos::RCP<BDFFnBase<Vector>>& bdf_fn)
    : plist_(plist), bdf_fn_(bdf_fn){};

  // SolverFnBase interface
  // ---------------------------
  // computes the non-linear functional r = F(u)
  virtual void
  Residual(const Teuchos::RCP<Vector>& u, const Teuchos::RCP<Vector>& r);

  // preconditioner application
  virtual int ApplyPreconditioner(const Teuchos::RCP<const Vector>& r,
                                  const Teuchos::RCP<Vector>& Pr);

  // preconditioner update
  virtual void UpdatePreconditioner(const Teuchos::RCP<const Vector>& u);

  // error norm
  virtual double ErrorNorm(const Teuchos::RCP<const Vector>& u,
                           const Teuchos::RCP<const Vector>& du);

  // Check the admissibility of an inner iterate (ensures preconditions for
  // F(u) to be defined).
  virtual bool IsAdmissible(const Teuchos::RCP<const Vector>& up);

  // Hack a correction for some reason.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(const Teuchos::RCP<const Vector>& res,
                   const Teuchos::RCP<const Vector>& u,
                   const Teuchos::RCP<Vector>& du);

  // parameter for continuation method
  virtual void UpdateContinuationParameter(double lambda);

  // bookkeeping for state
  virtual void ChangedSolution();

  // Other methods
  void SetTimes(double t_old, double t_new)
  {
    t_old_ = t_old;
    t_new_ = t_new;
    h_ = t_new - t_old;
  }

  void SetPreviousTimeSolution(const Teuchos::RCP<Vector>& u_old)
  {
    u_old_ = u_old;
  }

 protected:
  Teuchos::RCP<Teuchos::ParameterList> plist_;

  double t_new_;
  double t_old_;
  double h_;
  Teuchos::RCP<Vector> u_old_;

  Teuchos::RCP<BDFFnBase<Vector>> bdf_fn_;
};


// computes the non-linear functional r = F(u)
template <class Vector>
void
BDF1_SolverFnBase<Vector>::Residual(const Teuchos::RCP<Vector>& u,
                                    const Teuchos::RCP<Vector>& r)
{
  // std::cout << "BDF1_SolverFnBase::Residual u = " << Debug::get0(*u) << std::endl;
  bdf_fn_->FunctionalResidual(t_old_, t_new_, u_old_, u, r);
  // std::cout << "BDF1_SolverFnBase::Residual r = " << Debug::get0(*r) << std::endl;
}

// preconditioner application
template <class Vector>
int
BDF1_SolverFnBase<Vector>::ApplyPreconditioner(
  const Teuchos::RCP<const Vector>& r, const Teuchos::RCP<Vector>& Pr)
{
  return bdf_fn_->ApplyPreconditioner(r, Pr);
}

// preconditioner update
template <class Vector>
void
BDF1_SolverFnBase<Vector>::UpdatePreconditioner(
  const Teuchos::RCP<const Vector>& u)
{
  bdf_fn_->UpdatePreconditioner(t_new_, u, h_);
}

// error norm
template <class Vector>
double
BDF1_SolverFnBase<Vector>::ErrorNorm(const Teuchos::RCP<const Vector>& u,
                                     const Teuchos::RCP<const Vector>& du)
{
  return bdf_fn_->ErrorNorm(u, du);
}


// Check the admissibility of an inner iterate (ensures preconditions for
// F(u) to be defined).
template <class Vector>
bool
BDF1_SolverFnBase<Vector>::IsAdmissible(const Teuchos::RCP<const Vector>& up)
{
  return bdf_fn_->IsAdmissible(up);
}

// Hack a correction for some reason.
template <class Vector>
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
BDF1_SolverFnBase<Vector>::ModifyCorrection(
  const Teuchos::RCP<const Vector>& res, const Teuchos::RCP<const Vector>& u,
  const Teuchos::RCP<Vector>& du)
{
  return bdf_fn_->ModifyCorrection(h_, res, u, du);
}

template <class Vector>
void
BDF1_SolverFnBase<Vector>::UpdateContinuationParameter(double lambda)
{
  bdf_fn_->UpdateContinuationParameter(lambda);
}

// bookkeeping for state
template <class Vector>
void
BDF1_SolverFnBase<Vector>::ChangedSolution()
{
  bdf_fn_->ChangedSolution();
}

} // namespace Amanzi
#endif
