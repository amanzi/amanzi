/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Nonlinear Krylov Accelerator nonlinear solver.

/*!

Uses an accelerated nonlinear solver that is similar to a multidimensional
secant method.

 * `"max nka vectors`" ``[int]`` **10** Number of vectors used to span the
      secant space of previously explored directions.

 * `"nka vector tolerance`" ``[double]`` **0.05** Tolerance on when the dot
      product of two vectors suggests they are parallel.

 * `"nka lag iterations`" ``[int]`` **0** Number of iterations to wait to begin
        the NKA accelerator.  Note that no theoretical NKA work has been done
        with globalization, though experience suggests that it is robust to
        globalized corrections as long as the actual correction taken is given
        to the NKA algorithm.  Lagging NKA may be useful if the first few
        Jacobian corrections are always expected to be heavily modified via
        globalization.


 [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
      weighted moving finite element code I: in one dimension", SIAM J.
      Sci. Comput;, 19 (1998), pp. 728-765.  See section 9.
        
*/

#ifndef AMANZI_NKA_SOLVER_
#define AMANZI_NKA_SOLVER_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "SolverDefs.hh"
#include "SolverDefault.hh"
#include "NKA_Base.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector, class VectorSpace>
class SolverNKA : public SolverDefault<Vector, VectorSpace> {
 public:
  SolverNKA(Teuchos::ParameterList& plist)
      : SolverDefault<Vector,VectorSpace>(plist) {
    nka_dim_ = plist.get<int>("max nka vectors", 10);
    nka_dim_ = std::min<int>(nka_dim_, this->max_itrs_ - 1);
    nka_tol_ = plist.get<double>("nka vector tolerance", 0.05);
    nka_lag_ = plist.get<int>("nka lag iterations", 0);

    // update the verbose options
    this->vo_ = Teuchos::rcp(new VerboseObject(name(), plist));
  }

  SolverNKA(Teuchos::ParameterList& plist,
            const Teuchos::RCP<SolverFnBase<Vector>>& fn,
            const VectorSpace& map)
      : SolverNKA(plist)
  {
    Init(fn, map);
  }

  virtual void Init(const Teuchos::RCP<SolverFnBase<Vector>>& fn,
            const Teuchos::RCP<const VectorSpace>& map) override;

  virtual std::string name() const override { return "Solver::NKA"; }

 protected:  
  //
  // NKA's ModifyCorrection() calls the acceleration on the correction.
  virtual std::pair<MonitorStatus,double> ModifyCorrection_(Teuchos::RCP<Vector>& r,
          const Teuchos::RCP<Vector>& u, Teuchos::RCP<Vector>& du) override;

  virtual void Restart_() override { nka_->Restart(); }
  
 protected:
  Teuchos::RCP<NKA_Base<Vector, VectorSpace>> nka_;
  Teuchos::RCP<Vector> du_tmp_;

  double nka_tol_;
  int nka_dim_;
  int nka_lag_;
};


/* ******************************************************************
 * Public Init method.
 ****************************************************************** */
template <class Vector, class VectorSpace>
void
SolverNKA<Vector, VectorSpace>::Init(
  const Teuchos::RCP<SolverFnBase<Vector>>& fn,
  const Teuchos::RCP<const VectorSpace>& map)
{
  SolverDefault<Vector,VectorSpace>::Init(fn, map);

  // Allocate the NKA space
  nka_ = Teuchos::rcp(new NKA_Base<Vector, VectorSpace>(nka_dim_, nka_tol_, map));
  du_tmp_ = Teuchos::rcp(new Vector(map));
}


/* ******************************************************************
 * NKA's ModifyCorrection() calls the acceleration on the correction.
 ****************************************************************** */
template <class Vector, class VectorSpace>
std::pair<MonitorStatus,double>
SolverNKA<Vector, VectorSpace>::ModifyCorrection_(Teuchos::RCP<Vector>& r,
        const Teuchos::RCP<Vector>& u,
        Teuchos::RCP<Vector>& du)
{
  // Calculate the accelerated correction.
  if (this->num_itrs_ > nka_lag_) {
    // Calculate the NKA correction
    nka_->Correction(*du, *du_tmp_, du_tmp_.ptr());

    // Call the PK's modify
    if (this->modify_correction_) {
      bool hacked = this->fn_->ModifyCorrection(r, u, du_tmp_);
      if (hacked) {
        // If we had to hack things, it't not unlikely that the Jacobian
        // information is not representative of the local space. Take the hacked
        // correction, and restart NKA to start building a new Jacobian space.
        nka_->Restart();
      }
    }

    // Check the admissibility of the NKA iterate
    if (this->fn_->IsAdmissible(du_tmp_)) {
      du->assign(*du_tmp_); // NOTE: should be able to remove this if we are
                          // careful about RCPs and use swap instead? --etc
      return std::make_pair(MonitorStatus::CONTINUE, (double) -1.0);

    } else {
      // Yuck.  Restart NKA
      nka_->Restart();
    }
  }

  // Either we aren't using the NKA iteration due to lag, or because it was
  // bad.  Try the direct preconditioned correction.
  if (this->modify_correction_) {
    bool hacked = this->fn_->ModifyCorrection(r, u, du);
  }
  if (this->fn_->IsAdmissible(du)) return std::make_pair(MonitorStatus::CONTINUE, (double) -1.0);
  else return std::make_pair(MonitorStatus::INADMISSIBLE_SOLUTION, (double) -1.0);
}

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
