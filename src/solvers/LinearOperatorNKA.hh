/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Uses NKA method as a linear solver.

/*!

This is effectively equivalent to GMRES with a rolling restart, where vectors
fall off the end of the space.

Parameters:

* `"error tolerance`" ``[double]`` **1.e-6** Tolerance on which to declare success.

* `"maximum number of iterations`" ``[int]`` **100** Maximum iterations before declaring failure.

* `"overflow tolerance`" ``[double]`` **3.e50** Error above this value results in failure.

* `"convergence criterial`" ``[Array(string)]`` **"{relative rhs}"** A list of
  criteria, any of which can be applied.  Valid include:

  - `"relative rhs`" : measure error relative to the norm of the RHS vector
  - `"relative residual`" : measure error relative to the norm of the residual
  - `"absolute residual`" : measure error directly, norm of error
  - `"make one iteration`" : require at least one iteration to be performed before declaring success

* `"max nka vectors`" ``[int]`` **10** Size of the NKA space used to span the residual, conceptually equivalent to the size of the Krylov space.

* `"nka vector tolerance`" ``[double]`` **0.05** Vectors whose dot product are within this tolerance are considered parallel, and therefore the old vector is thrown out.


Calef et al. "Nonlinear Krylov acceleration applied to a discrete ordinates formulation of the k-eigenvalue problem." JCP 238 (2013): 188-209.

*/

#ifndef AMANZI_NKA_OPERATOR_HH_
#define AMANZI_NKA_OPERATOR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "VerboseObject.hh"

#include "LinearOperatorDefs.hh"
#include "LinearOperator.hh"
#include "NKA_Base.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class LinearOperatorNKA : public LinearOperator<Matrix, Vector, VectorSpace> {

 public:
  LinearOperatorNKA(const Teuchos::RCP<const Matrix>& m,
                    const Teuchos::RCP<const Matrix>& h) :
      LinearOperator<Matrix, Vector, VectorSpace>(m, h),
      tol_(1e-6),
      overflow_tol_(3.0e+50),  // mass of the Universe (J.Hopkins)
      max_itrs_(100),
      nka_dim_(10),
      criteria_(LIN_SOLVER_RELATIVE_RHS),
      initialized_(false) {};

  void Init(Teuchos::ParameterList& plist);
  void Init() { LinearOperator<Matrix,Vector,VectorSpace>::Init(); }

  int ApplyInverse(const Vector& v, Vector& hv) const {
    if (!initialized_) {
      Errors::Message msg("LinearOperatorNKA: has not been initialized.");
      Exceptions::amanzi_throw(msg);
    }
    int ierr = NKA_(v, hv, tol_, max_itrs_, criteria_);
    returned_code_ = ierr;
    //return (ierr > 0) ? 0 : 1;  // Positive ierr code means success.
    return returned_code_;
  }

  // mutators
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_criteria(int criteria) { criteria_ = criteria; }
  void add_criteria(int criteria) { criteria_ |= criteria; }

  // accessors
  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }
  int returned_code() { return returned_code_; }

 public:
  Teuchos::RCP<VerboseObject> vo_;

 private:
  int NKA_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;

  LinearOperatorNKA(const LinearOperatorNKA& other); // not implemented

 private:
  using LinearOperator<Matrix, Vector, VectorSpace>::m_;
  using LinearOperator<Matrix, Vector, VectorSpace>::h_;
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  int max_itrs_, nka_dim_, criteria_;
  double tol_, nka_tol_, overflow_tol_;
  mutable int num_itrs_, returned_code_;
  mutable double residual_;
  mutable bool initialized_;

  Teuchos::RCP<NKA_Base<Vector,VectorSpace> > nka_;
};



// Apply the inverse, x <-- A^-1 b
template<class Matrix, class Vector, class VectorSpace>
int LinearOperatorNKA<Matrix, Vector, VectorSpace>::NKA_(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  Teuchos::OSTab tab = vo_->getOSTab();

  //  AMANZI_ASSERT(f.Map().SameAs(m_->RangeMap()));
  //  AMANZI_ASSERT(x.Map().SameAs(m_->DomainMap()));
  nka_->Restart();

  residual_ = 0.0;
  num_itrs_ = 0;

  Teuchos::RCP<Vector> dx  = Teuchos::rcp(new Vector(x));
  Teuchos::RCP<Vector> dxp = Teuchos::rcp(new Vector(x));  // preconditioned correction
  Teuchos::RCP<Vector> r   = Teuchos::rcp(new Vector(x));

  double fnorm, xnorm;
  f.Norm2(&fnorm);
  if (fnorm == 0.0) {
    x.PutScalar(0.0);
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Converged, itr = " << num_itrs_ << ", ||r|| = " << residual_ << std::endl;
    return criteria;  // Zero solution satifies all criteria.
  }

  x.Norm2(&xnorm);

  int ierr = m_->Apply(x, *r);  // r = f - A * x
  AMANZI_ASSERT(!ierr);
  r->Update(1.0, f, -1.0);

  double rnorm0;
  r->Norm2(&rnorm0);
  residual_ = rnorm0;

  if (initialized_) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << num_itrs_ << " ||r||=" << residual_ << std::endl;
    }
  }
  if (criteria & LIN_SOLVER_RELATIVE_RHS) {
    if (rnorm0 < tol * fnorm) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Converged (relative ||RHS|| = " << fnorm << "), itr = " << num_itrs_ << " ||r|| = " << rnorm0 << std::endl;
      return LIN_SOLVER_RELATIVE_RHS;
    }
  }

  if (residual_ > overflow_tol_) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Diverged, ||r|| = " << residual_ << std::endl;
    return LIN_SOLVER_RESIDUAL_OVERFLOW;
  }

  bool done = false;
  while (!done) {
    ierr = h_->ApplyInverse(*r, *dxp);
    AMANZI_ASSERT(!ierr);

    nka_->Correction(*dxp, *dx);
    x.Update(1.0, *dx, 1.0);

    ierr = m_->Apply(x, *r);  // r = f - A * x
    AMANZI_ASSERT(!ierr);
    r->Update(1.0, f, -1.0);

    double rnorm;
    r->Norm2(&rnorm);
    residual_ = rnorm;

    num_itrs_++;

    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      *vo_->os() << num_itrs_ << " ||r||=" << residual_ << std::endl;
    }

    if (rnorm > overflow_tol_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Diverged, ||r|| = " << rnorm << std::endl;
      return LIN_SOLVER_RESIDUAL_OVERFLOW;
    }

    // Return the first criterion which is fulfilled.
    if (criteria & LIN_SOLVER_RELATIVE_RHS) {
      if (rnorm < tol * fnorm) {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (relative ||RHS|| = " << fnorm << "), itr = " << num_itrs_ << " ||r|| = " << rnorm << std::endl;
        return LIN_SOLVER_RELATIVE_RHS;
      }
    }

    if (criteria & LIN_SOLVER_RELATIVE_RESIDUAL) {
      if (rnorm < tol * rnorm0) {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (relative ||res|| = " << rnorm0 << "), itr = " << num_itrs_ << " ||r|| = " << rnorm << std::endl;
        return LIN_SOLVER_RELATIVE_RESIDUAL;
      }
    }
    if (criteria & LIN_SOLVER_ABSOLUTE_RESIDUAL) {
      if (rnorm < tol) {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (absolute), itr = " << num_itrs_ << " ||r|| = " << rnorm << std::endl;
        return LIN_SOLVER_ABSOLUTE_RESIDUAL;
      }
    }

    done = num_itrs_ > max_itrs;
  }

  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "Failed (" << num_itrs_ << " itrs) ||r|| = "
               << residual_ << std::endl;
  return LIN_SOLVER_MAX_ITERATIONS;
}


/* ******************************************************************
 * Initialization from a parameter list. Available parameters:
 * "error tolerance" [double] default = 1e-6
 * "maximum number of iterations" [int] default = 100
 * "convergence criteria" Array(string) default = "{relative rhs}"
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorNKA<Matrix, Vector, VectorSpace>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Solvers::NKA", plist));

  tol_ = plist.get<double>("error tolerance", 1e-6);
  max_itrs_ = plist.get<int>("maximum number of iterations", 100);
  overflow_tol_ = plist.get<double>("overflow tolerance", 3.0e+50);

  int criteria(0);
  if (plist.isParameter("convergence criteria")) {
    std::vector<std::string> names;
    names = plist.get<Teuchos::Array<std::string> > ("convergence criteria").toVector();

    for (int i = 0; i < names.size(); i++) {
      if (names[i] == "relative rhs") {
        criteria += LIN_SOLVER_RELATIVE_RHS;
      } else if (names[i] == "relative residual") {
        criteria += LIN_SOLVER_RELATIVE_RESIDUAL;
      } else if (names[i] == "absolute residual") {
        criteria += LIN_SOLVER_ABSOLUTE_RESIDUAL;
      } else if (names[i] == "make one iteration") {
        criteria += LIN_SOLVER_MAKE_ONE_ITERATION;
      } else {
	Errors::Message msg;
	msg << "LinearOperatorGMRES: \"convergence criteria\" type \"" << names[i] << "\" is not recognized.";
	Exceptions::amanzi_throw(msg);
      }
    }
  } else {
    criteria = LIN_SOLVER_RELATIVE_RHS;
  }
  set_criteria(criteria);

  // parameters for NKA
  nka_dim_ = plist.get<int>("max nka vectors", 10);
  nka_dim_ = std::min<int>(nka_dim_, max_itrs_);
  nka_tol_ = plist.get<double>("nka vector tolerance", 0.05);

  // NKA
  nka_ = Teuchos::rcp(new NKA_Base<Vector,VectorSpace>(nka_dim_, nka_tol_, m_->DomainMap()));
  nka_->Init(plist);

  initialized_ = true;
}

} // namespace
} // namespace

#endif
