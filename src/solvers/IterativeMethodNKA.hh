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

.. _iterative-method-nka-spec:
.. admonition:: iterative-method-nka-spec

    * `"error tolerance`" ``[double]`` **1.e-6** Tolerance on which to declare success.

    * `"maximum number of iterations`" ``[int]`` **100** Maximum iterations before declaring failure.

    * `"overflow tolerance`" ``[double]`` **3.e50** Error above this value results in failure.

    * `"convergence criterial`" ``[Array(string)]`` **{relative rhs}** A list of
      criteria, any of which can be applied.  Valid include:

      - `"relative rhs`" : measure error relative to the norm of the RHS vector
      - `"relative residual`" : measure error relative to the norm of the residual
      - `"absolute residual`" : measure error directly, norm of error
      - `"make one iteration`" : require at least one iteration to be performed before declaring success

    * `"max nka vectors`" ``[int]`` **10** Size of the NKA space used to span the residual, conceptually equivalent to the size of the Krylov space.

    * `"nka vector tolerance`" ``[double]`` **0.05** Vectors whose dot product are within this tolerance are considered parallel, and therefore the old vector is thrown out.


*/

#ifndef AMANZI_NKA_OPERATOR_HH_
#define AMANZI_NKA_OPERATOR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "VerboseObject.hh"

#include "InverseDefs.hh"
#include "Matrix.hh"
#include "NKA_Base.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Matrix,
          class Preconditioner = Matrix,
          class Vector = typename Matrix::Vector_t,
          class VectorSpace = typename Vector::VectorSpace_t>
class IterativeMethodNKA
  : public InverseIterativeMethod<Matrix, Preconditioner, Vector, VectorSpace> {
 private:
  using InvIt = InverseIterativeMethod<Matrix, Preconditioner, Vector, VectorSpace>;

 public:
  IterativeMethodNKA() : InvIt() {}

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override final;

  virtual int ApplyInverse(const Vector& v, Vector& hv) const override final
  {
    returned_code_ = NKA_(v, hv, tol_, max_itrs_, criteria_);
    if (returned_code_ <= 0) return 1;
    return 0;
  }

 protected:
  std::string MethodName_() const override final { return "NKA"; }

 private:
  int NKA_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;

  IterativeMethodNKA(const IterativeMethodNKA& other) = delete;

 private:
  using InvIt::m_;
  using InvIt::h_;
  using InvIt::vo_;
  using InvIt::num_itrs_;
  using InvIt::max_itrs_;
  using InvIt::tol_;
  using InvIt::residual_;
  using InvIt::returned_code_;
  using InvIt::CheckConvergence_;
  using InvIt::krylov_dim_;
  using InvIt::inited_;
  using InvIt::criteria_;
  using InvIt::rnorm0_;

  Teuchos::RCP<NKA_Base<Vector, VectorSpace>> nka_;
  int nka_dim_;
  double nka_tol_;
};


// Apply the inverse, x <-- A^-1 b
template <class Matrix, class Preconditioner, class Vector, class VectorSpace>
int
IterativeMethodNKA<Matrix, Preconditioner, Vector, VectorSpace>::NKA_(const Vector& f,
                                                                      Vector& x,
                                                                      double tol,
                                                                      int max_itrs,
                                                                      int criteria) const
{
  AMANZI_ASSERT(m_.get()); // set_matrices() called
  AMANZI_ASSERT(inited_);  // init called

  Teuchos::OSTab tab = vo_->getOSTab();

  nka_->Restart();

  residual_ = 0.0;
  num_itrs_ = 0;

  Teuchos::RCP<Vector> dx = Teuchos::rcp(new Vector(x));
  Teuchos::RCP<Vector> dxp = Teuchos::rcp(new Vector(x)); // preconditioned correction
  Teuchos::RCP<Vector> r = Teuchos::rcp(new Vector(x));

  double fnorm, xnorm;
  f.Norm2(&fnorm);
  x.Norm2(&xnorm);

  int ierr = m_->Apply(x, *r); // r = f - A * x
  AMANZI_ASSERT(!ierr);
  r->Update(1.0, f, -1.0);

  r->Norm2(&rnorm0_);
  residual_ = rnorm0_;

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << num_itrs_ << " ||r||=" << residual_ << std::endl;
  }

  if (!(criteria & LIN_SOLVER_MAKE_ONE_ITERATION)) {
    ierr = this->CheckConvergence_(rnorm0_, fnorm);
    if (ierr) return ierr;
  }

  bool done = false;
  while (!done) {
    ierr = h_->ApplyInverse(*r, *dxp);
    AMANZI_ASSERT(!ierr);

    nka_->Correction(*dxp, *dx);
    x.Update(1.0, *dx, 1.0);

    ierr = m_->Apply(x, *r); // r = f - A * x
    AMANZI_ASSERT(!ierr);
    r->Update(1.0, f, -1.0);

    r->Norm2(&residual_);

    num_itrs_++;

    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *vo_->os() << num_itrs_ << " ||r||=" << residual_ << std::endl;
    }
    ierr = this->CheckConvergence_(residual_, fnorm);
    if (ierr) return ierr;

    done = num_itrs_ > max_itrs;
  }

  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "Failed (" << num_itrs_ << " itrs) ||r|| = " << residual_ << std::endl;
  return LIN_SOLVER_MAX_ITERATIONS;
}


/* ******************************************************************
 * Initialization from a parameter list. Available parameters:
 * "error tolerance" [double] default = 1e-6
 * "maximum number of iterations" [int] default = 100
 * "convergence criteria" Array(string) default = "{relative rhs}"
 ****************************************************************** */
template <class Matrix, class Preconditioner, class Vector, class VectorSpace>
void
IterativeMethodNKA<Matrix, Preconditioner, Vector, VectorSpace>::set_inverse_parameters(
  Teuchos::ParameterList& plist)
{
  InvIt::set_inverse_parameters(plist);

  // parameters for NKA
  nka_dim_ = plist.get<int>("max nka vectors", krylov_dim_);
  nka_dim_ = std::min<int>(nka_dim_, max_itrs_);
  nka_tol_ = plist.get<double>("nka vector tolerance", 0.05);

  // NKA
  nka_ = Teuchos::rcp(new NKA_Base<Vector, VectorSpace>(nka_dim_, nka_tol_, m_->DomainMap()));
  nka_->Init(plist);

  inited_ = true;
}

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
