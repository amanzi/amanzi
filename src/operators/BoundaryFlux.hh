/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_BOUNDARY_FLUX_FUNC
#define AMANZI_BOUNDARY_FLUX_FUNC

#include "Brent.hh"


namespace Amanzi {

#define DEBUG_FLAG 0

/* ******************************************************************
* Nonlinear function F(x) = k(p0 - x) * [(T*dir) * (p - x) + g] - bc
* where k is the nonlinear function provided by the model
* and p0, T, dir, p, g, bc are constant parameters.
****************************************************************** */
template<class Model>
class BoundaryFluxFn {
 public:
  typedef double (Model::*NonlinFunc)(double p) const;

  BoundaryFluxFn(const double trans_f,
                 const double lambda,
                 double cell_p,
                 double bc_flux,
                 double g_flux,
                 int dir,
                 double patm,
                 Teuchos::RCP<const Model> model,
                 NonlinFunc nonlinfunc)
    : trans_f_(trans_f),
      lambda_(lambda),
      cell_p_(cell_p),
      bc_flux_(bc_flux),
      g_flux_(g_flux),
      dir_(dir),
      patm_(patm)
  {
    model_ = model;
    nonlinfunc_ = nonlinfunc;
  }

  double operator()(double face_p) const
  {
    lambda_ = face_p;
    double krel = ((*model_).*nonlinfunc_)(patm_ - lambda_);
    double flux = dir_ * trans_f_ * krel * (cell_p_ - lambda_);

    return flux + g_flux_ * krel - bc_flux_;
  }

 protected:
  double trans_f_;
  mutable double lambda_;
  int face_index_;
  double face_Mff_;
  double cell_p_;
  double bc_flux_;
  double g_flux_;
  int dir_;
  double patm_;
  Teuchos::RCP<const Model> model_;
  NonlinFunc nonlinfunc_;
};


/* ******************************************************************
* Bisection solver based on toms748 algorithm.
****************************************************************** */
template<class Model>
class BoundaryFaceSolver {
 public:
  typedef double (Model::*NonlinFunc)(double p) const;

  BoundaryFaceSolver(double trans_f,
                     double g_f,
                     double cell_val,
                     double lambda,
                     double bnd_flux,
                     int dir,
                     double patm,
                     double min_val,
                     double max_val,
                     double eps,
                     Teuchos::RCP<const Model> model,
                     NonlinFunc test_fun)
    : lambda_(lambda),
      cell_val_(cell_val),
      trans_f_(trans_f),
      g_f_(g_f),
      patm_(patm),
      bnd_flux_(bnd_flux),
      min_val_(min_val),
      max_val_(max_val),
      eps_(eps)
  {
    func_ = Teuchos::rcp(new BoundaryFluxFn<Model>(
      trans_f, lambda, cell_val, bnd_flux, g_f, dir, patm, model, test_fun));
  }

  double SolveBisection(double face_val,
                        double min_val,
                        double max_val,
                        BoundaryFluxFn<Model>& func,
                        int max_it,
                        int& actual_it);

  double FaceValue()
  {
    int max_it = 100;
    int actual_it(max_it);

    return SolveBisection(lambda_, min_val_, max_val_, *func_, max_it, actual_it);
  }

  double lambda_;                                     // initial guess
  double cell_val_, trans_f_, g_f_, patm_, bnd_flux_; // input data
  double min_val_, max_val_;                          // bracket values
  double eps_;
  Teuchos::RCP<BoundaryFluxFn<Model>> func_;
};


template<class Model>
double
BoundaryFaceSolver<Model>::SolveBisection(double face_val,
                                          double min_val,
                                          double max_val,
                                          BoundaryFluxFn<Model>& func,
                                          int max_it,
                                          int& actual_it)
{
  double res = func(face_val);
  double left(0.0), right(0.0), lres(0.0), rres(0.0);

#if DEBUG_FLAG
  std::cout << "STaRT interval: (" << min_val << " " << max_val << "),  fval=" << face_val << "\n";
#endif
  if (res > 0.0) {
    left = face_val;
    lres = res;
    right = std::max(face_val, max_val);
    rres = func(right);
#if DEBUG_FLAG
    std::cout << "---set right:      " << left << " (" << lres << "), " << right << " (" << rres
              << ")\n";
#endif
    while (rres > 0.0) {
      right += std::max(fabs(right), 100.0);
      rres = func(right);
#if DEBUG_FLAG
      std::cout << "   change right:   " << right << " (" << rres << ")\n";
#endif
    }
  } else {
    right = face_val;
    rres = res;
    left = std::min(min_val, face_val);
    lres = func(left);
#if DEBUG_FLAG
    std::cout << "---set left:       " << left << " (" << lres << "), " << right << " (" << rres
              << ")\n";
#endif
    while (lres < 0.0) {
      left -= std::max(fabs(left), 100.0);
      lres = func(left);
#if DEBUG_FLAG
      std::cout << "   change left:    " << left << " (" << lres << ")\n";
#endif
    }
  }
#if DEBUG_FLAG
  std::cout << "+++bracket (func): " << left << " (" << lres << "), " << right << " (" << rres
            << ")\n";
#endif

  if (std::fabs(right - left) <= eps_) {
    face_val = right;
    actual_it = 0;
  } else {
    actual_it = max_it;
    face_val = Utils::findRootBrent(func, left, right, eps_, &actual_it, eps_);
  }

#if DEBUG_FLAG
  std::cout << "   solution = " << face_val << " itrs = " << actual_it << "\n";
#endif
  return face_val;
}

} // namespace Amanzi

#endif
