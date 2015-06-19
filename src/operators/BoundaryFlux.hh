/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatskiy (dasvyat@lanl.gov)
*/

#ifndef AMANZI_BOUNDARYFLUX_FUNC
#define AMANZI_BOUNDARYFLUX_FUNC

#include <boost/math/tools/roots.hpp>


namespace Amanzi {

  //#define DEBUG_FLAG 1


template <class Model>
class BoundaryFluxFn {
 public:

  typedef double(Model::*NonlinFunc)(double p) const;

  BoundaryFluxFn(const double trans_f,
                 const double lambda,
                 double cell_p,
                 double bc_flux,
                 double g_flux,
                 int dir,
                 double patm,
                 Teuchos::RCP<const Model> model,
                 NonlinFunc nonlinfunc):
      trans_f_(trans_f), lambda_(lambda),
      cell_p_(cell_p), bc_flux_(bc_flux), g_flux_(g_flux),
      dir_(dir), patm_(patm) {
    model_ = model;
    nonlinfunc_ = nonlinfunc;
  }

  double operator()(double face_p) {
    lambda_ = face_p;
    double Krel;
    Krel = ((*model_).*nonlinfunc_)(patm_ - lambda_);

    return flux_() + g_flux_*Krel - bc_flux_;
  }

 protected:
  double flux_() {
    double s = 0.;
    double Krel;
    Krel = ((*model_).*nonlinfunc_)(patm_ - lambda_);

    s = dir_ * trans_f_ * Krel * (cell_p_ - lambda_);
    return s;
  }

 protected:
  double trans_f_;
  double lambda_;
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

struct Tol_ {
  Tol_(double eps) : eps_(eps) {};
  bool operator()(const double& a, const double& b) const {
    return std::abs(a - b) <= eps_;
  }
  double eps_;
};


template <class Model>
class BoundaryFaceSolver {
 public:

  typedef double(Model::*NonlinFunc)(double p) const;

  BoundaryFaceSolver(double trans_f, double g_f,
                     double cell_val, double lambda,
                     double bnd_flux, int dir, double patm, 
                     double min_val, double max_val, double eps,
                     Teuchos::RCP<const Model> model, NonlinFunc test_fun):
    trans_f_(trans_f), 
    g_f_(g_f), 
    cell_val_(cell_val), 
    lambda_(lambda), 
    patm_(patm), 
    min_val_(min_val),
    max_val_(max_val),
    eps_(eps),
    bnd_flux_(bnd_flux)
  {

    

    func_ = Teuchos::rcp(new BoundaryFluxFn<Model>(trans_f, lambda, cell_val, bnd_flux,
                                                   g_f, dir, patm, model, test_fun) );
    // nonlinfunc_ = test_fun;
    // std::cout<<((*model).*nonlinfunc_)(0)<<"\n";;


  };

  double SolveBisection(double face_val, 
                        double min_val,
                        double max_val,
                        BoundaryFluxFn<Model>& func,
                        Tol_& tol,
                        boost::uintmax_t max_it,
                        boost::uintmax_t& actual_it);

  double FaceValue(){
    // double eps = std::max(1.e-4 * std::abs(bnd_flux_), 1.e-8);
    // eps = 1e-8;
    Tol_ tol(eps_);
    boost::uintmax_t max_it = 100;
    boost::uintmax_t actual_it(max_it);
    
    
    return SolveBisection(lambda_, min_val_, max_val_, *func_, tol, max_it, actual_it);
  };

  double trans_f_;
  double lambda_;
  double cell_val_;
  double g_f_;
  double patm_;
  double bnd_flux_;
  double min_val_;
  double max_val_;
  double eps_;
  Teuchos::RCP<BoundaryFluxFn<Model> > func_;
  

};

template <class Model>
double BoundaryFaceSolver<Model>::SolveBisection(double face_val, 
                                                 double min_val, 
                                                 double max_val,
                                                 BoundaryFluxFn<Model>& func,
                                                 Tol_& tol,
                                                 boost::uintmax_t max_it,
                                                 boost::uintmax_t& actual_it){

  double res = func(face_val);
  double left = 0.;
  double right = 0.;
  double lres = 0.;
  double rres = 0.;

  if (res > 0.) {
    left = face_val;
    lres = res;
    right = std::max(face_val, max_val);
    rres = func(right);
    while (rres > 0.) {
      right += fabs(right);
      rres = func(right);
    }
  }
  else {
    right = face_val;
    rres = res;
#if DEBUG_FLAG
    std::cout << "RIGHT = " << right << ", " << rres << std::endl;
#endif
    left = std::min(min_val, face_val);
    lres = func(left);
    while (lres < 0.) {
#if DEBUG_FLAG
      std::cout << "LEFT = " << left << ", " << lres << std::endl;
#endif
      left -= fabs(left);
      lres = func(left);
    }
  }
#if DEBUG_FLAG
  std::cout << "   bracket (res): " << left << " (" << lres << "), "
            << right << " (" << rres << ")" << std::endl;
#endif

  std::pair<double,double> result;

  result = boost::math::tools::toms748_solve(func, left, right, lres, rres, tol, actual_it);
  // if (actual_it >= max_it) {
  //   std::cout << " Failed to converged in " << actual_it << " steps." << std::endl;
  //   return 3;
  // }
	  
  face_val = (result.first + result.second) / 2.;

#if DEBUG_FLAG
  std::cout << "face_val = "<<face_val<<"\n";
#endif
  return face_val;


}



}  // namespace Amanzi

#endif
