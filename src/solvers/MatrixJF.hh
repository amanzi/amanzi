/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  This decorator class wraps a nonlinear SolverFnBase as a Matrix
  class that approximates the linear problem at a point by doing a
  standard finite difference (i.e. Knoll and Keyes '04) to approximate
  the action of the Jacobian.  When used with, for instance, NKA, this
  residual evaluation reduces to GMRES on the linear system, giving as
  a solution the standard JFNK correction to the nonlinear problem.
*/

#ifndef AMANZI_JF_MATRIX_HH_
#define AMANZI_JF_MATRIX_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class MatrixJF {
 public:
  MatrixJF() {};  // default constructor for LinOp usage

  MatrixJF(Teuchos::ParameterList& plist,
           const Teuchos::RCP<SolverFnBase<Vector> > fn,
           const VectorSpace& map) :
      plist_(plist),
      fn_(fn) {
    map_ = Teuchos::rcp(new VectorSpace(map));
    r0_ = Teuchos::rcp(new Vector(*map_));
    // u0_ = Teuchos::rcp(new Vector(*map_));

    Init_();
  }

  // Space for the domain of the operator.
  const VectorSpace& DomainMap() const { return *map_; }

  // Space for the range of the operator.
  const VectorSpace& RangeMap() const { return *map_; }

  // Apply matrix, b <-- Ax, returns ierr = 0 if success, !0 otherwise
  int Apply(const Vector& x, Vector& b) const;

  // Apply the inverse, x <-- A^-1 b, returns ierr = 0 if success, !0 otherwise
  int ApplyInverse(const Vector& b, Vector& x) const;

  void set_linearization_point(const Teuchos::RCP<const Vector>& u);

 protected:
  double CalculateEpsilon_(const Vector& u, const Vector& x) const;
  void Init_();

 protected:
  Teuchos::RCP<const VectorSpace> map_;
  Teuchos::RCP<SolverFnBase<Vector> > fn_;
  Teuchos::ParameterList plist_;
  Teuchos::RCP<Vector> u0_;
  Teuchos::RCP<Vector> r0_;

  double eps_;
  std::string method_name_;
};


// Forward (Apply) operator
template<class Vector, class VectorSpace>
void MatrixJF<Vector,VectorSpace>::Init_() {
  eps_ = plist_.get<double>("finite difference epsilon", 1.0e-8);
  method_name_ = plist_.get<std::string>("method for epsilon", "Knoll-Keyes");
}


// Forward (Apply) operator
template<class Vector, class VectorSpace>
int MatrixJF<Vector,VectorSpace>::Apply(const Vector& x, Vector& b) const {
  Teuchos::RCP<Vector> r1 = Teuchos::rcp(new Vector(*map_));
 
  double eps = CalculateEpsilon_(*u0_, x);

  // std::cout<<"x \n";
  // x.Print(std::cout);
  // std::cout<<"u0\n";
  // u0_->Print(std::cout);

  // evaluate r1 = f(u0 + eps*x)
  u0_->Update(eps, x, 1.0);
  // std::cout<<"u1\n";
  // u0_->Print(std::cout);
  fn_->ChangedSolution();
  fn_->Residual(u0_, r1);


  // std::cout<<"r0\n";
  // r0_->Print(std::cout);
  // std::cout<<"r1\n";
  // r1->Print(std::cout);


  // evaluate Jx = (r1 - r0) / eps
  b = *r1;
  b.Update(-1.0/eps, *r0_, 1.0/eps);

  // std::cout<<"b\n";
  // b.Print(std::cout);

  // revert to old u0
  u0_->Update(-eps, x, 1.0);
  fn_->ChangedSolution();

  return 0;
}


// Forward (Apply) operator
template<class Vector, class VectorSpace>
int MatrixJF<Vector,VectorSpace>::ApplyInverse(const Vector& b, Vector& x) const {

  //std::cout<<"ApplyInverse\n";
  // ugliness in interfaces...
  Teuchos::RCP<const Vector> b_ptr = Teuchos::rcpFromRef(b);
  Teuchos::RCP<Vector> x_ptr = Teuchos::rcpFromRef(x);
  fn_->ApplyPreconditioner(b_ptr, x_ptr);
  return 0;
}


template<class Vector, class VectorSpace>
void MatrixJF<Vector,VectorSpace>::set_linearization_point(const Teuchos::RCP<const Vector>& u) {
  u0_ = Teuchos::rcp_const_cast<Vector>(u);
  fn_->Residual(u0_, r0_);
  //std::cout << "res( lin point ):" << std::endl;
  //r0_->Print(std::cout);
}


template<class Vector, class VectorSpace>
double MatrixJF<Vector,VectorSpace>::CalculateEpsilon_(const Vector& u, const Vector& x) const {
  // simple algorithm eqn 14 from Knoll and Keyes

  double eps = eps_;
  double typical_size_u = 100.;

  if (method_name_ == "Knoll-Keyes") {
    double xinf(0.0);
    x.NormInf(&xinf);
    
    if (xinf > 0) {
      double uinf(0.0);
      u.NormInf(&uinf);      
      eps = std::sqrt((1 + uinf) * 1.0e-12) / xinf;
    }
  }
  else if (method_name_ == "Knoll-Keyes L2") {
    double x_l2(0.0);
    x.Norm2(&x_l2);
    
    if (x_l2 > 0) {
      double u_l2(0.0);
      u.Norm2(&u_l2);      
      eps = std::sqrt((1 + u_l2) * 1.0e-12) / x_l2;
    }
  }
  else if (method_name_ == "Brown-Saad") {
    double x_l2(0.0);
    x.Norm2(&x_l2);
    double xinf(0.0);
    x.NormInf(&xinf);

    if (x_l2 > 0) {
      double alp(0.);
      u.Dot(x, &alp);
      // double sgn = (alp > 0) ? 1 : -1; 
      eps = (1e-12/x_l2)*std::max(fabs(alp), typical_size_u*xinf);
    }
  }              
  // else if (method_name_ == "simple") {
  //   double x_l2(0.0);
  //   x.Norm2(&x_l2);

  //   if (x_l2 > 0) {
  //     double b = 1e-13;
  //     double u_l1(0.);
  //     u.Norm1(&u_l1);

  //     int num = u.Size();
  //     //std::cout<<u_l1<<" "<<x_l2<<" "<<num<<"\n";    
  //     eps = (b*u_l1)/(num*x_l2) + b;
  //   }
  // }

  //eps = eps_;
  //eps = 1e-10;
  //std::cout<<"eps "<<eps<<"\n";
    
  return eps;
}

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif

