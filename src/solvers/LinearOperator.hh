/*
This is the Linear Solver component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
         Konstantin Lipnikov (lipnikov@lanl.gov)

Base class for linear solvers.
Usage: 
*/

#ifndef AMANZI_LINEAR_OPERATOR_HH_
#define AMANZI_LINEAR_OPERATOR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class LinearOperator : public Matrix {
 public:
  LinearOperator(const Teuchos::RCP<const Matrix>& m,
                 const Teuchos::RCP<const Matrix>& h) : 
      m_(m), h_(h) {};

  virtual ~LinearOperator() {};

  virtual void Init(Teuchos::ParameterList& plist) = 0;  

  virtual int Apply(const Vector& v, Vector& mv) const { return m_->Apply(v, mv); }
  virtual int ApplyInverse(const Vector& v, Vector& hv) const {
    int ierr = ApplyInverse_(v,hv);
    if (ierr > 0) return 0;
    else if (ierr < 0) return 1;
    else ASSERT(1);
  }

  virtual const VectorSpace& DomainMap() const { return m_->DomainMap(); }
  virtual const VectorSpace& RangeMap() const { return m_->RangeMap(); }

  double TrueResidual(const Vector& f, const Vector& v) const {
    Vector r(f);
    m_->Apply(v, r);  // r = f - M * x
    r.Update(1.0, f, -1.0);

    double true_residual;
    r.Norm2(&true_residual);
    return true_residual;
  }

  virtual double residual() = 0;
  virtual int num_itrs() = 0;
  virtual void add_criteria(int criteria) = 0;
  std::string name() { return name_; }
  void set_name(std::string name) { name_ = name; }

 protected:
  virtual int ApplyInverse_(const Vector& v, Vector& hv) const = 0;

 protected:
  Teuchos::RCP<const Matrix> m_;
  Teuchos::RCP<const Matrix> h_;
  std::string name_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi
 
#endif
               

