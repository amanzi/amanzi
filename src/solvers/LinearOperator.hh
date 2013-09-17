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

 
namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class LinearOperator : public Matrix {
 public:
  LinearOperator(Teuchos::RCP<const Matrix> m, Teuchos::RCP<const Matrix> h) : 
      m_(m), h_(h) {};
  ~LinearOperator() {};

  virtual void Init(Teuchos::ParameterList& plist) = 0;  

  void Apply(const Vector& v, Vector& mv) { m_->Apply(v, mv); }
  virtual int ApplyInverse(const Vector& v, Vector& hv) const = 0;

  Teuchos::RCP<const VectorSpace> domain() const { return m_->domain(); }
  Teuchos::RCP<const VectorSpace> range() const { return m_->range(); }

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
  std::string& name() { return name_; }

 protected:
  Teuchos::RCP<const Matrix> m_;
  Teuchos::RCP<const Matrix> h_;
  std::string name_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi
 
#endif
               

