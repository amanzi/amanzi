/*
 * NoxVector.hh
 *
 *  Created on: Mar 11, 2016
 *      Author: amklinv
 */

#ifndef SRC_SOLVERS_NOXVECTOR_HH_
#define SRC_SOLVERS_NOXVECTOR_HH_

#include <NOX_Abstract_Vector.H>

namespace Amanzi {

template<class Vector>
class NoxVector : public NOX::Abstract::Vector
{
public:
  NoxVector(const Teuchos::RCP<Vector>& vec) :
      vec_(vec) {}

  NoxVector(const NoxVector& other, NOX::CopyType type=NOX::DeepCopy);

  Teuchos::RCP<Vector> get_vector() { return vec_; }
  Teuchos::RCP<const Vector> get_vector() const { return vec_; }
  
  NOX::size_type length() const {
    return vec_->GlobalLength();
  }

  void print(std::ostream &stream) const {
    vec_->Print(stream);
  }
    
  NoxVector& init(double gamma) {
    vec_->PutScalar(gamma);
    return *this;
  }

  NoxVector& random(bool useSeed=false, int seed=1) {
    vec_->Random();
    return *this;
  }
    
  NoxVector& abs(const NoxVector& y) {
    vec_->Abs(*vec_);
    return *this;
  }

  NoxVector& operator=(const NoxVector& y) {
    *vec_ = *y.vec_;
    return *this;
  }

  NoxVector& reciprocal(const NoxVector& y) {
    vec_->Reciprocal(*vec_);
    return *this;
  }
  
  NoxVector& scale(double gamma) {
    vec_->Scale(gamma);
    return *this;
  }

  NoxVector& scale(const NoxVector& a) {
    vec_->Multiply(1., *vec_, *a.vec_, 0.);
    return *this;
  }    

  NoxVector& update(double alpha, const NoxVector& a, double gamma=0.0) {
    vec_->Update(alpha, *a.vec_, gamma);
    return *this;
  }

  NoxVector& update(double alpha, const NoxVector& a, double beta,
                    const NoxVector& b, double gamma=0.0) {
    vec_->Update(alpha, *a.vec_, beta, *b.vec_, gamma);
    return *this;
  }

  Teuchos::RCP<NoxVector> clone(NOX::CopyType type=NOX::DeepCopy) const {
    return Teuchos::rcp(new NoxVector<Vector>(*this, type));
  }

  double norm(NOX::Abstract::Vector::NormType type=NOX::Abstract::Vector::TwoNorm) const {
    double result = 0.;
    switch (type) {
      case NOX::Abstract::Vector::TwoNorm:
        vec_->Norm2(&result);
        break;
      case NOX::Abstract::Vector::OneNorm:
        vec_->Norm1(&result);
        break;
      case NOX::Abstract::Vector::MaxNorm:
        vec_->NormInf(&result);
        break;
      default:
        ASSERT(0);
        break;
    }
    return result;    
  }

  double norm(const NoxVector &weights) const {
    Vector v_tmp(*vec_);
    v_tmp->Multiply(1., *weights.vec_, *vec_, 0.);
    double result;
    v_tmp->Dot(*vec_, &result);
    return std::sqrt(result);
  }

  double innerProduct(const NoxVector& y) const {
    double result;
    vec_->Dot(*y.vec_, &result);
    return result;
  }

 private:
  Teuchos::RCP<Vector> vec_;
};


template<Vector>
NoxVector::NoxVector<Vector>(const NoxVector& other,
                             NOX::CopyType type=NOX::DeepCopy) {
  switch (type) {
    case NOX::DeepCopy:
      vec_ = Teuchos::rcp(new Vector(*other.vec_, INIT_MODE_COPY));
      break;
    case NOX::ShapeCopy:
      vec_ = Teuchos::rcp(new Vector(*other.vec_, INIT_MODE_NONE));
      break;
    default:
      vec_ = Teuchos::rcp(new Vector(*other.vec_, INIT_MODE_NONE));
      break;
  }
}


template<Vector>
NoxVector::NoxVector<Vector>(const NoxVector& other,
                             NOX::CopyType type=NOX::DeepCopy) {
  switch (type) {
    case NOX::DeepCopy:
      vec_ = Teuchos::rcp(new Vector(*other.vec_, INIT_MODE_COPY));
    case NOX::ShapeCopy:
      vec_ = Teuchos::rcp(new Vector(*other.vec_, INIT_MODE_NONE));
    default:
      vec_ = Teuchos::rcp(new Vector(*other.vec_, INIT_MODE_NONE));
  }
}

template<>
NoxVector::NoxVector<Epetra_Vector>(const NoxVector& other,
        NOX::CopyType type=NOX::DeepCopy) :
  vec_(other.vec_) {}


} // namespace

#endif /* SRC_SOLVERS_NOXVECTOR_HH_ */
