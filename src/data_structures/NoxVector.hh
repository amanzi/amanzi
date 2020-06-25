/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/*
 * NoxVector.hh
 *
 *  Created on: Mar 11, 2016
 *      Author: amklinv
 */

#ifndef AMANZI_NOXVECTOR_HH_
#define AMANZI_NOXVECTOR_HH_

#include "NOX_Abstract_Vector.H"

namespace Amanzi {

template <class VectorClass>
class NoxVector : public NOX::Abstract::Vector {
 public:
  NoxVector(const Teuchos::RCP<VectorClass>& vec) : vec_(vec){};

  NoxVector(const NoxVector& other, NOX::CopyType type = NOX::DeepCopy);

  Teuchos::RCP<VectorClass> get_vector() { return vec_; }
  Teuchos::RCP<const VectorClass> get_vector() const { return vec_; }

  NOX::size_type length() const { return vec_->getGlobalLength(); }

  void print(std::ostream& stream) const { vec_->Print(stream); }

  NOX::Abstract::Vector& init(double gamma)
  {
    vec_->putScalar(gamma);
    return *this;
  }

  NOX::Abstract::Vector& random(bool useSeed = false, int seed = 1)
  {
    vec_->random();
    return *this;
  }

  NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y)
  {
    vec_->abs(*vec_);
    return *this;
  }

  NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y)
  {
    *vec_ = *(dynamic_cast<const NoxVector&>(y).vec_);
    return *this;
  }

  NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y)
  {
    vec_->reciprocal(*vec_);
    return *this;
  }

  NOX::Abstract::Vector& scale(double gamma)
  {
    vec_->scale(gamma);
    return *this;
  }

  NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a)
  {
    vec_->elementWiseMultiply(
      1., *vec_, *(dynamic_cast<const NoxVector&>(a).vec_), 0.);
    return *this;
  }

  NOX::Abstract::Vector&
  update(double alpha, const NOX::Abstract::Vector& a, double gamma = 0.0)
  {
    vec_->update(alpha, *(dynamic_cast<const NoxVector&>(a).vec_), gamma);
    return *this;
  }

  NOX::Abstract::Vector&
  update(double alpha, const NOX::Abstract::Vector& a, double beta,
         const NOX::Abstract::Vector& b, double gamma = 0.0)
  {
    vec_->update(alpha,
                 *(dynamic_cast<const NoxVector&>(a).vec_),
                 beta,
                 *(dynamic_cast<const NoxVector&>(b).vec_),
                 gamma);
    return *this;
  }

  Teuchos::RCP<NOX::Abstract::Vector>
  clone(NOX::CopyType type = NOX::DeepCopy) const
  {
    return Teuchos::rcp(new NoxVector<VectorClass>(*this, type));
  }

  double norm(
    NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm) const
  {
    double result = 0.;
    switch (type) {
    case NOX::Abstract::Vector::TwoNorm:
      result = vec_->norm2();
      break;
    case NOX::Abstract::Vector::OneNorm:
      result = vec_->norm1();
      break;
    case NOX::Abstract::Vector::MaxNorm:
      result = vec_->normInf();
      break;
    default:
      assert(0);
      break;
    }
    return result;
  }

  double norm(const NOX::Abstract::Vector& weights) const
  {
    VectorClass v_tmp(*vec_);
    v_tmp.elementWiseMultiply(
      1., *(dynamic_cast<const NoxVector&>(weights).vec_), *vec_, 0.);
    double result;
    result = v_tmp.dot(*vec_);
    return std::sqrt(result);
  }

  double innerProduct(const NOX::Abstract::Vector& y) const
  {
    double result;
    result = vec_->dot(*(dynamic_cast<const NoxVector&>(y).vec_));
    return result;
  }

 private:
  Teuchos::RCP<VectorClass> vec_;
};

} // namespace Amanzi

#endif
