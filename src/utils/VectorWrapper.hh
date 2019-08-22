/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Wraps Tpetra Vector objects as Epetra Vectors.

/*
  NOTE: This is not really intended for general use, but is helpful in
  testing.
*/

#ifndef VECTOR_WRAPPER_HH_
#define VECTOR_WRAPPER_HH_

#include <numeric>
#include <type_traits>
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

namespace Amanzi {

template<class Vector, class Scalar=typename Vector::scalar_type>
class VectorWrapper {
 public:
  VectorWrapper(const Teuchos::RCP<Vector>& v) :
      v_(v) {}

  template<class OtherVector>
  VectorWrapper(const OtherVector& other) :
      v_(Teuchos::rcp(new Vector(*other.get()))) {}

  const Teuchos::RCP<Vector>& get() { return v_; }

  template<class Dummy=Vector>
  typename std::enable_if<!std::is_const<Dummy>::value, const Teuchos::RCP<const Vector>&>::type
  get() const { return (Teuchos::RCP<const Vector>) v_; }

  template<class Dummy=Vector>
  typename std::enable_if<std::is_const<Dummy>::value, const Teuchos::RCP<const Vector>&>::type
  get() const { return v_; }
  
  int MyLength() const {
    return v_->getLocalLength();
  }

  int NumVectors() const {
    return v_->getNumVectors();
  }
  
  int Norm2(Scalar* norm) const {
    Teuchos::Array<Scalar> norm_locs(NumVectors());
    v_->norm2(norm_locs);
    Scalar tmp = std::accumulate(norm_locs.begin(), norm_locs.end(), (Scalar) 0, [](Scalar x, Scalar y) { return x + y*y; } );
    *norm = sqrt(tmp);
    return 0;
  }  

  template<class Dummy=Vector>
  typename std::enable_if<!std::is_const<Dummy>::value, int>::type
  PutScalar(Scalar a) { v_->putScalar(a); return 0; }


  template<class OtherVector, class Dummy=Vector>
  typename std::enable_if<!std::is_const<Dummy>::value, int>::type
  Abs(const OtherVector& other) { v_->abs(*other.get()); return 0; }

  template<class Dummy=Vector>
  typename std::enable_if<!std::is_const<Dummy>::value, int>::type
  Scale(Scalar a) { v_->scale(a); return 0; }

  template<class OtherVector, class Dummy=Vector>
  typename std::enable_if<!std::is_const<Dummy>::value, int>::type
  Dot(const OtherVector& other, Scalar* result) {
    Teuchos::Array<Scalar> intermediate_result(NumVectors());
    v_->dot(*other.get(), intermediate_result);
    *result = std::accumulate(intermediate_result.begin(), intermediate_result.end(), 0);
    return 0;
  }

  template<class OtherVector, class Dummy=Vector>
  typename std::enable_if<!std::is_const<Dummy>::value, int>::type
  Update(const Scalar& alpha, const OtherVector& other, const Scalar& beta) {
    Scalar norm(0);
    Norm2(&norm);
    std::cout << "  in Update: |this| = " << norm;
    norm = 0;
    other.Norm2(&norm);
    std::cout << "  |other| = " << norm;
    v_->update(alpha, *other.get(), beta);

    norm = 0;
    Norm2(&norm);
    std::cout << "  |res| = " << norm << std::endl;
    return 0;
  }

  template<class OtherVector1, class OtherVector2, class Dummy=Vector>
  typename std::enable_if<!std::is_const<Dummy>::value, int>::type
  Update(const Scalar& alpha1, const OtherVector1& other1,
         const Scalar& alpha2, const OtherVector2& other2,
         const Scalar& beta) {
    v_->update(alpha1, *other1.get(), alpha2, *other2.get(), beta);
    return 0;
  }

  template<class OtherVector1, class OtherVector2, class Dummy=Vector>
  typename std::enable_if<!std::is_const<Dummy>::value, int>::type
  Multiply(const Scalar& alpha, const OtherVector1& other1, const OtherVector2& other2, const Scalar& beta) {
    v_->elementWiseMultiply(alpha, *other1.get(), *other2.get(), beta);
    return 0;
  }
  
  template<class OtherVector1, class OtherVector2, class Dummy=Vector>
  typename std::enable_if<!std::is_const<Dummy>::value, int>::type
  ReciprocalMultiply(const Scalar& alpha, const OtherVector1& other1, const OtherVector2& other2, const Scalar& beta) {
    VectorWrapper<Vector> other1_copy(other1);
    other1_copy.get()->reciprocal(*other1.get());
    Multiply(alpha, other1_copy, other2, beta);
    return 0;
  }


  
 protected:
  Teuchos::RCP<Vector> v_;

};


} // namespace Amanzi

#endif
