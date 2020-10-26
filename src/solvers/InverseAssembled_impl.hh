/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Base class for providing applyInverse() using assembled methods.

#pragma once

#include "SuperMap.hh"
#include "InverseHelpers.hh"
#include "Inverse.hh"
#include "InverseFactory.hh"

namespace Amanzi {
namespace AmanziSolvers {

namespace Impl {

//
// NOTE: In each case, there are two implementations -- one for our
// Amanzi-based vector/matrix implementations (e.g. Composite and Tree) and one
// for the default Trilinos ones.
//
// The latter is only used in testing.
//

//
// Move data from the vector to a flat vector.
//
template<class Vector>
void copyToSuperVector(const Teuchos::RCP<const Operators::SuperMap>& smap,
                       const Vector& v,
                       Teuchos::RCP<Vector_type>& sv) {
  AMANZI_ASSERT(smap.get());
  AMANZI_ASSERT(sv.get());
  Operators::copyToSuperVector(*smap, v, *sv);
}
template<>
void inline
copyToSuperVector<Vector_type>(const Teuchos::RCP<const Operators::SuperMap>& smap,
        const Vector_type& v,
        Teuchos::RCP<Vector_type>& sv) {
  sv = Teuchos::rcp(new Vector_type(v));
  *sv = v;
}

//
// Move data from the flat vector to the vector
//
template<class Vector>
void copyFromSuperVector(const Teuchos::RCP<const Operators::SuperMap>& smap,
                         const Vector_type& sv,
                         Vector& v) {
  AMANZI_ASSERT(smap.get());
  copyFromSuperVector(*smap, sv, v);
}
template<>
void inline
copyFromSuperVector<Vector_type>(const Teuchos::RCP<const Operators::SuperMap>& smap,
        const Vector_type& sv,
        Vector_type& v) {
  AMANZI_ASSERT(&sv == &v);
}

//
// Get a SuperMap if possible/needed.
//
template<class Operator>
typename std::enable_if<is_assembling<Operator>::value,
                        std::tuple<Teuchos::RCP<const Operators::SuperMap>,
                                   Teuchos::RCP<Vector_type>,
                                   Teuchos::RCP<Vector_type>>>::type
getSuperMap(Operator& m) {
  auto smap = m.getSuperMap();
  auto x = Teuchos::rcp(new Vector_type(smap->getMap()));
  auto y = Teuchos::rcp(new Vector_type(smap->getMap()));
  return std::make_tuple(smap, x, y);
}


template<class Operator>
typename std::enable_if<is_assembled<Operator>::value,
                        std::tuple<Teuchos::RCP<Operators::SuperMap>,
                                   Teuchos::RCP<Vector_type>,
                                   Teuchos::RCP<Vector_type>>>::type
getSuperMap(Operator& m) {
  auto x = Teuchos::rcp(new Vector_type(m.DomainMap()));
  auto y = Teuchos::rcp(new Vector_type(m.RangeMap()));
  return std::make_tuple(Teuchos::null, x, y);
}

//
// Assemble the matrix and get it.
//
template<class Operator>
typename std::enable_if<is_assembling<Operator>::value,
                        Teuchos::RCP<Matrix_type>>::type
assembleMatrix(const Teuchos::RCP<Operator>& m) {
  m->AssembleMatrix();
  return m->A();
}

template<class Operator>
typename std::enable_if<is_assembled<Operator>::value,
                        Teuchos::RCP<Matrix_type>>::type
assembleMatrix(const Teuchos::RCP<Operator>& m) {
  // already assembled matrix
  AMANZI_ASSERT(m->Filled());
  return m;
}

//
// Just get a matrix (may be unassembled)
//
template<class Operator>
typename std::enable_if<is_assembling<Operator>::value,
                        Teuchos::RCP<Matrix_type>>::type
getMatrix(const Teuchos::RCP<Operator>& h) {
  h->SymbolicAssembleMatrix();
  return h->A();
}

template<class Operator>
typename std::enable_if<is_assembled<Operator>::value,
                        Teuchos::RCP<Matrix_type>>::type
getMatrix(const Teuchos::RCP<Operator>& h) {
  return h;
}

} // namespace Impl


//
// Class for assembled inverse methods.
//

template<class Operator,
         class Preconditioner,
         class Vector,
         class VectorSpace>
void InverseAssembled<Operator,Preconditioner,Vector,VectorSpace>::set_inverse_parameters(
    Teuchos::ParameterList& plist)
{
  solver_ = createAssembledMethod<>(method_name_, plist);
}


template<class Operator,
         class Preconditioner,
         class Vector,
         class VectorSpace>
void InverseAssembled<Operator,Preconditioner,Vector,VectorSpace>::initializeInverse()
{
  AMANZI_ASSERT(h_.get()); // set_matrices was called
  AMANZI_ASSERT(solver_.get()); // set_inverse_parameters was called

  // note, this step is critical if the solver needs to destroy things
  // using the old matrix.
  solver_->set_matrix(Teuchos::null);

  Teuchos::RCP<Matrix_type> hA = Impl::getMatrix(h_);
  solver_->set_matrix(hA);
  solver_->initializeInverse();

  if (!updated_) {
    std::tie(smap_,Y_,X_) = Impl::getSuperMap(*m_);
    updated_ = true;
    computed_once_ = false;
  }
}


template<class Operator,
         class Preconditioner,
         class Vector,
         class VectorSpace>
void InverseAssembled<Operator,Preconditioner,Vector,VectorSpace>::computeInverse()
{
  AMANZI_ASSERT(updated_); // note, we could call Update here, but would prefer
                           // that developers to the right thing.
  Impl::assembleMatrix(h_);
  solver_->computeInverse();
  computed_once_ = true;
}


template<class Operator,
         class Preconditioner,
         class Vector,
         class VectorSpace>
int InverseAssembled<Operator,Preconditioner,Vector,VectorSpace>::applyInverse(const Vector& X, Vector& Y) const
{
  AMANZI_ASSERT(updated_);
  AMANZI_ASSERT(computed_once_); // see above, could also call Compute here, but
                            // would prefer that developers do the right thing,
                            // especially since they might change values
                            // between calls and forget to call compute.
  Impl::copyToSuperVector(smap_, X, X_);
  int returned_code = solver_->applyInverse(*X_, *Y_);
  Impl::copyFromSuperVector(smap_, *Y_, Y);
  return returned_code;
}

} // namespace AmanziSolvers
} // namespace Amanzi
