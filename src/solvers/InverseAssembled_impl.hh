/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Base class for providing ApplyInverse() using assembled methods.

#pragma once

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

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
void copyToSuperVector(const Teuchos::RCP<Operators::SuperMap>& smap,
                       const Vector& v,
                       Teuchos::RCP<Epetra_Vector>& sv) {
  AMANZI_ASSERT(smap.get());
  AMANZI_ASSERT(sv.get());
  Operators::copyToSuperVector(*smap, v, *sv);
}
template<>
void inline
copyToSuperVector<Epetra_Vector>(const Teuchos::RCP<Operators::SuperMap>& smap,
        const Epetra_Vector& v,
        Teuchos::RCP<Epetra_Vector>& sv) {
  sv = Teuchos::rcp(new Epetra_Vector(v));
  *sv = v;
}

//
// Move data from the flat vector to the vector
//
template<class Vector>
void copyFromSuperVector(const Teuchos::RCP<Operators::SuperMap>& smap,
                         const Epetra_Vector& sv,
                         Vector& v) {
  AMANZI_ASSERT(smap.get());
  copyFromSuperVector(*smap, sv, v);
}
template<>
void inline
copyFromSuperVector<Epetra_Vector>(const Teuchos::RCP<Operators::SuperMap>& smap,
        const Epetra_Vector& sv,
        Epetra_Vector& v) {
  AMANZI_ASSERT(&sv == &v);
}

//
// Get a SuperMap if possible/needed.
//
template<class Operator>
typename std::enable_if<is_assembling<Operator>::value,
                        std::tuple<Teuchos::RCP<Operators::SuperMap>,
                                   Teuchos::RCP<Epetra_Vector>,
                                   Teuchos::RCP<Epetra_Vector>>>::type
getSuperMap(Operator& m) {
  auto smap = m.getSuperMap();
  auto x = Teuchos::rcp(new Epetra_Vector(*smap->Map()));
  auto y = Teuchos::rcp(new Epetra_Vector(*smap->Map()));
  return std::make_tuple(smap, x, y);
}


template<class Operator>
typename std::enable_if<is_assembled<Operator>::value,
                        std::tuple<Teuchos::RCP<Operators::SuperMap>,
                                   Teuchos::RCP<Epetra_Vector>,
                                   Teuchos::RCP<Epetra_Vector>>>::type
getSuperMap(Operator& m) {
  auto x = Teuchos::rcp(new Epetra_Vector(m.DomainMap()));
  auto y = Teuchos::rcp(new Epetra_Vector(m.RangeMap()));
  return std::make_tuple(Teuchos::null, x, y);
}

//
// Assemble the matrix and get it.
//
template<class Operator>
typename std::enable_if<is_assembling<Operator>::value,
                        Teuchos::RCP<Epetra_CrsMatrix>>::type
assembleMatrix(const Teuchos::RCP<Operator>& m) {
  m->AssembleMatrix();
  return m->A();
}

template<class Operator>
typename std::enable_if<is_assembled<Operator>::value,
                        Teuchos::RCP<Epetra_CrsMatrix>>::type
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
                        Teuchos::RCP<Epetra_CrsMatrix>>::type
getMatrix(const Teuchos::RCP<Operator>& h) {
  h->SymbolicAssembleMatrix();
  return h->A();
}

template<class Operator>
typename std::enable_if<is_assembled<Operator>::value,
                        Teuchos::RCP<Epetra_CrsMatrix>>::type
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
void InverseAssembled<Operator,Preconditioner,Vector,VectorSpace>::set_parameters(
    Teuchos::ParameterList& plist)
{
  solver_ = createAssembledMethod<>(method_name_, plist);
}


template<class Operator,
         class Preconditioner,
         class Vector,
         class VectorSpace>
void InverseAssembled<Operator,Preconditioner,Vector,VectorSpace>::UpdateInverse()
{
  AMANZI_ASSERT(h_.get()); // set_matrices was called
  AMANZI_ASSERT(solver_.get()); // set_parameters was called

  Teuchos::RCP<Epetra_CrsMatrix> hA = Impl::getMatrix(h_);
  solver_->set_matrix(hA);
  solver_->UpdateInverse();

  if (!updated_) {
    std::tie(smap_,x_,y_) = Impl::getSuperMap(*m_);
    updated_ = true;
    computed_once_ = false;
  }
}


template<class Operator,
         class Preconditioner,
         class Vector,
         class VectorSpace>
void InverseAssembled<Operator,Preconditioner,Vector,VectorSpace>::ComputeInverse()
{
  AMANZI_ASSERT(updated_); // note, we could call Update here, but would prefer
                           // that developers to the right thing.
  Impl::assembleMatrix(h_);
  solver_->ComputeInverse();
  computed_once_ = true;
}


template<class Operator,
         class Preconditioner,
         class Vector,
         class VectorSpace>
int InverseAssembled<Operator,Preconditioner,Vector,VectorSpace>::ApplyInverse(const Vector& y, Vector& x) const
{
  AMANZI_ASSERT(updated_);
  AMANZI_ASSERT(computed_once_); // see above, could also call Compute here, but
                            // would prefer that developers do the right thing,
                            // especially since they might change values
                            // between calls and forget to call compute.
  Impl::copyToSuperVector(smap_, y, y_);
  int returned_code = solver_->ApplyInverse(*y_, *x_);
  Impl::copyFromSuperVector(smap_, *x_, x);
  return returned_code;
}

    
} // namespace AmanziSolvers
} // namespace Amanzi
