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
std::tuple<Teuchos::RCP<Operators::SuperMap>,
           Teuchos::RCP<Epetra_Vector>,
           Teuchos::RCP<Epetra_Vector>>
getSuperMap(Operator& m) {
  auto smap = m.getSuperMap();
  auto x = Teuchos::rcp(new Epetra_Vector(*smap->Map()));
  auto y = Teuchos::rcp(new Epetra_Vector(*smap->Map()));
  return std::make_tuple(smap, x, y);
}
template<>
std::tuple<Teuchos::RCP<Operators::SuperMap>,
           Teuchos::RCP<Epetra_Vector>,
           Teuchos::RCP<Epetra_Vector>>
inline
getSuperMap<Epetra_CrsMatrix>(Epetra_CrsMatrix& m) {
  auto x = Teuchos::rcp(new Epetra_Vector(m.DomainMap()));
  auto y = Teuchos::rcp(new Epetra_Vector(m.RangeMap()));
  return std::make_tuple(Teuchos::null, x, y);
}

//
// Assemble the matrix and get it.
//
template<class Operator>
Teuchos::RCP<Epetra_CrsMatrix> assembleMatrix(const Teuchos::RCP<Operator>& m) {
  m->AssembleMatrix();
  return m->A();
}
template<>
Teuchos::RCP<Epetra_CrsMatrix> assembleMatrix<Epetra_CrsMatrix>(
    const Teuchos::RCP<Epetra_CrsMatrix>& m) {
  // already assembled matrix
  AMANZI_ASSERT(m->Filled());
  return m;
}


} // namespace Impl


//
// Class for assembled inverse methods.
//
template<class Operator,
         class Preconditioner=Operator,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>
class InverseAssembled :
      public Inverse<Operator,Preconditioner,Vector,VectorSpace> {
 public:
  InverseAssembled() :
      updated_(false),
      Inverse<Operator,Preconditioner,Vector,VectorSpace>()
  {}

  virtual void InitInverse(Teuchos::ParameterList& plist) override final {
    plist_ = Teuchos::rcpFromRef(plist);
  }
  
  virtual void UpdateInverse() override final
  {
    AMANZI_ASSERT(m_.get() && h_.get()); // set_matrices was called
    if (!solver_.get()) {
      Teuchos::RCP<Epetra_CrsMatrix> hA = Impl::assembleMatrix(h_);
      //      solver_ = createInverse<Epetra_CrsMatrix,Epetra_CrsMatrix,Epetra_Vector,Epetra_Map>(*plist_, hA, hA);
    }
    solver_->UpdateInverse();

    if (!updated_) {
      std::tie(smap_,x_,y_) = Impl::getSuperMap(*m_);
      updated_ = true;
    }

  }

  virtual void ComputeInverse() override final
  {
    AMANZI_ASSERT(updated_);
    solver_->ComputeInverse();
  }

  virtual int ApplyInverse(const Vector& y, Vector& x) const override final
  {
    AMANZI_ASSERT(updated_);

    Impl::copyToSuperVector(smap_, y, y_);
    int returned_code = solver_->ApplyInverse(*y_, *x_);
    Impl::copyFromSuperVector(smap_, *x_, x);
    return returned_code;
  }

  virtual double residual() const override final {
    return solver_->residual();
  }

  virtual int num_itrs() const override final {
    return solver_->num_itrs();
  }

  virtual void add_criteria(int criteria) override final {
    return solver_->add_criteria(criteria);
  }

  virtual int returned_code() const override final {
    return solver_->returned_code();
  }

  virtual std::string returned_code_string() const override final {
    return solver_->returned_code_string();
  }

 protected:
  bool updated_;

  using Inverse<Operator,Preconditioner,Vector,VectorSpace>::m_;
  using Inverse<Operator,Preconditioner,Vector,VectorSpace>::h_;

  Teuchos::RCP<Inverse<Epetra_CrsMatrix,Epetra_CrsMatrix,Epetra_Vector,Epetra_Map>> solver_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<Operators::SuperMap> smap_;
  mutable Teuchos::RCP<Epetra_Vector> x_, y_;
};
  
    
} // namespace AmanziSolvers
} // namespace Amanzi
