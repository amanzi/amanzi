/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_PRECONDITIONER_IDENTITY_HH_
#define AMANZI_PRECONDITIONER_IDENTITY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template<class Matrix,class Vector>
class PreconditionerIdentity : public Preconditioner<Matrix,Vector> {
 public:
  PreconditionerIdentity() {};
  ~PreconditionerIdentity() {};

  void Init(const std::string& name, const Teuchos::ParameterList& list) override {};
  void Update(const Teuchos::RCP<const Matrix>& A) override {};
  void Destroy() override {};

  inline int ApplyInverse(const Vector& v, Vector& hv) const override;

  int returned_code() override { return 0; }
};


template<class Matrix,class Vector>
int
PreconditionerIdentity<Matrix,Vector>::ApplyInverse(const Vector& v, Vector& hv) const
{
  hv = v;
  return 0;
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi


// #ifdef HAVE_TPETRA_PRECONDITIONERS

#include "AmanziTypes.hh"
#include "AmanziVector.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template<>
int
PreconditionerIdentity<Matrix_type,Vector_type>::ApplyInverse(const Vector_type& v, Vector_type& hv) const
{
  hv.assign(v);
  return 0;
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi

// #endif // HAVE_TPETRA_PRECONDITIONERS


#endif
