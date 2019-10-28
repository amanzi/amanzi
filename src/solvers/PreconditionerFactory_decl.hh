/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziTypes.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template <class Matrix, class Vector>
class PreconditionerFactory {
 public:
  PreconditionerFactory(){};
  ~PreconditionerFactory(){};

  Teuchos::RCP<Preconditioner<Matrix, Vector>>
  Create(const std::string& name, const Teuchos::ParameterList& prec_list);

  Teuchos::RCP<Preconditioner<Matrix, Vector>>
  Create(Teuchos::ParameterList& prec_list);
};


//
// Specialization for Assembled matrix type, which is special in that it has more options!
template <>
class PreconditionerFactory<Matrix_type,Vector_type> {
 public:
  PreconditionerFactory(){};
  ~PreconditionerFactory(){};

  Teuchos::RCP<Preconditioner<Matrix_type,Vector_type>>
  Create(const std::string& name, const Teuchos::ParameterList& prec_list);

  Teuchos::RCP<Preconditioner<Matrix_type, Vector_type>>
  Create(Teuchos::ParameterList& prec_list);
};


} // namespace AmanziPreconditioners
} // namespace Amanzi
