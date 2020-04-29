/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Data_Initializers are functions that know how to initialize data.

/*
  Initialize from parameter list.
*/

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "FunctionFactory.hh"

#include "CompositeVector.hh"
#include "CompositeVectorFunction.hh"
//#include "Op.hh"
//#include "TensorVector.hh"
//#include "Operator.hh"
//#include "BCs.hh"
//#include "BoundaryFunction.hh"

#include "Checkpoint.hh"

namespace Amanzi {
namespace Data_Initializers {

//
// Initialize
// ======================================================================
// default initialize does nothing
template <typename T>
bool
Initialize(Teuchos::ParameterList& plist,
           const Teuchos::ParameterList& attrs, T& t)
{
  return true;
}


//
// Helper initializers
//

// Initialize from a single parameter value.
template <typename T>
bool
InitializePrimitiveByValue(Teuchos::ParameterList& plist,
                           const Teuchos::ParameterList& attrs, T& t)
{
  if (plist.isParameter("value")) {
    t = plist.template get<T>("value");
    return true;
  }
  return false;
}

// Initialize from a single parameter value.
template <template <typename> class V, typename T>
bool
InitializeVectorByValue(Teuchos::ParameterList& plist,
                        const Teuchos::ParameterList& attrs, V<T>& t)
{
  if (plist.isParameter("value")) {
    t.putScalar(plist.template get<T>("value"));
    return true;
  }
  return false;
}


//
// Specializations for simple data types
// ======================================================================
template <>
bool
Initialize<double>(Teuchos::ParameterList& plist,
                   const Teuchos::ParameterList& attrs, double& t);
template <>
bool
Initialize<int>(Teuchos::ParameterList& plist,
                const Teuchos::ParameterList& attrs, int& t);
template <>
bool
Initialize<std::string>(Teuchos::ParameterList& plist,
                        const Teuchos::ParameterList& attrs, std::string& t);

//
// Specializations for Vectors
// ======================================================================
template <typename Scalar>
bool
Initialize(Teuchos::ParameterList& plist, const Teuchos::ParameterList& attrs,
           Vector_type_<Scalar>& t)
{
  if (InitializeVectorByValue<Vector_type_, Scalar>(plist, attrs, t))
    return true;
  return false;
}

template <typename Scalar>
bool
Initialize(Teuchos::ParameterList& plist, const Teuchos::ParameterList& attrs,
           MultiVector_type_<Scalar>& t)
{
  if (InitializeVectorByValue<MultiVector_type_, Scalar>(plist, attrs, t))
    return true;
  return false;
}

template <>
bool
Initialize(Teuchos::ParameterList& plist, const Teuchos::ParameterList& attrs,
           CompositeVector_<int>& t);


//
// Specializations for Vectors
// ======================================================================
template <>
bool Initialize<CompositeVector>(
    Teuchos::ParameterList &plist, const Teuchos::ParameterList& attrs,
    CompositeVector &t);

} // namespace Data_Initializers
} // namespace Amanzi

