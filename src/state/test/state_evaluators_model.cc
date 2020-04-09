/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/*
  State

  Authors: Ethan Coon
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "evaluator/EvaluatorSecondaryMonotype.hh"
#include "evaluator/EvaluatorIndependent.hh"
#include "evaluator/EvaluatorPrimary.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;


/* ******************************************************************
 * Equation A = 2*B
 ****************************************************************** */
template <class DeviceType>
class AModel {
 public:
  AModel(OutputVector_type<DeviceType> A, InputVector_type<DeviceType> B,
         Teuchos::ParameterList& plist)
    : A_(A), B_(B)
  {}

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    A_(i) = 2 * B_(i);
  }

  class dAdB {};
  KOKKOS_INLINE_FUNCTION void operator()(dAdB, const int i) const
  {
    A_(i) = 2.0;
  }

 private:
  OutputVector_type<DeviceType> A_;
  InputVector_type<DeviceType> B_;
};
