/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "evaluator/EvaluatorSecondaryMonotype.hh"
#include "evaluator/EvaluatorPrimary.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

/*
  We will build the following dependencies tree:
    A -> {B, C, E, H}
    C -> {D, G}
    E -> {D, F}
    H -> F
    D -> G
    F -> G

  Primary fields are B=2 and G=3. The equations are
    A = 2*B + C*E*H = 6484
    C = 2*D + G     = 15
    E = D*F         = 36
    H = 2*F         = 12
    D = 2*G         = 6
    F = 2*G         = 6

  Derivatives are
    dA/dB = 2
    dA/dG = 8640

  WARNING: derivative of secondary field wrt to secondary field is
  not well defined. The code may throw an exception since
  intermediate derivatives are not saved.
*/

/* ******************************************************************
 * Equation A = 2*B + C*E*H
 ****************************************************************** */
// Device Type : GPUs or CPUs

template <class DeviceType>
class AModel {
 public:
  AModel(OutputVector_type<DeviceType> A, InputVector_type<DeviceType> B,
         InputVector_type<DeviceType> C, InputVector_type<DeviceType> E,
         InVector_type<DeviceType> H, Teuchos::ParameterList& plist)
    : A_(A), B_(B), C_(C), E_(E), H_(H)
  {}

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    A_(i) = 2 * B_(i) + C_(i) * E_(i) * H_(i);
  }

  class dAdB {};
  KOKKOS_INLINE_FUNCTION void operator()(dAdB, const int i) const
  {
    A_(i) = 2.0;
  }

  class dAdC {};
  KOKKOS_INLINE_FUNCTION void operator()(dAdC, const int i) const
  {
    A_(i) = E_(i) * H_(i);
  }
  class dAdE {};
  KOKKOS_INLINE_FUNCTION void operator()(dAdE, const int i) const
  {
    A_(i) = C_(i) * H_(i);
  }

  class dAdH {};
  KOKKOS_INLINE_FUNCTION void operator()(dAdH, const int i) const
  {
    A_(i) = C_(i) * E_(*);
  }

 private:
  OutputVector_type<DeviceType> A_;
  InputVector_type<DeviceType> B_;
  InputVector_type<DeviceType> C_;
  InputVector_type<DeviceType> E_;
  InputVector_type<DeviceType> H_;
};


/* ******************************************************************
 * Equation C = 2*D + G
 ****************************************************************** */


/* ******************************************************************
 * Equation D = 2*G
 ****************************************************************** */


/* ******************************************************************
 * Equation E = D*F
 ****************************************************************** */


/* ******************************************************************
 * Equation F = 2*G
 ****************************************************************** */


/* ******************************************************************
 * Equation H = 2*F
 ****************************************************************** */
