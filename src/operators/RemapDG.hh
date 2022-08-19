/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The helper advection-based base class for various remap methods. It
  provides support of time integration and calculation of various static
  and dynamic geometric quantities. The actual time-step loop could be
  implemented differently by an application.

  The integration is performed on the pseudo-time interval from 0 to 1. 
  The remap velocity u is constant, but since the integration is performed
  in the reference coordinate system associated with mesh0, the transformed
  velocity v is the time-dependent quantity. We call it as co-velocity,
  v = C^t u where C is the matrix of co-factors for the Jacobian matrix J.
  Recall that C = det(J) J^{-T}. The co-velocity is reprsented using a
  space-time polynomial.

  Input parameter list describes operators, limiters, and mesh maps,
  see native spec for more detail. 

  The design breaks implemetation into generic core capabilities (a helper
  class) and templated time integrator.
*/

#ifndef AMANZI_OPERATOR_REMAP_DG_HH_
#define AMANZI_OPERATOR_REMAP_DG_HH_

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "dbc.hh"
#include "Explicit_TI_RK.hh"

// Amanzi::Operators
#include "RemapDG_Helper.hh"

namespace Amanzi {
namespace Operators {

template<class Vector>
class RemapDG : public Explicit_TI::fnBase<Vector>,
                public RemapDG_Helper {
 public:
  RemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
          const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
          Teuchos::ParameterList& plist)
    : RemapDG_Helper(mesh0, mesh1, plist) {};
  ~RemapDG() {};

  // main members required by the time integration class
  // -- calculate functional f(t, u) where u is the conservative quantity
  virtual void FunctionalTimeDerivative(double t, const Vector& u, Vector& f) override;

  // -- limit solution at all steps of the RK scheme
  virtual void ModifySolution(double t, Vector& u) override;

  // change between conservative and non-conservative variable
  void ConservativeToNonConservative(double t, const Vector& u, Vector& v);
  void NonConservativeToConservative(double t, const Vector& u, Vector& v);
};

}  // namespace Operators
}  // namespace Amanzi

#endif
