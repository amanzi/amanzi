/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Markus Berndt

   Interface for the Prescribed Deformation PK.

   <ParameterList name="prescribed deformation">
      <Parameter name="PK model" type="string" value="Prescibed Mesh Deformation"/>
      <Parameter name="Deformation method" type="string" value="method name"/>
   </ParameterList>
   
   ------------------------------------------------------------------------- */

#include "prescribed_deformation.hh"

namespace Amanzi {
namespace Deform {

RegisteredPKFactory<PrescribedDeformation> PrescribedDeformation::reg_("prescribed deformation");


// -- Setup data
void PrescribedDeformation::setup(const Teuchos::Ptr<State>& S) {

}

// -- Initialize owned (dependent) variables.
void PrescribedDeformation::initialize(const Teuchos::Ptr<State>& S) {

}
  
// -- advance via one of a few methods
bool PrescribedDeformation::advance(double dt) {

}

} // namespace
} // namespace

