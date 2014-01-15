/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Implementation for the NullEnergy PK.  This PK simply provides a constant
   temperature, and is provided for testing with other PKs that depend upon an
   energy equation.  This could easily be provided by the state as an independent
   variable, but this is nice for testing the full hierarchy with a simple PK.

   Example usage:

   <ParameterList name="energy">
   <Parameter name="PK model" type="string" value="Constant Temperature"/>
   <Parameter name="Constant Temperature" type="double" value="290.0"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "test_snow_dist.hh"

namespace Amanzi {

void TestSnowDist::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);

  S->RequireField(key_, name_)->SetMesh(S->GetMesh("surface"))
          ->SetComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField("precipitation_snow")->SetMesh(S->GetMesh("surface"))
      ->AddComponent("cell", AmanziMesh::CELL,1);
  S->RequireFieldEvaluator("precipitation_snow");
};

// -- call your favorite
bool TestSnowDist::advance(double dt) {
  S_next_->GetFieldEvaluator("precipitation_snow")->HasFieldChanged(S_next_.ptr(), name_);
  S_next_->GetFieldData("snow_depth", name_)
      ->Update(dt, *S_next_->GetFieldData("precipitation_snow"), 1.);
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
};



} // namespace
