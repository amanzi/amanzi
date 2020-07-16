/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! ElevationEvaluatorColumn: evaluates the elevation (z-coordinate) and slope magnitude of a mesh.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

/*!
Evaluator type: `"subsidence`"

Evaluates the difference between the previous surface elevation and the new (subsidended) surface elevation.

* `"elevation key`" ``[string]`` **elevation** Name the elevation variable. [m]
* `"dynamic mesh`" ``[bool]`` **true** Lets the evaluator know that the elevation changes in time, and adds the `"deformation`".

Example:

.. code-block:: xml

  <ParameterList name="surface_star-elevation_init">
    <Parameter name="field evaluator type" type="string" value="initial elevation"/>
  </ParameterList>

*/

#ifndef AMANZI_FLOWRELATIONS_ELEVATION_INIT_EVALUATOR_
#define AMANZI_FLOWRELATIONS_ELEVATION_INIT_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class ElevationInitEvaluator : public SecondaryVariableFieldEvaluator {

public:
  explicit
  ElevationInitEvaluator(Teuchos::ParameterList& plist);

  ElevationInitEvaluator(const ElevationInitEvaluator& other);
  
  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new ElevationInitEvaluator(*this));
  }
  
protected:
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

protected:
  Key elev_key_;
  
private:
  static Utils::RegisteredFactory<FieldEvaluator,ElevationInitEvaluator> reg_;
  
};
  
} //namespace
} //namespace

#endif
