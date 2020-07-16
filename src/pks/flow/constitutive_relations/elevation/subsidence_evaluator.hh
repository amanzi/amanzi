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

  <ParameterList name="surface_star-subsidence">
    <Parameter name="field evaluator type" type="string" value="subsidence"/>
  </ParameterList>

*/

#ifndef AMANZI_FLOWRELATIONS_SUBSIDENCE_EVALUATOR_
#define AMANZI_FLOWRELATIONS_SUBSIDENCE_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class SubsidenceEvaluator : public SecondaryVariableFieldEvaluator {

public:
  explicit
  SubsidenceEvaluator(Teuchos::ParameterList& plist);

  SubsidenceEvaluator(const SubsidenceEvaluator& other);
  
  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new SubsidenceEvaluator(*this));
  }
  
protected:
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

protected:
  Key elev_key_, elev_init_key_;;
  
private:
  static Utils::RegisteredFactory<FieldEvaluator,SubsidenceEvaluator> reg_;
  
};
  
} //namespace
} //namespace

#endif
