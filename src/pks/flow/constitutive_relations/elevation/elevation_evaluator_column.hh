/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! ElevationEvaluatorColumn: evaluates the elevation (z-coordinate) and slope magnitude of a mesh.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
Evaluator type: `"elevation column`"

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h z|``

* `"elevation key`" ``[string]`` **elevation** Name the elevation variable. [m]
* `"slope magnitude key`" ``[string]`` **slope_magnitude** Name the elevation variable. [-]
* `"dynamic mesh`" ``[bool]`` **false** Lets the evaluator know that the elevation changes in time, and adds the `"deformation`" dependency.
* `"parent domain name`" ``[string]`` **DOMAIN** Domain name of the parent mesh, which is the 3D version of this domain.  Attempts to generate an intelligent default by stripping "surface" from this domain.

Example:

.. code-block:: xml

  <ParameterList name="elevation">
    <Parameter name="evaluator type" type="string" value="elevation column"/>
  </ParameterList>

*/

#ifndef AMANZI_FLOWRELATIONS_ELEVATION_EVALUATOR_COLUMN_
#define AMANZI_FLOWRELATIONS_ELEVATION_EVALUATOR_COLUMN_

#include "factory.hh"
#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

class ElevationEvaluatorColumn : public ElevationEvaluator {

 public:
  explicit
  ElevationEvaluatorColumn(Teuchos::ParameterList& plist);

  ElevationEvaluatorColumn(const ElevationEvaluatorColumn& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);

 private:
  static Utils::RegisteredFactory<FieldEvaluator,ElevationEvaluatorColumn> reg_;

  Key slope_key_, base_por_key_;
};

} //namespace
} //namespace

#endif
