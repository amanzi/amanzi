/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Field evaluator for specific molar enthalpy, h = u + p / rho.

.. code-block:: xml

  <ParameterList name="_ENERGY">  <!-- parent list -->
  <ParameterList name="enthalpy evaluator">
    <Parameter name="enthalpy key" type="string" value="enthalpy_liquid"/>
    <Parameter name="internal energy key" type="string" value="internal_energy_liquid"/>

    <Parameter name="include work term" type="bool" value="true"/>
    <Parameter name="pressure key" type="string" value="pressure"/>
    <Parameter name="molar density key" type="string" value="molar_density_liquid"/>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_ENTHALPY_EVALUATOR_HH_
#define AMANZI_ENTHALPY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class EnthalpyEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit EnthalpyEvaluator(Teuchos::ParameterList& plist);
  EnthalpyEvaluator(const EnthalpyEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  // Required for boundary conditions
  virtual double EvaluateFieldSingle(const Teuchos::Ptr<State>& S, int c, double T, double p);

 protected:
  Key pressure_key_, mol_density_key_, ie_key_;
  Tag tag_;
  bool include_work_, include_potential_;

 private:
  static Utils::RegisteredFactory<Evaluator, EnthalpyEvaluator> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
