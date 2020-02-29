/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for water content.

INTERFROST's comparison uses a very odd compressibility term that doesn't
quite fit into either compressible porosity or into a compressible density, so
it needs a special evaluator.

----------------------------------------------------------------------------- */


#ifndef AMANZI_INTERFROST_WATER_CONTENT_HH_
#define AMANZI_INTERFROST_WATER_CONTENT_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostWaterContent : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  InterfrostWaterContent(Teuchos::ParameterList& wc_plist);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  double beta_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,InterfrostWaterContent> reg_;
};

} // namespace
} // namespace
} // namespace

#endif
