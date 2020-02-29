/*
  Drainage rate.

  A simple model based on relaxation from current water content to a saturated water content.

          |
          | source
          V
         /   \
      I /     \
       V       |
  --Theta--    | T
       ^       |
       | D     | 
       V       V
  ----------------------

  This is the model for drainage D.  Note that this model allows for
  both drainage (positive, sink from the canopy/litter layer) and
  uptake (negative, source to the canopy/litter layer), the rewetting
  of litter by surface water.

  D = drainage - uptake.
  
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_DRAINAGE_EVALUATOR_HH_
#define AMANZI_RELATIONS_DRAINAGE_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class DrainageEvaluator : public SecondaryVariableFieldEvaluator {

 public:

  // constructor format for all derived classes
  explicit
  DrainageEvaluator(Teuchos::ParameterList& plist);

  DrainageEvaluator(const DrainageEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key ai_key_;
  Key wc_key_;
  Key pd_key_;
  Key source_key_;

  double tau_;
  double wc_sat_;
  double n_liq_;
  bool is_uptake_;

 private:
  static Amanzi::Utils::RegisteredFactory<FieldEvaluator,DrainageEvaluator> factory_;

};

} // namespace
} // namespace
} // namespace

#endif
