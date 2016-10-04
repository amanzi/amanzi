/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the porosity, given a small compressibility of rock.

<<<<<<< HEAD
=======
  Compressible grains are both physically realistic (based on bulk modulus)
  and a simple way to provide a non-elliptic, diagonal term for helping
  solvers to converge.


>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_EVALUATOR_HH_
#define AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_EVALUATOR_HH_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "compressible_porosity_model_partition.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class CompressiblePorosityEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit
  CompressiblePorosityEvaluator(Teuchos::ParameterList& plist);
  CompressiblePorosityEvaluator(const CompressiblePorosityEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<CompressiblePorosityModelPartition> get_Models() { return models_; }

protected:
  Key poro_key_;
  Key pres_key_;

  Teuchos::RCP<CompressiblePorosityModelPartition> models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,CompressiblePorosityEvaluator> fac_;



};

} // namespace
} // namespace
} // namespace

#endif
