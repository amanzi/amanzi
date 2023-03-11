/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Flow PK

  Model evaluator templated on a base class.
*/

#ifndef AMANZI_FLOW_MODEL_EVALUATOR_HH_
#define AMANZI_FLOW_MODEL_EVALUATOR_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

#include "ModelPartition.hh"

namespace Amanzi {
namespace Flow {

template <typename Model>
class ModelEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit ModelEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist) { Init_(); }
  ModelEvaluator(const ModelEvaluator<Model>& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new ModelEvaluator(*this));
  }

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override{};

 private:
  void Init_();

 private:
  Teuchos::RCP<ModelPartition<Model>> model_;
  std::vector<std::string> deps_;

  static Utils::RegisteredFactory<Evaluator, ModelEvaluator<Model>> reg_;
};


/* ******************************************************************
* Required member function
****************************************************************** */
template <typename Model>
void ModelEvaluator<Model>::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  if (!model_.get()) {
    auto mesh = results[0]->Map().Mesh();
    std::string model = plist_.get<std::string>("evaluator type");
    model_ = CreateModelPartition<Model>(mesh, plist_.sublist(model + " parameters"));
  }

  auto& data = *results[0]->ViewComponent("cell");
  int ncells = data.MyLength();

  if (deps_.size() > 0) {
    const auto& dep0 = *S.Get<CompositeVector>(deps_[0]).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) {
      data[0][c] = model_->second[(*model_->first)[c]]->Value(dep0[0][c]);
    }
  } else {
    for (int c = 0; c != ncells; ++c) {
      data[0][c] = model_->second[(*model_->first)[c]]->Value(0.0);
    }
  }
}


/* ******************************************************************
* Initialization.
****************************************************************** */
template <typename Model>
void ModelEvaluator<Model>::Init_()
{
  if (my_keys_.size() == 0) AMANZI_ASSERT(false);

  if (plist_.isParameter("dependencies")) {
    deps_ = plist_.get<Teuchos::Array<std::string>>("dependencies").toVector();
  }

  for (auto key : deps_) {
    dependencies_.insert(std::make_pair(key, Tags::DEFAULT));
  }
}

} // namespace Flow
} // namespace Amanzi

#endif
