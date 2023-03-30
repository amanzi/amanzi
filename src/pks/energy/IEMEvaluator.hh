/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  The internal energy model (IEM) evaluator simply calls the
  IEM with the correct arguments.
*/

#ifndef AMANZI_ENERGY_IEM_EVALUATOR_HH_
#define AMANZI_ENERGY_IEM_EVALUATOR_HH_

// Amanzi
#include "Key.hh"
#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

// Amanzi::Energy
#include "IEM.hh"

namespace Amanzi {
namespace Energy {

typedef std::vector<Teuchos::RCP<IEM>> IEMList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, IEMList> IEMPartition;

class IEMEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  explicit IEMEvaluator(Teuchos::ParameterList& plist);
  IEMEvaluator(const IEMEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  Teuchos::RCP<IEMPartition> iem_partition() { return iem_; }

  double EvaluateFieldSingle(int c, double T, double p);

 protected:
  void InitializeFromPlist_();

 private:
  void CreateIEMPartition_(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                           const Teuchos::ParameterList& plist);

  AmanziMesh::Entity_ID MyModel_(AmanziMesh::Entity_kind kind, AmanziMesh::Entity_ID id);

 protected:
  Key domain_;
  Key temperature_key_, pressure_key_;
  Tag tag_;
  Teuchos::RCP<IEMPartition> iem_;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  static Utils::RegisteredFactory<Evaluator, IEMEvaluator> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
