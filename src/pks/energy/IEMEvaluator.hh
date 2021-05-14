/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  The internal energy model (IEM) evaluator simply calls the
  IEM with the correct arguments.
*/

#ifndef AMANZI_ENERGY_IEM_EVALUATOR_HH_
#define AMANZI_ENERGY_IEM_EVALUATOR_HH_

// Amanzi
#include "Key.hh"
#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

// Amanzi::Energy
#include "IEM.hh"

namespace Amanzi {
namespace Energy {

typedef std::vector<Teuchos::RCP<IEM> > IEMList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, IEMList> IEMPartition;

class IEMEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  // constructor format for all derived classes
  explicit
  IEMEvaluator(Teuchos::ParameterList& plist);
  IEMEvaluator(const IEMEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& results);

  Teuchos::RCP<IEMPartition> iem_partition() { return iem_; }

  double EvaluateFieldSingle(int c, double T, double p);

 protected:
  void InitializeFromPlist_();

 private:
  void CreateIEMPartition_(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                           const Teuchos::ParameterList& plist);

  AmanziMesh::Entity_ID MyModel_(AmanziMesh::Entity_kind kind, AmanziMesh::Entity_ID id);

 protected:
  Key temperature_key_, pressure_key_;
  Key domain_;
  Teuchos::RCP<IEMPartition> iem_;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  static Utils::RegisteredFactory<FieldEvaluator,IEMEvaluator> factory_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
