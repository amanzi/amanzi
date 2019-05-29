/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Rel perm( pc ( sat ) ).

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_REL_PERM_EVALUATOR_
#define AMANZI_FLOWRELATIONS_REL_PERM_EVALUATOR_

#include "wrm.hh"
#include "wrm_partition.hh"
#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class RelPermEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  RelPermEvaluator(Teuchos::ParameterList& plist);

  RelPermEvaluator(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<WRMPartition>& wrms);

  RelPermEvaluator(const RelPermEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }

 protected:
  
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  void InitializeFromPlist_();

  Teuchos::RCP<WRMPartition> wrms_;
  Key sat_key_;
  Key dens_key_;
  Key visc_key_;
  Key surf_rel_perm_key_;

  bool is_dens_visc_;
  bool is_surf_;
  Key surf_domain_;
  
  double perm_scale_;
  double min_val_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,RelPermEvaluator> factory_;
};

} //namespace
} //namespace

#endif
