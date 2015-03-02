/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Reliative permeability as a function of capillary pressure, k=k(pc).
*/

#ifndef AMANZI_FLOW_REL_PERM_EVALUATOR_HH_
#define AMANZI_FLOW_REL_PERM_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"
#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class RelPermEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  RelPermEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrm);
  RelPermEvaluator(const RelPermEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  void InitializeFromPlist_();

  Teuchos::RCP<WRMPartition> wrm_;
  Key sat_key_;
  Key dens_key_;
  Key visc_key_;
  Key surf_rel_perm_key_;

  bool is_dens_visc_;
  bool is_surf_;
  Key surf_mesh_key_;
  
  double perm_scale_;
  double min_val_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
