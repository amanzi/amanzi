/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Evaluates a function to provide Dirichlet data on faces.

/*
This is moving toward a function approach, but is horribly non-functional.
Effectively this just passes all the work the stored class, which includes
boht data and functions.  Eventually the functions will migrate out of the
data structure into this class...
 */


#ifndef AMANZI_EVALUATOR_BOUNDARY_FUNCTION_HH_
#define AMANZI_EVALUATOR_BOUNDARY_FUNCTION_HH_

#include "EvaluatorIndependent.hh"
#include "Patch.hh"
#include "MeshFunction.hh"
#include "OperatorDefs.hh"

namespace Amanzi {

class Evaluator_BCsFunction
    : public EvaluatorIndependent<MultiPatch,MultiPatchSpace> {
 public:
  explicit Evaluator_BCsFunction(Teuchos::ParameterList& plist)
      : EvaluatorIndependent<MultiPatch,MultiPatchSpace>(plist) {}

  Evaluator_BCsFunction(const Evaluator_BCsFunction& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new Evaluator_BCsFunction(*this));
  }

  virtual void EnsureCompatibility(State& S) override
  {
    auto& fac = S.Require<MultiPatch,MultiPatchSpace>(my_key_, my_tag_, my_key_);

    if (fac.mesh.get() && !inited_) {
      EvaluatorIndependent<MultiPatch,MultiPatchSpace>::EnsureCompatibility(S);

      auto speclist_funcs = Functions::processListWithFunction(
          plist_.sublist(plist_.get("list name", "function")),
          plist_.get("function name", "function"));

      std::string bc_type = plist_.get<std::string>("boundary condition type");
      if (bc_type == "Dirichlet") fac.flag_type = Operators::OPERATOR_BC_DIRICHLET;
      else if (bc_type == "Neumann") fac.flag_type = Operators::OPERATOR_BC_NEUMANN;
      else {
        Errors::Message msg;
        msg << "Invalid boundary condition type \"" << bc_type << "\" Currently supported are \"Dirichlet\" and \"Neumann\"";
        throw(msg);
      }
      
      // merge the lists -- need only those that match this entity_kind, and
      // on this mesh
      for (int i=0; i!=speclist_funcs.first.size(); ++i) {
        const auto& patch = speclist_funcs.first[i];
        if (patch.entity_kind == fac.flag_entity) {
          fac.AddPatch(patch.region, patch.entity_kind, patch.n_dofs);
          funcs_.push_back(speclist_funcs.second[i]);
        }
      }
      
      inited_ = true;
    }      
  }

  virtual std::string name() const override { return "boundary condition function"; }

 protected:
  // ---------------------------------------------------------------------------
  // Does the actual work to update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override
  {
    auto& mp = S.GetW<MultiPatch>(my_key_, my_tag_, my_key_);
    Functions::computeMeshFunction(funcs_, S.Get<double>("time", my_tag_), mp);
  }

 protected:
  std::vector<Teuchos::RCP<const MultiFunction>> funcs_;
  

 private:
  static Utils::RegisteredFactory<Evaluator, Evaluator_BCsFunction> fac_;
};

} // namespace Amanzi

#endif
