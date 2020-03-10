/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "BCs.hh"
#include "Evaluator_BCs.hh"

namespace Amanzi {

Evaluator_BCs::Evaluator_BCs(Teuchos::ParameterList& plist)
  : EvaluatorSecondary(plist)
{
  AMANZI_ASSERT(my_keys_.size() == 1);
  for (auto sublist : plist.sublist("boundary functions")) {
    std::string sublist_name = Keys::cleanPListName(sublist.first);
    dependencies_.push_back(std::make_pair(sublist_name, my_keys_[0].second));
  }
}

void
Evaluator_BCs::EnsureCompatibility(State& S)
{
  auto& my_fac = S.Require<Operators::BCs, Operators::BCs_Factory>(
    my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);

  for (auto& dep : dependencies_) S.RequireEvaluator(dep.first, dep.second);

  // check plist for vis or checkpointing control
  EnsureCompatibility_Flags_(S);

  if (my_fac.mesh().get()) {
    for (const auto& dep : dependencies_) {
      auto& eval = S.RequireEvaluator(dep.first, dep.second);
      auto& fac = S.Require<Operators::BCs, Operators::BCs_Factory>(dep.first, dep.second);
      fac.set_mesh(my_fac.mesh());
      fac.set_entity_kind(my_fac.entity_kind());
      fac.set_type(my_fac.type());

      auto& bf_list = plist_.sublist("boundary functions").sublist(dep.first);
      fac.set_parameterlist(bf_list);
      eval.EnsureCompatibility(S);
    }
  }
}

void
Evaluator_BCs::Update_(State& S)
{
  auto& result = S.GetW<Operators::BCs>(
    my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);

  // overwrite with actual BCs
  {
    auto result_model = result.bc_model();
    auto result_val = result.bc_value();
    const AmanziMesh::Mesh* mesh = result.mesh().get();

    // initialize all boundary faces to Neumann, 0 flux
    Kokkos::parallel_for(
        "Evaluator_BCs::Init",
        result_model.extent(0),
        KOKKOS_LAMBDA(const int f) {
          AmanziMesh::Entity_ID_View cells;
          result.mesh()->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          if (cells.size() == 1) {
            result_model(f) = Operators::OPERATOR_BC_NEUMANN;
            result_val(f) = 0.0;
          }
        });

    // loop over dependencies and accumulate them
    for (const auto& dep : dependencies_) {
      const auto& i_bcs = S.Get<Operators::BCs>(dep.first, dep.second);
      auto i_model = i_bcs.bc_model();
      auto i_val = i_bcs.bc_value();

      Kokkos::parallel_for(
          "Evaluator_BCs::Combine",
          result_model.extent(0),
          KOKKOS_LAMBDA(const int& j) {
            if (i_model(j)) {
              result_model(j) = i_model(j);
              result_val(j) = i_val(j);
            }
          });
    }
  }
}


// Evaluator_BCsFunction::Evaluator_BCsFunction(Teuchos::ParameterList& plist) :
//     EvaluatorIndependent<BCs,BCs_Factory>(plist) {
//   std::string bc_model = plist.get<std::string>("boundary condition type");
//   if (bc_model == "Dirichlet") {
//     bc_model_ = Operators::OPERATOR_BC_DIRICHLET;
//   } else if (bc_model == "Neumann") {
//     bc_model_ = Operators::OPERATOR_BC_NEUMANN;
//   } else {
//     Errors::Message msg;
//     msg << "Invalid boundary condition type: \"" << bc_model << "\", must be Dirichlet or Neumann.";
//     throw(msg);
//   }
// }

// void
// Evaluator_BCsFunction::Update_(State& S) override
// {
//   auto& bc = S.GetW<BCs>(my_key_, my_tag_, my_key_);

//   if (!computed_once_) {
//     // Create the function.
//     AMANZI_ASSERT(plist_.isSublist("function"));
//     func_ = Functions::CreateCompositeVectorFunction(plist_.sublist("function"), bc.space());
//   }

//   // NOTE: EvaluatorIndependentFunctions own their own data.
//   time_ = S.time(my_tag_);
//   func_->Compute(time_, bc.value());
  
// }




} // namespace Amanzi
