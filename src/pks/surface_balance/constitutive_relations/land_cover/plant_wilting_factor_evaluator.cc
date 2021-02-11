/*
  The plant wilting factor evaluator is an algebraic evaluator of a given model.
Wilting factor.

Beta, or the water availability factor, or the plant wilting factor.

Beta =  (p_closed - p) / (p_closed - p_open)

where p is the capillary pressure or soil mafic potential, and closed
and open indicate the values at which stomates are fully open or fully
closed (the wilting point).

  
  Generated via evaluator_generator.
*/

#include "plant_wilting_factor_evaluator.hh"
#include "plant_wilting_factor_model.hh"

namespace Amanzi {
namespace LandCover {
namespace Relations {

// Constructor from ParameterList
PlantWiltingFactorEvaluator::PlantWiltingFactorEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("plant_wilting_factor parameters");
  for (auto p : sublist) {
    if (!sublist.isSublist(p.first)) {
      Errors::Message message("PlantWiltingFactorEvaluator: expected list of models.");
      Exceptions::amanzi_throw(message);
    }
    auto& model_plist = sublist.sublist(p.first);
    models_.emplace_back(std::make_pair(model_plist.get<std::string>("region"),
            Teuchos::rcp(new PlantWiltingFactorModel(model_plist))));
  }
  InitializeFromPlist_();
}

// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
PlantWiltingFactorEvaluator::Clone() const
{
  return Teuchos::rcp(new PlantWiltingFactorEvaluator(*this));
}


// Initialize by setting up dependencies
void
PlantWiltingFactorEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomainPrefix(my_key_);

  // - pull Keys from plist
  // dependency: capillary_pressure_gas_liq
  pc_key_ = plist_.get<std::string>("capillary pressure gas liq key",
          domain_name+"capillary_pressure_gas_liq");
  dependencies_.insert(pc_key_);
}


void
PlantWiltingFactorEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> pc = S->GetFieldData(pc_key_);

  for (auto region_model : models_) {
    const Epetra_MultiVector& pc_v = *pc->ViewComponent("cell", false);
    Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

    AmanziMesh::Entity_ID_List cells;
    pc->Mesh()->get_set_entities(region_model.first, AmanziMesh::CELL,
            AmanziMesh::Parallel_type::OWNED, &cells);
    for (auto c : cells) {
      result_v[0][c] = region_model.second->PlantWiltingFactor(pc_v[0][c]);
    }
  }
}


void
PlantWiltingFactorEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> pc = S->GetFieldData(pc_key_);

  if (wrt_key == pc_key_) {
  for (auto region_model : models_) {
    const Epetra_MultiVector& pc_v = *pc->ViewComponent("cell", false);
    Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

    AmanziMesh::Entity_ID_List cells;
    pc->Mesh()->get_set_entities(region_model.first, AmanziMesh::CELL,
            AmanziMesh::Parallel_type::OWNED, &cells);
    for (auto c : cells) {
      result_v[0][c] = region_model.second->DPlantWiltingFactorDCapillaryPressureGasLiq(pc_v[0][c]);
    }
  }
  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
