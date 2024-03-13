/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#include "Evaluator.hh"
#include "PDE_CouplingFlux.hh"
#include "State.hh"
#include "UniqueLocalIndex.hh"

#include "FractureInsertion_Helper.hh"

namespace Amanzi {

/* *******************************************************************
* Populate advective coupling fluxes
******************************************************************* */
void
UpdateEnthalpyCouplingFluxes(
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_matrix,
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_fracture,
  const std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux>>& adv_coupling_ops,
  bool flag)
{
  Key temperature_key_m("temperature"), temperature_key_f("fracture-temperature");
  Key enthalpy_key_m("enthalpy"), enthalpy_key_f("fracture-enthalpy");

  S->Get<CompositeVector>(enthalpy_key_m).ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>(temperature_key_m).ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("molar_flow_rate").ScatterMasterToGhosted("face");

  // extract enthalpy fields
  S->GetEvaluator(enthalpy_key_m).Update(*S, enthalpy_key_m);
  const auto& H_m = *S->Get<CompositeVector>(enthalpy_key_m).ViewComponent("cell", true);
  const auto& T_m = *S->Get<CompositeVector>(temperature_key_m).ViewComponent("cell", true);

  S->GetEvaluator(enthalpy_key_f).Update(*S, enthalpy_key_f);
  const auto& H_f = *S->Get<CompositeVector>(enthalpy_key_f).ViewComponent("cell", true);
  const auto& T_f = *S->Get<CompositeVector>(temperature_key_f).ViewComponent("cell", true);

  Teuchos::RCP<Epetra_MultiVector> dHdT_m, dHdT_f;
  if (!flag) {
    auto tmp1 = S->GetDerivative<CompositeVector>(
      enthalpy_key_m, Tags::DEFAULT, temperature_key_m, Tags::DEFAULT);
    tmp1.ScatterMasterToGhosted("cell");
    dHdT_m = tmp1.ViewComponent("cell");

    auto tmp2 = S->GetDerivative<CompositeVector>(
      enthalpy_key_f, Tags::DEFAULT, temperature_key_f, Tags::DEFAULT);
    tmp2.ScatterMasterToGhosted("cell");
    dHdT_f = tmp2.ViewComponent("cell");
  }

  // update coupling terms for advection
  int ncells_owned_f =
    mesh_fracture->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  auto values1 = std::make_shared<std::vector<double>>(2 * ncells_owned_f, 0.0);
  auto values2 = std::make_shared<std::vector<double>>(2 * ncells_owned_f, 0.0);

  int np(0), dir, shift;
  AmanziMesh::Entity_ID_List cells;
  const auto& flux = *S->Get<CompositeVector>("molar_flow_rate").ViewComponent("face", true);
  const auto& mmap = flux.Map();
  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    int first = mmap.FirstPointInElement(f);

    const auto& cells = mesh_matrix->getFaceCells(f);
    int ncells = cells.size();
    mesh_matrix->getFaceNormal(f, cells[0], &dir);
    shift = Operators::UniqueIndexFaceToCells(*mesh_matrix, f, cells[0]);

    for (int k = 0; k < ncells; ++k) {
      double tmp = flux[0][first + shift] * dir;

      // since we multiply by temperature, the model for the flux is
      // q (H / T) * T for both matrix and preconditioner
      if (tmp > 0) {
        int c1 = cells[k];
        double factor = (flag) ? H_m[0][c1] / T_m[0][c1] : (*dHdT_m)[0][c1];
        (*values1)[np] = tmp * factor;
      } else {
        double factor = (flag) ? H_f[0][c] / T_f[0][c] : (*dHdT_f)[0][c];
        (*values2)[np] = -tmp * factor;
      }

      dir = -dir;
      shift = 1 - shift;
      np++;
    }
  }

  // setup coupling operators with new data
  adv_coupling_ops[0]->Setup(values1, 1.0);
  adv_coupling_ops[0]->UpdateMatrices(Teuchos::null, Teuchos::null);

  adv_coupling_ops[1]->Setup(values2, -1.0);
  adv_coupling_ops[1]->UpdateMatrices(Teuchos::null, Teuchos::null);

  adv_coupling_ops[2]->Setup(values1, -1.0);
  adv_coupling_ops[2]->UpdateMatrices(Teuchos::null, Teuchos::null);

  adv_coupling_ops[3]->Setup(values2, 1.0);
  adv_coupling_ops[3]->UpdateMatrices(Teuchos::null, Teuchos::null);
}

} // namespace Amanzi
