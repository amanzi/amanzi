/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#include "Teuchos_RCP.hpp"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "InverseFactory.hh"
#include "PDE_Accumulation.hh"
#include "PDE_DiffusionFactory.hh"
#include "MFD3D_Diffusion.hh"
#include "nlfv.hh"
#include "Tensor.hh"

#include "TransportDefs.hh"
#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Calculate diffusion tensor and add it to the dispersion tensor.
******************************************************************* */
void
Transport_PK::CalculateDiffusionTensorEffective_(double mdl,
                                                 double mdg,
                                                 double kH,
                                                 const Epetra_MultiVector& porosity,
                                                 const Epetra_MultiVector& saturation)
{
  if (D_.size() == 0) {
    D_.resize(ncells_owned);
    for (int c = 0; c < ncells_owned; c++) D_[c].Init(dim, 1);
  }

  for (int mb = 0; mb < mat_properties_.size(); mb++) {
    Teuchos::RCP<MaterialProperties> spec = mat_properties_[mb];

    std::vector<AmanziMesh::Entity_ID> block;
    for (int r = 0; r < (spec->regions).size(); r++) {
      std::string region = (spec->regions)[r];
      mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

      AmanziMesh::Entity_ID_List::iterator c;
      for (c = block.begin(); c != block.end(); c++) {
        double sl = saturation[0][*c];
        D_[*c](0, 0) =
          (mdl * sl * spec->tau[0] + kH * mdg * (1.0 - sl) * spec->tau[1]) * porosity[0][*c];
      }
    }
  }
}


/* ******************************************************************
* One-phase solver based on effective diffusion
****************************************************************** */
void
Transport_PK::DiffusionSolverEffective(Epetra_MultiVector& tcc_next, double t_old, double t_new)
{
  double dt_MPC = t_new - t_old;

  const auto& wc = S_->Get<CompositeVector>(wc_key_, Tags::DEFAULT);
  const auto& sat_c =
    *S_->Get<CompositeVector>(saturation_liquid_key_, Tags::DEFAULT).ViewComponent("cell");

  Teuchos::ParameterList& op_list =
    tp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");

  // default boundary conditions (none inside domain and Neumann on its boundary)
  auto bc_dummy =
    Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  std::vector<int>& bc_model = bc_dummy->bc_model();
  std::vector<double>& bc_value = bc_dummy->bc_value();

  // create the Dispersion and Accumulation operators
  Operators::PDE_DiffusionFactory opfactory;
  Teuchos::RCP<Operators::PDE_Diffusion> op1 = opfactory.Create(op_list, mesh_, bc_dummy);
  op1->SetBCs(bc_dummy, bc_dummy);
  auto op = op1->global_operator();
  auto op2 = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op));

  // Create the preconditioner and solver.
  auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(dispersion_preconditioner,
                                                                *preconditioner_list_,
                                                                dispersion_solver,
                                                                *linear_solver_list_,
                                                                true);
  inv_list.setName(dispersion_preconditioner);
  op->set_inverse_parameters(inv_list);
  op->InitializeInverse();

  const CompositeVectorSpace& cvs = op1->global_operator()->DomainMap();
  CompositeVector sol(cvs);

  int il, ig, phase;
  double sl, total, mdl, mdg;

  // Disperse and diffuse aqueous components
  for (int i = 0; i < num_aqueous; i++) {
    ComputeBCs_(bc_model, bc_value, i);

    ig = num_aqueous + i;
    il = air_water_map_[i];
    for (int c = 0; c < ncells_owned; ++c) {
      sl = sat_c[0][c];
      total = tcc_next[il][c] * sl + tcc_next[ig][c] * (1.0 - sl);

      tcc_next[il][c] = total / sl;
      tcc_next[ig][c] = 0.0;
    }

    FindDiffusionValue(component_names_[il], &mdl, &phase);
    FindDiffusionValue(component_names_[ig], &mdg, &phase);
    CalculateDiffusionTensorEffective_(mdl, mdg, kH_[i], *transport_phi, sat_c);

    // set the initial guess
    auto& sol_cell = *sol.ViewComponent("cell");
    for (int c = 0; c < ncells_owned; c++) { sol_cell[0][c] = tcc_next[i][c]; }
    if (sol.HasComponent("face")) { sol.ViewComponent("face")->PutScalar(0.0); }

    op->Init();
    Teuchos::RCP<std::vector<WhetStone::Tensor>> Dptr = Teuchos::rcpFromRef(D_);
    op1->Setup(Dptr, Teuchos::null, Teuchos::null);
    op1->UpdateMatrices(Teuchos::null, Teuchos::null);

    // add accumulation term
    op2->AddAccumulationDelta(sol, wc, wc, dt_MPC, "cell");

    op1->ApplyBCs(true, true, true);
    op->ComputeInverse();

    CompositeVector& rhs = *op->rhs();
    int ierr = op->ApplyInverse(rhs, sol);

    if (ierr < 0) {
      Errors::Message msg("TransportExplicit_PK solver failed with message: \"");
      msg << op->returned_code_string() << "\"";
      Exceptions::amanzi_throw(msg);
    }

    for (int c = 0; c < ncells_owned; c++) { tcc_next[i][c] = sol_cell[0][c]; }
  }
}

} // namespace Transport
} // namespace Amanzi
