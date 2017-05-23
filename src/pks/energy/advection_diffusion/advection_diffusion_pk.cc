/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "advection_factory.hh"
#include "energy_bc_factory.hh"
#include "CompositeVectorSpace.hh"

#include "OperatorAdvection.hh"
#include "OperatorDiffusion.hh"
#include "OperatorDiffusionFactory.hh"
#include "advection_diffusion.hh"

namespace Amanzi {
namespace Energy {

void AdvectionDiffusion::Setup(const Teuchos::Ptr<State>& S) {
  PK_PhysicalBDF_Default::Setup(S);

  // require fields
  Teuchos::RCP<CompositeVectorSpace> factory;

  // -- temperature -- on cells and faces
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  Teuchos::ParameterList mfd_plist = plist_->sublist("diffusion");
  std::string disc = mfd_plist.get<std::string>("discretization primary");
  bool tpfa = (disc == std::string("mfd: two-point flux approximation"))
      || (disc == std::string("fv: default"));

  factory = S->RequireField("temperature", "energy");
  factory->SetMesh(S->GetMesh());
  factory->SetGhosted(true);
  if (tpfa) {
    factory->SetComponent("cell", AmanziMesh::CELL, 1);
  } else {
    factory->SetComponents(names2, locations2, num_dofs2);
  }

  // -- thermal conductivity -- just cells
  factory = S->RequireField("thermal_conductivity", "energy");
  factory->SetMesh(S->GetMesh());
  factory->SetGhosted(true);
  factory->SetComponent("face", AmanziMesh::FACE, 1);

  // -- independent variables (not owned by this pk)
  factory = S->RequireField("mass_flux");
  factory->SetMesh(S->GetMesh());
  factory->SetGhosted(true);
  factory->SetComponent("face", AmanziMesh::FACE, 1);

  factory = S->RequireField("cell_volume");
  factory->SetMesh(S->GetMesh());
  factory->SetGhosted(true);
  factory->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("cell_volume");

  // boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(S->GetMesh(), bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_flux_ = bc_factory.CreateTotalFlux();

  int nfaces = S->GetMesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::OPERATOR_BC_NONE);
  bc_values_.resize(nfaces, 0.0);
  std::vector<double> mixed;

  bc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_markers_, bc_values_, mixed));

  // operator for the diffusion terms
  Operators::OperatorDiffusionFactory opfactory;
  AmanziGeometry::Point g;

  matrix_diff_ = opfactory.Create(mfd_plist, mesh_, bc_);
  // matrix_diff_->SetTensorCoefficient(Teuchos::null);

  // operator for advection terms
  Teuchos::ParameterList advect_plist = plist_->sublist("advection");
  matrix_adv_ = Teuchos::rcp(new Operators::OperatorAdvection(advect_plist, mesh_));
  
  // preconditioner
  Teuchos::ParameterList mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  preconditioner_diff_ = opfactory.Create(mfd_pc_plist, mesh_, bc_);
  // preconditioner_diff_->SetTensorCoefficient(Teuchos::null);
  // preconditioner_ = preconditioner_diff_->global_operator();
  // preconditioner_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, preconditioner_));

  // implicit_advection_ = plist_->get<bool>("implicit advection", false);
  // if (implicit_advection_) {
  //   preconditioner_adv_ = Teuchos::rcp(new Operators::OperatorAdvection(advect_plist, preconditioner_));
  // }

  // preconditioner_->SymbolicAssembleMatrix();
  // preconditioner_->InitPreconditioner("preconditioner", mfd_pc_plist);
  
  
};


// -- Initialize owned (dependent) variables.
void AdvectionDiffusion::Initialize(const Teuchos::Ptr<State>& S) {
  PK_PhysicalBDF_Default::Initialize(S);

  double time = S->time();
  bc_temperature_->Compute(time);
  bc_flux_->Compute(time);
  UpdateBoundaryConditions_();
};


// Evaluate BCs
void AdvectionDiffusion::UpdateBoundaryConditions_() {
  // mark all as null
  for (unsigned int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::OPERATOR_BC_NONE;
    bc_values_[n] = 0.0;
  }

  // Dirichlet bcs
  Functions::BoundaryFunction::Iterator bc;
  for (bc=bc_temperature_->begin(); bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  // flux bcs
  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_values_[f] = bc->second;
  }

  // mark all remaining boundary conditions as zero flux conditions
  AmanziMesh::Entity_ID_List cells;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_markers_[f] == Operators::OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_values_[f] = 0.0;
      }
    }
  }

};


// Add a boundary marker to owned faces.
void AdvectionDiffusion::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature) {
  Epetra_MultiVector& temp_f = *temperature->ViewComponent("face",true);
  int nfaces = temperature->size("face",false);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) {
      temp_f[0][f] = bc_values_[f];
    }
  }
};

} // namespace Energy
} // namespace Amanzi
