/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "gmv_mesh.hh"

#include "flow.hh"

namespace Amanzi {
namespace Flow {

// constructor
Flow::Flow(Teuchos::ParameterList& flow_plist, const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution) :
    flow_plist_(flow_plist) {

  // data layouts for fields
  std::vector<AmanziMesh::Entity_kind> locations(1);
  std::vector<std::string> names(1);
  locations[0] = AmanziMesh::CELL;
  names[0] = "cell";

  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  // require fields
  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  std::vector< std::vector<std::string> > subfield_names(2);
  subfield_names[0].resize(1); subfield_names[0][0] = "pressure";
  subfield_names[1].resize(1); subfield_names[1][0] = "pressure_lambda";
  S->RequireField("pressure", "flow", names2, locations2, 1, true);
  S->GetRecord("pressure","flow")->set_io_vis(true);
  Teuchos::RCP<CompositeVector> pressure = S->GetFieldData("pressure", "flow");
  pressure->set_subfield_names(subfield_names);
  solution->set_data(pressure);
  solution_ = solution;

  // -- secondary variables
  S->RequireField("darcy_flux", "flow", AmanziMesh::FACE, 1, true);

  // -- secondary variables -- flow deals with EOS
  S->RequireField("density_liquid", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("viscosity_liquid", "flow", AmanziMesh::CELL, 1, true);
  
  // -- parameters
  S->RequireField("permeability", "flow", AmanziMesh::CELL, 1, true);

  // -- parameters and variables provided by other PKs/state
  S->RequireField("cell_volume", AmanziMesh::CELL, 1, true);
  S->RequireField("temperature", AmanziMesh::CELL, 1, true);

  // boundary conditions
  Teuchos::ParameterList bc_plist = flow_plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->mesh(), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_head_ = bc_factory.CreateHead();
  bc_flux_ = bc_factory.CreateFlux();

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = flow_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->mesh()));
  matrix_->SetSymmetryProperty(true);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner
  // NOTE: may want to allow these to be the same/different?
  Teuchos::ParameterList mfd_pc_plist = flow_plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->mesh()));
  preconditioner_->SetSymmetryProperty(true);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  Teuchos::ParameterList mfd_pc_ml_plist = mfd_pc_plist.sublist("ML Parameters");
  preconditioner_->InitMLPreconditioner(mfd_pc_ml_plist);
};
  

/* ******************************************************************
* Add a boundary marker to used faces.                                          
****************************************************************** */
void Flow::UpdateBoundaryConditions_() {
  for (int n=0; n<bc_markers.size(); n++) {
    bc_markers[n] = MFD_BC_NULL;
    bc_values[n] = 0.0;
  }

  BoundaryFunction::Iterator bc;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers[f] = MFD_BC_DIRICHLET;
    bc_values[f] = bc->second;
  }

  for (bc=bc_head_->begin(); bc!=bc_head_->end(); ++bc) {
    int f = bc->first;
    bc_markers[f] = MFD_BC_DIRICHLET;
    bc_values[f] = bc->second;
  }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers[f] = MFD_BC_NEUMANN;
    bc_values[f] = bc->second;
  }
};


/* ******************************************************************
* Add a boundary marker to owned faces.                                          
****************************************************************** */
void Flow::ApplyBoundaryConditions_(const State& S,
        const Teuchos::RCP<CompositeVector>& temperature) {
  int nfaces = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers[f] == MFD_BC_DIRICHLET) {
      pressure_faces[f] = bc_values[f];
    }
  }
}


/* ******************************************************************
* Routine updates elemental discretization matrices and must be 
* called before applying boundary conditions and global assembling.                                             
****************************************************************** */
void Flow::AddGravityFluxesToOperator_(const Teuchos::RCP<const State>& S,
        const std::vector<WhetStone::Tensor>& K, const CompositeVector& Krel,
        const Teuchos::RCP<MatrixMFD>& matrix) {

  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("density_liquid");
  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");
  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = g_vec[i];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<Epetra_SerialDenseVector>& Ff_cells = matrix->Ff_cells();
  std::vector<double>& Fc_cells = matrix->Fc_cells();

  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Epetra_SerialDenseVector& Ff = Ff_cells[c];
    double& Fc = Fc_cells[c];

    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = S->mesh()->face_normal(f);

      double outward_flux = (((*rho)("cell", 0, c) * K[c] * gravity) * normal)
                                        * dirs[n] * Krel("face",0,f);
      Ff[n] += outward_flux;
      Fc -= outward_flux;  // Nonzero-sum contribution when flag_upwind = false.
    }
  }
};


/* ******************************************************************
* Updates global Darcy vector calculated by a discretization method.                                             
****************************************************************** */
void Flow::AddGravityFluxesToVector_(Teuchos::RCP<State>& S,
        std::vector<WhetStone::Tensor>& K, const Epetra_Vector& Krel,
        Teuchos::RCP<CompositeVector>& darcy_mass_flux) {

  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("density_liquid");
  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");
  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = g_vec[i];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int f_used = S_->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int f_owned = S_->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  std::vector<bool> done(f_used, false);

  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = S->mesh()->face_normal(f);

      if (f<f_owned && !done[f]) {
        (*darcy_mass_flux)(f) += (((*rho)("cell",0,c) * K[c] * gravity) * normal)
                                        * Krel_faces[f];
        flag[f] = true;
      }
    }
  }
};


/* ****************************************************************
* DEBUG: creating GMV file 
**************************************************************** */
void Flow::WriteGMVfile_(Teuchos::RCP<State> S) const {
  Teuchos::RCP<AmanziMesh::Mesh> mesh = S->mesh();

  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(*S_->GetFieldData("pressure")->ViewComponent("pressure",false),0, "pressure");
  GMV::write_cell_data(*S_->GetFieldData("saturation")->ViewComponent("saturation",false),0, "saturation");
  GMV::write_cell_data(*S_->GetFieldData("darcy_velocity")->ViewComponent("darcy_velocity",false), 0, "velocity_h");
  GMV::write_cell_data(*S_->GetFieldData("darcy_velocity")->ViewComponent("darcy_velocity",false), dim-1, "velocity_v");
  GMV::close_data_file();
};

}  // namespace AmanziFlow
}  // namespace Amanzi

