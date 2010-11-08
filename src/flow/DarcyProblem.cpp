#include "DarcyProblem.hpp"
#include "DarcyMatvec.hpp"

DarcyProblem::DarcyProblem(Teuchos::RCP<Mesh_maps_base> &mesh, Teuchos::RCP<FlowBC> &bc) : mesh_(mesh), bc_(bc)
{
  // Create the combined cell/face DoF map.
  dof_map_ = create_dof_map_(CellMap(), FaceMap());

  face_importer_ = new Epetra_Import(FaceMap(true),FaceMap(false));

  // Create the MimeticHexLocal objects.
  init_mimetic_disc_(*mesh, MD);
  md_ = new MimeticHex(mesh); // evolving replacement for mimetic_hex

  // Create the matvec operator for the system
  matvec_ = new DarcyMatvec(this);

  // Create the diffusion matrix (structure only, no values)
  D_ = Teuchos::rcp<DiffusionMatrix>(create_diff_matrix_(mesh, bc));

  // Create the preconditioner (structure only, no values)
  precon_ = new DiffusionPrecon(D_, Map());

  // Create the RHS vector.
  rhs_ = new Epetra_Vector(Map());

  // DEFINE DEFAULTS FOR PROBLEM PARAMETERS
  rho_ = 1.0;
  mu_  = 1.0;
  k_.resize(CellMap().NumMyElements(), 1.0);
  for (int i = 0; i < 3; ++i) g_[i] = 0.0; // no gravity
}


DarcyProblem::~DarcyProblem()
{
  delete dof_map_;
  delete rhs_;
  delete matvec_;
  delete precon_;
}


DiffusionMatrix* DarcyProblem::create_diff_matrix_(Teuchos::RCP<Mesh_maps_base> &mesh, Teuchos::RCP<FlowBC> &bc) const
{
  // Generate the list of all Dirichlet-type faces.
  std::vector<int> dir_faces;
  for (int j = 0; j < bc->NumBC(); ++j) {
    FlowBC::bc_spec& BC = (*bc)[j];
    if (BC.Type == FlowBC::PRESSURE_CONSTANT) {
      dir_faces.reserve(dir_faces.size() + BC.Faces.size());
      for (int i = 0; i < BC.Faces.size(); ++i)
        dir_faces.push_back(BC.Faces[i]);
    }
  }
  return new DiffusionMatrix(mesh, dir_faces);
}


void DarcyProblem::init_mimetic_disc_(Mesh_maps_base &mesh, std::vector<MimeticHexLocal> &MD) const
{
  // Local storage for the 8 vertex coordinates of a hexahedral cell.
  double x[8][3];
  double *xBegin = &x[0][0];  // begin iterator
  double *xEnd = xBegin+24;   // end iterator

  MD.resize(mesh.cell_map(false).NumMyElements());
  for (int j = 0; j < MD.size(); ++j) {
    mesh.cell_to_coordinates((unsigned int) j, xBegin, xEnd);
    MD[j].update(x);
  }
}


void DarcyProblem::Assemble()
{
  // Fill the diffusion matrix with values.
  std::vector<double> K(k_);
  for (int j = 0; j < K.size(); ++j) K[j] = rho_ * K[j] / mu_;
  D_->Compute(K);

  // Compute the face Schur complement of the diffusion matrix.
  D_->ComputeFaceSchur(); //D_->Print();

  // Compute the preconditioner from the newly computed diffusion matrix and Schur complement.
  precon_->Compute();

  // Compute the RHS.
  ComputeF(*rhs_, *rhs_); // bad form, I know -- input and output are aliased
  (*rhs_).Scale(-1.0);
}


void DarcyProblem::SetFluidDensity(double rho)
{
  /// should verify rho > 0
  rho_ = rho;
}


void DarcyProblem::SetFluidViscosity(double mu)
{
  /// should verify mu > 0
  mu_ = mu;
}


void DarcyProblem::SetGravity(double g[3])
{
  for (int i = 0; i < 3; ++i) g_[i] = g[i];
}


void DarcyProblem::SetPermeability(double k)
{
  /// should verify k > 0
  for (int i = 0; i < k_.size(); ++i) k_[i] = k;
}


void DarcyProblem::SetPermeability(std::vector<double> &k)
{
  /// should verify k has the expected size
  /// should verify k values are all > 0
  for (int i = 0; i < k_.size(); ++i) k_[i] = k[i];
}


void DarcyProblem::ComputeF(const Epetra_Vector &X, Epetra_Vector &F)
{
  // The cell and face-based DoF are packed together into the X and F Epetra
  // vectors: cell-based DoF in the first part, followed by the face-based DoF.
  // In addition, only the owned DoF belong to the vectors.

  // Create views into the cell and face segments of X and F
  Epetra_Vector &Pcell = *CreateCellView(X);
  Epetra_Vector &Pface_own = *CreateFaceView(X);

  Epetra_Vector &Fcell = *CreateCellView(F);
  Epetra_Vector &Fface_own = *CreateFaceView(F);

  // Create face vectors that include ghosts.
  Epetra_Vector Pface(FaceMap(true));
  Epetra_Vector Fface(FaceMap(true));

  // Populate the face pressure vector from the input.
  Pface.Import(Pface_own, *face_importer_, Insert);

  // Apply initial BC fixups to PFACE.
  apply_BC_initial_(Pface); // modifies used values

  int cface[6];
  double aux1[6], aux2[6];

  Fface.PutScalar(0.0);
  for (int j = 0; j < Pcell.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces((unsigned int) j, (unsigned int*) cface, (unsigned int*) cface+6);
    // Gather the local face pressures int AUX1.
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    // Compute the local value of the diffusion operator.
    MD[j].diff_op((rho_ * k_[j] / mu_), Pcell[j], aux1, Fcell[j], aux2);
    // Scatter the local face result into FFACE.
    for (int k = 0; k < 6; ++k) Fface[cface[k]] += aux2[k];
  }
  Fface_own.Export(Fface, *face_importer_, Add);

  // Apply final BC fixups to FFACE.
  apply_BC_final_(Fface_own); // only modifies owned values

  delete &Pcell, &Pface_own, &Fcell, &Fface_own;
}


// BC fixups for F computation: initial pass.
void DarcyProblem::apply_BC_initial_(Epetra_Vector &Pface)
{
  for (int j = 0; j < (*bc_).NumBC(); ++j) {
    FlowBC::bc_spec& BC = (*bc_)[j];
    switch (BC.Type) {
      case FlowBC::PRESSURE_CONSTANT:
        for (int i = 0; i < BC.Faces.size(); ++i) {
          int n = BC.Faces[i];
          BC.Aux[i] = Pface[n] - BC.Value;
          Pface[n] = BC.Value;
        }
        break;
    }
  }
}


// BC fixups for F computation: final pass.
void DarcyProblem::apply_BC_final_(Epetra_Vector &Fface)
{
  for (int j = 0; j < (*bc_).NumBC(); ++j) {
    FlowBC::bc_spec& BC = (*bc_)[j];
    switch (BC.Type) {
      case FlowBC::PRESSURE_CONSTANT:
        for (int i = 0; i < BC.Faces.size(); ++i) {
          int n = BC.Faces[i];
          Fface[n] = BC.Aux[i];
        }
        break;
// NOT IMPLEMENTED YET -- MISSING DATA
      case FlowBC::DARCY_CONSTANT:
        for (int i = 0; i < BC.Faces.size(); ++i) {
          int n = BC.Faces[i];
          Fface[n] += rho_ * BC.Value * md_->face_area_[n]; // NB: assumes face oriented out of domain
        }
        break;
      case FlowBC::NO_FLOW:
        // The do-nothing boundary condition.
        break;
    }
  }
}


Epetra_Map* DarcyProblem::create_dof_map_(const Epetra_Map &cell_map, const Epetra_Map &face_map) const
{
  // Create the combined cell/face DoF map
  int ncell_tot = cell_map.NumGlobalElements();
  int ndof_tot = ncell_tot + face_map.NumGlobalElements();
  int ncell = cell_map.NumMyElements();
  int ndof = ncell + face_map.NumMyElements();
  int *gids = new int[ndof];
  cell_map.MyGlobalElements(&(gids[0]));
  face_map.MyGlobalElements(&(gids[ncell]));
  for (int i = ncell; i < ndof; ++i) gids[i] += ncell_tot;
  Epetra_Map *map = new Epetra_Map(ndof_tot, ndof, gids, 0, cell_map.Comm());
  delete [] gids;
  return map;
}


Epetra_Vector* DarcyProblem::CreateCellView(const Epetra_Vector &X)
{
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  return new Epetra_Vector(View, CellMap(), data);
}


Epetra_Vector* DarcyProblem::CreateFaceView(const Epetra_Vector &X)
{
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  int ncell = CellMap().NumMyElements();
  return new Epetra_Vector(View, FaceMap(), data+ncell);
}
