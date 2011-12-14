#include "Epetra_IntVector.h"

#include "DarcyProblem.hh"
#include "boundary-function.hh"
#include "DarcyMatvec.hpp"
#include "flow-bc-factory.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi
{

DarcyProblem::DarcyProblem(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
			   Teuchos::ParameterList &list) : mesh_(mesh)
{
  // Create the combined cell/face DoF map.
  dof_map_ = create_dof_map_(CellMap(), FaceMap());

  face_importer_ = Teuchos::rcp(new Epetra_Import(FaceMap(true),FaceMap(false)));
  cell_importer_ = Teuchos::rcp(new Epetra_Import(CellMap(true),CellMap(false)));

  // Create the MimeticHexLocal objects.
  init_mimetic_disc_(mesh_, MD_);
  md_ = new MimeticHex(mesh); // evolving replacement for mimetic_hex

  // Create the matvec operator for the system
  matvec_ = new DarcyMatvec(this);
};

void  DarcyProblem::InitializeProblem(Teuchos::ParameterList& list) {

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(list.sublist("boundary conditions",true));
  FlowBCFactory bc_factory(mesh_, bc_list);
  bc_press_ = Teuchos::rcp(bc_factory.CreatePressure());
  bc_head_  = Teuchos::rcp(bc_factory.CreateStaticHead(0.0, rho_, gravity_));
  bc_flux_  = Teuchos::rcp(bc_factory.CreateMassFlux());
  validate_boundary_conditions();

  // Evaluate the BC at dummy time 0.
  //TODO: introduce some method to have a time-parametrized, steady Darcy problem.
  bc_press_->Compute(0.0);
  bc_head_->Compute(0.0);
  bc_flux_->Compute(0.0);

  // Create the diffusion matrix (structure only, no values)
  D_ = Teuchos::rcp<DiffusionMatrix>(create_diff_matrix_(mesh_));

  // Create the preconditioner (structure only, no values)
  Teuchos::ParameterList diffprecon_plist = list.sublist("Diffusion Preconditioner");
  precon_ = new DiffusionPrecon(D_, diffprecon_plist, Map());

  // Create the RHS vector.
  rhs_ = new Epetra_Vector(Map());

  // DEFINE DEFAULTS FOR PROBLEM PARAMETERS
  k_.resize(CellMap(true).NumMyElements(), 1.0);
}


DarcyProblem::~DarcyProblem()
{
  delete md_;
  delete precon_;

  delete rhs_;
  delete matvec_;
}


void DarcyProblem::validate_boundary_conditions() const
{
  // Create sets of the face indices belonging to each BC type.
  std::set<int> press_faces, head_faces, flux_faces;
  BoundaryFunction::Iterator bc;
  for (bc = bc_press_->begin(); bc != bc_press_->end(); ++bc) press_faces.insert(bc->first);
  for (bc = bc_head_->begin();  bc != bc_head_->end();  ++bc) head_faces.insert(bc->first);
  for (bc = bc_flux_->begin();  bc != bc_flux_->end();  ++bc) flux_faces.insert(bc->first);
  
  std::set<int> overlap;
  std::set<int>::iterator overlap_end;
  int local_overlap, global_overlap;
  
  // Check for overlap between pressure and static head BC.
  std::set_intersection(press_faces.begin(), press_faces.end(),
                         head_faces.begin(),  head_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  Comm().SumAll(&local_overlap, &global_overlap, 1); //TODO: this will over count ghost faces
  if (global_overlap != 0) {
    Errors::Message m;
    std::stringstream s;
    s << global_overlap;
    m << "DarcyProblem: \"static head\" BC overlap \"dirichlet\" BC on "
      << s.str().c_str() << " faces";
    Exceptions::amanzi_throw(m);
  }
  
  // Check for overlap between pressure and flux BC.
  overlap.clear();
  std::set_intersection(press_faces.begin(), press_faces.end(),
                         flux_faces.begin(),  flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  Comm().SumAll(&local_overlap, &global_overlap, 1); //TODO: this will over count ghost faces
  if (global_overlap != 0) {
    Errors::Message m;
    std::stringstream s;
    s << global_overlap;
    m << "DarcyProblem: \"flux\" BC overlap \"dirichlet\" BC on "
      << s.str().c_str() << " faces";
    Exceptions::amanzi_throw(m);
  }
  
  // Check for overlap between static head and flux BC.
  overlap.clear();
  std::set_intersection(head_faces.begin(), head_faces.end(),
                        flux_faces.begin(), flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  Comm().SumAll(&local_overlap, &global_overlap, 1); //TODO: this will over count ghost faces
  if (global_overlap != 0) {
    Errors::Message m;
    std::stringstream s;
    s << global_overlap;
    m << "DarcyProblem: \"flux\" BC overlap \"static head\" BC on "
      << s.str().c_str() << " faces";
    Exceptions::amanzi_throw(m);
  }
  
  //TODO: Verify that a BC has been applied to every boundary face.
  //      Right now faces without BC are considered no-mass-flux.
}


DiffusionMatrix* DarcyProblem::create_diff_matrix_(Teuchos::RCP<AmanziMesh::Mesh>& mesh) const {
  // Generate the list of all Dirichlet-type faces.
  // The provided list should include USED BC faces.
  std::vector<int> dir_faces;
  for (BoundaryFunction::Iterator i = bc_press_->begin(); i != bc_press_->end(); ++i)
    dir_faces.push_back(i->first);
  for (BoundaryFunction::Iterator i = bc_head_->begin();  i != bc_head_->end();  ++i)
    dir_faces.push_back(i->first);
  return new DiffusionMatrix(mesh, dir_faces);
}


void DarcyProblem::init_mimetic_disc_(Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                                      std::vector<MimeticHexLocal>& MD) const {
  // Local storage for the 8 vertex coordinates of a hexahedral cell.
  double x[8][3];
  double *xBegin = &x[0][0];  // begin iterator
  double *xEnd = xBegin+24;   // end iterator

  MD.resize(mesh->cell_map(true).NumMyElements());
  for (int j = 0; j < MD.size(); ++j) {
    mesh->cell_to_coordinates((unsigned int) j, xBegin, xEnd);
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
  D_->ComputeFaceSchur();

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


void DarcyProblem::SetGravity(const double g[3])
{
  if (g[0] != 0. || g[1] != 0.) {
    Errors::Message message("RichardsProblem currently requires gravity g = g_z");
    Exceptions::amanzi_throw(message);
  } else {
    g_[0] = 0.;
    g_[1] = 0.;
    g_[2] = g[2];
    gravity_ = -g[2];
  }
};

void DarcyProblem::SetPermeability(double k)
{
  /// should verify k > 0
  for (int i = 0; i < k_.size(); ++i) k_[i] = k;
}


void DarcyProblem::SetPermeability(const Epetra_Vector &k)
{
  /// should verify k.Map() is CellMap()
  /// should verify k values are all > 0
  Epetra_Vector k_ovl(mesh_->cell_map(true));
  Epetra_Import importer(mesh_->cell_map(true), mesh_->cell_map(false));

  k_ovl.Import(k, importer, Insert);

  for (int i = 0; i < k_.size(); ++i) k_[i] = k_ovl[i];
}


void DarcyProblem::ComputeF(const Epetra_Vector &X, Epetra_Vector &F)
{
  // The cell and face-based DoF are packed together into the X and F Epetra
  // vectors: cell-based DoF in the first part, followed by the face-based DoF.
  // In addition, only the owned DoF belong to the vectors.

  // Create views into the cell and face segments of X and F
  Epetra_Vector &Pcell_own = *CreateCellView(X);
  Epetra_Vector &Pface_own = *CreateFaceView(X);

  Epetra_Vector &Fcell_own = *CreateCellView(F);
  Epetra_Vector &Fface_own = *CreateFaceView(F);
  
  // Create input cell and face pressure vectors that include ghosts.
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);
  
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);

  // Precompute the Dirichlet-type BC residuals for later use.
  // Impose the Dirichlet boundary data on the face pressure vector.
  //bc_press_->Compute(time);
  std::map<int,double> Fpress;
  for (BoundaryFunction::Iterator bc = bc_press_->begin(); bc != bc_press_->end(); ++bc) {
    Fpress[bc->first] = Pface[bc->first] - bc->second;
    Pface[bc->first] = bc->second;
  }
  //bc_head_->Compute(time);
  std::map<int,double> Fhead;
  for (BoundaryFunction::Iterator bc = bc_head_->begin(); bc != bc_head_->end(); ++bc) {
    Fhead[bc->first] = Pface[bc->first] - bc->second;
    Pface[bc->first] = bc->second;
  }

  // Create cell and face result vectors that include ghosts.
  Epetra_Vector Fcell(CellMap(true));
  Epetra_Vector Fface(FaceMap(true));

  int cface[6];
  double aux1[6], aux2[6], gflux[6];

  Fface.PutScalar(0.0);
  for (int j = 0; j < Pcell.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces((unsigned int) j, (unsigned int*) cface, (unsigned int*) cface+6);
    // Gather the local face pressures int AUX1.
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    // Compute the local value of the diffusion operator.
    double K = (rho_ * k_[j] / mu_);
    MD_[j].diff_op(K, Pcell[j], aux1, Fcell[j], aux2);
    // Gravity contribution
    MD_[j].GravityFlux(g_, gflux);
    for (int k = 0; k < 6; ++k) aux2[k] -= rho_ * K * gflux[k];
    // Scatter the local face result into FFACE.
    for (int k = 0; k < 6; ++k) Fface[cface[k]] += aux2[k];
  }

  // Dirichlet-type condition residuals; overwrite with the pre-computed values.
  std::map<int,double>::const_iterator i;
  for (i = Fpress.begin(); i != Fpress.end(); ++i) Fface[i->first] = i->second;
  for (i = Fhead.begin();  i != Fhead.end();  ++i) Fface[i->first] = i->second;
  
  // Mass flux BC contribution.
  //bc_flux_->Compute(time);
  for (BoundaryFunction::Iterator bc = bc_flux_->begin(); bc != bc_flux_->end(); ++bc)
    Fface[bc->first] += bc->second * md_->face_area_[bc->first];
  
  // Copy owned part of result into the output vectors.
  for (int j = 0; j < Fcell_own.MyLength(); ++j) Fcell_own[j] = Fcell[j];
  for (int j = 0; j < Fface_own.MyLength(); ++j) Fface_own[j] = Fface[j];

  delete &Pcell_own, &Pface_own, &Fcell_own, &Fface_own;
}


Teuchos::RCP<Epetra_Map> DarcyProblem::create_dof_map_(const Epetra_Map &cell_map, const Epetra_Map &face_map) const
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
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(ndof_tot, ndof, gids, 0, cell_map.Comm()));
  delete [] gids;
  return map;
}


Epetra_Vector* DarcyProblem::CreateCellView(const Epetra_Vector &X) const
{
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  return new Epetra_Vector(View, CellMap(), data);
}


Epetra_Vector* DarcyProblem::CreateFaceView(const Epetra_Vector &X) const
{
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  int ncell = CellMap().NumMyElements();
  return new Epetra_Vector(View, FaceMap(), data+ncell);
}


void DarcyProblem::DeriveDarcyVelocity(const Epetra_Vector &X, Epetra_MultiVector &Q) const
{
  // should verify X.Map() is Map()
  // should verify Q.Map() is CellMap(false)
  // should verify Q.NumVectors() is 3

  // Create views into the cell and face segments of X and F
  Epetra_Vector &Pcell = *CreateCellView(X);
  Epetra_Vector &Pface_own = *CreateFaceView(X);

  // Create face vectors that include ghosts.
  Epetra_Vector Pface(FaceMap(true));

  // Populate the face pressure vector from the input.
  Pface.Import(Pface_own, *face_importer_, Insert);

  int cface[6];
  double aux1[6], aux2[6], aux3[3], gflux[6], dummy;

  for (int j = 0; j < CellMap(false).NumMyElements(); ++j) {

  }
  for (int j = 0; j < Pcell.MyLength(); ++j) {
    mesh_->cell_to_faces((unsigned int) j, (unsigned int*) cface, (unsigned int*) cface+6);
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    double K = (k_[j] / mu_);
    MD_[j].diff_op(K, Pcell[j], aux1, dummy, aux2);
    MD_[j].GravityFlux(g_, gflux);
    for (int k = 0; k < 6; ++k) aux2[k] = rho_ * K * gflux[k] - aux2[k];
    MD_[j].CellFluxVector(aux2, aux3);
    Q[0][j] = aux3[0];
    Q[1][j] = aux3[1];
    Q[2][j] = aux3[2];
  }

  delete &Pcell, &Pface_own;
}


void DarcyProblem::DeriveDarcyFlux(const Epetra_Vector &P, Epetra_Vector &F, double &l2_error) const
{
  /// should verify P.Map() is Map()
  /// should verify F.Map() is FaceMap()

  int fdirs[6];
  unsigned int cface[6];
  double aux1[6], aux2[6], gflux[6], dummy;

  // Create a view into the cell pressure segment of P.
  Epetra_Vector &Pcell_own = *CreateCellView(P);

  // Create a copy of the face pressure segment of P that includes ghosts.
  Epetra_Vector &Pface_own = *CreateFaceView(P);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);

  // Create face flux and face count vectors that include ghosts.
  Epetra_Vector Fface(FaceMap(true)); // fills with 0
  Epetra_Vector error(FaceMap(true)); // fills with 0
  Epetra_IntVector count(FaceMap(true)); // fills with 0

  // Process-local assembly of the face mass fluxes.
  for (unsigned int j = 0; j < Pcell_own.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces(j, cface, cface+6);
    // Gather the local face pressures int AUX1.
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    // Compute the local value of the diffusion operator.
    double K = (rho_ * k_[j] / mu_);
    MD_[j].diff_op(K, Pcell_own[j], aux1, dummy, aux2);
    // Gravity contribution
    MD_[j].GravityFlux(g_, gflux);
    for (int k = 0; k < 6; ++k) aux2[k] = rho_ * K * gflux[k] - aux2[k];
    mesh_->cell_to_face_dirs(j, fdirs, fdirs+6);
    // Scatter the local face result into FFACE.
    for (int k = 0; k < 6; ++k) {
      Fface[cface[k]] += fdirs[k] * aux2[k];
      error[cface[k]] += aux2[k]; // sums to the flux discrepancy
      count[cface[k]]++;
    }
  }

  // Global assembly of face mass fluxes into the return vector.
  F.Export(Fface, *face_importer_, Add);

  // Create an owned face count vector that overlays the count vector with ghosts.
  int *count_data;
  count.ExtractView(&count_data);
  Epetra_IntVector count_own(View, FaceMap(false), count_data);

  // Final global assembly of the face count vector.
  count_own.Export(count, *face_importer_, Add);

  // Correct the double counting of fluxes on interior faces and convert
  // mass flux to Darcy flux by dividing by the constant fluid density.
  for (int j = 0; j < F.MyLength(); ++j)
    F[j] = F[j] / (rho_ * count[j]);

  // Create an owned face error vector that overlays the error vector with ghosts.
  double *error_data;
  error.ExtractView(&error_data);
  Epetra_Vector error_own(View, FaceMap(false), error_data);

  // Final global assembly of the flux discrepancy vector.
  error_own.Export(error, *face_importer_, Add);

  // Set the flux discrepancy error to 0 on boundary faces where there was only one flux computed.
  for (int j = 0; j < F.MyLength(); ++j)
    if (count[j] == 1) error_own[j] = 0.0;

  // Compute the norm of the flux discrepancy.
  error_own.Norm2(&l2_error);

  delete &Pcell_own, &Pface_own;
}

} // close namespace Amanzi
