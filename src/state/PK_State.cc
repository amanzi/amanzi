
#include "PK_State.hh"


namespace Amanzi {

PK_State::PK_State(std::string name, Teuchos::RCP<AmanziMesh::Mesh> mesh) :
    name_(name), mesh_(mesh), ghosted_(false) {
  S_ = Teuchos::rcp(new State());
  S_->RegisterDomainMesh(mesh);
}

PK_State::PK_State(std::string name, Teuchos::RCP<State> S) :
    name_(name), S_(S), ghosted_(false) {
  mesh_ = S_->GetMesh();
}

PK_State::PK_State(std::string name, State& S) :
    name_(name), ghosted_(false) {
  S_ = Teuchos::rcpFromRef(S);
  mesh_ = S_->GetMesh();
}

PK_State::PK_State(PK_State& other, StateConstructMode mode) :
    name_(other.name_),
    mesh_(other.mesh_),
    ghosted_(other.ghosted_) {
  S_ = Teuchos::rcp(new State(*other.S_, mode));
}


/* *******************************************************************
* Copy cell-based data from master to ghost positions.              
* WARNING: vector v must contain ghost cells.              
******************************************************************* */
void PK_State::CopyMasterCell2GhostCell(Epetra_Vector& v)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->cell_map(true);
  Epetra_Import importer(target_cmap, source_cmap);

  double* vdata;
  v.ExtractView(&vdata);
  Epetra_Vector vv(View, source_cmap, vdata);

  v.Import(vv, importer, Insert);
#endif
}



/* *******************************************************************
* Copy cell-based data from master to ghost positions.              
* WARNING: MultiVector v must contain ghost cells.              
******************************************************************* */
void PK_State::CopyMasterMultiCell2GhostMultiCell(Epetra_MultiVector& v)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->cell_map(true);
  Epetra_Import importer(target_cmap, source_cmap);

  double** vdata;
  v.ExtractView(&vdata);
  Epetra_MultiVector vv(View, source_cmap, vdata, v.NumVectors());

  v.Import(vv, importer, Insert);
#endif
}


/* *******************************************************************
* Routine imports a short multivector to a parallel overlapping multivector.                
******************************************************************* */
void PK_State::CopyMasterMultiCell2GhostMultiCell(const Epetra_MultiVector& v, 
                                                         Epetra_MultiVector& vv,
                                                         int parallel_comm)
{
#ifdef HAVE_MPI
  if (parallel_comm == 1) {
    const Epetra_BlockMap& source_cmap = mesh_->cell_map(false);
    const Epetra_BlockMap& target_cmap = mesh_->cell_map(true);
    Epetra_Import importer(target_cmap, source_cmap);

    vv.Import(v, importer, Insert);
  } else {
    int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    int num_vectors = v.NumVectors();
    for (int c = 0; c < ncells_owned; c++) {
      for (int i = 0; i < num_vectors; i++) vv[i][c] = v[i][c];
    }
  }
#else
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int num_vectors = v.NumVectors();
  for (int c = 0; c < ncells_owned; c++) {
    for (int i = 0; i < num_vectors; i++) vv[i][c] = v[i][c];
  }
#endif
}

/* *******************************************************************
* Copy face-based data from master to ghost positions.              
* WARNING: vector v must contain ghost cells.              
******************************************************************* */
void PK_State::CopyMasterFace2GhostFace(Epetra_Vector& v)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->face_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->face_map(true);
  Epetra_Import importer(target_cmap, source_cmap);

  double* vdata;
  v.ExtractView(&vdata);
  Epetra_Vector vv(View, source_cmap, vdata);

  v.Import(vv, importer, Insert);
#endif
}


/* *******************************************************************
* Transfers face-based data from ghost to master positions and 
* performs the operation 'mode' there. 
* WARNING: Vector v must contain ghost faces.              
******************************************************************* */
void PK_State::CombineGhostFace2MasterFace(Epetra_Vector& v, Epetra_CombineMode mode)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_fmap = mesh_->face_map(false);
  const Epetra_BlockMap& target_fmap = mesh_->face_map(true);
  Epetra_Import importer(target_fmap, source_fmap);

  double* vdata;
  v.ExtractView(&vdata);
  Epetra_Vector vv(View, source_fmap, vdata);

  vv.Export(v, importer, mode);
#endif
}


/* *******************************************************************
* Copy face-based data from master to ghost positions.              
* WARNING: vector vghost must contain ghost cells.              
******************************************************************* */
void PK_State::CopyMasterCell2GhostCell(const Epetra_Vector& v, Epetra_Vector& vghost)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->cell_map(true);
  Epetra_Import importer(target_cmap, source_cmap);
 
  vghost.Import(v, importer, Insert);
#else
  vghost = v;
#endif
}


/* *******************************************************************
* Transfers cell-based data from ghost to master positions and 
* performs the operation 'mode' there. 
* WARNING: Vector v must contain ghost cells.              
******************************************************************* */
void PK_State::CombineGhostCell2MasterCell(Epetra_Vector& v, Epetra_CombineMode mode)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->cell_map(true);
  Epetra_Import importer(target_cmap, source_cmap);

  double* vdata;
  v.ExtractView(&vdata);
  Epetra_Vector vv(View, source_cmap, vdata);

  vv.Export(v, importer, mode);
#endif
}


/* *******************************************************************
* Copy face-based data from master to ghost positions.              
* WARNING: vector vhost must contain ghost cells.              
******************************************************************* */
void PK_State::CopyMasterFace2GhostFace(const Epetra_Vector& v, Epetra_Vector& vghost)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->face_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->face_map(true);
  Epetra_Import importer(target_cmap, source_cmap);
 
  vghost.Import(v, importer, Insert);
#else
  vghost = v;
#endif
}


/* *******************************************************************
* Lp norm of the vector v1.    
******************************************************************* */
double PK_State::normLpCell(const Epetra_Vector& v1, double p)
{
  int ncells = (mesh_->cell_map(false)).NumMyElements();

  double Lp_norm, Lp = 0.0;
  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->cell_volume(c);
    Lp += volume * pow(v1[c], p);
  }
  v1.Comm().MaxAll(&Lp, &Lp_norm, 1);

  return pow(Lp_norm, 1.0/p);
}


/* *******************************************************************
* Lp norm of the component-wise product v1 .* v2             
******************************************************************* */
double PK_State::normLpCell(const Epetra_Vector& v1, const Epetra_Vector& v2, double p)
{
  int ncells = (mesh_->cell_map(false)).NumMyElements();

  double Lp_norm, Lp = 0.0;
  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->cell_volume(c);
    Lp += volume * pow(v1[c] * v2[c], p);
  }
  v1.Comm().MaxAll(&Lp, &Lp_norm, 1);

  return pow(Lp_norm, 1.0/p);
}


/* *******************************************************************
* Extract cells from a supervector             
******************************************************************* */
Epetra_Vector* PK_State::CreateCellView(const Epetra_Vector& u) const
{
  double* data;
  u.ExtractView(&data);
  return new Epetra_Vector(View, mesh_->cell_map(false), data);
}


/* *******************************************************************
* Extract faces from a supervector             
******************************************************************* */
Epetra_Vector* PK_State::CreateFaceView(const Epetra_Vector& u) const
{
  double* data;
  u.ExtractView(&data);
  int ncells = (mesh_->cell_map(false)).NumMyElements();
  return new Epetra_Vector(View, mesh_->face_map(false), data+ncells);
}


/* *******************************************************************
* Calculate minimum values in a multivector using master cells.            
******************************************************************* */
void PK_State::MinValueMasterCells(Epetra_MultiVector& v, double* vmin)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& cmap = mesh_->cell_map(false);

  double** vdata;
  v.ExtractView(&vdata);
  Epetra_MultiVector vv(View, cmap, vdata, v.NumVectors());
  
  vv.MinValue(vmin);
#else
  v.MinValues(vmin);
#endif
}


/* *******************************************************************
* Calculate maximum values in a multivector using master cells.             
******************************************************************* */
void PK_State::MaxValueMasterCells(Epetra_MultiVector& v, double* vmax)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& cmap = mesh_->cell_map(false);

  double** vdata;
  v.ExtractView(&vdata);
  Epetra_MultiVector vv(View, cmap, vdata, v.NumVectors());
  
  vv.MaxValue(vmax);
#else
  v.MaxValues(vmax);
#endif
}


} // namespace
