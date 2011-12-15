/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

//namespace Amanzi {

MimeticDiscretization::MimeticDiscretization(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
        Teuchos::ParameterList& disc_plist) : mesh_(mesh) {
  init_mimetic_disc_(mesh_, MD_);
  md_ = new MimeticHex(mesh_); // evolving replacement for mimetic_hex
};

void PermafrostProblem::init_mimetic_disc_(Teuchos::RCP<AmanziMesh::Mesh> &mesh,
        std::vector<MimeticHexLocal> &MD) const {
  // Local storage for the 8 vertex coordinates of a hexahedral cell.
  double x[8][3];
  double *xBegin = &x[0][0];  // begin iterator
  double *xEnd = xBegin+24;   // end iterator

  MD.resize(mesh->cell_map(true).NumMyElements());
  for (int j = 0; j < MD.size(); ++j) {
    mesh->cell_to_coordinates((unsigned int) j, xBegin, xEnd);
    MD[j].update(x);
  }
};

void MimeticDiscretization::CalcCellVolumes(Teuchos::RCP<Epetra_Vector> &cell_volumes) const {
  ASSERT(cell_volumes->Map().SameAs(CellMap()));
  // cell volumes must be non-ghosted
  for (int n=0; n<(CellMap()).NumMyElements(); n++) {
    (*cell_volumes)[n] = md_->Volume(n);
  }
};
// } // namespace
