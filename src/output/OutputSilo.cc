/*
  Output

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Silo implementation of an Output object.
*/

#include <locale>
#include <iomanip>

#include "errors.hh"
#include "Point.hh"
#include "Mesh.hh"
#include "OutputSilo.hh"

namespace Amanzi {

OutputSilo::OutputSilo(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                       bool is_vis,
                       bool is_dynamic)
    : mesh_(mesh),
      fid_(NULL),
      count_(0)
{
  if (mesh_->vis_mesh().get_comm()->NumProc() > 1) {
    Errors::Message msg("OutputSilo does not yet support parallel runs.");
    Exceptions::amanzi_throw(msg);  
  }
  if (mesh_->vis_mesh().space_dimension() != 3 || mesh_->vis_mesh().manifold_dimension() != 3) {
    Errors::Message msg("OutputSilo is untested on non-3D meshes.");
    Exceptions::amanzi_throw(msg);  
  }
  Init_(plist);
}


// destructor must release file resource on non-finalized
OutputSilo::~OutputSilo() {
  if (fid_) {
    CloseFile_();
  }
}


// open and close files
void
OutputSilo::InitializeCycle(double time, int cycle, const std::string& tag) {
  // check not open
  AMANZI_ASSERT(fid_ == NULL);
  if (fid_) {
    CloseFile_();
  }

  // open file
  std::stringstream filename;
  filename << filenamebase_ << "_"
           << std::setfill('0') << std::setw(sigfigs_) 
           << count_;
  std::string fname = filename.str() + ".silo";
  fid_ = DBCreate(fname.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);

  // write mesh
  // -- coordinate names
  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"x-coords";
  coordnames[1] = (char*)"y-coords";
  coordnames[2] = (char*)"z-coords";

  // -- nodal coordinates
  int nnodes = mesh_->vis_mesh().num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  std::vector<double> x(nnodes), y(nnodes), z(nnodes);

  AmanziGeometry::Point xyz;
  for (int i=0; i!=nnodes; ++i) {
    mesh_->vis_mesh().node_get_coordinates(i, &xyz);
    x[i] = xyz[0];
    y[i] = xyz[1];
    z[i] = xyz.dim() > 2 ? xyz[2] : 0.0;
  }
  double *coords[3];
  coords[0] = &x[0];
  coords[1] = &y[0];
  coords[2] = &z[0];

  // -- DBOptions
  DBoptlist* optlist = DBMakeOptlist(10);
  std::string str = "zonelist";
  std::vector<char> writable(str.begin(), str.end());
  writable.push_back('\0');

  DBAddOption(optlist, DBOPT_PHZONELIST, (void*)&writable[0]);
  DBAddOption(optlist, DBOPT_CYCLE, &cycle);
  DBAddOption(optlist, DBOPT_DTIME, &time);

  // -- write the base mesh UCD object
  int ncells = mesh_->vis_mesh().num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ierr = DBPutUcdmesh(fid_, "mesh", mesh_->vis_mesh().space_dimension(),
                          (char const* const*)coordnames, coords,
                          nnodes, ncells, 0, 0, DB_DOUBLE, optlist);
  AMANZI_ASSERT(!ierr);

  // -- Construct the silo face-node info.
  // We rely on the mesh having the faces nodes arranged counter-clockwise
  // around the face.  This should be satisfied by AmanziMesh.
  int nfaces = mesh_->vis_mesh().num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  std::vector<int> face_node_counts(nfaces);
  std::vector<char> ext_faces(nfaces, 0x0);
  std::vector<int> face_node_list;

  for (int f=0; f!=nfaces; ++f) {
    AmanziMesh::Entity_ID_List fnodes;
    mesh_->vis_mesh().face_get_nodes(f, &fnodes);
    face_node_counts[f] = fnodes.size();
    face_node_list.insert(face_node_list.end(), fnodes.begin(), fnodes.end());

    AmanziMesh::Entity_ID_List fcells;
    mesh_->vis_mesh().face_get_cells(f, AmanziMesh::Parallel_type::ALL, &fcells);
    if (fcells.size() == 1) {
      ext_faces[f] = 0x1;
    }
  }

  // -- Construct the silo cell-face info
  std::vector<int> cell_face_counts(ncells);
  std::vector<int> cell_face_list;
  
  for (int c=0; c!=ncells; ++c) {
    AmanziMesh::Entity_ID_List cfaces;
    std::vector<int> dirs;
    mesh_->vis_mesh().cell_get_faces_and_dirs(c, &cfaces, &dirs, false);
    for (int i=0; i!=cfaces.size(); ++i) {
      if (dirs[i] < 0) cfaces[i] = ~cfaces[i];
    }
    
    cell_face_counts[c] = cfaces.size();

    cell_face_list.insert(cell_face_list.end(), cfaces.begin(), cfaces.end());
  }
  
  ierr |= DBPutPHZonelist(fid_, "zonelist", nfaces, &face_node_counts[0],
                          face_node_list.size(), &face_node_list[0], &ext_faces[0],
                          ncells, &cell_face_counts[0],
                          cell_face_list.size(), &cell_face_list[0],
                          0, 0, ncells-1, optlist);
  AMANZI_ASSERT(!ierr);
                          
  // -- clean up (could be done on finalize if needed? --etc)
  DBFreeOptlist(optlist);

}

void
OutputSilo::FinalizeCycle() {
  if (fid_) {
    CloseFile_();
    count_++;
  }
}


// write data to file
void
OutputSilo::WriteVector(const Epetra_Vector& vec, const std::string& name,
                        const AmanziMesh::Entity_kind& kind) const {
  if (kind == AmanziMesh::CELL) {
    int ierr = DBPutUcdvar1(fid_, FixName_(name).c_str(), "mesh",
                            (void*)&vec[0], vec.MyLength(), NULL, 0,
                            DB_DOUBLE, DB_ZONECENT, NULL);
    AMANZI_ASSERT(!ierr);
  } else if (kind == AmanziMesh::NODE) {
    int ierr = DBPutUcdvar1(fid_, FixName_(name).c_str(), "mesh",
                            (void*)&vec[0], vec.MyLength(), NULL, 0,
                            DB_DOUBLE, DB_NODECENT, NULL);
    AMANZI_ASSERT(!ierr);
  } else {
    Errors::Message msg("OutputSilo only knows how to write CELL and NODE based quantities.");
    Exceptions::amanzi_throw(msg);
  }    
}


void
OutputSilo::WriteMultiVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names,
                             const AmanziMesh::Entity_kind& kind) const {
  AMANZI_ASSERT(vec.NumVectors() == names.size());
  for (int i=0; i!=vec.NumVectors(); ++i) {
    WriteVector(*vec(i), names[i], kind);
  }
}

// can we template this?
void
OutputSilo::WriteAttribute(const double& val, const std::string& name) const {
  AMANZI_ASSERT(0);
}

void
OutputSilo::WriteAttribute(const int& val, const std::string& name) const {
  AMANZI_ASSERT(0);
}


void
OutputSilo::WriteAttribute(const std::string& val, const std::string& name) const {
  AMANZI_ASSERT(0);
}


void
OutputSilo::Init_(Teuchos::ParameterList& plist) {
  filenamebase_ = plist.get<std::string>("file name base", "amanzi_vis");
  sigfigs_ = plist.get<int>("file name counter digits", 5);
}


void
OutputSilo::CloseFile_() {
  DBClose(fid_);
  fid_ = NULL;
}

std::string
OutputSilo::FixName_(const std::string& s) const {
  int n = s.size(), wp = 0;
  std::vector<char> result(n);
  for (int i=0; i<n; ++i) {
    if (!std::isalnum(s[i])) {
      result[wp++] = '_';
    } else {
      result[wp++] = s[i];
    }
  }
  return std::string(&result[0], &result[wp]);
}

void
OutputSilo::ReadThrowsError_() const {
  Errors::Message msg("OutputSilo does not yet support read -- only write.");
  Exceptions::amanzi_throw(msg);  
}

} // namespace Amanzi
