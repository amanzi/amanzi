/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   Debugging object for writing debug cells using VerboseObject.

   ------------------------------------------------------------------------- */
#include "errors.hh"
#include "dbc.hh"
#include "CompositeVector.hh"

#include "Debugger.hh"

namespace Amanzi {

// Constructor
Debugger::Debugger(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                   std::string name,
                   Teuchos::ParameterList& plist,
                   Teuchos::EVerbosityLevel verb_level)
  : verb_level_(verb_level),
    plist_(plist),
    mesh_(mesh),
    formatter_(plist.get<int>("column width", 15),
               plist.get<int>("precision", 7),
               plist.get<int>("header width", 20),
               plist.get<int>("cell number width", 5)),
    name_(name)
{
  vo_ = Teuchos::rcp(new VerboseObject(name, plist));

  std::vector<AmanziMesh::Entity_ID> vcells;

  // cells to debug
  if (plist_.isParameter("debug cells")) {
    auto dcs = plist_.get<Teuchos::Array<int>>("debug cells");
    vcells.insert(vcells.end(), dcs.begin(), dcs.end());
  }

  // faces to debug
  if (plist_.isParameter("debug faces")) {
    auto dfs = plist_.get<Teuchos::Array<int>>("debug faces");
    const auto& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE, true);
    const auto& cell_map = mesh->getMap(AmanziMesh::Entity_kind::CELL, false);

    for (const auto& f : dfs) {
      AmanziMesh::Entity_ID lf = face_map.LID(f);
      if (lf >= 0) {
        // debug the neighboring cells
        auto fcells = mesh->getFaceCells(lf);
        for (const auto& c : fcells)
          if (c < cell_map.NumMyElements()) vcells.emplace_back(cell_map.GID(c));
      }
    }
  }
  AmanziMesh::Entity_ID_View cells;
  vectorToView(cells, vcells);
  set_cells(cells);
}


void
Debugger::set_cells(const AmanziMesh::Entity_ID_View& dc)
{
  dc_.clear();
  dcvo_.clear();

  // make sure they are unique (they may not be because of add_cells()
  // implementation decisions!
  std::set<AmanziMesh::Entity_ID> dc_set(dc.begin(), dc.end());
  Kokkos::resize(dc_gid_, dc_set.size());

  const auto& cell_map = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  int my_pid = mesh_->getComm()->MyPID();
  int size = 0;
  for (const auto& c : dc_set) {
    AmanziMesh::Entity_ID lc = cell_map.LID(c);
    if (lc >= 0) {
      // include the LID
      dc_.emplace_back(lc);
      dc_gid_[size++] = c;

      // make a verbose object for each case
      Teuchos::ParameterList vo_plist;
      vo_plist.sublist("verbose object");
      vo_plist.sublist("verbose object") = plist_.sublist("verbose object");
      vo_plist.sublist("verbose object").set("write on rank", my_pid);
      vo_plist.sublist("verbose object").set("show rank", false);
      dcvo_.emplace_back(Teuchos::rcp(new VerboseObject(mesh_->getComm(), name_, vo_plist)));
    }
  }
}


void
Debugger::add_cells(const AmanziMesh::Entity_ID_View& dc)
{
  AmanziMesh::Entity_ID_View dc_new = get_cells();
  dc_new.insert(dc_new.end(), dc.begin(), dc.end());
  set_cells(dc_new);
}


const AmanziMesh::Entity_ID_View&
Debugger::get_cells() const
{
  return dc_gid_;
}


// Write cell + face info
void
Debugger::WriteCellInfo(bool include_faces)
{
  Teuchos::OSTab tab1 = vo_->getOSTab();
  if (vo_->os_OK(verb_level_)) { *vo_->os() << "Debug Cells Information:" << std::endl; }

  for (int i = 0; i != dc_.size(); ++i) {
    AmanziMesh::Entity_ID c0 = dc_[i];
    AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
    Teuchos::OSTab tab = dcvo_[i]->getOSTab();

    AmanziGeometry::Point c0_centroid = mesh_->getCellCentroid(c0);
    if (dcvo_[i]->os_OK(verb_level_)) {
      *dcvo_[i]->os() << "Cell c(" << c0_gid << ") centroid = " << c0_centroid << std::endl;

      if (include_faces) {
        auto [fnums0, dirs] = mesh_->getCellFacesAndDirections(c0);

        if (dcvo_[i]->os_OK(verb_level_)) {
          for (unsigned int n = 0; n != fnums0.size(); ++n) {
            AmanziMesh::Entity_ID f_gid =
              mesh_->getMap(AmanziMesh::Entity_kind::FACE, true).GID(fnums0[n]);
            AmanziGeometry::Point f_centroid = mesh_->getFaceCentroid(fnums0[n]);
            *dcvo_[i]->os() << "  neighbor face(" << f_gid << ") [dir=" << dirs[n]
                            << "] centroid = " << f_centroid << std::endl;
          }
        }
      }
    }
  }
}


// Write a vector individually.
void
Debugger::WriteVector(const std::string& vname,
                      const Teuchos::Ptr<const CompositeVector>& vec,
                      bool include_faces)
{
  int n_vecs = 0;
  Teuchos::RCP<const Epetra_MultiVector> vec_c;
  if (vec->HasComponent("cell")) {
    vec_c = vec->ViewComponent("cell", false);
    n_vecs = vec_c->NumVectors();
  }

  Teuchos::RCP<const Epetra_MultiVector> vec_f;
  Teuchos::RCP<const Epetra_MultiVector> vec_bf;
  if (vec->HasComponent("face")) {
    vec_f = vec->ViewComponent("face", true);
    n_vecs = vec_f->NumVectors();
  } else if (vec->HasComponent("boundary_face")) {
    vec_bf = vec->ViewComponent("boundary_face", true);
    n_vecs = vec_bf->NumVectors();
  }

  for (int i = 0; i != dc_.size(); ++i) {
    for (int j = 0; j != n_vecs; ++j) {
      AmanziMesh::Entity_ID c0 = dc_[i];
      AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
      Teuchos::OSTab tab = dcvo_[i]->getOSTab();

      if (dcvo_[i]->os_OK(verb_level_)) {
        *dcvo_[i]->os() << formatter_.formatHeader(vname, c0_gid);

        if (vec_c != Teuchos::null) *dcvo_[i]->os() << " " << formatter_.format((*vec_c)[j][c0]);

        if (include_faces && vec_f != Teuchos::null) {
          auto [fnums0, dirs] = mesh_->getCellFacesAndDirections(c0);

          for (unsigned int n = 0; n != fnums0.size(); ++n)
            *dcvo_[i]->os() << " " << formatter_.format((*vec_f)[j][fnums0[n]]);
        }

        if (include_faces && vec_bf != Teuchos::null) {
          auto [fnums0, dirs] = mesh_->getCellFacesAndDirections(c0);

          for (unsigned int n = 0; n != fnums0.size(); ++n) {
            AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(*mesh_, fnums0[n]);
            if (bf >= 0) *dcvo_[i]->os() << " " << formatter_.format((*vec_bf)[j][bf]);
          }
        }
        *dcvo_[i]->os() << std::endl;
      }
    }
  }
}


// Write a vector individually.
void
Debugger::WriteCellVector(const std::string& name, const Epetra_MultiVector& vec)
{
  int n_vecs = vec.NumVectors();
  for (int i = 0; i != dc_.size(); ++i) {
    for (int j = 0; j != n_vecs; ++j) {
      AmanziMesh::Entity_ID c0 = dc_[i];
      AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
      Teuchos::OSTab tab = dcvo_[i]->getOSTab();

      if (dcvo_[i]->os_OK(verb_level_)) {
        *dcvo_[i]->os() << formatter_.formatHeader(name, c0_gid) << formatter_.format(vec[j][c0])
                        << std::endl;
      }
    }
  }
}


// Write list of vectors.
void
Debugger::WriteVectors(const std::vector<std::string>& names,
                       const std::vector<Teuchos::Ptr<const CompositeVector>>& vecs,
                       bool include_faces)
{
  AMANZI_ASSERT(names.size() == vecs.size());

  for (int i = 0; i != dc_.size(); ++i) {
    AmanziMesh::Entity_ID c0 = dc_[i];
    AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
    Teuchos::OSTab tab = dcvo_[i]->getOSTab();

    if (dcvo_[i]->os_OK(verb_level_)) {
      for (int lcv = 0; lcv != names.size(); ++lcv) {
        std::string name = names[lcv];
        Teuchos::Ptr<const CompositeVector> vec = vecs[lcv];

        int n_vec = 1;
        Teuchos::RCP<const Epetra_MultiVector> vec_c;
        if (vec->HasComponent("cell")) {
          vec_c = vec->ViewComponent("cell", false);
          n_vec = vec_c->NumVectors();
        }

        Teuchos::RCP<const Epetra_MultiVector> vec_f;
        Teuchos::RCP<const Epetra_MultiVector> vec_bf;
        if (vec->HasComponent("face")) {
          vec_f = vec->ViewComponent("face", false);
          n_vec = vec_f->NumVectors();
        } else if (vec->HasComponent("boundary_face")) {
          vec_bf = vec->ViewComponent("boundary_face", false);
          n_vec = vec_bf->NumVectors();
        }

        for (int j = 0; j != n_vec; ++j) {
          *dcvo_[i]->os() << formatter_.formatHeader(name, c0_gid);
          if (vec_c != Teuchos::null) *dcvo_[i]->os() << " " << formatter_.format((*vec_c)[j][c0]);

          if (include_faces && vec_f != Teuchos::null) {
            auto [fnums0, dirs] = mesh_->getCellFacesAndDirections(c0);

            for (unsigned int n = 0; n != fnums0.size(); ++n)
              *dcvo_[i]->os() << " " << formatter_.format((*vec_f)[j][fnums0[n]]);
          }

          if (include_faces && vec_bf != Teuchos::null) {
            auto [fnums0, dirs] = mesh_->getCellFacesAndDirections(c0);

            for (unsigned int n = 0; n != fnums0.size(); ++n) {
              AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(*mesh_, fnums0[n]);
              if (bf >= 0) *dcvo_[i]->os() << " " << formatter_.format((*vec_bf)[j][bf]);
            }
          }

          *dcvo_[i]->os() << std::endl;
        }
      }
    }
  }
}


// Write boundary condition data.
void
Debugger::WriteBoundaryConditions(const std::vector<int>& flag, const std::vector<double>& data)
{
  // bcs use 3 extra characters
  int width = formatter_.getWidth();
  int precision = formatter_.getPrecision();
  formatter_.setWidth(width - 2); // note, this may also change precision

  for (int i = 0; i != dc_.size(); ++i) {
    AmanziMesh::Entity_ID c0 = dc_[i];
    AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
    Teuchos::OSTab tab = dcvo_[i]->getOSTab();

    if (dcvo_[i]->os_OK(verb_level_)) {
      *dcvo_[i]->os() << formatter_.formatHeader("BCs", c0_gid);
      auto [fnums0, dirs] = mesh_->getCellFacesAndDirections(c0);

      for (unsigned int n = 0; n != fnums0.size(); ++n)
        *dcvo_[i]->os() << flag[fnums0[n]] << "(" << formatter_.format(data[fnums0[n]]) << ")";
      *dcvo_[i]->os() << std::endl;
    }
  }

  formatter_.setWidth(width);
  formatter_.setPrecision(precision);
}


// call MPI_Comm_Barrier to sync between writing steps
void
Debugger::Barrier()
{
  mesh_->getComm()->Barrier();
}


// write a line of ----
void
Debugger::WriteDivider()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(verb_level_))
    *vo_->os() << "------------------------------------------------------------------------"
               << std::endl;
}


// Reverse order... get the VerboseObject for Entity_
Teuchos::RCP<VerboseObject>
Debugger::GetVerboseObject(AmanziMesh::Entity_ID id, int rank)
{
  std::vector<AmanziMesh::Entity_ID>::iterator loc = std::find(dc_.begin(), dc_.end(), id);
  if (loc == dc_.end()) {
    return Teuchos::null;
  } else {
    int index = loc - dc_.begin();
    return dcvo_[index];
  }
}

} // namespace Amanzi
