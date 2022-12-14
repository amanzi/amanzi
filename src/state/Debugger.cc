/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include <boost/format.hpp>

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
    width_(15),
    header_width_(20),
    precision_(10),
    cellnum_width_(5),
    decimal_width_(7),
    name_(name)
{
  vo_ = Teuchos::rcp(new VerboseObject(name, plist));

  AmanziMesh::Entity_ID_List cells;

  // cells to debug
  if (plist_.isParameter("debug cells")) {
    auto dcs = plist_.get<Teuchos::Array<int>>("debug cells");
    cells.insert(cells.end(), dcs.begin(), dcs.end());
  }

  // faces to debug
  if (plist_.isParameter("debug faces")) {
    auto dfs = plist_.get<Teuchos::Array<int>>("debug faces");
    const auto& face_map = mesh->face_map(true);
    const auto& cell_map = mesh->cell_map(false);

    for (const auto& f : dfs) {
      AmanziMesh::Entity_ID lf = face_map.LID(f);
      if (lf >= 0) {
        // debug the neighboring cells
        AmanziMesh::Entity_ID_List fcells;
        mesh->face_get_cells(lf, AmanziMesh::Parallel_type::OWNED, &fcells);
        for (const auto& c : fcells) cells.emplace_back(cell_map.GID(c));
      }
    }
  }

  set_cells(cells);

  // formatting
  cellnum_width_ = plist_.get<int>("cell number width", cellnum_width_);
  decimal_width_ = plist_.get<int>("decimal width", decimal_width_);
  header_width_ = plist_.get<int>("header width", header_width_);
  width_ = plist_.get<int>("column width", width_);
  precision_ = plist_.get<int>("precision", precision_);
}


const AmanziMesh::Entity_ID_List&
Debugger::get_cells() const
{
  return dc_gid_;
}


void
Debugger::set_cells(const AmanziMesh::Entity_ID_List& dc)
{
  dc_.clear();
  dc_gid_.clear();
  dcvo_.clear();

  // make sure they are unique (they may not be because of add_cells()
  // implementation decisions!
  std::set<AmanziMesh::Entity_ID> dc_set(dc.begin(), dc.end());

  const auto& cell_map = mesh_->cell_map(false);
  int my_pid = mesh_->get_comm()->MyPID();
  for (const auto& c : dc_set) {
    AmanziMesh::Entity_ID lc = cell_map.LID(c);
    if (lc >= 0) {
      // include the LID
      dc_.emplace_back(lc);
      dc_gid_.emplace_back(c);

      // make a verbose object for each case
      Teuchos::ParameterList vo_plist;
      vo_plist.sublist("verbose object");
      vo_plist.sublist("verbose object") = plist_.sublist("verbose object");
      vo_plist.sublist("verbose object").set("write on rank", my_pid);
      vo_plist.sublist("verbose object").set("show rank", false);
      dcvo_.emplace_back(Teuchos::rcp(new VerboseObject(mesh_->get_comm(), name_, vo_plist)));
    }
  }
}


void
Debugger::add_cells(const AmanziMesh::Entity_ID_List& dc)
{
  AmanziMesh::Entity_ID_List dc_new = get_cells();
  dc_new.insert(dc_new.end(), dc.begin(), dc.end());
  set_cells(dc_new);
}


std::string
Debugger::Format_(double dat)
{
  std::stringstream datastream;
  if (dat == 0.) {
    std::stringstream formatstream1;
    formatstream1 << boost::format("%%%ds") % (decimal_width_);
    std::stringstream formatstream2;
    formatstream2 << boost::format("%%%ds") % (decimal_width_ - 1);

    datastream << boost::format(formatstream1.str()) % std::string("") << "0."
               << boost::format(formatstream2.str()) % std::string("");
  } else {
    int mag = std::floor(std::log10(std::abs(dat)));
    if (mag < decimal_width_ && mag > -2) { // fixed format
      std::stringstream formatstream1;
      formatstream1 << boost::format("%%%ds") % (decimal_width_ - mag - 1);

      std::stringstream formatstream2;
      //formatstream2 << boost::format("%%%d.%df") % (width_ - (decimal_width_-mag-1)) % decimal_width_;
      if (mag < 0) {
        formatstream2 << boost::format("%%%df") % (width_ - (decimal_width_ - mag - 2));
      } else {
        formatstream2 << boost::format("%%%df") % (width_ - (decimal_width_ - mag - 1));
      }

      datastream << boost::format(formatstream1.str()) % std::string("")
                 << boost::format(formatstream2.str()) % dat;

    } else { // sci format
      std::stringstream formatstream1;
      formatstream1 << boost::format("%%%ds") % (decimal_width_ - 5);

      std::stringstream formatstream2;
      //formatstream2 << boost::format("%%%d.%df") % (width_ - (decimal_width_-mag-1)) % decimal_width_;
      formatstream2 << boost::format("%%%de") % (width_ - (decimal_width_ - 5));

      datastream << boost::format(formatstream1.str()) % std::string("")
                 << boost::format(formatstream2.str()) % dat;
    }
  }
  return datastream.str();
}


// note: vector name argument taken by value intentionally
std::string
Debugger::FormatHeader_(std::string header_prefix, int c)
{
  int header_prefix_width = header_width_ - cellnum_width_ - 4;
  if (header_prefix.size() > header_prefix_width) {
    header_prefix.erase(header_prefix_width);
  } else if (header_prefix.size() < header_prefix_width) {
    header_prefix.append(header_prefix_width - header_prefix.size(), ' ');
  }

  std::stringstream headerstream;
  headerstream << header_prefix << "(" << std::setw(cellnum_width_) << std::right << c << "): ";
  return headerstream.str();
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

    AmanziGeometry::Point c0_centroid = mesh_->cell_centroid(c0);
    if (dcvo_[i]->os_OK(verb_level_)) {
      *dcvo_[i]->os() << "Cell c(" << c0_gid << ") centroid = " << c0_centroid << std::endl;

      if (include_faces) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(c0, &fnums0, &dirs);

        if (dcvo_[i]->os_OK(verb_level_)) {
          for (unsigned int n = 0; n != fnums0.size(); ++n) {
            AmanziMesh::Entity_ID f_gid = mesh_->face_map(true).GID(fnums0[n]);
            AmanziGeometry::Point f_centroid = mesh_->face_centroid(fnums0[n]);
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
  int nfaces_valid = 0;
  if (vec->HasComponent("face")) {
    vec_f = vec->ViewComponent("face", true);
    nfaces_valid = vec_f->MyLength();
    n_vecs = vec_f->NumVectors();
  }

  Teuchos::RCP<const Epetra_MultiVector> vec_bf;
  int nbfaces_valid = 0;
  if (vec->HasComponent("boundary_face")) {
    vec_bf = vec->ViewComponent("boundary_face", true);
    nbfaces_valid = vec_bf->MyLength();
    n_vecs = vec_bf->NumVectors();
  }

  for (int i = 0; i != dc_.size(); ++i) {
    for (int j = 0; j != n_vecs; ++j) {
      AmanziMesh::Entity_ID c0 = dc_[i];
      AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
      Teuchos::OSTab tab = dcvo_[i]->getOSTab();

      if (dcvo_[i]->os_OK(verb_level_)) {
        *dcvo_[i]->os() << FormatHeader_(vname, c0_gid);

        if (vec_c != Teuchos::null) *dcvo_[i]->os() << Format_((*vec_c)[j][c0]);

        if (include_faces && vec_f != Teuchos::null) {
          AmanziMesh::Entity_ID_List fnums0;
          std::vector<int> dirs;
          mesh_->cell_get_faces_and_dirs(c0, &fnums0, &dirs);

          for (unsigned int n = 0; n != fnums0.size(); ++n)
            if (fnums0[n] < nfaces_valid) *dcvo_[i]->os() << " " << Format_((*vec_f)[j][fnums0[n]]);
        }

        if (include_faces && vec_bf != Teuchos::null) {
          AmanziMesh::Entity_ID_List fnums0;
          std::vector<int> dirs;
          mesh_->cell_get_faces_and_dirs(c0, &fnums0, &dirs);

          for (unsigned int n = 0; n != fnums0.size(); ++n) {
            AmanziMesh::Entity_ID f = fnums0[n];
            AmanziMesh::Entity_ID bf =
              mesh_->exterior_face_map(true).LID(mesh_->face_map(true).GID(f));
            if (bf >= 0 && bf < nbfaces_valid) *dcvo_[i]->os() << " " << Format_((*vec_bf)[j][bf]);
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
        *dcvo_[i]->os() << FormatHeader_(name, c0_gid) << Format_(vec[j][c0]) << std::endl;
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

  std::stringstream formatstream;
  formatstream << "%_" << width_ << "." << precision_ << "g";
  std::string format = formatstream.str();

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
        if (vec->HasComponent("face")) {
          vec_f = vec->ViewComponent("face", false);
          n_vec = vec_f->NumVectors();
        }

        Teuchos::RCP<const Epetra_MultiVector> vec_bf;
        if (vec->HasComponent("boundary_face")) {
          vec_bf = vec->ViewComponent("boundary_face", false);
          n_vec = vec_bf->NumVectors();
        }

        for (int j = 0; j != n_vec; ++j) {
          *dcvo_[i]->os() << FormatHeader_(name, c0_gid);
          if (vec_c != Teuchos::null) *dcvo_[i]->os() << Format_((*vec_c)[j][c0]);

          if (include_faces && vec_f != Teuchos::null) {
            AmanziMesh::Entity_ID_List fnums0;
            std::vector<int> dirs;
            mesh_->cell_get_faces_and_dirs(c0, &fnums0, &dirs);

            for (unsigned int n = 0; n != fnums0.size(); ++n)
              *dcvo_[i]->os() << " " << Format_((*vec_f)[j][fnums0[n]]);
          }

          if (include_faces && vec_bf != Teuchos::null) {
            AmanziMesh::Entity_ID_List fnums0;
            std::vector<int> dirs;
            mesh_->cell_get_faces_and_dirs(c0, &fnums0, &dirs);

            for (unsigned int n = 0; n != fnums0.size(); ++n) {
              AmanziMesh::Entity_ID f = fnums0[n];
              AmanziMesh::Entity_ID bf =
                mesh_->exterior_face_map(true).LID(mesh_->face_map(true).GID(f));
              if (bf >= 0) *dcvo_[i]->os() << " " << Format_((*vec_bf)[j][bf]);
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
  std::stringstream formatstream;
  formatstream << "%_" << width_ << "." << precision_ << "g";
  std::string format = formatstream.str();

  for (int i = 0; i != dc_.size(); ++i) {
    AmanziMesh::Entity_ID c0 = dc_[i];
    AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
    Teuchos::OSTab tab = dcvo_[i]->getOSTab();

    if (dcvo_[i]->os_OK(verb_level_)) {
      *dcvo_[i]->os() << FormatHeader_("BCs", c0_gid);
      AmanziMesh::Entity_ID_List fnums0;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(c0, &fnums0, &dirs);

      for (unsigned int n = 0; n != fnums0.size(); ++n)
        *dcvo_[i]->os() << " " << flag[fnums0[n]] << "(" << Format_(data[fnums0[n]]) << ")";
      *dcvo_[i]->os() << std::endl;
    }
  }
}


// call MPI_Comm_Barrier to sync between writing steps
void
Debugger::Barrier()
{
  mesh_->get_comm()->Barrier();
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
