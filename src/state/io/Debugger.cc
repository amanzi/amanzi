/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include <boost/format.hpp>

#include "CompositeVector.hh"
#include "dbc.hh"
#include "errors.hh"

#include "Debugger.hh"

namespace Amanzi {

// Constructor
Debugger::Debugger(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                   std::string name, Teuchos::ParameterList& plist,
                   Teuchos::EVerbosityLevel verb_level)
  : mesh_(mesh),
    verb_level_(verb_level),
    precision_(10),
    width_(15),
    header_width_(20),
    cellnum_width_(5),
    decimal_width_(7)
{
  vo_ = Teuchos::rcp(new VerboseObject(name, plist));

  // cells to debug
  if (plist.isParameter("debug cells")) {
    Teuchos::Array<int> dcs = plist.get<Teuchos::Array<int>>("debug cells");
    for (Teuchos::Array<int>::const_iterator c = dcs.begin(); c != dcs.end();
         ++c) {
      AmanziMesh::Entity_ID lc = mesh->cell_map(false)->getLocalElement(*c);
      if (lc >= 0) {
        // include the LID
        dc_.push_back(lc);
        dc_gid_.push_back(*c);

        // make a verbose object for each case
        Teuchos::ParameterList vo_plist;
        vo_plist.sublist("verbose object");
        vo_plist.sublist("verbose object") = plist.sublist("verbose object");
        vo_plist.sublist("verbose object")
          .set("write on rank", mesh->get_comm()->getRank());
        dcvo_.push_back(
          Teuchos::rcp(new VerboseObject(*mesh_->get_comm(), name, vo_plist)));
      }
    }
  }

  // faces to debug
  if (plist.isParameter("debug faces")) {
    Teuchos::Array<int> dfs = plist.get<Teuchos::Array<int>>("debug faces");
    for (Teuchos::Array<int>::const_iterator f = dfs.begin(); f != dfs.end();
         ++f) {
      AmanziMesh::Entity_ID lf = mesh->face_map(true)->getLocalElement(*f);
      if (lf >= 0) {
        // debug the neighboring cells
        AmanziMesh::Entity_ID_View cells;
        mesh->face_get_cells(lf, AmanziMesh::Parallel_type::OWNED, cells);

        for (LO i = 0; i != cells.extent(0); ++i) {
          auto lc = cells[i];
          // include the LID
          dc_.push_back(lc);
          dc_gid_.push_back(mesh->cell_map(false)->getGlobalElement(lc));

          // make a verbose object for each case
          Teuchos::ParameterList vo_plist;
          vo_plist.sublist("verbose object");
          vo_plist.sublist("verbose object") = plist.sublist("verbose object");
          vo_plist.sublist("verbose object")
            .set("write on rank", mesh->get_comm()->getRank());
          dcvo_.push_back(Teuchos::rcp(
            new VerboseObject(*mesh_->get_comm(), name, vo_plist)));
        }
      }
    }
  }

  // formatting
  cellnum_width_ = plist.get<int>("cell number width", cellnum_width_);
  decimal_width_ = plist.get<int>("decimal width", decimal_width_);
  header_width_ = plist.get<int>("header width", header_width_);
  width_ = plist.get<int>("column width", width_);
  precision_ = plist.get<int>("precision", precision_);
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
      // formatstream2 << boost::format("%%%d.%df") % (width_ -
      // (decimal_width_-mag-1)) % decimal_width_;
      if (mag < 0) {
        formatstream2 << boost::format("%%%df") %
                           (width_ - (decimal_width_ - mag - 2));
      } else {
        formatstream2 << boost::format("%%%df") %
                           (width_ - (decimal_width_ - mag - 1));
      }

      datastream << boost::format(formatstream1.str()) % std::string("")
                 << boost::format(formatstream2.str()) % dat;

    } else { // sci format
      std::stringstream formatstream1;
      formatstream1 << boost::format("%%%ds") % (decimal_width_ - 5);

      std::stringstream formatstream2;
      // formatstream2 << boost::format("%%%d.%df") % (width_ -
      // (decimal_width_-mag-1)) % decimal_width_;
      formatstream2 << boost::format("%%%de") % (width_ - (decimal_width_ - 5));

      datastream << boost::format(formatstream1.str()) % std::string("")
                 << boost::format(formatstream2.str()) % dat;
    }
  }
  return datastream.str();
}

std::string
Debugger::FormatHeader_(std::string name, int c)
{
  std::string header_prefix(name);
  int header_prefix_width = header_width_ - cellnum_width_ - 4;
  if (header_prefix.size() > header_prefix_width) {
    header_prefix.erase(header_prefix_width);
  } else if (header_prefix.size() < header_prefix_width) {
    header_prefix.append(header_prefix_width - header_prefix.size(), ' ');
  }

  std::stringstream headerstream;
  headerstream << header_prefix << "(" << std::setw(cellnum_width_)
               << std::right << c << "): ";
  return headerstream.str();
}

// Write cell + face info
void
Debugger::WriteCellInfo(bool include_faces)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(verb_level_)) {
    *vo_->os() << "Debug Cells Information:" << std::endl;
  }

  for (int i = 0; i != dc_.size(); ++i) {
    AmanziMesh::Entity_ID c0 = dc_[i];
    AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
    Teuchos::OSTab itab = dcvo_[i]->getOSTab();

    AmanziGeometry::Point c0_centroid = mesh_->cell_centroid(c0);
    if (dcvo_[i]->os_OK(verb_level_)) {
      *dcvo_[i]->os() << "Cell c(" << c0_gid << ") centroid = " << c0_centroid
                      << std::endl;

      if (include_faces) {
        AmanziMesh::Entity_ID_View fnums0;
        Kokkos::View<int*> dirs;
        mesh_->cell_get_faces_and_dirs(c0, fnums0, dirs);

        if (dcvo_[i]->os_OK(verb_level_)) {
          for (unsigned int n = 0; n != fnums0.size(); ++n) {
            AmanziMesh::Entity_ID f_gid =
              mesh_->face_map(true)->getGlobalElement(fnums0[n]);
            AmanziGeometry::Point f_centroid = mesh_->face_centroid(fnums0[n]);
            *dcvo_[i]->os()
              << "  neighbor face(" << f_gid << ") [dir=" << dirs[n]
              << "] centroid = " << f_centroid << std::endl;
          }
        }
      }
    }
  }
}

// Write a vector individually.
void
Debugger::WriteVector(const std::string& name,
                      const CompositeVector& vec,
                      bool include_faces)
{
  CompositeVector::cMultiVectorView_type<DefaultHost> vec_c;
  if (vec.HasComponent("cell")) vec_c = vec.ViewComponent("cell", false);

  CompositeVector::cMultiVectorView_type<DefaultHost> vec_f;
  if (vec.HasComponent("face")) vec_f = vec.ViewComponent("face", true);

  for (int i = 0; i != dc_.size(); ++i) {
    AmanziMesh::Entity_ID c0 = dc_[i];
    AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
    Teuchos::OSTab tab = dcvo_[i]->getOSTab();

    if (dcvo_[i]->os_OK(verb_level_)) {
      *dcvo_[i]->os() << FormatHeader_(name, c0_gid);

      if (vec_c.extent(0)) *dcvo_[i]->os() << Format_(vec_c(c0, 0));

      if (include_faces && vec_f.extent(0)) {
        AmanziMesh::Entity_ID_View fnums0;
        Kokkos::View<int*> dirs;
        mesh_->cell_get_faces_and_dirs(c0, fnums0, dirs);

        for (unsigned int n = 0; n != fnums0.size(); ++n)
          if (fnums0[n] < vec_f.extent(0))
            *dcvo_[i]->os() << " " << Format_(vec_f(fnums0[n], 0));
      }
      *dcvo_[i]->os() << std::endl;
    }
  }
}

// Write list of vectors.
void
Debugger::WriteVectors(
  const std::vector<std::string>& names,
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

        CompositeVector::cMultiVectorView_type<DefaultHost> vec_c;
        if (vec->HasComponent("cell"))
          vec_c = vec->ViewComponent("cell", false);

        CompositeVector::cMultiVectorView_type<DefaultHost> vec_f;
        if (vec->HasComponent("face")) vec_f = vec->ViewComponent("face", true);

        *dcvo_[i]->os() << FormatHeader_(name, c0_gid);
        if (vec_c.extent(0)) *dcvo_[i]->os() << Format_(vec_c(c0, 0));

        if (include_faces && vec_f(c0, 0)) {
          AmanziMesh::Entity_ID_View fnums0;
          Kokkos::View<int*> dirs;
          mesh_->cell_get_faces_and_dirs(c0, fnums0, dirs);

          for (unsigned int n = 0; n != fnums0.size(); ++n)
            *dcvo_[i]->os() << " " << Format_(vec_f(fnums0[n], 0));
        }
        *dcvo_[i]->os() << std::endl;
      }
    }
  }
}

// // Write boundary condition data.
// void Debugger::WriteBoundaryConditions(const std::vector<int> &flag,
//                                        const std::vector<double> &data) {
//   std::stringstream formatstream;
//   formatstream << "%_" << width_ << "." << precision_ << "g";
//   std::string format = formatstream.str();

//   for (int i = 0; i != dc_.size(); ++i) {
//     AmanziMesh::Entity_ID c0 = dc_[i];
//     AmanziMesh::Entity_ID c0_gid = dc_gid_[i];
//     Teuchos::OSTab tab = dcvo_[i]->getOSTab();

//     if (dcvo_[i]->os_OK(verb_level_)) {
//       *dcvo_[i]->os() << FormatHeader_("BCs", c0_gid);
//       AmanziMesh::Entity_ID_View fnums0;
//       std::vector<int> dirs;
//       mesh_->cell_get_faces_and_dirs(c0, &fnums0, dirs);

//       for (unsigned int n = 0; n != fnums0.size(); ++n)
//         *dcvo_[i]->os() << " " << flag[fnums0[n]] << "("
//                         << Format_(data[fnums0[n]]) << ")";
//       *dcvo_[i]->os() << std::endl;
//     }
//   }
// }

// call MPI_Comm_Barrier to sync between writing steps
void
Debugger::Barrier()
{
  mesh_->get_comm()->barrier();
}

// write a line of ----
void
Debugger::WriteDivider()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(verb_level_))
    *vo_->os() << "------------------------------------------------------------"
                  "------------"
               << std::endl;
}

// Reverse order... get the VerboseObject for Entity_
Teuchos::RCP<VerboseObject>
Debugger::GetVerboseObject(AmanziMesh::Entity_ID id, int rank)
{
  std::vector<AmanziMesh::Entity_ID>::iterator loc =
    std::find(dc_.begin(), dc_.end(), id);
  if (loc == dc_.end()) {
    return Teuchos::null;
  } else {
    int index = loc - dc_.begin();
    return dcvo_[index];
  }
}

} // namespace Amanzi
