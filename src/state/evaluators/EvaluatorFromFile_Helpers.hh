/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Shared machinery for "from file" independent evaluators.
/*
  State

  TimeSeriesInterpolatorCV wraps TimeSeriesInterpolator.hh's
  paging/windowing/interpolation machinery to populate a CompositeVector.
  It is shared by EvaluatorIndependentFromFile and
  EvaluatorIndependentTensorFromFile.

  space is assumed to have exactly one component, used to size the page
  and to derive the on-disk entity-location string Input::readRange's
  attrs key expects.

  checkpoint_file (file layout: no time dimension on disk) and
  temporally_variable (semantics: does the value actually change in time)
  are independent -- a time-dependent file may hold a field that is
  constant in time (temporally_variable false, checkpoint_file false), but
  a checkpoint_file is always read as if constant in time.
*/

#pragma once

#include <limits>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Tpetra_Access.hpp"

#include "errors.hh"
#include "exceptions.hh"

#include "FileReaderUtils.hh"
#include "IOFileFormat.hh"
#include "Input.hh"
#include "InputFactory.hh"
#include "Mesh.hh"
#include "CompositeSpace.hh"
#include "CompositeVector.hh"
#include "TimeSeriesInterpolator.hh"
#include "ViewUtils.hh"

namespace Amanzi {
namespace EvaluatorFromFile_Helpers {

class TimeSeriesInterpolatorCV {
 public:
  using Page_type = Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace>;

  TimeSeriesInterpolatorCV(std::string filename,
                           std::string varname,
                           const Teuchos::RCP<const CompositeSpace>& space,
                           bool temporally_variable,
                           int page_size = 2,
                           Utils::IOFileFormat_kind file_format = Utils::IOFileFormat_kind::HDF5,
                           bool checkpoint_file = false)
    : filename_(std::move(filename)),
      varname_(std::move(varname)),
      space_(space),
      page_size_(page_size),
      file_format_(file_format)
  {
    AMANZI_ASSERT(page_size_ >= 2);
    AMANZI_ASSERT(space_->size() == 1);
    AMANZI_ASSERT(!checkpoint_file || !temporally_variable);

    compname_ = *space_->begin();
    ndofs_ = space_->getNumVectors(compname_);
    num_entities_ = space_->getComponentMap(compname_, false)->getLocalNumElements();
    locname_ = AmanziMesh::to_string(space_->getLocation(compname_));

    Teuchos::ParameterList file_plist("time series interpolator cv");
    file_plist.set<std::string>("file name", filename_);
    file_plist.set<std::string>("file format", to_string(file_format_));
    file_ = Teuchos::RCP<Input>(
      InputFactory::createForFile(file_plist, space_->getComm(), space_->getMesh()).release());

    if (checkpoint_file) {
      constant_ =
        Teuchos::rcp(new MultiVector_type(space_->getComponentMap(compname_, false), ndofs_));
      Teuchos::ParameterList attrs(varname_ + "." + locname_);
      attrs.set("always write subfield dof", true);
      file_->read(attrs, *constant_);
      return;
    }

    std::vector<double> times;
    if (temporally_variable) {
      // NOTE: all ranks read this.  Could do on rank 0 and broadcast.
      Kokkos::View<double*, Kokkos::HostSpace> times_view;
      readLocalArray(file_format_, space_->getComm(), filename_, "time", times_view);
      times = asVector(Kokkos::MeshView<double*, Kokkos::HostSpace>(times_view));

      for (std::size_t j = 1; j < times.size(); ++j) {
        if (times[j] <= times[j - 1]) {
          Errors::Message m;
          m << "TimeSeriesInterpolatorCV: times values in \"" << filename_
            << "\" are not strictly increasing";
          Exceptions::amanzi_throw(m);
        }
      }
    } else {
      times.push_back(std::numeric_limits<double>::max());
    }

    Page_type page("TimeSeriesInterpolatorCV::page_", page_size_, ndofs_, num_entities_);
    Page_type result("TimeSeriesInterpolatorCV::result_", 1, ndofs_, num_entities_);

    // does not capture this, so the closure stays valid regardless of
    // this object's own lifetime management.
    Input* file_ptr = file_.get();
    std::string attrs_name = varname_ + "." + locname_;
    int ndofs = ndofs_;
    int num_entities = num_entities_;
    interp_ = std::make_unique<TimeSeriesInterpolator<Page_type>>(
      times,
      page_size_,
      [file_ptr, page, attrs_name, ndofs, num_entities](int start, int count, int dest_block) {
        Teuchos::ParameterList attrs(attrs_name);
        Page_type dest(&page(dest_block, 0, 0), count, ndofs, num_entities);
        file_ptr->readRange(attrs, start, count, dest);
      },
      page,
      result);
  }

  // dest must share this object's space.
  void interpolate(double t, CompositeVector& dest)
  {
    auto dvec = dest.viewComponent<MemSpace_kind::HOST>(compname_, false);
    if (interp_ == nullptr) {
      // checkpoint_file case -- constant_ was read once at construction.
      auto cvec = constant_->getLocalViewHost(Tpetra::Access::ReadOnly);
      for (int j = 0; j != ndofs_; ++j) {
        for (int i = 0; i != num_entities_; ++i) { dvec(i, j) = cvec(i, j); }
      }
      return;
    }

    const auto& result = interp_->interpolate(t);
    for (int j = 0; j != ndofs_; ++j) {
      for (int i = 0; i != num_entities_; ++i) { dvec(i, j) = result(0, j, i); }
    }
  }

  Teuchos::RCP<const CompositeSpace> getSpace() const { return space_; }

 private:
  std::string filename_, varname_, compname_, locname_;
  Teuchos::RCP<const CompositeSpace> space_;
  int page_size_;
  Utils::IOFileFormat_kind file_format_;
  int ndofs_, num_entities_;

  Teuchos::RCP<Input> file_;
  std::unique_ptr<TimeSeriesInterpolator<Page_type>> interp_;
  // set instead of interp_ when checkpoint_file is true.
  Teuchos::RCP<MultiVector_type> constant_;
};

} // namespace EvaluatorFromFile_Helpers
} // namespace Amanzi
