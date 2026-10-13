/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "errors.hh"
#include "exceptions.hh"
#include <memory>

#include "ViewUtils.hh"

#include "FunctionBilinearAndTime.hh"

namespace Amanzi {

FunctionBilinearAndTime::FunctionBilinearAndTime(const std::string& filename,
                                                 const std::string& time_header,
                                                 const std::string& row_header,
                                                 const std::string& row_coordinate,
                                                 const std::string& col_header,
                                                 const std::string& col_coordinate,
                                                 const std::string& val_header,
                                                 Form_kind form,
                                                 int page_size,
                                                 Utils::IOFileFormat_kind file_format,
                                                 const Comm_ptr_type& comm)
  : filename_(filename),
    time_header_(time_header),
    row_header_(row_header),
    col_header_(col_header),
    val_header_(val_header),
    row_coordinate_(row_coordinate),
    col_coordinate_(col_coordinate),
    form_(form),
    page_size_(page_size),
    file_format_(file_format),
    comm_(comm)
{
  if (row_coordinate_ == "x")
    row_index_ = 1;
  else if (row_coordinate_ == "y")
    row_index_ = 2;
  else if (row_coordinate_ == "z")
    row_index_ = 3;

  if (col_coordinate_ == "x")
    col_index_ = 1;
  else if (col_coordinate_ == "y")
    col_index_ = 2;
  else if (col_coordinate_ == "z")
    col_index_ = 3;

  initialize_();
}


FunctionBilinearAndTime::FunctionBilinearAndTime(const FunctionBilinearAndTime& other)
  : filename_(other.filename_),
    time_header_(other.time_header_),
    row_header_(other.row_header_),
    col_header_(other.col_header_),
    val_header_(other.val_header_),
    row_coordinate_(other.row_coordinate_),
    col_coordinate_(other.col_coordinate_),
    row_index_(other.row_index_),
    col_index_(other.col_index_),
    form_(other.form_),
    page_size_(other.page_size_),
    file_format_(other.file_format_),
    comm_(other.comm_)
{
  // file_, interp_, and their associated views cannot be cheaply copied
  // (file_ is an open FileReader; interp_'s load_page closure captures
  // this) -- rebuild from the stored construction parameters instead,
  // exactly as a fresh construction would.
  initialize_();
}


// (Re-)builds file_, row_/col_, and interp_ from the stored construction
// parameters.  Reads the row/column grids and the time series once; the
// first page is loaded by TimeSeriesInterpolator's own constructor.
void
FunctionBilinearAndTime::initialize_()
{
  file_ = createFileReader(file_format_, comm_, filename_);

  Kokkos::View<double*, Kokkos::HostSpace> row_h, col_h, times_h;
  file_->read1D(row_header_, row_h);
  file_->read1D(col_header_, col_h);
  file_->read1D(time_header_, times_h);

  nx_ = row_h.extent(0);
  ny_ = col_h.extent(0);
  row_ = asDualView(row_h);
  col_ = asDualView(col_h);

  std::vector<double> times(times_h.extent(0));
  for (std::size_t i = 0; i != times_h.extent(0); ++i) times[i] = times_h(i);

  Kokkos::resize(page_, page_size_, nx_, ny_);
  Kokkos::resize(result_, 1, nx_, ny_);
  Kokkos::resize(result_d_, nx_, ny_);

  auto form = (form_ == Form_kind::CONSTANT) ?
                TimeSeriesInterpolator<Page_type>::Form_kind::CONSTANT :
                TimeSeriesInterpolator<Page_type>::Form_kind::LINEAR;

  // captures file_'s raw pointer, page_ (via Kokkos' reference-counted
  // View copy -- the lambda's page and this->page_ share the same
  // underlying buffer), and a copy of val_header_ -- never this itself,
  // so the closure stays valid across the copy constructor's
  // reassignment of file_/page_ (each FunctionBilinearAndTime owns its
  // own interp_, never shared).
  Utils::FileReader* file_ptr = file_.get();
  Page_type page = page_;
  std::string val_header = val_header_;
  interp_ = std::make_unique<TimeSeriesInterpolator<Page_type>>(
    times,
    page_size_,
    [file_ptr, page, val_header](int start, int count, int dest_block) {
      Page_type dest(&page(dest_block, 0, 0), count, page.extent(1), page.extent(2));
      file_ptr->readPagedChunk(
        val_header,
        Utils::FileReader::TimeRange{ (std::size_t)start, (std::size_t)(start + count) },
        dest);
    },
    page_,
    result_,
    form);
}


double
FunctionBilinearAndTime::operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& x) const
{
  const auto& result = interp_->interpolate(x(0, 0));
  Kokkos::View<const double**, Kokkos::LayoutRight, Kokkos::HostSpace> result_slice(
    &result(0, 0, 0), nx_, ny_);
  auto f = Impl::FunctionBilinearFunctor(
    row_.view_host(), col_.view_host(), row_index_, col_index_, result_slice, x);
  return f(0);
}


void
FunctionBilinearAndTime::apply(
  const Kokkos::View<const double**>& in,
  Kokkos::View<double*>& out,
  const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
{
  auto t_view = Kokkos::create_mirror_view_and_copy(Amanzi::DefaultHostMemorySpace(),
                                                    Kokkos::subview(in, 0, 0));
  const auto& result = interp_->interpolate(t_view());
  Kokkos::View<const double**, Kokkos::LayoutRight, Kokkos::HostSpace> result_h(
    &result(0, 0, 0), nx_, ny_);
  Kokkos::deep_copy(result_d_, result_h);

  auto f = Impl::FunctionBilinearFunctor(
    row_.view_device(), col_.view_device(), row_index_, col_index_, result_d_, in);

  if (ids) {
    auto ids_loc = *ids;
    Kokkos::parallel_for(
      "FunctionBilinearAndTime::apply1", ids_loc.extent(0), KOKKOS_LAMBDA(const int i) {
        out(ids_loc(i)) = f(i);
      });
  } else {
    Kokkos::parallel_for(
      "FunctionBilinearAndTime::apply2", in.extent(1), KOKKOS_LAMBDA(const int i) {
        out(i) = f(i);
      });
  }
}

} // namespace Amanzi
