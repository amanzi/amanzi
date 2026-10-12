/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Reads and writes NetCDF files for extruded-column (MPAS-style) meshes.
/*
  Template definitions for FileNetCDF -- see FileNetCDF.hh for the class
  documentation.
*/

#include "Tpetra_Access.hpp"
#include "Tpetra_Vector.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi {

namespace NetCDFDetails {

// a compile-time map from type to NetCDF's nc_type
template <typename T>
struct NC_TypeMap;
template <>
struct NC_TypeMap<int> {
  static const nc_type type = NC_INT;
};
template <>
struct NC_TypeMap<float> {
  static const nc_type type = NC_FLOAT;
};
template <>
struct NC_TypeMap<double> {
  static const nc_type type = NC_DOUBLE;
};
template <>
struct NC_TypeMap<bool> {
  static const nc_type type = NC_BYTE;
};

} // namespace NetCDFDetails


//
// Attributes -- scalar.  Always global (NC_GLOBAL).
//
template <typename scalar_type>
void
FileNetCDF::writeAttribute(const std::string& attr_name, scalar_type value)
{
  int ierr = nc_put_att(
    ncid_, NC_GLOBAL, attr_name.c_str(), NetCDFDetails::NC_TypeMap<scalar_type>::type, 1, &value);
  checkErr_(ierr, "writeAttribute(\"" + attr_name + "\")");
}


template <typename scalar_type>
scalar_type
FileNetCDF::readAttribute(const std::string& attr_name)
{
  scalar_type value;
  int ierr = nc_get_att(ncid_, NC_GLOBAL, attr_name.c_str(), &value);
  checkErr_(ierr, "readAttribute(\"" + attr_name + "\")");
  return value;
}


// Extending an UNLIMITED dimension needs NC_COLLECTIVE and every rank to
// call in; rank 0 passes the real count, every other rank passes 0.
template <typename view_type>
std::enable_if_t<view_type::rank == 1>
FileNetCDF::writeAttribute(const std::string& attr_name, view_type vals, TimeRange time_chunk)
{
  using scalar_type = typename view_type::non_const_value_type;
  const std::string& varname = attr_name;

  auto it = varids_.find(varname);
  int varid;
  if (it != varids_.end()) {
    varid = it->second;
  } else {
    int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
    if (ierr != NC_NOERR) {
      int dimid = getOrCreateDim_("time", 0, true);
      ierr = nc_def_var(
        ncid_, varname.c_str(), NetCDFDetails::NC_TypeMap<scalar_type>::type, 1, &dimid, &varid);
      checkErr_(ierr, "writeAttribute def_var(\"" + varname + "\")");
    }
    // NC_VAR_PAR_ACCESS is a property of this open file handle, not
    // persisted on disk, so it must be set whenever this instance first
    // sees the variable.
    ierr = nc_var_par_access(ncid_, varid, NC_COLLECTIVE);
    checkErr_(ierr, "writeAttribute par_access(\"" + varname + "\")");
    varids_[varname] = varid;
  }

  size_t start = time_chunk.first;
  size_t count = (comm_->getRank() == 0) ? (time_chunk.second - time_chunk.first) : 0;
  AMANZI_ASSERT(comm_->getRank() != 0 || vals.extent(0) == count);
  int ierr;
  if constexpr (std::is_same_v<scalar_type, double>) {
    ierr = nc_put_vara_double(ncid_, varid, &start, &count, vals.data());
  } else {
    ierr = nc_put_vara_int(ncid_, varid, &start, &count, vals.data());
  }
  checkErr_(ierr, "writeAttribute put_vara(\"" + varname + "\")");
}


template <typename view_type>
std::enable_if_t<view_type::rank == 1>
FileNetCDF::readAttribute(const std::string& attr_name, view_type vals, TimeRange time_chunk)
{
  using scalar_type = typename view_type::non_const_value_type;
  const std::string& varname = attr_name;
  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  checkErr_(ierr, "readAttribute inq_varid(\"" + varname + "\")");

  size_t start = time_chunk.first;
  size_t count = time_chunk.second - time_chunk.first;
  AMANZI_ASSERT(vals.extent(0) == count);
  if constexpr (std::is_same_v<scalar_type, double>) {
    ierr = nc_get_vara_double(ncid_, varid, &start, &count, vals.data());
  } else {
    ierr = nc_get_vara_int(ncid_, varid, &start, &count, vals.data());
  }
  checkErr_(ierr, "readAttribute get_vara(\"" + varname + "\")");
}


template <typename scalar_type>
void
FileNetCDF::writeAttribute(const std::string& attr_name, scalar_type value, std::size_t time_index)
{
  Kokkos::View<const scalar_type*, Kokkos::LayoutRight, Kokkos::HostSpace> v(&value, 1);
  writeAttribute(attr_name, v, TimeRange{ time_index, time_index + 1 });
}


template <typename scalar_type>
scalar_type
FileNetCDF::readAttribute(const std::string& attr_name, std::size_t time_index)
{
  scalar_type value;
  Kokkos::View<scalar_type*, Kokkos::LayoutRight, Kokkos::HostSpace> v(&value, 1);
  readAttribute(attr_name, v, TimeRange{ time_index, time_index + 1 });
  return value;
}


namespace NetCDFDetails {

// Builds the on-disk DimSpec list (outermost first) and the start/count
// arrays (same order) for one readDistributedView/writeDistributedView
// call.  view_extent(k) returns the k-th extent of view (the actual view
// passed by the caller, time axis included if it has one); rank is that
// view's rank, and dim_names has exactly `rank` entries naming every one
// of its axes, in view order -- including the time axis's name, when
// has_time is true, as dim_names[0].  The caller names every axis
// explicitly; this function has no opinion on what any axis represents.
// has_time/time_start/time_count describe that leading time
// (unlimited/record) axis, which is set directly from them rather than
// view_extent/decomposed_axis (time is never the decomposed axis).
template <typename ExtentFn>
void
buildDimsAndSlab(const Comm_ptr_type& comm,
                 int rank,
                 std::size_t decomposed_axis,
                 const std::vector<std::string>& dim_names,
                 ExtentFn view_extent,
                 bool has_time,
                 std::size_t time_start,
                 std::size_t time_count,
                 std::vector<FileNetCDF::DimSpec>& dims,
                 std::vector<size_t>& start,
                 std::vector<size_t>& count)
{
  AMANZI_ASSERT((int)dim_names.size() == rank);
  dims.clear();
  start.clear();
  count.clear();

  int k0 = 0;
  if (has_time) {
    dims.push_back({ dim_names[0], 0, true });
    start.push_back(time_start);
    count.push_back(time_count);
    k0 = 1;
  }

  for (int k = k0; k != rank; ++k) {
    std::size_t extent = view_extent(k);
    if (static_cast<std::size_t>(k) == decomposed_axis) {
      long long local = static_cast<long long>(extent);
      long long global = 0, inclusive_scan = 0;
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &local, &global);
      Teuchos::scan(*comm, Teuchos::REDUCE_SUM, 1, &local, &inclusive_scan);
      std::size_t offset = static_cast<std::size_t>(inclusive_scan - local);
      dims.push_back({ dim_names[k], static_cast<std::size_t>(global), false });
      start.push_back(offset);
      count.push_back(extent);
    } else {
      dims.push_back({ dim_names[k], extent, false });
      start.push_back(0);
      count.push_back(extent);
    }
  }
}

// Checks varid's actual on-disk dimension lengths against dims' expected
// lengths, throwing Errors::Message on any mismatch.
inline void
checkOnDiskDims(int ncid,
                int varid,
                const std::string& varname,
                const std::vector<FileNetCDF::DimSpec>& dims)
{
  int ndims;
  int ierr = nc_inq_varndims(ncid, varid, &ndims);
  if (ierr != NC_NOERR || ndims != static_cast<int>(dims.size())) {
    Errors::Message message;
    message << "FileNetCDF \"" << varname << "\" expected rank " << dims.size()
            << ", on-disk rank is " << ndims << ".";
    throw(message);
  }

  std::vector<int> dimids(ndims);
  ierr = nc_inq_vardimid(ncid, varid, dimids.data());
  if (ierr != NC_NOERR) {
    Errors::Message message;
    message << "FileNetCDF \"" << varname << "\" failed to query dimension ids.";
    throw(message);
  }

  for (int d = 0; d != ndims; ++d) {
    std::size_t on_disk_len;
    ierr = nc_inq_dimlen(ncid, dimids[d], &on_disk_len);
    if (ierr != NC_NOERR) {
      Errors::Message message;
      message << "FileNetCDF \"" << varname << "\" failed to query dimension " << d << "'s length.";
      throw(message);
    }
    if (dims[d].unlimited) continue; // "time" -- not bounds-checked here
    if (on_disk_len != dims[d].length) {
      Errors::Message message;
      message << "FileNetCDF \"" << varname << "\" on-disk dimension " << d << " (" << on_disk_len
              << ") does not match requested global length (" << dims[d].length << ").";
      throw(message);
    }
  }
}

} // namespace NetCDFDetails


template <typename view_type>
std::enable_if_t<!(view_type::rank == 2 &&
                   std::is_same_v<typename view_type::array_layout, Kokkos::LayoutLeft>)>
FileNetCDF::writeDistributedView(const std::string& varname,
                                 view_type view,
                                 std::size_t decomposed_axis,
                                 const std::vector<std::string>& dim_names)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  std::vector<DimSpec> dims;
  std::vector<size_t> start, count;
  NetCDFDetails::buildDimsAndSlab(
    comm_,
    view_type::rank,
    decomposed_axis,
    dim_names,
    [&](int k) { return view.extent(k); },
    /*has_time=*/false,
    0,
    0,
    dims,
    start,
    count);

  std::vector<int> dimids;
  int varid =
    getOrCreateVariable_(varname, dims, dimids, NetCDFDetails::NC_TypeMap<scalar_type>::type);

  int ierr;
  if constexpr (std::is_same_v<scalar_type, double>) {
    ierr = nc_put_vara_double(ncid_, varid, start.data(), count.data(), view.data());
  } else {
    ierr = nc_put_vara_int(ncid_, varid, start.data(), count.data(), view.data());
  }
  checkErr_(ierr, "writeDistributedView(\"" + varname + "\")");
}


template <typename view_type>
std::enable_if_t<!(view_type::rank == 2 &&
                   std::is_same_v<typename view_type::array_layout, Kokkos::LayoutLeft>)>
FileNetCDF::readDistributedView(const std::string& varname,
                                view_type view,
                                std::size_t decomposed_axis,
                                const std::vector<std::string>& dim_names)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  checkErr_(ierr, "readDistributedView inq_varid(\"" + varname + "\")");

  std::vector<DimSpec> dims;
  std::vector<size_t> start, count;
  NetCDFDetails::buildDimsAndSlab(
    comm_,
    view_type::rank,
    decomposed_axis,
    dim_names,
    [&](int k) { return view.extent(k); },
    /*has_time=*/false,
    0,
    0,
    dims,
    start,
    count);
  NetCDFDetails::checkOnDiskDims(ncid_, varid, varname, dims);

  if constexpr (std::is_same_v<scalar_type, double>) {
    ierr = nc_get_vara_double(ncid_, varid, start.data(), count.data(), view.data());
  } else {
    ierr = nc_get_vara_int(ncid_, varid, start.data(), count.data(), view.data());
  }
  checkErr_(ierr, "readDistributedView(\"" + varname + "\")");
}


//
// Rank-2, LayoutLeft, decomposed_axis fixed at 0 -- see FileNetCDF.hh.
// dim_names[0] names the on-disk NDOFS axis, dim_names[1] the decomposed
// N axis -- that on-disk ([NDOFS, N]) order, not the view's own order.
//
template <typename view_type>
std::enable_if_t<view_type::rank == 2 &&
                 std::is_same_v<typename view_type::array_layout, Kokkos::LayoutLeft>>
FileNetCDF::writeDistributedView(const std::string& varname,
                                 view_type view,
                                 const std::vector<std::string>& dim_names)
{
  AMANZI_ASSERT(dim_names.size() == 2);
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;

  long long local = static_cast<long long>(view.extent(0));
  long long global = 0, inclusive_scan = 0;
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &local, &global);
  Teuchos::scan(*comm_, Teuchos::REDUCE_SUM, 1, &local, &inclusive_scan);
  std::size_t offset = static_cast<std::size_t>(inclusive_scan - local);

  std::vector<DimSpec> dims{ { dim_names[0], (std::size_t)view.extent(1), false },
                             { dim_names[1], static_cast<std::size_t>(global), false } };
  std::vector<size_t> start{ 0, offset };
  std::vector<size_t> count{ (size_t)view.extent(1), (size_t)view.extent(0) };

  std::vector<int> dimids;
  int varid =
    getOrCreateVariable_(varname, dims, dimids, NetCDFDetails::NC_TypeMap<scalar_type>::type);

  int ierr;
  if constexpr (std::is_same_v<scalar_type, double>) {
    ierr = nc_put_vara_double(ncid_, varid, start.data(), count.data(), view.data());
  } else {
    ierr = nc_put_vara_int(ncid_, varid, start.data(), count.data(), view.data());
  }
  checkErr_(ierr, "writeDistributedView(\"" + varname + "\")");
}


template <typename view_type>
std::enable_if_t<view_type::rank == 2 &&
                 std::is_same_v<typename view_type::array_layout, Kokkos::LayoutLeft>>
FileNetCDF::readDistributedView(const std::string& varname,
                                view_type view,
                                const std::vector<std::string>& dim_names)
{
  AMANZI_ASSERT(dim_names.size() == 2);
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;

  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  checkErr_(ierr, "readDistributedView inq_varid(\"" + varname + "\")");

  long long local = static_cast<long long>(view.extent(0));
  long long global = 0, inclusive_scan = 0;
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &local, &global);
  Teuchos::scan(*comm_, Teuchos::REDUCE_SUM, 1, &local, &inclusive_scan);
  std::size_t offset = static_cast<std::size_t>(inclusive_scan - local);

  std::vector<DimSpec> dims{ { dim_names[0], (std::size_t)view.extent(1), false },
                             { dim_names[1], static_cast<std::size_t>(global), false } };
  NetCDFDetails::checkOnDiskDims(ncid_, varid, varname, dims);

  std::vector<size_t> start{ 0, offset };
  std::vector<size_t> count{ (size_t)view.extent(1), (size_t)view.extent(0) };

  if constexpr (std::is_same_v<scalar_type, double>) {
    ierr = nc_get_vara_double(ncid_, varid, start.data(), count.data(), view.data());
  } else {
    ierr = nc_get_vara_int(ncid_, varid, start.data(), count.data(), view.data());
  }
  checkErr_(ierr, "readDistributedView(\"" + varname + "\")");
}


template <typename view_type>
void
FileNetCDF::writeDistributedView(const std::string& varname,
                                 view_type view,
                                 std::size_t decomposed_axis,
                                 const std::vector<std::string>& dim_names,
                                 TimeRange time_chunk)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  std::size_t ntimes = time_chunk.second - time_chunk.first;
  AMANZI_ASSERT(view.extent(0) == ntimes);

  // view's axis 0 is time; decomposed_axis and dim_names are both
  // relative to the full view (time axis included, dim_names[0] names it).
  std::vector<DimSpec> dims;
  std::vector<size_t> start, count;
  NetCDFDetails::buildDimsAndSlab(
    comm_,
    view_type::rank,
    decomposed_axis,
    dim_names,
    [&](int k) { return view.extent(k); },
    /*has_time=*/true,
    time_chunk.first,
    ntimes,
    dims,
    start,
    count);

  std::vector<int> dimids;
  int varid =
    getOrCreateVariable_(varname, dims, dimids, NetCDFDetails::NC_TypeMap<scalar_type>::type);

  int ierr;
  if constexpr (std::is_same_v<scalar_type, double>) {
    ierr = nc_put_vara_double(ncid_, varid, start.data(), count.data(), view.data());
  } else {
    ierr = nc_put_vara_int(ncid_, varid, start.data(), count.data(), view.data());
  }
  checkErr_(ierr, "writeDistributedView(\"" + varname + "\")");
}


template <typename view_type>
void
FileNetCDF::readDistributedView(const std::string& varname,
                                view_type view,
                                std::size_t decomposed_axis,
                                const std::vector<std::string>& dim_names,
                                TimeRange time_chunk)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  std::size_t ntimes = time_chunk.second - time_chunk.first;
  AMANZI_ASSERT(view.extent(0) == ntimes);

  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  checkErr_(ierr, "readDistributedView inq_varid(\"" + varname + "\")");

  // See writeDistributedView's TimeRange overload above re: decomposed_axis.
  std::vector<DimSpec> dims;
  std::vector<size_t> start, count;
  NetCDFDetails::buildDimsAndSlab(
    comm_,
    view_type::rank,
    decomposed_axis,
    dim_names,
    [&](int k) { return view.extent(k); },
    /*has_time=*/true,
    time_chunk.first,
    ntimes,
    dims,
    start,
    count);
  NetCDFDetails::checkOnDiskDims(ncid_, varid, varname, dims);

  if constexpr (std::is_same_v<scalar_type, double>) {
    ierr = nc_get_vara_double(ncid_, varid, start.data(), count.data(), view.data());
  } else {
    ierr = nc_get_vara_int(ncid_, varid, start.data(), count.data(), view.data());
  }
  checkErr_(ierr, "readDistributedView(\"" + varname + "\")");
}


// Folds a single time_index to a TimeRange of length 1, via an unmanaged
// view with an extra leading extent-1 axis over the same memory.
// dim_names is relative to that folded (one-more-axis) view, same as it
// is for the TimeRange overload below -- dim_names[0] names the time
// axis, unchanged by the fold.
template <typename view_type>
void
FileNetCDF::writeDistributedView(const std::string& varname,
                                 view_type view,
                                 std::size_t decomposed_axis,
                                 const std::vector<std::string>& dim_names,
                                 std::size_t time_index)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  if constexpr (view_type::rank == 1) {
    Kokkos::View<const scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace> v(
      view.data(), 1, view.extent(0));
    writeDistributedView(
      varname, v, decomposed_axis + 1, dim_names, TimeRange{ time_index, time_index + 1 });
  } else if constexpr (view_type::rank == 2) {
    Kokkos::View<const scalar_type***, Kokkos::LayoutRight, Kokkos::HostSpace> v(
      view.data(), 1, view.extent(0), view.extent(1));
    writeDistributedView(
      varname, v, decomposed_axis + 1, dim_names, TimeRange{ time_index, time_index + 1 });
  } else {
    static_assert(view_type::rank == 3);
    Kokkos::View<const scalar_type****, Kokkos::LayoutRight, Kokkos::HostSpace> v(
      view.data(), 1, view.extent(0), view.extent(1), view.extent(2));
    writeDistributedView(
      varname, v, decomposed_axis + 1, dim_names, TimeRange{ time_index, time_index + 1 });
  }
}


template <typename view_type>
void
FileNetCDF::readDistributedView(const std::string& varname,
                                view_type view,
                                std::size_t decomposed_axis,
                                const std::vector<std::string>& dim_names,
                                std::size_t time_index)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  if constexpr (view_type::rank == 1) {
    Kokkos::View<scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace> v(
      view.data(), 1, view.extent(0));
    readDistributedView(
      varname, v, decomposed_axis + 1, dim_names, TimeRange{ time_index, time_index + 1 });
  } else if constexpr (view_type::rank == 2) {
    Kokkos::View<scalar_type***, Kokkos::LayoutRight, Kokkos::HostSpace> v(
      view.data(), 1, view.extent(0), view.extent(1));
    readDistributedView(
      varname, v, decomposed_axis + 1, dim_names, TimeRange{ time_index, time_index + 1 });
  } else {
    static_assert(view_type::rank == 3);
    Kokkos::View<scalar_type****, Kokkos::LayoutRight, Kokkos::HostSpace> v(
      view.data(), 1, view.extent(0), view.extent(1), view.extent(2));
    readDistributedView(
      varname, v, decomposed_axis + 1, dim_names, TimeRange{ time_index, time_index + 1 });
  }
}


template <typename scalar_type>
int
FileNetCDF::getOrCreateLocalVariable_(const std::string& varname,
                                      const std::vector<DimSpec>& dims,
                                      std::vector<int>& dimids_out,
                                      int par_access)
{
  dimids_out.clear();
  for (const auto& dim : dims)
    dimids_out.push_back(getOrCreateDim_(dim.name, dim.length, dim.unlimited));

  auto it = varids_.find(varname);
  if (it != varids_.end()) return it->second;

  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  if (ierr != NC_NOERR) {
    ierr = nc_def_var(ncid_,
                      varname.c_str(),
                      NetCDFDetails::NC_TypeMap<scalar_type>::type,
                      (int)dimids_out.size(),
                      dimids_out.data(),
                      &varid);
    checkErr_(ierr, "getOrCreateLocalVariable_(\"" + varname + "\")");
  }

  ierr = nc_var_par_access(ncid_, varid, par_access);
  checkErr_(ierr, "getOrCreateLocalVariable_ par_access(\"" + varname + "\")");

  varids_[varname] = varid;
  return varid;
}


namespace NetCDFDetails {

template <typename scalar_type>
int
putVara_(int ncid, int varid, const size_t* start, const size_t* count, const scalar_type* buf)
{
  if constexpr (std::is_same_v<scalar_type, double>) {
    return nc_put_vara_double(ncid, varid, start, count, buf);
  } else if constexpr (std::is_same_v<scalar_type, int>) {
    return nc_put_vara_int(ncid, varid, start, count, buf);
  } else {
    static_assert(std::is_same_v<scalar_type, double> || std::is_same_v<scalar_type, int>,
                  "FileNetCDF::LocalView supports only double or int scalar_type.");
  }
}

template <typename scalar_type>
int
getVara_(int ncid, int varid, const size_t* start, const size_t* count, scalar_type* buf)
{
  if constexpr (std::is_same_v<scalar_type, double>) {
    return nc_get_vara_double(ncid, varid, start, count, buf);
  } else if constexpr (std::is_same_v<scalar_type, int>) {
    return nc_get_vara_int(ncid, varid, start, count, buf);
  } else {
    static_assert(std::is_same_v<scalar_type, double> || std::is_same_v<scalar_type, int>,
                  "FileNetCDF::LocalView supports only double or int scalar_type.");
  }
}

template <typename view_type>
std::enable_if_t<view_type::rank == 1>
resizeToDims_(view_type& view, const std::vector<std::size_t>& dims)
{
  Kokkos::resize(view, dims[0]);
}

template <typename view_type>
std::enable_if_t<view_type::rank == 2>
resizeToDims_(view_type& view, const std::vector<std::size_t>& dims)
{
  Kokkos::resize(view, dims[0], dims[1]);
}

template <typename view_type>
std::enable_if_t<view_type::rank == 3>
resizeToDims_(view_type& view, const std::vector<std::size_t>& dims)
{
  Kokkos::resize(view, dims[0], dims[1], dims[2]);
}

} // namespace NetCDFDetails


//
// LocalView -- scalar, 1D, or 2D.
//
template <typename scalar_type>
std::enable_if_t<std::is_arithmetic_v<scalar_type>>
FileNetCDF::writeLocalView(const std::string& varname, scalar_type value)
{
  std::vector<DimSpec> dims{ { varname + "_dim0", 1, false } };
  std::vector<int> dimids;
  int varid = getOrCreateLocalVariable_<scalar_type>(varname, dims, dimids, NC_COLLECTIVE);
  size_t start[] = { 0 };
  size_t count[] = { (comm_->getRank() == 0) ? (size_t)1 : (size_t)0 };
  int ierr = NetCDFDetails::putVara_(ncid_, varid, start, count, &value);
  checkErr_(ierr, "writeLocalView put_vara(\"" + varname + "\")");
}


template <typename scalar_type>
std::enable_if_t<std::is_arithmetic_v<scalar_type>, scalar_type>
FileNetCDF::readLocalView(const std::string& varname)
{
  scalar_type value{};
  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  checkErr_(ierr, "readLocalView inq_varid(\"" + varname + "\")");
  size_t start[] = { 0 };
  size_t count[] = { 1 };
  ierr = NetCDFDetails::getVara_(ncid_, varid, start, count, &value);
  checkErr_(ierr, "readLocalView get_vara(\"" + varname + "\")");
  return value;
}


template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2)>
FileNetCDF::writeLocalView(const std::string& varname, view_type view)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  bool have_data = (comm_->getRank() == 0);

  std::vector<DimSpec> dims;
  std::vector<size_t> start, count;
  for (int d = 0; d != view_type::rank; ++d) {
    dims.push_back({ varname + "_dim" + std::to_string(d), (std::size_t)view.extent(d), false });
    start.push_back(0);
    count.push_back(have_data ? view.extent(d) : 0);
  }
  std::vector<int> dimids;
  int varid = getOrCreateLocalVariable_<scalar_type>(varname, dims, dimids, NC_COLLECTIVE);
  int ierr = NetCDFDetails::putVara_(ncid_, varid, start.data(), count.data(), view.data());
  checkErr_(ierr, "writeLocalView put_vara(\"" + varname + "\")");
}


template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2)>
FileNetCDF::readLocalView(const std::string& varname, view_type& view)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;

  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  checkErr_(ierr, "readLocalView inq_varid(\"" + varname + "\")");

  int ndims;
  ierr = nc_inq_varndims(ncid_, varid, &ndims);
  checkErr_(ierr, "readLocalView inq_varndims(\"" + varname + "\")");
  if (ndims != view_type::rank) {
    Errors::Message message;
    message << "FileNetCDF::readLocalView \"" << varname << "\" expected rank " << view_type::rank
            << ", on-disk rank is " << ndims << ".";
    throw(message);
  }

  std::vector<int> dimids(ndims);
  ierr = nc_inq_vardimid(ncid_, varid, dimids.data());
  checkErr_(ierr, "readLocalView inq_vardimid(\"" + varname + "\")");
  std::vector<std::size_t> on_disk_dims(ndims);
  for (int d = 0; d != ndims; ++d) {
    ierr = nc_inq_dimlen(ncid_, dimids[d], &on_disk_dims[d]);
    checkErr_(ierr, "readLocalView inq_dimlen(\"" + varname + "\")");
  }

  if (view.size() == 0) {
    NetCDFDetails::resizeToDims_(view, on_disk_dims);
  } else {
    for (int d = 0; d != ndims; ++d) {
      if (on_disk_dims[d] != static_cast<std::size_t>(view.extent(d))) {
        Errors::Message message;
        message << "FileNetCDF::readLocalView \"" << varname << "\" on-disk dimension " << d << " ("
                << on_disk_dims[d] << ") does not match view's dimension ("
                << static_cast<std::size_t>(view.extent(d)) << ").";
        throw(message);
      }
    }
  }

  std::vector<size_t> start(ndims, 0), count(ndims);
  for (int d = 0; d != ndims; ++d) count[d] = view.extent(d);
  ierr = NetCDFDetails::getVara_(ncid_, varid, start.data(), count.data(), view.data());
  checkErr_(ierr, "readLocalView get_vara(\"" + varname + "\")");
}


//
// LocalView, paged by a single time_index -- on-disk array has one more
// (leading) dimension than value.
//
template <typename scalar_type>
std::enable_if_t<std::is_arithmetic_v<scalar_type>>
FileNetCDF::writeLocalView(const std::string& varname, scalar_type value, std::size_t time_index)
{
  std::vector<DimSpec> dims{ { varname + "_time", 0, true } };
  std::vector<int> dimids;
  int varid = getOrCreateLocalVariable_<scalar_type>(varname, dims, dimids, NC_COLLECTIVE);

  size_t start[] = { static_cast<size_t>(time_index) };
  size_t count[] = { (comm_->getRank() == 0) ? (size_t)1 : (size_t)0 };
  int ierr = NetCDFDetails::putVara_(ncid_, varid, start, count, &value);
  checkErr_(ierr, "writeLocalView put_vara(\"" + varname + "\")");
}


template <typename scalar_type>
std::enable_if_t<std::is_arithmetic_v<scalar_type>, scalar_type>
FileNetCDF::readLocalView(const std::string& varname, std::size_t time_index)
{
  scalar_type value{};
  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  checkErr_(ierr, "readLocalView inq_varid(\"" + varname + "\")");
  size_t start[] = { static_cast<size_t>(time_index) };
  size_t count[] = { 1 };
  ierr = NetCDFDetails::getVara_(ncid_, varid, start, count, &value);
  checkErr_(ierr, "readLocalView get_vara(\"" + varname + "\")");
  return value;
}


// On-disk array has rank view_type::rank + 1 (dimension 0 is the time
// axis, which view never represents -- only dims 1..ndims-1 correspond to
// view's own dims 0..rank-1).
template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2)>
FileNetCDF::writeLocalView(const std::string& varname, view_type view, std::size_t time_index)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;

  std::vector<DimSpec> dims{ { varname + "_time", 0, true } };
  for (int d = 0; d != view_type::rank; ++d)
    dims.push_back({ varname + "_dim" + std::to_string(d), (std::size_t)view.extent(d), false });
  std::vector<int> dimids;
  int varid = getOrCreateLocalVariable_<scalar_type>(varname, dims, dimids, NC_COLLECTIVE);

  std::vector<size_t> start(view_type::rank + 1, 0), count(view_type::rank + 1);
  start[0] = static_cast<size_t>(time_index);
  bool have_data = (comm_->getRank() == 0);
  count[0] = have_data ? 1 : 0;
  for (int d = 0; d != view_type::rank; ++d) count[d + 1] = have_data ? view.extent(d) : 0;
  int ierr = NetCDFDetails::putVara_(ncid_, varid, start.data(), count.data(), view.data());
  checkErr_(ierr, "writeLocalView put_vara(\"" + varname + "\")");
}


template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2)>
FileNetCDF::readLocalView(const std::string& varname, view_type& view, std::size_t time_index)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  constexpr int ndims = view_type::rank + 1;

  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  checkErr_(ierr, "readLocalView inq_varid(\"" + varname + "\")");

  int on_disk_ndims;
  ierr = nc_inq_varndims(ncid_, varid, &on_disk_ndims);
  checkErr_(ierr, "readLocalView inq_varndims(\"" + varname + "\")");
  if (on_disk_ndims != ndims) {
    Errors::Message message;
    message << "FileNetCDF::readLocalView \"" << varname << "\" expected rank " << ndims
            << ", on-disk rank is " << on_disk_ndims << ".";
    throw(message);
  }

  std::vector<int> dimids(ndims);
  ierr = nc_inq_vardimid(ncid_, varid, dimids.data());
  checkErr_(ierr, "readLocalView inq_vardimid(\"" + varname + "\")");
  std::vector<std::size_t> on_disk_dims(ndims);
  for (int d = 0; d != ndims; ++d) {
    ierr = nc_inq_dimlen(ncid_, dimids[d], &on_disk_dims[d]);
    checkErr_(ierr, "readLocalView inq_dimlen(\"" + varname + "\")");
  }

  if (view.size() == 0) {
    std::vector<std::size_t> trailing_dims(on_disk_dims.begin() + 1, on_disk_dims.end());
    NetCDFDetails::resizeToDims_(view, trailing_dims);
  } else {
    for (int d = 1; d != ndims; ++d) {
      if (on_disk_dims[d] != static_cast<std::size_t>(view.extent(d - 1))) {
        Errors::Message message;
        message << "FileNetCDF::readLocalView \"" << varname << "\" on-disk dimension " << d << " ("
                << on_disk_dims[d] << ") does not match view's dimension ("
                << static_cast<std::size_t>(view.extent(d - 1)) << ").";
        throw(message);
      }
    }
  }

  std::vector<size_t> start(ndims, 0), count(ndims);
  start[0] = static_cast<size_t>(time_index);
  count[0] = 1;
  for (int d = 1; d != ndims; ++d) count[d] = view.extent(d - 1);
  ierr = NetCDFDetails::getVara_(ncid_, varid, start.data(), count.data(), view.data());
  checkErr_(ierr, "readLocalView get_vara(\"" + varname + "\")");
}


//
// LocalView, paged by a time_chunk range -- on-disk array has the same
// rank as value, value's 0th dimension being the page length.
//
template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2 || view_type::rank == 3)>
FileNetCDF::writeLocalView(const std::string& varname, view_type view, TimeRange time_chunk)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  constexpr int ndims = view_type::rank;
  std::size_t ntimes = time_chunk.second - time_chunk.first;
  bool have_data = (comm_->getRank() == 0);
  AMANZI_ASSERT(!have_data || view.extent(0) == ntimes);

  std::vector<DimSpec> dims{ { varname + "_time", 0, true } };
  for (int d = 1; d != ndims; ++d)
    dims.push_back({ varname + "_dim" + std::to_string(d), (std::size_t)view.extent(d), false });
  std::vector<int> dimids;
  int varid = getOrCreateLocalVariable_<scalar_type>(varname, dims, dimids, NC_COLLECTIVE);

  std::vector<size_t> start(ndims, 0), count(ndims);
  start[0] = static_cast<size_t>(time_chunk.first);
  count[0] = have_data ? ntimes : 0;
  for (int d = 1; d != ndims; ++d) count[d] = have_data ? view.extent(d) : 0;
  int ierr = NetCDFDetails::putVara_(ncid_, varid, start.data(), count.data(), view.data());
  checkErr_(ierr, "writeLocalView put_vara(\"" + varname + "\")");
}


template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2 || view_type::rank == 3)>
FileNetCDF::readLocalView(const std::string& varname, view_type& view, TimeRange time_chunk)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  constexpr int ndims = view_type::rank;
  std::size_t ntimes = time_chunk.second - time_chunk.first;

  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  checkErr_(ierr, "readLocalView inq_varid(\"" + varname + "\")");

  int on_disk_ndims;
  ierr = nc_inq_varndims(ncid_, varid, &on_disk_ndims);
  checkErr_(ierr, "readLocalView inq_varndims(\"" + varname + "\")");
  if (on_disk_ndims != ndims) {
    Errors::Message message;
    message << "FileNetCDF::readLocalView \"" << varname << "\" expected rank " << ndims
            << ", on-disk rank is " << on_disk_ndims << ".";
    throw(message);
  }

  std::vector<int> dimids(ndims);
  ierr = nc_inq_vardimid(ncid_, varid, dimids.data());
  checkErr_(ierr, "readLocalView inq_vardimid(\"" + varname + "\")");
  std::vector<std::size_t> on_disk_dims(ndims);
  for (int d = 0; d != ndims; ++d) {
    ierr = nc_inq_dimlen(ncid_, dimids[d], &on_disk_dims[d]);
    checkErr_(ierr, "readLocalView inq_dimlen(\"" + varname + "\")");
  }

  if (view.size() == 0) {
    std::vector<std::size_t> resize_dims(on_disk_dims);
    resize_dims[0] = ntimes;
    NetCDFDetails::resizeToDims_(view, resize_dims);
  } else {
    AMANZI_ASSERT(view.extent(0) == ntimes);
    for (int d = 1; d != ndims; ++d) {
      if (on_disk_dims[d] != static_cast<std::size_t>(view.extent(d))) {
        Errors::Message message;
        message << "FileNetCDF::readLocalView \"" << varname << "\" on-disk dimension " << d << " ("
                << on_disk_dims[d] << ") does not match view's dimension ("
                << static_cast<std::size_t>(view.extent(d)) << ").";
        throw(message);
      }
    }
  }

  std::vector<size_t> start(ndims, 0), count(ndims);
  start[0] = static_cast<size_t>(time_chunk.first);
  count[0] = ntimes;
  for (int d = 1; d != ndims; ++d) count[d] = view.extent(d);
  ierr = NetCDFDetails::getVara_(ncid_, varid, start.data(), count.data(), view.data());
  checkErr_(ierr, "readLocalView get_vara(\"" + varname + "\")");
}


//
// Vector -- delegates to the 1D writeDistributedView/readDistributedView.
//
template <typename scalar_type>
void
FileNetCDF::writeVector(const std::string& varname,
                        const Vector_type_<scalar_type>& vec,
                        const std::vector<std::string>& dim_names)
{
  auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);
  writeDistributedView(varname, Kokkos::subview(vecv, Kokkos::ALL(), 0), 0, dim_names);
}


template <typename scalar_type>
void
FileNetCDF::readVector(const std::string& varname,
                       Vector_type_<scalar_type>& vec,
                       const std::vector<std::string>& dim_names)
{
  auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
  readDistributedView(varname, Kokkos::subview(vecv, Kokkos::ALL(), 0), 0, dim_names);
}


//
// STD::STRING -- declared here, defined in FileNetCDF.cc.
//
template <>
void
FileNetCDF::writeAttribute<std::string>(const std::string& attr_name, std::string value);

template <>
std::string
FileNetCDF::readAttribute<std::string>(const std::string& attr_name);

} // namespace Amanzi
