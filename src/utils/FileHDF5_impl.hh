/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! reads and writes HDF5 files.
/*
  reads and writes HDF5 files via parallelIO library.
*/

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Access.hpp"
#include "Tpetra_Vector.hpp"

#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace Utils {

namespace IODetails {

// horrendous class to deal with const issues in ascemio
class DangerousString {
 public:
  DangerousString(std::string str) : str_(std::move(str)) {}
  DangerousString(const char* chars) : DangerousString(std::string(chars)) {}
  char* c_str() { return const_cast<char*>(str_.c_str()); }

 private:
  std::string str_;
};

// a compile-time map from type to ASCEMIO::datatype_t
template <typename T>
struct PIO_DatatypeMap;
template <>
struct PIO_DatatypeMap<int> {
  static const datatype_t type = PIO_INTEGER;
};
template <>
struct PIO_DatatypeMap<float> {
  static const datatype_t type = PIO_FLOAT;
};
template <>
struct PIO_DatatypeMap<double> {
  static const datatype_t type = PIO_DOUBLE;
};
template <>
struct PIO_DatatypeMap<long> {
  static const datatype_t type = PIO_LONG;
};
template <>
struct PIO_DatatypeMap<bool> {
  static const datatype_t type = PIO_BYTE;
};

// a compile-time map from type to HDF5's native hid_t, for the plain
// (non-ASCEMIO) HDF5 calls LocalView uses.
template <typename T>
struct H5_DatatypeMap;
template <>
struct H5_DatatypeMap<int> {
  static hid_t type() { return H5T_NATIVE_INT; }
};
template <>
struct H5_DatatypeMap<double> {
  static hid_t type() { return H5T_NATIVE_DOUBLE; }
};


template <typename scalar_type>
void
writeLocalDataset_(hid_t file_id,
                   const std::string& full_path,
                   const scalar_type* buf,
                   const std::vector<hsize_t>& dims)
{
  if (H5Lexists(file_id, full_path.c_str(), H5P_DEFAULT) <= 0) {
    hid_t filespace = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
    hid_t dset = H5Dcreate2(file_id,
                            full_path.c_str(),
                            H5_DatatypeMap<scalar_type>::type(),
                            filespace,
                            H5P_DEFAULT,
                            H5P_DEFAULT,
                            H5P_DEFAULT);
    AMANZI_ASSERT(dset >= 0);
    H5Dwrite(dset, H5_DatatypeMap<scalar_type>::type(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    H5Dclose(dset);
    H5Sclose(filespace);
    return;
  }

  hid_t dset = H5Dopen2(file_id, full_path.c_str(), H5P_DEFAULT);
  AMANZI_ASSERT(dset >= 0);
  hid_t filespace = H5Dget_space(dset);
  int ndims = H5Sget_simple_extent_ndims(filespace);
  if (ndims != static_cast<int>(dims.size())) {
    H5Sclose(filespace);
    H5Dclose(dset);
    Errors::Message message;
    message << "FileHDF5::writeLocalView \"" << full_path << "\" expected rank " << dims.size()
            << ", on-disk rank is " << ndims << ".";
    throw(message);
  }
  std::vector<hsize_t> on_disk_dims(ndims);
  H5Sget_simple_extent_dims(filespace, on_disk_dims.data(), nullptr);
  for (int d = 0; d != ndims; ++d) {
    if (on_disk_dims[d] != dims[d]) {
      H5Sclose(filespace);
      H5Dclose(dset);
      Errors::Message message;
      message << "FileHDF5::writeLocalView \"" << full_path << "\" on-disk dimension " << d << " ("
              << static_cast<std::size_t>(on_disk_dims[d]) << ") does not match view's dimension ("
              << static_cast<std::size_t>(dims[d]) << ").";
      throw(message);
    }
  }
  H5Dwrite(dset, H5_DatatypeMap<scalar_type>::type(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  H5Sclose(filespace);
  H5Dclose(dset);
}

template <typename view_type>
std::enable_if_t<view_type::rank == 1>
resizeToDims_(view_type& view, const std::vector<hsize_t>& dims)
{
  Kokkos::resize(view, dims[0]);
}

template <typename view_type>
std::enable_if_t<view_type::rank == 2>
resizeToDims_(view_type& view, const std::vector<hsize_t>& dims)
{
  Kokkos::resize(view, dims[0], dims[1]);
}

template <typename view_type>
std::enable_if_t<view_type::rank == 3>
resizeToDims_(view_type& view, const std::vector<hsize_t>& dims)
{
  Kokkos::resize(view, dims[0], dims[1], dims[2]);
}

template <typename scalar_type>
void
readLocalScalarDataset_(hid_t file_id, const std::string& full_path, scalar_type* buf)
{
  hid_t dset = H5Dopen2(file_id, full_path.c_str(), H5P_DEFAULT);
  AMANZI_ASSERT(dset >= 0);
  H5Dread(dset, H5_DatatypeMap<scalar_type>::type(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  H5Dclose(dset);
}

// view may be default-constructed (size 0), in which case it is resized
// to the on-disk shape -- either all ranks must pass an empty view, or
// all must pass one already sized to match.
template <typename view_type>
void
readLocalDataset_(hid_t file_id, const std::string& full_path, view_type& view)
{
  using scalar_type = typename view_type::non_const_value_type;
  constexpr int ndims = view_type::rank;

  hid_t dset = H5Dopen2(file_id, full_path.c_str(), H5P_DEFAULT);
  AMANZI_ASSERT(dset >= 0);
  hid_t filespace = H5Dget_space(dset);
  int on_disk_ndims = H5Sget_simple_extent_ndims(filespace);
  if (on_disk_ndims != ndims) {
    H5Sclose(filespace);
    H5Dclose(dset);
    Errors::Message message;
    message << "FileHDF5::readLocalView \"" << full_path << "\" expected rank " << ndims
            << ", on-disk rank is " << on_disk_ndims << ".";
    throw(message);
  }
  std::vector<hsize_t> on_disk_dims(ndims);
  H5Sget_simple_extent_dims(filespace, on_disk_dims.data(), nullptr);

  if (view.size() == 0) {
    resizeToDims_(view, on_disk_dims);
  } else {
    for (int d = 0; d != ndims; ++d) {
      if (on_disk_dims[d] != static_cast<hsize_t>(view.extent(d))) {
        H5Sclose(filespace);
        H5Dclose(dset);
        Errors::Message message;
        message << "FileHDF5::readLocalView \"" << full_path << "\" on-disk dimension " << d << " ("
                << static_cast<std::size_t>(on_disk_dims[d])
                << ") does not match view's dimension (" << static_cast<std::size_t>(view.extent(d))
                << ").";
        throw(message);
      }
    }
  }
  H5Dread(dset, H5_DatatypeMap<scalar_type>::type(), H5S_ALL, H5S_ALL, H5P_DEFAULT, view.data());
  H5Sclose(filespace);
  H5Dclose(dset);
}


template <typename view_type>
void
writeLocalViewPage_(hid_t file_id,
                    const std::string& var_path,
                    view_type view,
                    std::size_t time_index)
{
  using scalar_type = typename view_type::non_const_value_type;
  constexpr int ndims = view_type::rank + 1;

  hid_t dset;
  if (H5Lexists(file_id, var_path.c_str(), H5P_DEFAULT) <= 0) {
    std::vector<hsize_t> dims(ndims), maxdims(ndims);
    dims[0] = 1;
    maxdims[0] = H5S_UNLIMITED;
    for (int d = 1; d != ndims; ++d) {
      dims[d] = view.extent(d - 1);
      maxdims[d] = view.extent(d - 1);
    }
    hid_t filespace = H5Screate_simple(ndims, dims.data(), maxdims.data());
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, ndims, dims.data());
    dset = H5Dcreate2(file_id,
                      var_path.c_str(),
                      H5_DatatypeMap<scalar_type>::type(),
                      filespace,
                      H5P_DEFAULT,
                      plist,
                      H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(filespace);
  } else {
    dset = H5Dopen2(file_id, var_path.c_str(), H5P_DEFAULT);
  }
  AMANZI_ASSERT(dset >= 0);

  hsize_t needed = static_cast<hsize_t>(time_index) + 1;
  hid_t filespace = H5Dget_space(dset);
  std::vector<hsize_t> cur_dims(ndims);
  H5Sget_simple_extent_dims(filespace, cur_dims.data(), nullptr);
  H5Sclose(filespace);
  if (cur_dims[0] < needed) {
    std::vector<hsize_t> new_dims(cur_dims);
    new_dims[0] = needed;
    H5Dset_extent(dset, new_dims.data());
  }

  filespace = H5Dget_space(dset);
  std::vector<hsize_t> start(ndims, 0), count(ndims);
  start[0] = static_cast<hsize_t>(time_index);
  count[0] = 1;
  for (int d = 1; d != ndims; ++d) count[d] = view.extent(d - 1);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);
  hid_t memspace = H5Screate_simple(ndims, count.data(), nullptr);
  H5Dwrite(
    dset, H5_DatatypeMap<scalar_type>::type(), memspace, filespace, H5P_DEFAULT, view.data());
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset);
}

// On-disk array has rank view_type::rank + 1 (dimension 0 is the time
// axis, which view never represents -- only dims 1..ndims-1 correspond to
// view's own dims 0..rank-1).  If view is empty, it is resized here from
// the on-disk trailing dims (1..ndims-1), not dimension 0.
template <typename view_type>
void
readLocalViewPage_(hid_t file_id,
                   const std::string& var_path,
                   view_type& view,
                   std::size_t time_index)
{
  using scalar_type = typename view_type::non_const_value_type;
  constexpr int ndims = view_type::rank + 1;

  hid_t dset = H5Dopen2(file_id, var_path.c_str(), H5P_DEFAULT);
  AMANZI_ASSERT(dset >= 0);
  hid_t filespace = H5Dget_space(dset);
  int on_disk_ndims = H5Sget_simple_extent_ndims(filespace);
  if (on_disk_ndims != ndims) {
    H5Sclose(filespace);
    H5Dclose(dset);
    Errors::Message message;
    message << "FileHDF5::readLocalView \"" << var_path << "\" expected rank " << ndims
            << ", on-disk rank is " << on_disk_ndims << ".";
    throw(message);
  }
  std::vector<hsize_t> on_disk_dims(ndims);
  H5Sget_simple_extent_dims(filespace, on_disk_dims.data(), nullptr);

  if (view.size() == 0) {
    std::vector<hsize_t> trailing_dims(on_disk_dims.begin() + 1, on_disk_dims.end());
    resizeToDims_(view, trailing_dims);
  } else {
    for (int d = 1; d != ndims; ++d) {
      if (on_disk_dims[d] != static_cast<hsize_t>(view.extent(d - 1))) {
        H5Sclose(filespace);
        H5Dclose(dset);
        Errors::Message message;
        message << "FileHDF5::readLocalView \"" << var_path << "\" on-disk dimension " << d << " ("
                << static_cast<std::size_t>(on_disk_dims[d])
                << ") does not match view's dimension ("
                << static_cast<std::size_t>(view.extent(d - 1)) << ").";
        throw(message);
      }
    }
  }

  std::vector<hsize_t> start(ndims, 0), count(ndims);
  start[0] = static_cast<hsize_t>(time_index);
  count[0] = 1;
  for (int d = 1; d != ndims; ++d) count[d] = view.extent(d - 1);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);
  hid_t memspace = H5Screate_simple(ndims, count.data(), nullptr);
  H5Dread(dset, H5_DatatypeMap<scalar_type>::type(), memspace, filespace, H5P_DEFAULT, view.data());
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset);
}


template <typename view_type>
void
writeLocalViewChunk_(hid_t file_id,
                     const std::string& var_path,
                     view_type view,
                     FileHDF5::TimeRange time_chunk)
{
  using scalar_type = typename view_type::non_const_value_type;
  constexpr int ndims = view_type::rank;
  std::size_t ntimes = time_chunk.second - time_chunk.first;
  AMANZI_ASSERT(view.extent(0) == ntimes);

  hid_t dset;
  if (H5Lexists(file_id, var_path.c_str(), H5P_DEFAULT) <= 0) {
    std::vector<hsize_t> dims(ndims), maxdims(ndims);
    dims[0] = time_chunk.second;
    maxdims[0] = H5S_UNLIMITED;
    for (int d = 1; d != ndims; ++d) {
      dims[d] = view.extent(d);
      maxdims[d] = view.extent(d);
    }
    hid_t filespace = H5Screate_simple(ndims, dims.data(), maxdims.data());
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    std::vector<hsize_t> chunkdims(dims);
    chunkdims[0] = std::max<hsize_t>(ntimes, 1);
    H5Pset_chunk(plist, ndims, chunkdims.data());
    dset = H5Dcreate2(file_id,
                      var_path.c_str(),
                      H5_DatatypeMap<scalar_type>::type(),
                      filespace,
                      H5P_DEFAULT,
                      plist,
                      H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(filespace);
  } else {
    dset = H5Dopen2(file_id, var_path.c_str(), H5P_DEFAULT);
  }
  AMANZI_ASSERT(dset >= 0);

  hsize_t needed = static_cast<hsize_t>(time_chunk.second);
  hid_t filespace = H5Dget_space(dset);
  std::vector<hsize_t> cur_dims(ndims);
  H5Sget_simple_extent_dims(filespace, cur_dims.data(), nullptr);
  H5Sclose(filespace);
  if (cur_dims[0] < needed) {
    std::vector<hsize_t> new_dims(cur_dims);
    new_dims[0] = needed;
    H5Dset_extent(dset, new_dims.data());
  }

  filespace = H5Dget_space(dset);
  std::vector<hsize_t> start(ndims, 0), count(ndims);
  start[0] = static_cast<hsize_t>(time_chunk.first);
  count[0] = ntimes;
  for (int d = 1; d != ndims; ++d) count[d] = view.extent(d);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);
  hid_t memspace = H5Screate_simple(ndims, count.data(), nullptr);
  H5Dwrite(
    dset, H5_DatatypeMap<scalar_type>::type(), memspace, filespace, H5P_DEFAULT, view.data());
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset);
}

// On-disk array has the same rank as view, dimension 0 being the full
// time extent written so far; the read only ever touches the time_chunk
// sub-range of it, so view's dimension 0 is ntimes = time_chunk's length,
// not on_disk_dims[0].  If view is empty, it is resized to [ntimes,
// on_disk_dims[1], ..., on_disk_dims[rank-1]].
template <typename view_type>
void
readLocalViewChunk_(hid_t file_id,
                    const std::string& var_path,
                    view_type& view,
                    FileHDF5::TimeRange time_chunk)
{
  using scalar_type = typename view_type::non_const_value_type;
  constexpr int ndims = view_type::rank;
  std::size_t ntimes = time_chunk.second - time_chunk.first;

  hid_t dset = H5Dopen2(file_id, var_path.c_str(), H5P_DEFAULT);
  AMANZI_ASSERT(dset >= 0);
  hid_t filespace = H5Dget_space(dset);
  int on_disk_ndims = H5Sget_simple_extent_ndims(filespace);
  if (on_disk_ndims != ndims) {
    H5Sclose(filespace);
    H5Dclose(dset);
    Errors::Message message;
    message << "FileHDF5::readLocalView \"" << var_path << "\" expected rank " << ndims
            << ", on-disk rank is " << on_disk_ndims << ".";
    throw(message);
  }
  std::vector<hsize_t> on_disk_dims(ndims);
  H5Sget_simple_extent_dims(filespace, on_disk_dims.data(), nullptr);

  if (view.size() == 0) {
    std::vector<hsize_t> resize_dims(on_disk_dims);
    resize_dims[0] = ntimes;
    resizeToDims_(view, resize_dims);
  } else {
    AMANZI_ASSERT(view.extent(0) == ntimes);
    for (int d = 1; d != ndims; ++d) {
      if (on_disk_dims[d] != static_cast<hsize_t>(view.extent(d))) {
        H5Sclose(filespace);
        H5Dclose(dset);
        Errors::Message message;
        message << "FileHDF5::readLocalView \"" << var_path << "\" on-disk dimension " << d << " ("
                << static_cast<std::size_t>(on_disk_dims[d])
                << ") does not match view's dimension (" << static_cast<std::size_t>(view.extent(d))
                << ").";
        throw(message);
      }
    }
  }

  std::vector<hsize_t> start(ndims, 0), count(ndims);
  start[0] = static_cast<hsize_t>(time_chunk.first);
  count[0] = ntimes;
  for (int d = 1; d != ndims; ++d) count[d] = view.extent(d);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);
  hid_t memspace = H5Screate_simple(ndims, count.data(), nullptr);
  H5Dread(dset, H5_DatatypeMap<scalar_type>::type(), memspace, filespace, H5P_DEFAULT, view.data());
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset);
}

} // namespace IODetails


template <typename Fn>
void
FileHDF5::runOnRank0AndBroadcastErrors_(Fn&& fn)
{
  int failed = 0;
  std::string what;
  if (comm_->getRank() == 0) {
    try {
      fn();
    } catch (Errors::Message& e) {
      failed = 1;
      what = e.what();
    }
  }

  Teuchos::broadcast(*comm_, 0, 1, &failed);
  if (failed) {
    int len = static_cast<int>(what.size());
    Teuchos::broadcast(*comm_, 0, 1, &len);
    what.resize(len);
    if (len > 0) Teuchos::broadcast(*comm_, 0, len, &what[0]);
    throw Errors::Message(what);
  }
}


//
// Attributes -- scalar, at a group given by path.
//
template <typename scalar_type>
void
FileHDF5::writeAttribute(const std::string& name, const std::string& path, scalar_type value)
{
  ensureGroupExists_(path);
  IODetails::DangerousString attr_name(name);
  IODetails::DangerousString h5path(path);
  parallelIO_write_simple_attr(attr_name.c_str(),
                               (void*)&value,
                               IODetails::PIO_DatatypeMap<scalar_type>::type,
                               data_file_,
                               h5path.c_str(),
                               &IOgroup_);
}


template <typename scalar_type>
scalar_type
FileHDF5::readAttribute(const std::string& name, const std::string& path)
{
  IODetails::DangerousString attr_name(name);
  IODetails::DangerousString h5path(path);

  scalar_type* loc_value;
  parallelIO_read_simple_attr(attr_name.c_str(),
                              reinterpret_cast<void**>(&loc_value),
                              IODetails::PIO_DatatypeMap<scalar_type>::type,
                              data_file_,
                              h5path.c_str(),
                              &IOgroup_);
  scalar_type value = *loc_value;
  free(loc_value);
  return value;
}


//
// DistributedView -- 1D, decomposed across ranks.
//
template <typename view_type>
std::enable_if_t<view_type::rank == 1>
FileHDF5::writeDistributedView(const std::string& name, const std::string& path, view_type view)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  std::string var_path = joinPath_(path, name);
  ensureGroupExists_(path);

  GO globallength = view.extent(0);
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &globallength, &globallength);

  int globaldims[] = { static_cast<int>(globallength), 1 };
  int localdims[] = { static_cast<int>(view.extent(0)), 1 };

  IODetails::DangerousString full_h5path(var_path);
  parallelIO_write_dataset((void*)view.data(),
                           IODetails::PIO_DatatypeMap<scalar_type>::type,
                           2,
                           globaldims,
                           localdims,
                           data_file_,
                           full_h5path.c_str(),
                           &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
}


template <typename view_type>
std::enable_if_t<view_type::rank == 1>
FileHDF5::readDistributedView(const std::string& name, const std::string& path, view_type view)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  std::string var_path = joinPath_(path, name);
  IODetails::DangerousString full_h5path(var_path);

  int ndims;
  parallelIO_get_dataset_ndims(&ndims, data_file_, full_h5path.c_str(), &IOgroup_);
  if (ndims < 0) {
    Errors::Message message("FileHDF5::readDistributedView data has negative dimension.");
    throw(message);
  }

  std::vector<int> globaldims(ndims), localdims(ndims);
  parallelIO_get_dataset_dims(globaldims.data(), data_file_, full_h5path.c_str(), &IOgroup_);

  GO globallength = view.extent(0);
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &globallength, &globallength);
  if (globaldims[0] != static_cast<int>(globallength)) {
    Errors::Message message;
    message << "FileHDF5::readDistributedView \"" << var_path << "\" on-disk global length ("
            << globaldims[0] << ") does not match view's global length (" << globallength << ").";
    throw(message);
  }

  localdims[0] = static_cast<int>(view.extent(0));
  for (int d = 1; d != ndims; ++d) localdims[d] = globaldims[d];

  parallelIO_read_dataset((void*)view.data(),
                          IODetails::PIO_DatatypeMap<scalar_type>::type,
                          ndims,
                          globaldims.data(),
                          localdims.data(),
                          data_file_,
                          full_h5path.c_str(),
                          &IOgroup_,
                          NONUNIFORM_CONTIGUOUS_READ);
}


//
// Tpetra::Vector, delegates to the 1D DistributedView.
//
template <typename scalar_type>
void
FileHDF5::writeVector(const std::string& name,
                      const std::string& path,
                      const Vector_type_<scalar_type>& vec)
{
  auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);
  writeDistributedView(name, path, Kokkos::subview(vecv, Kokkos::ALL(), 0));
}


template <typename scalar_type>
void
FileHDF5::readVector(const std::string& name,
                     const std::string& path,
                     Vector_type_<scalar_type>& vec)
{
  auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
  readDistributedView(name, path, Kokkos::subview(vecv, Kokkos::ALL(), 0));
}


//
// DistributedView -- 2D [N, NDOFS], decomposed along dimension 0 (N);
// NDOFS is the fastest-varying dimension.
//
template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 view_type::rank == 2>
FileHDF5::writeDistributedView(const std::string& name, const std::string& path, view_type view)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  std::string var_path = joinPath_(path, name);
  ensureGroupExists_(path);

  GO globallength = view.extent(0);
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &globallength, &globallength);

  int globaldims[] = { static_cast<int>(globallength), static_cast<int>(view.extent(1)) };
  int localdims[] = { static_cast<int>(view.extent(0)), static_cast<int>(view.extent(1)) };

  IODetails::DangerousString full_h5path(var_path);
  parallelIO_write_dataset((void*)view.data(),
                           IODetails::PIO_DatatypeMap<scalar_type>::type,
                           2,
                           globaldims,
                           localdims,
                           data_file_,
                           full_h5path.c_str(),
                           &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
}


template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 view_type::rank == 2>
FileHDF5::readDistributedView(const std::string& name, const std::string& path, view_type view)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  std::string var_path = joinPath_(path, name);
  IODetails::DangerousString full_h5path(var_path);

  int ndims;
  parallelIO_get_dataset_ndims(&ndims, data_file_, full_h5path.c_str(), &IOgroup_);
  if (ndims != 2) {
    Errors::Message message;
    message << "FileHDF5::readDistributedView \"" << var_path << "\" expected rank 2, got " << ndims
            << ".";
    throw(message);
  }

  std::vector<int> globaldims(ndims), localdims(ndims);
  parallelIO_get_dataset_dims(globaldims.data(), data_file_, full_h5path.c_str(), &IOgroup_);

  GO globallength = view.extent(0);
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &globallength, &globallength);
  if (globaldims[0] != static_cast<int>(globallength) ||
      globaldims[1] != static_cast<int>(view.extent(1))) {
    Errors::Message message;
    message << "FileHDF5::readDistributedView \"" << var_path << "\" on-disk shape ("
            << globaldims[0] << ", " << globaldims[1] << ") does not match view's shape ("
            << globallength << ", " << view.extent(1) << ").";
    throw(message);
  }

  localdims[0] = static_cast<int>(view.extent(0));
  localdims[1] = static_cast<int>(view.extent(1));

  parallelIO_read_dataset((void*)view.data(),
                          IODetails::PIO_DatatypeMap<scalar_type>::type,
                          ndims,
                          globaldims.data(),
                          localdims.data(),
                          data_file_,
                          full_h5path.c_str(),
                          &IOgroup_,
                          NONUNIFORM_CONTIGUOUS_READ);
}


template <typename view_type>
void
FileHDF5::resizeIfEmptyAndBroadcast_(view_type& view, bool was_empty)
{
  if (was_empty) {
    std::vector<hsize_t> dims(view_type::rank, 0);
    if (comm_->getRank() == 0)
      for (int d = 0; d != view_type::rank; ++d) dims[d] = view.extent(d);
    Teuchos::broadcast(*comm_, 0, static_cast<int>(view_type::rank), dims.data());
    if (comm_->getRank() != 0) IODetails::resizeToDims_(view, dims);
  }
  Teuchos::broadcast(*comm_, 0, static_cast<int>(view.span()), view.data());
}


//
// LocalView -- scalar, 1D, or 2D.
//
template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2)>
FileHDF5::writeLocalView(const std::string& name, const std::string& path, view_type view)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  using scalar_type = typename view_type::non_const_value_type;
  ensureGroupExists_(path);
  std::string var_path = joinPath_(path, name);

  runOnRank0AndBroadcastErrors_([&]() {
    std::vector<hsize_t> dims(view_type::rank);
    for (int d = 0; d != view_type::rank; ++d) dims[d] = view.extent(d);
    IODetails::writeLocalDataset_(rawFileId_(), var_path, view.data(), dims);
  });
}


template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2)>
FileHDF5::readLocalView(const std::string& name, const std::string& path, view_type& view)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  std::string var_path = joinPath_(path, name);

  bool was_empty = (view.size() == 0);
  runOnRank0AndBroadcastErrors_(
    [&]() { IODetails::readLocalDataset_(rawFileId_(), var_path, view); });
  resizeIfEmptyAndBroadcast_(view, was_empty);
}


//
// LocalView, scalar -- no view, just a single value.
//
template <typename scalar_type>
std::enable_if_t<std::is_arithmetic_v<scalar_type>>
FileHDF5::writeLocalView(const std::string& name, const std::string& path, scalar_type value)
{
  ensureGroupExists_(path);
  std::string var_path = joinPath_(path, name);
  runOnRank0AndBroadcastErrors_([&]() {
    IODetails::writeLocalDataset_(rawFileId_(), var_path, &value, std::vector<hsize_t>{ 1 });
  });
}


template <typename scalar_type>
std::enable_if_t<std::is_arithmetic_v<scalar_type>, scalar_type>
FileHDF5::readLocalView(const std::string& name, const std::string& path)
{
  std::string var_path = joinPath_(path, name);
  scalar_type value{};
  runOnRank0AndBroadcastErrors_(
    [&]() { IODetails::readLocalScalarDataset_(rawFileId_(), var_path, &value); });
  Teuchos::broadcast(*comm_, 0, 1, &value);
  return value;
}


//
// LocalView, paged by a single time_index -- on-disk array has one more
// (leading) dimension than value.
//
template <typename scalar_type>
std::enable_if_t<std::is_arithmetic_v<scalar_type>>
FileHDF5::writeLocalView(const std::string& name,
                         const std::string& path,
                         scalar_type value,
                         std::size_t time_index)
{
  ensureGroupExists_(path);
  std::string var_path = joinPath_(path, name);
  if (comm_->getRank() == 0) {
    hid_t file_id = rawFileId_();
    hid_t dset;
    if (H5Lexists(file_id, var_path.c_str(), H5P_DEFAULT) <= 0) {
      hsize_t dims[] = { 1 };
      hsize_t maxdims[] = { H5S_UNLIMITED };
      hid_t filespace = H5Screate_simple(1, dims, maxdims);
      hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
      hsize_t chunkdims[] = { 1 };
      H5Pset_chunk(plist, 1, chunkdims);
      dset = H5Dcreate2(file_id,
                        var_path.c_str(),
                        IODetails::H5_DatatypeMap<scalar_type>::type(),
                        filespace,
                        H5P_DEFAULT,
                        plist,
                        H5P_DEFAULT);
      H5Pclose(plist);
      H5Sclose(filespace);
    } else {
      dset = H5Dopen2(file_id, var_path.c_str(), H5P_DEFAULT);
    }
    AMANZI_ASSERT(dset >= 0);

    hsize_t needed = static_cast<hsize_t>(time_index) + 1;
    hid_t filespace = H5Dget_space(dset);
    hsize_t cur_dims[1];
    H5Sget_simple_extent_dims(filespace, cur_dims, nullptr);
    H5Sclose(filespace);
    if (cur_dims[0] < needed) {
      hsize_t new_dims[] = { needed };
      H5Dset_extent(dset, new_dims);
    }

    filespace = H5Dget_space(dset);
    hsize_t start[] = { static_cast<hsize_t>(time_index) };
    hsize_t count[] = { 1 };
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);
    hid_t memspace = H5Screate_simple(1, count, nullptr);
    H5Dwrite(dset,
             IODetails::H5_DatatypeMap<scalar_type>::type(),
             memspace,
             filespace,
             H5P_DEFAULT,
             &value);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset);
  }
}


template <typename scalar_type>
std::enable_if_t<std::is_arithmetic_v<scalar_type>, scalar_type>
FileHDF5::readLocalView(const std::string& name, const std::string& path, std::size_t time_index)
{
  std::string var_path = joinPath_(path, name);
  scalar_type value{};
  if (comm_->getRank() == 0) {
    hid_t file_id = rawFileId_();
    hid_t dset = H5Dopen2(file_id, var_path.c_str(), H5P_DEFAULT);
    AMANZI_ASSERT(dset >= 0);
    hid_t filespace = H5Dget_space(dset);
    hsize_t start[] = { static_cast<hsize_t>(time_index) };
    hsize_t count[] = { 1 };
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);
    hid_t memspace = H5Screate_simple(1, count, nullptr);
    H5Dread(dset,
            IODetails::H5_DatatypeMap<scalar_type>::type(),
            memspace,
            filespace,
            H5P_DEFAULT,
            &value);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset);
  }
  Teuchos::broadcast(*comm_, 0, 1, &value);
  return value;
}


template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2)>
FileHDF5::writeLocalView(const std::string& name,
                         const std::string& path,
                         view_type view,
                         std::size_t time_index)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  ensureGroupExists_(path);
  std::string var_path = joinPath_(path, name);
  runOnRank0AndBroadcastErrors_(
    [&]() { IODetails::writeLocalViewPage_(rawFileId_(), var_path, view, time_index); });
}


template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2)>
FileHDF5::readLocalView(const std::string& name,
                        const std::string& path,
                        view_type& view,
                        std::size_t time_index)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  std::string var_path = joinPath_(path, name);

  bool was_empty = (view.size() == 0);
  runOnRank0AndBroadcastErrors_(
    [&]() { IODetails::readLocalViewPage_(rawFileId_(), var_path, view, time_index); });
  resizeIfEmptyAndBroadcast_(view, was_empty);
}


//
// LocalView, paged by a time_chunk range -- on-disk array has the same
// rank as value, value's 0th dimension being the page length.
//
template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2 || view_type::rank == 3)>
FileHDF5::writeLocalView(const std::string& name,
                         const std::string& path,
                         view_type view,
                         TimeRange time_chunk)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  ensureGroupExists_(path);
  std::string var_path = joinPath_(path, name);
  runOnRank0AndBroadcastErrors_(
    [&]() { IODetails::writeLocalViewChunk_(rawFileId_(), var_path, view, time_chunk); });
}


template <typename view_type>
std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                 (view_type::rank == 1 || view_type::rank == 2 || view_type::rank == 3)>
FileHDF5::readLocalView(const std::string& name,
                        const std::string& path,
                        view_type& view,
                        TimeRange time_chunk)
{
  AMANZI_ASSERT(view.span_is_contiguous());
  std::string var_path = joinPath_(path, name);

  bool was_empty = (view.size() == 0);
  runOnRank0AndBroadcastErrors_(
    [&]() { IODetails::readLocalViewChunk_(rawFileId_(), var_path, view, time_chunk); });
  resizeIfEmptyAndBroadcast_(view, was_empty);
}


//
// STD::STRING
//
template <>
void
FileHDF5::writeAttribute<std::string>(const std::string& name,
                                      const std::string& path,
                                      std::string value);

template <>
std::string
FileHDF5::readAttribute<std::string>(const std::string& name, const std::string& path);

} // namespace Utils
} // namespace Amanzi
