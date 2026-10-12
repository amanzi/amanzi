/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! reads and writes HDF5 files.
/*
  Three families of read/write, by decomposition:
    - Attribute: scalar, at a group.
    - DistributedView: 1D or 2D [N, NDOFS], MPI-decomposed along N.
      Vector overloads (1D only) delegate here.
    - LocalView: scalar, 1D, or 2D, NOT decomposed -- identical on every
      rank.  Also supports paging the leading dimension of the on-disk
      array, by time_index (on-disk array has one more dimension than
      value) or time_chunk (same rank, value's dim 0 is the page length).

  Every view argument must be contiguous; LocalViews and the 2D
  DistributedView must additionally be LayoutRight.
*/

#ifndef AMANZI_FILE_HDF5_HH_
#define AMANZI_FILE_HDF5_HH_

#include <string>
#include <type_traits>
#include <utility>
#include <vector>

extern "C"
{
#include "hdf5.h"
#include "parallelIO.h"
};

#include "Tpetra_ConfigDefs.hpp"
#include "AmanziTypes.hh"

namespace Amanzi {

class FileHDF5 {
 public:
  // [begin, end) range of indices, for paged LocalView reads/writes.
  using TimeRange = Kokkos::pair<std::size_t, std::size_t>;

  FileHDF5(const Comm_ptr_type& comm, const std::string& filename, file_mode_t mode);
  ~FileHDF5();

  // not copyable or movable -- owns a parallelIO/HDF5 file handle with no
  // defined moved-from state; a shallow copy would double-close on
  // destruction.
  FileHDF5(const FileHDF5&) = delete;
  FileHDF5& operator=(const FileHDF5&) = delete;
  FileHDF5(FileHDF5&&) = delete;
  FileHDF5& operator=(FileHDF5&&) = delete;

  void openFile(file_mode_t mode = FILE_READWRITE);
  void closeFile();

  //
  // Attributes -- scalar, at a group given by path (default root).
  //
  template <typename scalar_type>
  void writeAttribute(const std::string& name, const std::string& path, scalar_type value);

  template <typename scalar_type>
  scalar_type readAttribute(const std::string& name, const std::string& path = "/");

  //
  // LocalView -- scalar, 1D, or 2D, NOT decomposed.
  //
  template <typename scalar_type>
  std::enable_if_t<std::is_arithmetic_v<scalar_type>>
  writeLocalView(const std::string& name, const std::string& path, scalar_type value);

  template <typename scalar_type>
  std::enable_if_t<std::is_arithmetic_v<scalar_type>, scalar_type>
  readLocalView(const std::string& name, const std::string& path);

  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2)>
  writeLocalView(const std::string& name, const std::string& path, view_type view);

  // view may be default-constructed (size 0), in which case it is resized
  // to fit the on-disk shape -- either all ranks must pass an empty view,
  // or all must pass one already sized to match.
  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2)>
  readLocalView(const std::string& name, const std::string& path, view_type& view);

  //
  // LocalView, paged by a single time_index -- on-disk array has one
  // more dimension than value.
  //
  template <typename scalar_type>
  std::enable_if_t<std::is_arithmetic_v<scalar_type>> writeLocalView(const std::string& name,
                                                                     const std::string& path,
                                                                     scalar_type value,
                                                                     std::size_t time_index);

  template <typename scalar_type>
  std::enable_if_t<std::is_arithmetic_v<scalar_type>, scalar_type>
  readLocalView(const std::string& name, const std::string& path, std::size_t time_index);

  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2)>
  writeLocalView(const std::string& name,
                 const std::string& path,
                 view_type view,
                 std::size_t time_index);

  // See the unpaged readLocalView's comment above re: empty views.
  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2)>
  readLocalView(const std::string& name,
                const std::string& path,
                view_type& view,
                std::size_t time_index);

  //
  // LocalView, paged by a time_chunk range -- on-disk array has the
  // same rank as value, value's 0th dimension being the page length.
  //
  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2 || view_type::rank == 3)>
  writeLocalView(const std::string& name,
                 const std::string& path,
                 view_type view,
                 TimeRange time_chunk);

  // See the unpaged readLocalView's comment above re: empty views.
  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2 || view_type::rank == 3)>
  readLocalView(const std::string& name,
                const std::string& path,
                view_type& view,
                TimeRange time_chunk);

  //
  // DistributedView -- 1D, decomposed across ranks.
  //
  template <typename view_type>
  std::enable_if_t<view_type::rank == 1>
  writeDistributedView(const std::string& name, const std::string& path, view_type view);

  template <typename view_type>
  std::enable_if_t<view_type::rank == 1>
  readDistributedView(const std::string& name, const std::string& path, view_type view);

  //
  // Tpetra::Vector, delegates to the 1D DistributedView.
  //
  template <typename scalar_type>
  void writeVector(const std::string& name,
                   const std::string& path,
                   const Vector_type_<scalar_type>& vec);

  template <typename scalar_type>
  void readVector(const std::string& name, const std::string& path, Vector_type_<scalar_type>& vec);

  //
  // DistributedView -- 2D [N, NDOFS], decomposed along dimension 0 (N);
  // NDOFS is the fastest-varying dimension.  View must be LayoutRight
  // (not the transpose Tpetra::MultiVector uses internally).
  //
  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   view_type::rank == 2>
  writeDistributedView(const std::string& name, const std::string& path, view_type view);

  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   view_type::rank == 2>
  readDistributedView(const std::string& name, const std::string& path, view_type view);

 private:
  // Runs fn() on rank 0 only, then broadcasts whether it threw
  // Errors::Message so every rank either throws together or proceeds
  // together -- otherwise a rank-0-only throw strands every other rank at
  // the broadcast that follows a successful read.
  template <typename Fn>
  void runOnRank0AndBroadcastErrors_(Fn&& fn);

  // Shared tail of every LocalView read: if view arrived empty, broadcasts
  // rank 0's now-known shape and resizes every other rank's view to
  // match; then broadcasts the data itself.  Called after rank 0 has
  // already read view's data via runOnRank0AndBroadcastErrors_.
  template <typename view_type>
  void resizeIfEmptyAndBroadcast_(view_type& view, bool was_empty);

  void ensureGroupExists_(const std::string& h5path);

  static std::string joinPath_(const std::string& path, const std::string& name);

  void createGroup_(const std::string& h5path);
  void closeGroup_();
  bool hasGroup_(const std::string& h5path) const;

  hid_t rawFileId_() const;

  Comm_ptr_type comm_;

  iogroup_conf_t IOconfig_;
  iogroup_t IOgroup_;
  hid_t data_file_;
  std::vector<int> open_groups_;
  std::vector<std::string> all_groups_;

  std::string filename_;
};


} // namespace Amanzi


#include "FileHDF5_impl.hh"


#endif
