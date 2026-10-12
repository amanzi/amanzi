/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Reads and writes NetCDF files for extruded-column (MPAS-style) meshes.
/*
  FileNetCDF has no opinion on what a view's axes represent -- the caller
  picks dim_names (see InputOutputVerticallyStructuredNetCDF.hh for the
  columnar convention this is designed to support) and decomposed_axis
  (which axis of the view is MPI-decomposed); FileNetCDF just reads or
  writes that shape.

  All views must be contiguous and Kokkos::LayoutRight, except for one
  rank-2, decomposed_axis=0, LayoutLeft overload -- exactly a
  Tpetra::MultiVector's local view.

  Attributes are always global (NC_GLOBAL); there is no group hierarchy.
*/

#pragma once

#include <map>
#include <string>
#include <vector>

extern "C"
{
#include "netcdf.h"
#include "netcdf_par.h"
#include "parallelIO.h"
};

#include "Tpetra_ConfigDefs.hpp"
#include "AmanziTypes.hh"

namespace Amanzi {
class FileNetCDF {
 public:
  // [begin, end) range of time indices, for chunked time-dependent reads/writes.
  using TimeRange = Kokkos::pair<std::size_t, std::size_t>;

  FileNetCDF(const Comm_ptr_type& comm, const std::string& filename, file_mode_t mode);
  ~FileNetCDF();

  // not copyable or movable -- owns a NetCDF file handle with no defined
  // moved-from state; a shallow copy would double-close on destruction.
  FileNetCDF(const FileNetCDF&) = delete;
  FileNetCDF& operator=(const FileNetCDF&) = delete;
  FileNetCDF(FileNetCDF&&) = delete;
  FileNetCDF& operator=(FileNetCDF&&) = delete;

  void openFile(file_mode_t mode = FILE_READWRITE);
  void closeFile();

  //
  // Attributes -- scalar.  Always global (NC_GLOBAL).
  //
  template <typename scalar_type>
  void writeAttribute(const std::string& attr_name, scalar_type value);

  template <typename scalar_type>
  scalar_type readAttribute(const std::string& attr_name);

  //
  // Attributes -- time-dependent (scalar_type = double or int), chunk or
  // single index.  Always a plain global [NTIMES] variable.
  //
  template <typename view_type>
  std::enable_if_t<view_type::rank == 1>
  writeAttribute(const std::string& attr_name, view_type vals, TimeRange time_chunk);

  template <typename view_type>
  std::enable_if_t<view_type::rank == 1>
  readAttribute(const std::string& attr_name, view_type vals, TimeRange time_chunk);

  template <typename scalar_type>
  void writeAttribute(const std::string& attr_name, scalar_type value, std::size_t time_index);

  template <typename scalar_type>
  scalar_type readAttribute(const std::string& attr_name, std::size_t time_index);

  //
  // writeDistributedView/readDistributedView: see the class-level comment
  // for the full shape/decomposed_axis table.  view must be, except for
  // the one LayoutLeft overload immediately below, Kokkos::LayoutRight.
  //
  // dim_names has exactly view_type::rank entries, one per view axis in
  // the view's own order (never including the separate "time" axis, which
  // is always named "time" and handled by the TimeRange/time_index
  // overloads below) -- every on-disk dimension is named explicitly by
  // the caller, so two variables only ever share a dimension when the
  // caller says they mean the same axis.  FileNetCDF itself has no
  // opinion on what any axis represents.
  template <typename view_type>
  std::enable_if_t<!(view_type::rank == 2 &&
                     std::is_same_v<typename view_type::array_layout, Kokkos::LayoutLeft>)>
  writeDistributedView(const std::string& varname,
                       view_type view,
                       std::size_t decomposed_axis,
                       const std::vector<std::string>& dim_names);
  template <typename view_type>
  std::enable_if_t<!(view_type::rank == 2 &&
                     std::is_same_v<typename view_type::array_layout, Kokkos::LayoutLeft>)>
  readDistributedView(const std::string& varname,
                      view_type view,
                      std::size_t decomposed_axis,
                      const std::vector<std::string>& dim_names);

  // Rank-2, LayoutLeft, decomposed_axis fixed at 0 -- exactly a
  // Tpetra::MultiVector's local view.  The on-disk variable's dimensions
  // are declared in the buffer's memory order ([NDOFS, N], N innermost),
  // so a file written this way is only byte-compatible with a LayoutLeft
  // reader of the same shape.  dim_names has exactly 2 entries, in that
  // same on-disk ([NDOFS, N]) order.
  template <typename view_type>
  std::enable_if_t<view_type::rank == 2 &&
                   std::is_same_v<typename view_type::array_layout, Kokkos::LayoutLeft>>
  writeDistributedView(const std::string& varname,
                       view_type view,
                       const std::vector<std::string>& dim_names);
  template <typename view_type>
  std::enable_if_t<view_type::rank == 2 &&
                   std::is_same_v<typename view_type::array_layout, Kokkos::LayoutLeft>>
  readDistributedView(const std::string& varname,
                      view_type view,
                      const std::vector<std::string>& dim_names);

  template <typename view_type>
  void writeDistributedView(const std::string& varname,
                            view_type view,
                            std::size_t decomposed_axis,
                            const std::vector<std::string>& dim_names,
                            TimeRange time_chunk);
  template <typename view_type>
  void readDistributedView(const std::string& varname,
                           view_type view,
                           std::size_t decomposed_axis,
                           const std::vector<std::string>& dim_names,
                           TimeRange time_chunk);

  template <typename view_type>
  void writeDistributedView(const std::string& varname,
                            view_type view,
                            std::size_t decomposed_axis,
                            const std::vector<std::string>& dim_names,
                            std::size_t time_index);
  template <typename view_type>
  void readDistributedView(const std::string& varname,
                           view_type view,
                           std::size_t decomposed_axis,
                           const std::vector<std::string>& dim_names,
                           std::size_t time_index);

  //
  // LocalView -- scalar, 1D, or 2D, identical on every rank, NOT
  // decomposed.  Views must be Kokkos::LayoutRight.
  //
  template <typename scalar_type>
  std::enable_if_t<std::is_arithmetic_v<scalar_type>>
  writeLocalView(const std::string& varname, scalar_type value);

  template <typename scalar_type>
  std::enable_if_t<std::is_arithmetic_v<scalar_type>, scalar_type>
  readLocalView(const std::string& varname);

  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2)>
  writeLocalView(const std::string& varname, view_type view);

  // view may be default-constructed (size 0), in which case it is resized
  // to fit the on-disk shape -- either all ranks must pass an empty view,
  // or all must pass one already sized to match.
  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2)>
  readLocalView(const std::string& varname, view_type& view);

  //
  // LocalView, paged by a single time_index -- on-disk array has one more
  // (leading) dimension than value.
  //
  template <typename scalar_type>
  std::enable_if_t<std::is_arithmetic_v<scalar_type>>
  writeLocalView(const std::string& varname, scalar_type value, std::size_t time_index);

  template <typename scalar_type>
  std::enable_if_t<std::is_arithmetic_v<scalar_type>, scalar_type>
  readLocalView(const std::string& varname, std::size_t time_index);

  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2)>
  writeLocalView(const std::string& varname, view_type view, std::size_t time_index);

  // See the unpaged readLocalView's comment above re: empty views.
  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2)>
  readLocalView(const std::string& varname, view_type& view, std::size_t time_index);

  //
  // LocalView, paged by a time_chunk range -- on-disk array has the same
  // rank as value, value's 0th dimension being the page length.
  //
  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2 || view_type::rank == 3)>
  writeLocalView(const std::string& varname, view_type view, TimeRange time_chunk);

  // See the unpaged readLocalView's comment above re: empty views.
  template <typename view_type>
  std::enable_if_t<std::is_same_v<typename view_type::array_layout, Kokkos::LayoutRight> &&
                   (view_type::rank == 1 || view_type::rank == 2 || view_type::rank == 3)>
  readLocalView(const std::string& varname, view_type& view, TimeRange time_chunk);

  //
  // Vector -- delegates to the 1D writeDistributedView/readDistributedView
  // (decomposed_axis=0).  dim_names has exactly 1 entry, naming that axis.
  //
  template <typename scalar_type>
  void writeVector(const std::string& varname,
                   const Vector_type_<scalar_type>& vec,
                   const std::vector<std::string>& dim_names);

  template <typename scalar_type>
  void readVector(const std::string& varname,
                  Vector_type_<scalar_type>& vec,
                  const std::vector<std::string>& dim_names);

  // Describes one on-disk dimension for getOrCreateVariable_.  length is
  // the dimension's full global length.
  struct DimSpec {
    std::string name;
    std::size_t length;
    bool unlimited = false;
  };

 private:
  // Ensures a variable of the given shape exists, creating it (with
  // on-disk type nc_type, and its dimensions, if not already defined) if
  // not.  dims is the on-disk dimension list in order, outermost first.
  // dimids_out is populated with the resulting dimension ids, same order
  // as dims.  nc_type is ignored if the variable already exists.
  int getOrCreateVariable_(const std::string& varname,
                           const std::vector<DimSpec>& dims,
                           std::vector<int>& dimids_out,
                           nc_type type);

  // Same as getOrCreateVariable_, but for a LocalView variable, with
  // par_access (NC_INDEPENDENT or NC_COLLECTIVE) and the on-disk type
  // picked from scalar_type.
  template <typename scalar_type>
  int getOrCreateLocalVariable_(const std::string& varname,
                                const std::vector<DimSpec>& dims,
                                std::vector<int>& dimids_out,
                                int par_access);

  int getOrCreateDim_(const std::string& dimname, size_t len, bool unlimited = false);

  void checkErr_(int ierr, const std::string& what) const;

  Comm_ptr_type comm_;
  MPI_Comm mpi_comm_;
  std::string filename_;
  file_mode_t mode_;

  int ncid_;
  bool file_open_;

  // cached dimension ids, by name
  std::map<std::string, int> dimids_;
  // cached variable ids, by name
  std::map<std::string, int> varids_;
};

} // namespace Amanzi

#include "FileNetCDF_impl.hh"
