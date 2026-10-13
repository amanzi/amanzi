/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Shared helpers for parsing "tensor type" and writing per-cell DoFs into a
//! TensorVector.
/*
  State

  Small, free-function helpers shared by tensor-valued independent
  evaluators (e.g. EvaluatorIndependentTensorFromFile).

  - parseTensorType() maps a "tensor type" string and spatial dimension to a
    tensor rank and the number of degrees of freedom (columns) needed.
  - copyVectorToTensorVector() takes a host view of a CompositeVector
    component, laid out as [local id, dof], and writes the DoFs into the
    corresponding entries of a TensorVector, following the same DoF-to-entry
    convention as the "tensor type" options above.
*/

#pragma once

#include <string>

#include "errors.hh"
#include "exceptions.hh"

#include "TensorVector.hh"

namespace Amanzi {
namespace EvaluatorTensorFromFile_Helpers {

//
// Maps a "tensor type" string and spatial dimension to a tensor rank and the
// number of degrees of freedom (columns) needed to represent it:
//
//   "scalar"                  -> rank 1, isotropic                (1 dof)
//   "horizontal and vertical" -> rank 2, diagonal, x=y != z        (2 dof)
//   "diagonal"                -> rank 2, diagonal                  (dimension dof)
//   "full symmetric"          -> rank 2, symmetric                 (3 dof in 2D, 6 in 3D)
//   "full"                    -> rank 2, general                   (dimension^2 dof)
//
// error_prefix is used to prefix the exception message on an invalid
// tensor_type (typically the name of the calling evaluator class).
inline void
parseTensorType(const std::string& tensor_type,
                int dimension,
                const std::string& error_prefix,
                int& rank,
                int& num_funcs)
{
  AMANZI_ASSERT(dimension > 0);

  if (tensor_type == "scalar") {
    rank = 1;
    num_funcs = 1;
  } else if (tensor_type == "horizontal and vertical") {
    rank = 2;
    num_funcs = 2;
  } else if (tensor_type == "diagonal") {
    rank = 2;
    num_funcs = dimension;
  } else if (tensor_type == "full symmetric") {
    rank = 2;
    num_funcs = dimension == 2 ? 3 : 6;
  } else if (tensor_type == "full") {
    rank = 2;
    num_funcs = dimension * dimension;
  } else {
    Errors::Message msg;
    msg << error_prefix << ": invalid parameter \"tensor type\" = \"" << tensor_type
        << "\", must be one of: \"scalar\", \"horizontal and vertical\", \"diagonal\", "
           "\"full symmetric\", or \"full\"";
    Exceptions::amanzi_throw(msg);
  }
}


//
// Writes a host view of a single CompositeVector component (indexed [local
// id, dof]) into the (row,col) entries of the cells [j, j+v.extent(0)) of a
// TensorVector, following the DoF ordering implied by parseTensorType()
// above:
//
//   ndofs == 1                    -> isotropic:            (0,0)
//   ndofs == 2 && space_dim == 3  -> horiz/vert:            (0,0)=(1,1), (2,2)
//   ndofs == space_dim            -> diagonal:              (d,d) for d in [0,space_dim)
//   ndofs == 3 (2D) or 6 (3D)     -> full symmetric:        diagonal + symmetric off-diagonal
//   ndofs == dim*dim              -> full (nonsymmetric):   all (row,col), row-major
//
// Entries are written via TensorVector::at_host(), so the caller must call
// tv.data.update_entries_device() once all writes (across all components)
// are complete.
template <class HostView_type>
void
copyVectorToTensorVector(const HostView_type& v, int j, int space_dim, TensorVector& tv)
{
  const int ni = v.extent(0);
  const int ndofs = v.extent(1);

  if (ndofs == 1) {
    // isotropic
    for (int i = 0; i != ni; ++i) { tv.at_host(i + j)(0, 0) = v(i, 0); }

  } else if (ndofs == 2 && space_dim == 3) {
    // horizontal and vertical
    for (int i = 0; i != ni; ++i) {
      auto Ti = tv.at_host(i + j);
      Ti(0, 0) = v(i, 0);
      Ti(1, 1) = v(i, 0);
      Ti(2, 2) = v(i, 1);
    }

  } else if (ndofs >= space_dim) {
    // diagonal entries
    for (int d = 0; d != space_dim; ++d) {
      for (int i = 0; i != ni; ++i) { tv.at_host(i + j)(d, d) = v(i, d); }
    }

    if (ndofs > space_dim) {
      // off-diagonal entries -- full symmetric or full (nonsymmetric)
      if (ndofs == 3) {
        // 2D full symmetric: xy & yx
        for (int i = 0; i != ni; ++i) {
          auto Ti = tv.at_host(i + j);
          Ti(0, 1) = Ti(1, 0) = v(i, 2);
        }

      } else if (ndofs == 6) {
        // 3D full symmetric: xy & yx, xz & zx, yz & zy
        for (int i = 0; i != ni; ++i) {
          auto Ti = tv.at_host(i + j);
          Ti(0, 1) = Ti(1, 0) = v(i, 3);
          Ti(0, 2) = Ti(2, 0) = v(i, 4);
          Ti(1, 2) = Ti(2, 1) = v(i, 5);
        }

      } else if (ndofs == 9) {
        // 3D full (nonsymmetric) tensor, row-major
        for (int i = 0; i != ni; ++i) {
          auto Ti = tv.at_host(i + j);
          Ti(0, 0) = v(i, 0);
          Ti(0, 1) = v(i, 1);
          Ti(0, 2) = v(i, 2);
          Ti(1, 0) = v(i, 3);
          Ti(1, 1) = v(i, 4);
          Ti(1, 2) = v(i, 5);
          Ti(2, 0) = v(i, 6);
          Ti(2, 1) = v(i, 7);
          Ti(2, 2) = v(i, 8);
        }
      } else if (ndofs == 4) {
        // 2D full (nonsymmetric) tensor, row-major
        for (int i = 0; i != ni; ++i) {
          auto Ti = tv.at_host(i + j);
          Ti(0, 0) = v(i, 0);
          Ti(0, 1) = v(i, 1);
          Ti(1, 0) = v(i, 2);
          Ti(1, 1) = v(i, 3);
        }
      } else {
        AMANZI_ASSERT(0);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace EvaluatorTensorFromFile_Helpers
} // namespace Amanzi
