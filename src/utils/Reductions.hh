/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
#pragma once

#include <limits>
#include "Epetra_Vector.h"
#include "dbc.hh"
#include "AmanziComm.hh"

namespace Amanzi {
namespace Reductions {

// Compute pairs of max value + location
struct MaxLoc {
  double value;
  int gid;
};

inline MaxLoc createEmptyMaxLoc()
{
  return { std::numeric_limits<double>::lowest(), -1 };
}

inline std::ostream&
operator<<(std::ostream& os, const MaxLoc& ml) {
  os << "[" << ml.gid << "] " << ml.value;
  return os;
}

// Compute {min value, location}, {max value, location}
using MinMaxLoc = std::array<MaxLoc,2>;

inline MinMaxLoc createEmptyMinMaxLoc()
{
  return { MaxLoc{ std::numeric_limits<double>::max(), -1 },
           MaxLoc{ std::numeric_limits<double>::lowest(), -1 } };
}

inline std::ostream&
operator<<(std::ostream& os, const MinMaxLoc& mml) {
  os << mml[0] << ", " << mml[1];
  return os;
}


inline MaxLoc
reduceAllMaxLoc(const Comm_type& comm, const MaxLoc& local)
{
  MpiComm_type const* mpi_comm = dynamic_cast<const MpiComm_type*>(&comm);
  const MPI_Comm& mpi_comm_raw = mpi_comm->Comm();
  MaxLoc global;
  int ierr = MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpi_comm_raw);
  AMANZI_ASSERT(!ierr);
  return global;
}


inline MinMaxLoc
reduceAllMinMaxLoc(const Comm_type& comm, const MinMaxLoc& local)
{
  MpiComm_type const* mpi_comm = dynamic_cast<const MpiComm_type*>(&comm);
  const MPI_Comm& mpi_comm_raw = mpi_comm->Comm();

  // negate the min value to allow a single Allreduce
  MinMaxLoc global(local);
  global[0].value = -global[0].value;
  int ierr = MPI_Allreduce(MPI_IN_PLACE, &global, 2, MPI_DOUBLE_INT, MPI_MAXLOC, mpi_comm_raw);
  AMANZI_ASSERT(!ierr);
  global[0].value = -global[0].value;
  return global;
}


inline MaxLoc
reduceAllMaxLoc(const Epetra_Vector& vec)
{
  MaxLoc local{ std::numeric_limits<double>::lowest(), -1 };
  for (int i = 0; i != vec.MyLength(); ++i) {
    if (vec[i] > local.value) {
      local.value = vec[i];
      local.gid = vec.Map().GID(i);
    }
  }
  return reduceAllMaxLoc(vec.Comm(), local);
}


inline MinMaxLoc
reduceAllMinMaxLoc(const Epetra_Vector& vec)
{
  MinMaxLoc local = createEmptyMinMaxLoc();

  for (int i = 0; i != vec.MyLength(); ++i) {
    if (vec[i] < local[0].value) {
      local[0].value = vec[i];
      local[0].gid = vec.Map().GID(i);
    }
    if (vec[i] > local[1].value) {
      local[1].value = vec[i];
      local[1].gid = vec.Map().GID(i);
    }
  }

  return reduceAllMinMaxLoc(vec.Comm(), local);
}


} // namespace Reductions
} // namespace Amanzi
