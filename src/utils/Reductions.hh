/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
#pragma once

namespace Amanzi {
namespace Reductions {

// Compute pairs of value + location
struct ValLoc {
  double value;
  AmanziMesh::Entity_GID gid;
};


inline int
reduceAllMaxLoc(const Comm_type& comm, const ValLoc& local, ValLoc& global)
{
  MpiComm_type const* mpi_comm = dynamic_cast<const MpiComm_type*>(&comm);
  const MPI_Comm& mpi_comm_raw = mpi_comm->Comm();
  return MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpi_comm_raw);
}

inline ValLoc
reduceAllMaxLoc(const Epetra_Vector& vec)
{
  ValLoc local{ 0., 0 };
  for (int i = 0; i != vec.MyLength(); ++i) {
    if (vec[i] > local.value) {
      local.value = vec[i];
      local.gid = vec.Map().GID(i);
    }
  }
  ValLoc global{ 0., 0 };
  int ierr = reduceAllMaxLoc(vec.Comm(), local, global);
  AMANZI_ASSERT(!ierr);
  return global;
}


} // namespace Reductions
} // namespace Amanzi
