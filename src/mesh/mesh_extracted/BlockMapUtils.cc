/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Mesh

*/

// TPLs
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"

// Amanzi
#include "BlockMapUtils.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

/* ******************************************************************
*  Convert discontinuous map to continuous map
*
*  Input:
*  mesh - pointer to Amanzi mesh
*  parent_maps - pair of continuous master and ghosted maps
*  subset_maps - pair of discontinuous master and ghosted maps
*                this must be a subset of parent_maps
*
*  Returns:
*  pair of continuous master and ghosted maps that corresponds to
*  subset maps
****************************************************************** */
std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map>>
CreateContinuousMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                     const std::pair<Teuchos::RCP<const Epetra_BlockMap>,
                                     Teuchos::RCP<const Epetra_BlockMap>>& parent_maps,
                     const std::pair<Teuchos::RCP<const Epetra_BlockMap>,
                                     Teuchos::RCP<const Epetra_BlockMap>>& subset_maps)
{
  int n_owned = subset_maps.first->NumMyElements();

  const auto& comm = subset_maps.first->Comm();
  auto continuous_map = Teuchos::rcp(new Epetra_Map(-1, n_owned, 0, comm));

  int n_all = subset_maps.second->NumMyElements();
  int n_ghost = n_all - n_owned;
  int m_ghost = std::max(n_ghost, 1);
  std::vector<int> gl_id(m_ghost), pr_id(m_ghost), lc_id(m_ghost);

  int total_proc = mesh->get_comm()->NumProc();
  int my_pid = mesh->get_comm()->MyPID();
  std::vector<int> min_global_id(total_proc, 0), tmp(total_proc, 0);

  tmp[my_pid] = continuous_map->GID(0);

  mesh->get_comm()->SumAll(tmp.data(), min_global_id.data(), total_proc);

  for (int n = n_owned; n < subset_maps.second->NumMyElements(); n++) {
    int f = parent_maps.second->LID(subset_maps.second->GID(n));
    gl_id[n - n_owned] = parent_maps.second->GID(f);
  }

  subset_maps.first->RemoteIDList(n_ghost, gl_id.data(), pr_id.data(), lc_id.data());

  int n_all_new = n_owned;
  for (int i = 0; i < n_ghost; i++) {
    if (pr_id[i] >= 0) n_all_new++;
  }

  int m_all_new = std::max(n_all_new, 1);
  std::vector<int> global_id_ghosted(m_all_new);
  for (int i = 0; i < n_owned; i++) { global_id_ghosted[i] = continuous_map->GID(i); }

  int k = n_owned;
  for (int i = 0; i < n_ghost; i++) {
    if (pr_id[i] >= 0) {
      int proc_id = pr_id[i];
      global_id_ghosted[k++] = min_global_id[proc_id] + lc_id[i];
    }
  }

  Teuchos::RCP<Epetra_Map> continuous_map_ghosted =
    Teuchos::rcp(new Epetra_Map(-1, n_all_new, global_id_ghosted.data(), 0, comm));

  return std::make_pair(continuous_map, continuous_map_ghosted);
}

} // namespace AmanziMesh
} // namespace Amanzi
