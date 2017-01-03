/*
  Author: Ethan Coon (ecoon@lanl.gov)

  Takes non-contiguous data structure spaces (CompositeVector, TreeVector) and
  converts them into a single map.
*/


#include "CompositeVectorSpace.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector_Utils.hh"
#include "SuperMap.hh"

namespace Amanzi {
namespace Operators {

SuperMap::SuperMap(const Epetra_MpiComm& comm,
                   const std::vector<std::string>& compnames,
                   const std::vector<int>& dofnums,
                   const std::vector<Teuchos::RCP<const Epetra_Map> >& maps,
                   const std::vector<Teuchos::RCP<const Epetra_Map> >& ghosted_maps) :
    compnames_(compnames)
{
  ASSERT(compnames.size() == dofnums.size());
  ASSERT(compnames.size() == maps.size());
  ASSERT(compnames.size() == ghosted_maps.size());

  int offset = 0;

  // fill the offsets
  for (int i=0; i!=compnames.size(); ++i) {
    std::string compname = compnames[i];
    num_dofs_[compname] = dofnums[i];
    counts_[compname] = maps[i]->NumMyElements();
    ghosted_counts_[compname] = ghosted_maps[i]->NumMyElements() - counts_[compname];
    offsets_[compname] = offset;

    offset += dofnums[i]*counts_[compnames[i]];
  }
  int n_local = offset;

  // fill the ghosted offsets
  for (int i=0; i!=compnames.size(); ++i) {
    std::string compname = compnames[i];
    ghosted_offsets_[compname] = offset;
    offset += dofnums[i]*ghosted_counts_[compname];
  }
  int n_local_ghosted = offset;

  // populate the map GIDs
  int global_offset = 0;
  int n(0);
  std::vector<int> gids(offset, -1);
  for (int i=0; i!=compnames.size(); ++i) {
    for (int j=0; j!=counts_[compnames[i]]; ++j) {
      int base = global_offset + ghosted_maps[i]->GID(j)*dofnums[i];
      for (int dof=0; dof!=dofnums[i]; ++dof) gids[n++] = base+dof;
    }
    global_offset += maps[i]->NumGlobalElements() * dofnums[i];
  }
  int n_global = global_offset;

  int n_global_ghosted = 0;
  global_offset = 0;
  for (int i=0; i!=compnames.size(); ++i) {
    for (int j=counts_[compnames[i]]; j!=(counts_[compnames[i]]+ghosted_counts_[compnames[i]]); ++j) {
      int base = global_offset + ghosted_maps[i]->GID(j)*dofnums[i];
      for (int dof=0; dof!=dofnums[i]; ++dof) gids[n++] = base+dof;
    }
    global_offset += maps[i]->NumGlobalElements() * dofnums[i];
    n_global_ghosted += ghosted_maps[i]->NumGlobalElements() * dofnums[i];
  }
  
  // create the maps
  map_ = Teuchos::rcp(new Epetra_Map(n_global, n_local, &gids[0], 0, comm));
  ghosted_map_ = Teuchos::rcp(new Epetra_Map(n_global_ghosted, n_local_ghosted, &gids[0], 0, comm));
}


// Copy constructor specifically does not copy indices and ghosted_indices.
// These will be lazily recreated if needed.
SuperMap::SuperMap(const SuperMap& other) :
    compnames_(other.compnames_),
    offsets_(other.offsets_),
    num_dofs_(other.num_dofs_),
    counts_(other.counts_),
    ghosted_offsets_(other.ghosted_offsets_),
    ghosted_counts_(other.ghosted_counts_),
    map_(other.map_),
    ghosted_map_(other.ghosted_map_) {};


bool SuperMap::HasComponent(const std::string& key) const {
  std::map<std::string,int>::const_iterator lb = offsets_.lower_bound(key);
  if (lb != offsets_.end() && !(offsets_.key_comp()(key, lb->first))) {
    return true;
  } else {
    return false;
  }
}
  

const std::vector<int>&
SuperMap::Indices(const std::string& compname, int dofnum) const {
  if (indices_.count(compname)) {
    if (indices_[compname].count(dofnum)) {
      return indices_[compname][dofnum];
    }
  }
  return CreateIndices_(compname, dofnum, false);
}


const std::vector<int>&
SuperMap::GhostIndices(const std::string& compname, int dofnum) const {
  if (ghosted_indices_.count(compname)) {
    if (ghosted_indices_[compname].count(dofnum)) {
      return ghosted_indices_[compname][dofnum];
    }
  }
  return CreateIndices_(compname, dofnum, true);
}


const std::vector<int>&
SuperMap::CreateIndices_(const std::string& compname, int dofnum, bool ghosted) const {
  if (ghosted) {
    if (ghosted_indices_.count(compname) == 0) {
      ghosted_indices_[compname];
    }

    // create the vector
    int nentities_owned = counts_.at(compname);
    int nentities = nentities_owned + ghosted_counts_.at(compname);

    std::vector<int> indices(nentities, -1);
    int offset = offsets_.at(compname);
    int num_dof = num_dofs_.at(compname);
    for (int i=0; i!=nentities_owned; ++i) {
      indices[i] = offset + dofnum + i*num_dof;
    }

    int ghosted_offset = ghosted_offsets_.at(compname);
    for (int i=nentities_owned; i!=nentities; ++i) {
      indices[i] = ghosted_offset + dofnum + (i-nentities_owned)*num_dof;
    }

    // assign
    ghosted_indices_[compname][dofnum] = indices;
    return ghosted_indices_[compname][dofnum];

  } else {
    if (indices_.count(compname) == 0) {
      indices_[compname];
    }

    // create the vector
    int nentities = counts_.at(compname);

    std::vector<int> indices(nentities, -1);
    int offset = offsets_.at(compname);
    int num_dof = num_dofs_.at(compname);
    for (int i=0; i!=nentities; ++i) {
      indices[i] = offset + dofnum + i*num_dof;
    }

    // assign
    indices_[compname][dofnum] = indices;

    return indices_[compname][dofnum];
  }
}


std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> >
getMaps(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_kind location) {
  switch(location) {
    case AmanziMesh::CELL:
      return std::make_pair(Teuchos::rcpFromRef(mesh.cell_map(false)),
                            Teuchos::rcpFromRef(mesh.cell_map(true)));

    case AmanziMesh::FACE:
      return std::make_pair(Teuchos::rcpFromRef(mesh.face_map(false)),
                            Teuchos::rcpFromRef(mesh.face_map(true)));

    case AmanziMesh::EDGE:
      return std::make_pair(Teuchos::rcpFromRef(mesh.edge_map(false)),
                            Teuchos::rcpFromRef(mesh.edge_map(true)));

    case AmanziMesh::NODE:
      return std::make_pair(Teuchos::rcpFromRef(mesh.node_map(false)),
                            Teuchos::rcpFromRef(mesh.node_map(true)));

    case AmanziMesh::BOUNDARY_FACE:
      return std::make_pair(Teuchos::rcpFromRef(mesh.exterior_face_map(false)),
                            Teuchos::rcpFromRef(mesh.exterior_face_map(false)));
    default:
      ASSERT(false);
      return std::make_pair(Teuchos::null, Teuchos::null);
  }
}


// Nonmember contructors/factories
Teuchos::RCP<SuperMap> createSuperMap(const CompositeVectorSpace& cv) {
  std::vector<std::string> names;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  std::vector<Teuchos::RCP<const Epetra_Map> > ghost_maps;

  for (CompositeVectorSpace::name_iterator it=cv.begin();
       it!=cv.end(); ++it) {
    names.push_back(*it);
    dofnums.push_back(cv.NumVectors(*it));

    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cv.Mesh(), cv.Location(*it));
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }

  return Teuchos::rcp(new SuperMap(cv.Comm(), names, dofnums, maps, ghost_maps));
}


Teuchos::RCP<SuperMap> createSuperMap(const TreeVectorSpace& tv) {
  std::vector<std::string> names;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  std::vector<Teuchos::RCP<const Epetra_Map> > ghost_maps;

  if (tv.Data() != Teuchos::null) {
    // TV with only a CV inside
    return createSuperMap(*tv.Data());
  } else {
    // multiple children
    // grab the leaf nodes
    std::vector<Teuchos::RCP<const TreeVectorSpace> > tvs =
        collectTreeVectorLeaves_const<TreeVectorSpace>(tv);

    // loop over nodes, finding unique component names on unique meshes
    for (std::vector<Teuchos::RCP<const TreeVectorSpace> >::const_iterator it=tvs.begin();
         it!=tvs.end(); ++it) {
      for (CompositeVectorSpace::name_iterator compname=(*it)->Data()->begin();
           compname!=(*it)->Data()->end(); ++compname) {
        int index = std::find(names.begin(), names.end(), *compname) - names.begin();
        if (index == names.size()) {
          names.push_back(*compname);
          dofnums.push_back(1);

          std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
              getMaps(*(*it)->Data()->Mesh(), (*it)->Data()->Location(*compname));
          maps.push_back(meshmaps.first);
          ghost_maps.push_back(meshmaps.second);
        } else {
          dofnums[index]++;
        }
      }
    }
    return Teuchos::rcp(new SuperMap(tvs[0]->Data()->Comm(), names, dofnums, maps, ghost_maps));
  }
}

} // namespace Operators
} // namespace Amanzi
