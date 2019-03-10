/*
  Data Structures

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Takes non-contiguous data structure spaces (CompositeVector, TreeVector) 
  and converts them into a single map.
*/

#include "CompositeVectorSpace.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector_Utils.hh"
#include "SuperMap.hh"
#include "Epetra_Vector.h"

namespace Amanzi {
namespace Operators {

SuperMap::SuperMap(const Comm_ptr_type& comm,
                   const std::vector<std::string>& compnames,
                   const std::vector<int>& dofnums,
                   const std::vector<Teuchos::RCP<const Epetra_BlockMap> >& maps,
                   const std::vector<Teuchos::RCP<const Epetra_BlockMap> >& ghosted_maps) :
    compnames_(compnames)
{
  AMANZI_ASSERT(compnames.size() == dofnums.size());
  AMANZI_ASSERT(compnames.size() == maps.size());
  AMANZI_ASSERT(compnames.size() == ghosted_maps.size());

  int MyPID = comm.MyPID();
  int NumProc = comm.NumProc();

  int ncomps = compnames.size();
    
  // fill the offsets
  int offset = 0;
  for (int i=0; i!=ncomps; ++i) {
    std::string compname = compnames[i];
    num_dofs_[compname] = dofnums[i];

    counts_[compname] = maps[i]->NumMyPoints();
    ghosted_counts_[compname] = ghosted_maps[i]->NumMyPoints() - counts_[compname];

    offsets_[compname] = offset;
    offset += dofnums[i] * counts_[compname];

    indices_map_[compname] = maps[i];
    ghosted_indices_map_[compname] = ghosted_maps[i];
  }
  int n_local = offset;
                
  // fill the ghosted offsets
  for (int i=0; i!=ncomps; ++i) {
    std::string compname = compnames[i];
    ghosted_offsets_[compname] = offset;
    offset += dofnums[i] * ghosted_counts_[compname];
  }
  int n_local_ghosted = offset;

  std::vector<int> tmp(NumProc * ncomps, 0);
  std::vector<int> tmp_ghosted(NumProc * ncomps, 0);  
  std::vector<int> counts_all((NumProc + 1) * ncomps, 0);
  std::vector<int> counts_ghosted_all((NumProc + 1) * ncomps, 0);  

  for (int i=0; i!=ncomps; ++i) {
    tmp[MyPID * ncomps + i] = counts_[compnames[i]];
    tmp_ghosted[MyPID * ncomps + i] = counts_[compnames[i]] + ghosted_counts_[compnames[i]];
  }
      
  comm.SumAll(&tmp[0], &counts_all[0], NumProc * ncomps);
  comm.SumAll(&tmp_ghosted[0], &counts_ghosted_all[0], NumProc * ncomps);

  for (int i=0; i!=ncomps; ++i) {
    int tmp_val = 0;
    int sum_val = 0;
    int sum_ght_val = 0;    

    for (int j=0; j<=NumProc; j++) {
      tmp_val = counts_all[j * ncomps + i];
      counts_all[j * ncomps + i] = sum_val;
      sum_val += tmp_val;

      tmp_val = counts_ghosted_all[j * ncomps + i];
      counts_ghosted_all[j * ncomps + i] = sum_ght_val;
      sum_ght_val += tmp_val;
    }
  }
  
  std::vector<Teuchos::RCP<Epetra_Vector> > global_ids_owned(ncomps);
  std::vector<Teuchos::RCP<Epetra_Vector> > global_ids_ght(ncomps);
  std::vector<Teuchos::RCP<Epetra_Import> > importer(ncomps);

  int n(0);
  for (int i=0; i!=ncomps; ++i) {
    global_ids_owned[i] = Teuchos::rcp(new Epetra_Vector(*(maps[i])));
    global_ids_ght[i] = Teuchos::rcp(new Epetra_Vector(*(ghosted_maps[i])));
    importer[i] = Teuchos::rcp(new Epetra_Import(*(ghosted_maps[i]), *(maps[i])));

    for (int j=0; j!=counts_[compnames[i]]; ++j) {
      int id = counts_all[MyPID * ncomps + i] + j;
      (*global_ids_owned[i])[j] = id;
    }

    int ierr = global_ids_ght[i]->Import(*(global_ids_owned[i]), *(importer[i]), Insert);
  }
  
  // populate the map GIDs
  n = 0;
  int global_offset = 0;

  std::vector<int> gids(offset, -1);
  for (int i=0; i!=ncomps; ++i) {
    for (int j=0; j!=counts_[compnames[i]]; ++j) {
      int base = global_offset + (*global_ids_owned[i])[j] * dofnums[i];
      for (int dof=0; dof!=dofnums[i]; ++dof) {        
        gids[n++] = base + dof;
      }
    }
    global_offset += maps[i]->NumGlobalPoints() * dofnums[i];
  }
  int n_global = global_offset;

  global_offset = 0;
  int n_global_ghosted = 0;
  for (int i=0; i!=ncomps; ++i) {
    int j0 = counts_[compnames[i]];
    for (int j=j0; j!=(j0 + ghosted_counts_[compnames[i]]); ++j) {
      int base = global_offset + (*global_ids_ght[i])[j] * dofnums[i];
      for (int dof=0; dof!=dofnums[i]; ++dof) {
        gids[n++] = base+dof;
      }
    }
    global_offset += maps[i]->NumGlobalPoints() * dofnums[i];
    n_global_ghosted += ghosted_maps[i]->NumGlobalPoints() * dofnums[i];
  }

  // create the maps
  map_ = Teuchos::rcp(new Epetra_Map(n_global, n_local, &gids[0], 0, *comm)); 
  ghosted_map_ = Teuchos::rcp(new Epetra_Map(n_global_ghosted, n_local_ghosted, &gids[0], 0, *comm));
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


std::pair<int, Teuchos::RCP<std::vector<int> > >
SuperMap::BlockIndices() const {
  auto block_indices = Teuchos::rcp(new std::vector<int>(Map()->NumMyElements()));
  int block_id = 0;
  for (const auto& comp : *this) {
    int ndofs = NumDofs(comp);
    for (int d=0; d!=ndofs; ++d) {
      const auto& inds = Indices(comp, d);
      for (int i=0; i!=inds.size(); ++i) (*block_indices)[inds[i]] = block_id;
      block_id++;
    }
  }
  return std::make_pair(block_id, block_indices);
}


const std::vector<int>&
SuperMap::CreateIndices_(const std::string& compname, int dofnum, bool ghosted) const {
  if (ghosted) {
    if (ghosted_indices_.count(compname) == 0) {
      ghosted_indices_[compname];
    }

    // create the vector
    int num_dof = num_dofs_.at(compname);
    int nentities_owned = counts_.at(compname);
    int nentities = nentities_owned + ghosted_counts_.at(compname);

    std::vector<int> indices(nentities, -1);
    int offset = offsets_.at(compname);

    for (int i=0; i!=nentities_owned; ++i) {
      indices[i] = offset + dofnum + i*num_dof;
    }
      
    int MyPID = BaseGhostedMap(compname)->Comm().MyPID();
    int ghosted_offset = ghosted_offsets_.at(compname);

    for (int i=nentities_owned; i!=nentities; ++i) {
      indices[i] = ghosted_offset + dofnum + (i-nentities_owned)*num_dof;
      //if (MyPID==1) std::cout << "ghost "<<i<<" "<<indices[i] <<"\n";
    }

    // assign
    ghosted_indices_[compname][dofnum] = indices;
    return ghosted_indices_[compname][dofnum];

  } else {
    if (indices_.count(compname) == 0) {
      indices_[compname];
    }
    
    int num_dof = num_dofs_.at(compname);
    // create the vector
    int nentities = counts_.at(compname);
    
    std::vector<int> indices(nentities, -1);
    int offset = offsets_.at(compname);

    for (int i=0; i!=nentities; ++i) {
      indices[i] = offset + dofnum + i*num_dof;
    }

    // assign
    indices_[compname][dofnum] = indices;

    return indices_[compname][dofnum];
  }
}


// Nonmember contructors/factories
Teuchos::RCP<SuperMap> createSuperMap(const CompositeVectorSpace& cv) {
  std::vector<std::string> names;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > maps;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > ghost_maps;

  for (auto it=cv.begin(); it!=cv.end(); ++it) {
    names.push_back(*it);
    dofnums.push_back(cv.NumVectors(*it));

    maps.push_back(cv.Map(*it, false));
    ghost_maps.push_back(cv.Map(*it, true));
  }

  return Teuchos::rcp(new SuperMap(cv.Comm(), names, dofnums, maps, ghost_maps));
}


Teuchos::RCP<SuperMap> createSuperMap(const TreeVectorSpace& tv) {
  std::vector<std::string> names;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > maps;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > ghost_maps;

  if (tv.Data() != Teuchos::null) {
    // TV with only a CV inside
    return createSuperMap(*tv.Data());
  } else {
    // multiple children
    // grab the leaf nodes
    std::vector<Teuchos::RCP<const TreeVectorSpace> > tvs =
        collectTreeVectorLeaves_const<TreeVectorSpace>(tv);

    // loop over nodes, finding unique component names on unique meshes
    for (auto it=tvs.begin(); it!=tvs.end(); ++it) {
      for (CompositeVectorSpace::name_iterator compname=(*it)->Data()->begin();
           compname!=(*it)->Data()->end(); ++compname) {
        int index = std::find(names.begin(), names.end(), *compname) - names.begin();
        if (index == names.size()) {
          names.push_back(*compname);
          dofnums.push_back(1);

          maps.push_back((*it)->Data()->Map(*compname, false));
          ghost_maps.push_back((*it)->Data()->Map(*compname, true));
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
