/*
  Mesh

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Lipnikov Konstantin (lipnikov@lanl.gov)
           Rao Garimella (rao@lanl.gov)

  Implementation of algorithms independent of the actual mesh 
  framework.
*/

#include <iterator>
#include <map>
#include <string>
#include <vector>

// Amanzi
#include "errors.hh"
#include "Point.hh"
#include "Region.hh"
#include "RegionLabeledSet.hh"

// Amanzi::AmanziMesh
#include "MeshDefs.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

//---------------------------------------------------------
// Is there a set with this id and entity type
//---------------------------------------------------------
bool Mesh::valid_set_id(Set_ID id, Entity_kind kind) const
{
  if (!geometric_model_.get()) return false;
  Teuchos::RCP<const AmanziGeometry::Region> rgn;
  try {
    rgn = geometric_model_->FindRegion(id);
  } catch (...) {
    return false;
  }
  return valid_set_name(rgn->name(), kind);
}


//---------------------------------------------------------
// Is there a set with this name and entity type
//---------------------------------------------------------
bool Mesh::valid_set_name(std::string name, Entity_kind kind) const
{
  if (!geometric_model_.get()) {
    Errors::Message mesg("Mesh sets not enabled because mesh was created without reference to a geometric model");
    Exceptions::amanzi_throw(mesg);
  }

  Teuchos::RCP<const AmanziGeometry::Region> rgn;
  try {
    rgn = geometric_model_->FindRegion(name);
  } catch (...) {
    return false;
  }

  unsigned int rdim = rgn->manifold_dimension();

  // For regions of type Color Function, the dimension
  // parameter is not guaranteed to be correct
  if (rgn->type() == AmanziGeometry::COLORFUNCTION) return true;

  // For regions of type Labeled set, extract some more info and verify
  if (rgn->type() == AmanziGeometry::LABELEDSET) {
    auto lsrgn = Teuchos::rcp_dynamic_cast<const AmanziGeometry::RegionLabeledSet>(rgn);
    std::string entity_type = lsrgn->entity_str();

    if (parent() == Teuchos::null) {
      if ((kind == CELL && entity_type == "CELL") ||
          (kind == FACE && entity_type == "FACE") ||
          (kind == EDGE && entity_type == "EDGE") ||
          (kind == NODE && entity_type == "NODE"))
        return true;
    } else {
      if ((kind == CELL && entity_type == "FACE") ||
          (kind == CELL && entity_type == "CELL") ||
          (kind == FACE && entity_type == "FACE") ||
          (kind == NODE && entity_type == "NODE")) {
        // guaranteed we can call, though it may be empty.
        return true;
      }
    } 
    return false;
  }

  // If we are looking for a cell set the region has to be
  // of the same topological dimension as the cells or it
  // has to be a point region
  if (kind == CELL && (rdim >= manifold_dim_ || rdim == 1 || rdim == 0)) return true;
  // If we are looking for a side set, the region has to be
  // one topological dimension less than the cells
  if (kind == FACE && (rdim >= manifold_dim_-1 || rdim == 0)) return true;

  // If we are looking for a node set, the region can be of any
  // dimension upto the spatial dimension of the domain
  if (kind == NODE) return true;

  return false;
}


//---------------------------------------------------------
// TBW
//---------------------------------------------------------
void Mesh::get_set_entities(const std::string& setname,
                            const Entity_kind kind,
                            const Parallel_type ptype,
                            Entity_ID_List *entids) const
{
  std::vector<double> vofs;
  get_set_entities_and_vofs(setname, kind, ptype, entids, &vofs);
}


//---------------------------------------------------------
// TBW
//---------------------------------------------------------
void Mesh::get_set_entities_box_vofs_(
    Teuchos::RCP<const AmanziGeometry::Region> region,
    const Entity_kind kind, 
    const Parallel_type ptype, 
    std::vector<Entity_ID>* setents,
    std::vector<double>* volume_fractions) const
{
  AMANZI_ASSERT(volume_fractions != NULL);

  setents->clear();
  volume_fractions->clear();

  switch (kind) {      
  case CELL:
  {
    std::string name = region->name() + "_cell";

    if (region_ids.find(name) != region_ids.end()) {
      *setents = region_ids[name];
      *volume_fractions = region_vofs[name];
    } else {
      int ncells = num_entities(CELL, ptype);  
      double volume;

      Entity_ID_List faces, cnodes, fnodes;
      std::vector<int> dirs;
      std::vector<AmanziGeometry::Point> polytope_nodes;
      std::vector<std::vector<int> > polytope_faces;

      for (int c = 0; c < ncells; ++c) {
        cell_get_coordinates(c, &polytope_nodes);

        if (space_dimension() == 3) { 
          cell_get_nodes(c, &cnodes);
          cell_get_faces_and_dirs(c, &faces, &dirs);
          int nfaces = faces.size();

          polytope_faces.clear();
          polytope_faces.resize(nfaces);
          for (int n = 0; n < nfaces; ++n) {
            face_get_nodes(faces[n], &fnodes);
            int nnodes = fnodes.size();

            for (int i = 0; i < nnodes; ++i) {
              int j = (dirs[n] > 0) ? i : nnodes - i - 1; 
              int pos = std::distance(cnodes.begin(), std::find(cnodes.begin(), cnodes.end(), fnodes[j]));
              polytope_faces[n].push_back(pos);
            }
          }
        }

        if ((volume = region->intersect(polytope_nodes, polytope_faces)) > 0.0) {
          setents->push_back(c);
          if (region->type()==AmanziGeometry::LINE_SEGMENT) volume_fractions->push_back(volume);
          else volume_fractions->push_back(volume /cell_volume(c));
        }
      }

      region_ids[name] = *setents;
      region_vofs[name] = *volume_fractions;
    }
    break;
  }

  case FACE:
  {
    std::string name = region->name() + "_face";

    if (region_ids.find(name) != region_ids.end()) {
      *setents = region_ids[name];
      *volume_fractions = region_vofs[name];
    } else {
      int nfaces = num_entities(FACE, ptype);
      double area;
        
      std::vector<AmanziGeometry::Point> polygon;

      for (int f = 0; f < nfaces; ++f) {
        face_get_coordinates(f, &polygon);
        if ((area = region->intersect(polygon)) > 0.0) {
          setents->push_back(f);
          volume_fractions->push_back(area / face_area(f));
        }
      }

      region_ids[name] = *setents;
      region_vofs[name] = *volume_fractions;
    }
    break;
  }

  case EDGE:
  {
    int nedges = num_entities(EDGE, ptype);
        
    for (int e = 0; e < nedges; ++e) {
      if (region->inside(edge_centroid(e))) {
        setents->push_back(e);
      }
    }
    break;
  }

  case NODE:
  {
    int nnodes = num_entities(NODE, ptype);
        
    AmanziGeometry::Point xv(space_dimension());

    for (int v = 0; v < nnodes; ++v) {
      node_get_coordinates(v, &xv);
      if (region->inside(xv)) {
        setents->push_back(v);
      }
    }
    break;
  }

  default:
    break;
  }

  // Check if no processor got any mesh entities
  int nents, nents_tmp = setents->size();
  get_comm()->SumAll(&nents_tmp, &nents, 1);
  if (nents == 0) {
    Errors::Message msg;
    msg << "Could not retrieve any mesh entities for set \"" << region->name() << "\".\n";
    Exceptions::amanzi_throw(msg);
  }
}


//---------------------------------------------------------
// Generic implemnetation of set routines.
//---------------------------------------------------------
unsigned int Mesh::get_set_size(const std::string& setname, 
                                const Entity_kind kind, 
                                const Parallel_type ptype) const 
{
  Entity_ID_List setents;
  std::vector<double> vofs;

  get_set_entities_and_vofs(setname, kind, ptype, &setents, &vofs);
  
  return setents.size();
}

}  // namespace AmanziMesh
}  // namespace Amanzi
