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

#include <string>
#include <vector>

#include "errors.hh"
#include "Point.hh"
#include "Region.hh"

#include "Mesh.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

void Mesh::get_set_entities_box_vof(
    Teuchos::RCP<const AmanziGeometry::Region> region,
    const Entity_kind kind, 
    const Parallel_type ptype, 
    std::vector<Entity_ID>* setents,
    std::vector<double>* volume_fractions) const
{
  /*
  std::cout << "region=" << region << std::endl;
  std::cout << "kind=" << kind << std::endl;
  std::cout << "ptype=" << ptype << std::endl;
  std::cout << "data=" << setents << std::endl;
  */

  setents->clear();

  int ncells, nfaces, nedges, nnodes;
  AmanziGeometry::Point xv(space_dimension());

  switch (kind) {      
  case CELL:
    ncells = num_entities(CELL, ptype);  
    for (int c = 0; c < ncells; ++c) {
      if (region->inside(cell_centroid(c))) {
        setents->push_back(c);
      }
    }
    break;

  case FACE:
    nfaces = num_entities(FACE, ptype);
        
    for (int f = 0; f < nfaces; ++f) {
      if (region->inside(face_centroid(f))) {
        // std::cout << "found:" << f << "  xf=" << face_centroid(f) << std::endl;
        setents->push_back(f);
      }
    }
    break;

  case EDGE:
    nedges = num_entities(EDGE, ptype);
        
    for (int e = 0; e < nedges; ++e) {
      if (region->inside(edge_centroid(e))) {
        setents->push_back(e);
      }
    }
    break;

  case NODE:
    nnodes = num_entities(NODE, ptype);
        
    for (int v = 0; v < nnodes; ++v) {
      node_get_coordinates(v, &xv);
      if (region->inside(xv)) {
        setents->push_back(v);
      }
    }
    break;

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

}  // namespace AmanziMesh
}  // namespace Amanzi
