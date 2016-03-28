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

#include "Mesh.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

void Mesh::get_set_entities_box_vof(
    Teuchos::RCP<const AmanziGeometry::Region> region,
    const Entity_kind kind, 
    const Parallel_type ptype, 
    std::vector<Entity_ID>* setents) const
{
  /*
  std::cout << "region=" << region << std::endl;
  std::cout << "kind=" << kind << std::endl;
  std::cout << "ptype=" << ptype << std::endl;
  std::cout << "data=" << setents << std::endl;
  */

  setents->clear();

  int ncells, nfaces;
  switch (kind) {      
  case CELL:
    ncells = num_entities(CELL, USED);              
    for (int c = 0; c < ncells; ++c) {
      if (region->inside(cell_centroid(c))) {
        setents->push_back(c);
      }
    }
    break;

  case FACE:
    nfaces = num_entities(FACE, USED);
        
    for (int f = 0; f < nfaces; ++f) {
      if (region->inside(face_centroid(f))) {
        // std::cout << "found:" << f << "  xf=" << face_centroid(f) << std::endl;
        setents->push_back(f);
      }
    }
    break;

  default:
    break;
  }
}

}  // namespace AmanziMesh
}  // namespace Amanzi
