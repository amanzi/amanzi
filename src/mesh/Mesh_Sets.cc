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


  
  return setents.size();
}

}  // namespace AmanziMesh
}  // namespace Amanzi
