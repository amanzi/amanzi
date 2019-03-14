/*
  This is the data structures component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATORS_SUPER_MAP_SURFACE_HH_
#define AMANZI_OPERATORS_SUPER_MAP_SURFACE_HH_


#include "SuperMapLumped.hh"

/*

This is a class takes a SuperMapLumped created on a subsurface mesh and allows it to
be used by Operators defined on a surface mesh to assemble into the subsurface
parent mesh's unknowns.

*/

namespace Amanzi {
namespace Operators {

class SuperMapLumpedSurface : public SuperMapLumped {
  
 public:
  // Constructor
  SuperMapLumpedSurface(const SuperMapLumped& map,
                  const Teuchos::RCP<const AmanziMesh::Mesh>& surf_mesh);

 protected:
  virtual const std::vector<int>& CreateIndices_(const std::string& compname, int dofnum, bool ghosted) const;
  
 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;
};

} // namespace Operators
} // namespace Amanzi

#endif
