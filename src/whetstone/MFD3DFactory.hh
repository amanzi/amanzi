/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Factory of mimetic methods.
*/

#ifndef AMANZI_MFD3D_FACTORY_HH_
#define AMANZI_MFD3D_FACOTRY_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "MFD3D.hh"
#include "MFD3D_Diffusion.hh"
#include "MFD3D_Generalized_Diffusion.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3DFactory {
 public:
  explicit MFD3DFactory() {};
  ~MFD3DFactory() {};

  // geometry methods
  Teuchos::RCP<MFD3D> Create(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             const std::string& physics,
                             const std::string& schema_unique_name);
};


/* ******************************************************************
* Implementation of factory
****************************************************************** */
Teuchos::RCP<MFD3D> MFD3DFactory::Create(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    const std::string& physics, const std::string& schema_unique_name)
{
  if (physics == "diffusion") {
    Teuchos::RCP<MFD3D_Diffusion> mfd = Teuchos::rcp(new MFD3D_Diffusion(mesh));
    return mfd;
  }
  else if (physics == "diffusion generalized") {
    Teuchos::RCP<MFD3D_Generalized_Diffusion> mfd = 
        Teuchos::rcp(new MFD3D_Generalized_Diffusion(mesh));
    return mfd;
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif

