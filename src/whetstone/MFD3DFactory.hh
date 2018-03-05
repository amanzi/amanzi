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
#define AMANZI_MFD3D_FACTORY_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "BilinearForm.hh"
#include "DG_Modal.hh"
#include "MFD3D.hh"
#include "MFD3D_BernardiRaugel.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Diffusion.hh"
#include "MFD3D_Generalized_Diffusion.hh"
#include "MFD3D_Lagrange.hh"
#include "MFD3D_LagrangeSerendipity.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3DFactory {
 public:
  explicit MFD3DFactory() {};
  ~MFD3DFactory() {};

  // select numerical scheme using its name and order 
  Teuchos::RCP<BilinearForm> Create(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                    const std::string& method, int method_order);
};


/* ******************************************************************
* Implementation of factory
****************************************************************** */
inline
Teuchos::RCP<BilinearForm> MFD3DFactory::Create(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    const std::string& method, int method_order)
{
  if (method == "diffusion") {
    Teuchos::RCP<MFD3D_Diffusion> mfd = Teuchos::rcp(new MFD3D_Diffusion(mesh));
    return mfd;
  }
  else if (method == "diffusion generalized") {
    Teuchos::RCP<MFD3D_Generalized_Diffusion> mfd = 
        Teuchos::rcp(new MFD3D_Generalized_Diffusion(mesh));
    return mfd;
  }
  else if (method == "BernardiRaugel") {
    Teuchos::RCP<MFD3D_BernardiRaugel> mfd = Teuchos::rcp(new MFD3D_BernardiRaugel(mesh));
    return mfd;
  } 
  else if (method == "CrouzeixRaviart") {
    Teuchos::RCP<MFD3D_CrouzeixRaviart> mfd = Teuchos::rcp(new MFD3D_CrouzeixRaviart(mesh));
    mfd->set_order(method_order);
    return mfd;
  } 
  else if (method == "Lagrange") {
    Teuchos::RCP<MFD3D_Lagrange> mfd = Teuchos::rcp(new MFD3D_Lagrange(mesh));
    mfd->set_order(method_order);
    return mfd;
  } 
  else if (method == "Lagrange serendipity") {
    Teuchos::RCP<MFD3D_LagrangeSerendipity> mfd = Teuchos::rcp(new MFD3D_LagrangeSerendipity(mesh));
    mfd->set_order(method_order);
    return mfd;
  } 
  else if (method == "dg modal") {
    Teuchos::RCP<DG_Modal> mfd = Teuchos::rcp(new DG_Modal(mesh));
    mfd->set_order(method_order);
    return mfd;
  }
  else if (method == "dg modal simple") {
    Teuchos::RCP<DG_Modal> mfd = Teuchos::rcp(new DG_Modal(mesh));
    mfd->set_order(method_order);
    mfd->set_basis(TAYLOR_BASIS_NATURAL);
    return mfd;
  }

  return Teuchos::null;
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif

