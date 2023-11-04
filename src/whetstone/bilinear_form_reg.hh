/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

*/

#include "DG_Modal.hh"
#include "MFD3D_BernardiRaugel.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_CrouzeixRaviartAnyOrder.hh"
#include "MFD3D_CrouzeixRaviartSerendipity.hh"
#include "MFD3D_Diffusion.hh"
#include "MFD3D_Diffusion_CurvedFace.hh"
#include "MFD3D_GeneralizedDiffusion.hh"
#include "MFD3D_Electromagnetics.hh"
#include "MFD3D_Elasticity.hh"
#include "MFD3D_ElasticityWeakSymmetry.hh"
#include "MFD3D_Lagrange.hh"
#include "MFD3D_LagrangeAnyOrder.hh"
#include "MFD3D_LagrangeSerendipity.hh"

namespace Amanzi {
namespace WhetStone {

RegisteredFactory<MFD3D_BernardiRaugel> MFD3D_BernardiRaugel::reg_("BernardiRaugel");
RegisteredFactory<MFD3D_CrouzeixRaviart> MFD3D_CrouzeixRaviart::reg_("CrouzeixRaviart");
RegisteredFactory<MFD3D_CrouzeixRaviartAnyOrder>
  MFD3D_CrouzeixRaviartAnyOrder::reg_("CrouzeixRaviart high order");
RegisteredFactory<MFD3D_CrouzeixRaviartSerendipity>
  MFD3D_CrouzeixRaviartSerendipity::reg_("CrouzeixRaviart serendipity");

RegisteredFactory<MFD3D_Diffusion> MFD3D_Diffusion::reg_("diffusion");
RegisteredFactory<MFD3D_GeneralizedDiffusion>
  MFD3D_GeneralizedDiffusion::reg_("diffusion generalized");
RegisteredFactory<MFD3D_Diffusion_CurvedFace>
  MFD3D_Diffusion_CurvedFace::reg_("diffusion curved face");

RegisteredFactory<MFD3D_Elasticity> MFD3D_Elasticity::reg_("elasticity");
RegisteredFactory<MFD3D_ElasticityWeakSymmetry>
  MFD3D_ElasticityWeakSymmetry::reg_("elasticity weak symmetry");
RegisteredFactory<MFD3D_Electromagnetics> MFD3D_Electromagnetics::reg_("electromagnetics");
RegisteredFactory<MFD3D_Lagrange> MFD3D_Lagrange::reg_("Lagrange");
RegisteredFactory<MFD3D_LagrangeAnyOrder> MFD3D_LagrangeAnyOrder::reg_("Lagrange high order");
RegisteredFactory<MFD3D_LagrangeSerendipity>
  MFD3D_LagrangeSerendipity::reg_("Lagrange serendipity");

RegisteredFactory<DG_Modal> DG_Modal::reg_("dg modal");

} // namespace WhetStone
} // namespace Amanzi
