/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Registration of MPC PKs.
*/

#include "PK_MPCStrong.hh"
#include "PK_MPCSubcycled.hh"

#include "FlowEnergy_PK.hh"
#include "FlowEnergyMatrixFracture_PK.hh"
#include "FlowMatrixFracture_PK.hh"
#include "FlowReactiveTransport_PK.hh"

#include "ChemistryMatrixFracture_PK.hh"
#include "MultiphaseMatrixFracture_PK.hh"
#include "ReactiveTransport_PK.hh"
#include "ReactiveTransportMatrixFracture_PK.hh"
#include "ShallowWaterTransport_PK.hh"
#include "SurfaceSubsurface_PK.hh"
#include "TransportMatrixFracture_PK.hh"
#include "TransportMatrixFractureImplicit_PK.hh"

namespace Amanzi {

template <>
RegisteredPKFactory<PK_MPCStrong<PK_BDF>> PK_MPCStrong<PK_BDF>::reg_("mpc strong");
RegisteredPKFactory<PK_MPCSubcycled> PK_MPCSubcycled::reg_("mpc subcycled");
RegisteredPKFactory<PK_MPCWeak> PK_MPCWeak::reg_("mpc weak");

RegisteredPKFactory<FlowEnergy_PK> FlowEnergy_PK::reg_("thermal flow");
RegisteredPKFactory<FlowEnergyMatrixFracture_PK>
  FlowEnergyMatrixFracture_PK::reg_("thermal flow matrix fracture");
RegisteredPKFactory<FlowReactiveTransport_PK>
  FlowReactiveTransport_PK::reg_("flow reactive transport");
RegisteredPKFactory<FlowMatrixFracture_PK> FlowMatrixFracture_PK::reg_("darcy matrix fracture");
RegisteredPKFactory<MultiphaseMatrixFracture_PK>
  MultiphaseMatrixFracture_PK::reg_("multiphase matrix fracture");

RegisteredPKFactory<ReactiveTransport_PK> ReactiveTransport_PK::reg_("reactive transport");

// integrated matrix-fracture models
RegisteredPKFactory<TransportMatrixFracture_PK>
  TransportMatrixFracture_PK::reg_("transport matrix fracture");
RegisteredPKFactory<ChemistryMatrixFracture_PK>
  ChemistryMatrixFracture_PK::reg_("chemistry matrix fracture");

RegisteredPKFactory<TransportMatrixFractureImplicit_PK>
  TransportMatrixFractureImplicit_PK::reg_("transport matrix fracture implicit");
RegisteredPKFactory<ReactiveTransportMatrixFracture_PK>
  ReactiveTransportMatrixFracture_PK::reg_("reactive transport matrix fracture");

// integrated surface-subsurface models
RegisteredPKFactory<SurfaceSubsurface_PK> SurfaceSubsurface_PK::reg_("surface subsurface");
RegisteredPKFactory<ShallowWaterTransport_PK>
  ShallowWaterTransport_PK::reg_("shallow water transport");

} // namespace Amanzi
