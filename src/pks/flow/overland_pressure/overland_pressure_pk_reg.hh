/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Author: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "overland_pressure.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory_ATS<OverlandPressureFlow> OverlandPressureFlow::reg_("overland flow, pressure basis");


} // namespace
} // namespace

