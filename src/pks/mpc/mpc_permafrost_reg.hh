/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

MPC for the Coupled Permafrost model.  This MPC sits at the top of the
subtree:

                    MPCPermafrost
                     /          \
                    /            \
                   /              \
         surf/subsurf            surf/subsurf
           water                   energy
         /      \                  /      \
        /        \                /        \
    flow/        flow/         energy/     energy/
  permafrost  icy_overland    threephase    surface_ice

------------------------------------------------------------------------- */
#include "mpc_permafrost.hh"
#include "mpc_permafrost2.hh"
#include "mpc_permafrost3.hh"

namespace Amanzi {

  //RegisteredPKFactory<MPCPermafrost> MPCPermafrost::reg_("permafrost model");
  //RegisteredPKFactory<MPCPermafrost2> MPCPermafrost2::reg_("new permafrost model");
RegisteredPKFactory<MPCPermafrost3> MPCPermafrost3::reg_("new permafrost model no SC");

} // namespace
