/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived WeakMPC class.  Provides only the advance()
method missing from MPC.hh.  In weak coupling, we simply loop over the
sub-PKs, calling their advance() methods and returning failure if any fail.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "weak_mpc.hh"

namespace Amanzi {

RegisteredPKFactory<WeakMPC> WeakMPC::reg_("weak MPC");

} // namespace
/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived StrongMPC class.  Is both a PK and a Model
Evalulator, providing needed methods for BDF time integration of the coupled
system.

Completely automated and generic to any sub PKs, this uses a block diagonal
preconditioner.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */


#include "strong_mpc.hh"

namespace Amanzi {

template<>
RegisteredPKFactory<StrongMPC<PKBDFBase> > StrongMPC<PKBDFBase>::reg_("strong MPC");

} // namespace
#include "mpc_subsurface.hh"
Amanzi::RegisteredPKFactory<Amanzi::MPCSubsurface> Amanzi::MPCSubsurface::reg_("subsurface permafrost");

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
#include "mpc_permafrost4.hh"

namespace Amanzi {

RegisteredPKFactory<MPCPermafrost4> MPCPermafrost4::reg_("permafrost model");

} // namespace
#include "mpc_coupled_water.hh"

namespace Amanzi {

RegisteredPKFactory<MPCCoupledWater> MPCCoupledWater::reg_("coupled water");

}
#include "weak_mpc_semi_coupled.hh"

namespace Amanzi {

RegisteredPKFactory<WeakMPCSemiCoupled> WeakMPCSemiCoupled::reg_("weak MPC semi coupled");

}
