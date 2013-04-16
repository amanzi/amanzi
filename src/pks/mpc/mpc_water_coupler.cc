/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

To be used with either MPCSurfaceSubsurfaceDirichletCoupler or
MPCSurfaceSubsurfaceFluxCoupler.

------------------------------------------------------------------------- */

#include "matrix_mfd_tpfa.hh"
#include "pk_physical_bdf_base.hh"

#include "mpc_surface_subsurface_dirichlet_coupler.hh"
#include "mpc_surface_subsurface_flux_coupler.hh"
#include "mpc_water_coupler.hh"


namespace Amanzi {

template<>
RegisteredPKFactory< MPCWaterCoupler<MPCSurfaceSubsurfaceDirichletCoupler> >
MPCWaterCoupler<MPCSurfaceSubsurfaceDirichletCoupler>::reg_("surface-subsurface Dirichlet water coupler");


template<>
RegisteredPKFactory< MPCWaterCoupler<MPCSurfaceSubsurfaceFluxCoupler> >
MPCWaterCoupler<MPCSurfaceSubsurfaceFluxCoupler>::reg_("surface-subsurface flux water coupler");

}
