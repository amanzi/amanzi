/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

To be used with either MPCSurfaceSubsurfaceDirichletCoupler or
MPCSurfaceSubsurfaceFluxCoupler.

------------------------------------------------------------------------- */

#ifndef PKS_MPC_WATER_FLUX_COUPLER_HH_
#define PKS_MPC_WATER_FLUX_COUPLER_HH_


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

#include "TreeVector.hh"
#include "CompositeVector.hh"
#include "pk_factory.hh"
#include "pk_default_base.hh"
#include "FieldEvaluator.hh"
#include "MatrixMFD_Surf.hh"
#include "mpc_surface_subsurface_flux_coupler.hh"
#include "mpc_water_coupler.hh"

namespace Amanzi {

class MPCWaterFluxCoupler : public MPCWaterCoupler<MPCSurfaceSubsurfaceFluxCoupler> {
 public:

  MPCWaterFluxCoupler(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                      Teuchos::ParameterList& FElist,
                      const Teuchos::RCP<TreeVector>& soln) :
      MPCWaterCoupler<MPCSurfaceSubsurfaceFluxCoupler>(plist, FElist, soln),
      PKDefaultBase(plist, FElist, soln)
  {}

  // evaluate the flux
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // virtual bool PreconPostprocess_(Teuchos::RCP<const TreeVector> u,
  //         Teuchos::RCP<TreeVector> Pu);
  
 private:
  static RegisteredPKFactory<MPCWaterFluxCoupler> reg_;


};

} // namespapce

#endif

