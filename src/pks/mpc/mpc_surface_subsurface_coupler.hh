/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

This class should never be instantiated -- it just provides a base class for
multiple coupler types.
------------------------------------------------------------------------- */

#ifndef PKS_MPC_SURFACE_SUBSURFACE_COUPLER_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_COUPLER_HH_

#include "pk_physical_bdf_base.hh"
#include "strong_mpc.hh"

namespace Amanzi {

namespace Operators {
  class MatrixMFD_Surf;
  class MatrixMFD_TPFA;
}

class MPCSurfaceSubsurfaceCoupler : public StrongMPC<PKPhysicalBDFBase> {

 public:
  MPCSurfaceSubsurfaceCoupler(const Teuchos::RCP<Teuchos::ParameterList>& plist,
          Teuchos::ParameterList& FElist,
          const Teuchos::RCP<TreeVector>& soln);

  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);

  Teuchos::RCP<Operators::MatrixMFD_Surf> coupled_preconditioner() {
    return mfd_preconditioner_; }

  // Accessor for debugger, for use by coupling MPCs
  Teuchos::RCP<Debugger> domain_debugger() { return domain_db_; }

  // Accessor for debugger, for use by coupling MPCs
  Teuchos::RCP<Debugger> surface_debugger() { return surf_db_; }


 protected:
  bool modify_predictor_copy_surf_to_subsurf_(double h, Teuchos::RCP<TreeVector> up);
  bool modify_predictor_copy_subsurf_to_surf_(double h, Teuchos::RCP<TreeVector> up);



 protected:
  // mesh info
  Key domain_mesh_key_;
  Key surf_mesh_key_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;

  // pk info
  int domain_pk_index_;
  Key domain_pk_name_;
  int surf_pk_index_;
  Key surf_pk_name_;
  Teuchos::RCP<PKPhysicalBDFBase> surf_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> domain_pk_;

  // PC info
  Teuchos::RCP<Operators::MatrixMFD_Surf> mfd_preconditioner_;
  Teuchos::RCP<Operators::MatrixMFD_TPFA> surf_preconditioner_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;

};

} // namespace


#endif
