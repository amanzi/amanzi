/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Markus Berndt

   Interface for the Volumetric Deformation PK.

   <ParameterList name="volumetric deformation">
   <Parameter name="PK model" type="string" value="Prescibed Mesh Deformation"/>
   <Parameter name="Deformation method" type="string" value="method name"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#ifndef PKS_VOLUMETRIC_DEFORMATION_HH_
#define PKS_VOLUMETRIC_DEFORMATION_HH_

#include "CompositeMatrix.hh"
#include "CompositeVectorFunction.hh"
#include "Function.hh"

#include "pk_factory.hh"
#include "pk_default_base.hh"
#include "pk_physical_base.hh"
#include "MatrixVolumetricDeformation.hh"

namespace Amanzi {
namespace Deform {

class VolumetricDeformation : public PKPhysicalBase {

 public:

  VolumetricDeformation(Teuchos::Ptr<State> S, const Teuchos::RCP<Teuchos::ParameterList>& plist,
                        Teuchos::ParameterList& FElist,
                        const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~VolumetricDeformation() {}

  // ConstantTemperature is a PK
  // -- Setup data
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // -- advance via one of a few methods
  virtual bool advance(double dt);

  virtual double get_dt() {
    return dt_;
  }

 private:
  Key poro_key_;

  // strategy for calculating nodal deformation from cell volume change
  enum DeformStrategy {
    DEFORM_STRATEGY_GLOBAL_OPTIMIZATION,
    DEFORM_STRATEGY_MSTK,
    DEFORM_STRATEGY_AVERAGE
  };
  DeformStrategy strategy_;

  // fixed regions (bottom faces/nodes)
  std::vector<std::string> fixed_regions_;
  std::string fixed_region_type_;

  // function describing d(cv)/dT
  enum DeformMode {
    DEFORM_MODE_THAW_FRONT,
    DEFORM_MODE_DVDT,
    DEFORM_MODE_SATURATION
  };
  DeformMode deform_mode_;
  std::string deform_region_;
  double  deform_value_;
  Teuchos::RCP<Function> thaw_front_func_;
  double min_vol_frac_, min_S_liq_;
  Teuchos::RCP<Functions::CompositeVectorFunction> deform_func_;
  double dt_;

  // meshes
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf3d_mesh_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_nc_;
  Teuchos::RCP<AmanziMesh::Mesh> surf_mesh_nc_;
  Teuchos::RCP<AmanziMesh::Mesh> surf3d_mesh_nc_;

  // operator
  bool global_solve_;
  Teuchos::RCP<CompositeMatrix> operator_;
  Teuchos::RCP<Operators::MatrixVolumetricDeformation> def_matrix_;

  // factory registration
  static RegisteredPKFactory<VolumetricDeformation> reg_;
};

} // namespace
} // namespace

#endif
