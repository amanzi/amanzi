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

#include "PK.hh"
#include "PK_Factory.hh"
#include "pk_physical_default.hh"
#include "MatrixVolumetricDeformation.hh"

namespace Amanzi {
namespace Deform {

class VolumetricDeformation : public PK_Physical_Default {

 public:

  VolumetricDeformation(Teuchos::ParameterList& pk_tree,
                        const Teuchos::RCP<Teuchos::ParameterList>& glist,
                        const Teuchos::RCP<State>& S,
                        const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~VolumetricDeformation() {}

  // ConstantTemperature is a PK
  // -- Setup data
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {}

  // -- advance via one of a few methods
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  virtual double get_dt() {
    return dt_max_;
  }

  virtual void set_dt(double dt) {
    dt_ = dt;
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
    DEFORM_MODE_SATURATION,
    DEFORM_MODE_STRUCTURAL
  };
  DeformMode deform_mode_;
  double overpressured_limit_;

  std::string deform_region_;

  // DEFORM_MODE_THAW_FRONT
  Teuchos::RCP<Function> thaw_front_func_;

  // DEFORM_MODE_DVDT
  Teuchos::RCP<Functions::CompositeVectorFunction> deform_func_;

  // DEFORM_MODE_SATURATION
  double min_vol_frac_, min_S_liq_;

  // DEFORM_MODE_STRUCTURAL
  double time_scale_, structural_vol_frac_;

  double dt_, dt_max_;

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
