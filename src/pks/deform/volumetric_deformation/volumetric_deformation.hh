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

#include "composite_matrix.hh"
#include "composite_vector_function.hh"

#include "pk_factory.hh"
#include "pk_default_base.hh"
#include "pk_physical_base.hh"
#include "MatrixVolumetricDeformation.hh"

namespace Amanzi {
namespace Deform {

class VolumetricDeformation : public PKPhysicalBase {

 public:

  VolumetricDeformation(Teuchos::ParameterList& plist,
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
    return 1.0e99;
  }

 private:
  Key poro_key_;

  // function describing d(cv)/dT
  Teuchos::RCP<Functions::CompositeVectorFunction> deform_func_;

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
