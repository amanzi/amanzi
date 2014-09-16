/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Markus Berndt

   Interface for the Prescribed Deformation PK.

   <ParameterList name="prescribed deformation">
   <Parameter name="PK model" type="string" value="Prescibed Mesh Deformation"/>
   <Parameter name="Deformation method" type="string" value="method name"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#ifndef PKS_PRESCRIBED_VOLUMETRIC_DEFORMATION_HH_
#define PKS_PRESCRIBED_VOLUMETRIC_DEFORMATION_HH_

#include "pk_factory.hh"
#include "pk_default_base.hh"
#include "pk_physical_base.hh"

namespace Amanzi {
namespace Deform {

class PrescribedVolumetricDeformation : public PKPhysicalBase {

 public:

  PrescribedVolumetricDeformation(const Teuchos::RCP<Teuchos::ParameterList>& plist,
          Teuchos::ParameterList& FElist,
          const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~PrescribedVolumetricDeformation() {}

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
    if (S_next_->time() < deform_time_end_) 
      return deform_time_end_ - S_next_->time();
    else 
      return 1.0e99;
  }

 private:

  Key poro_key_;

  std::vector<std::string>  bottom_surface_;
  std::string deform_region_;
  double deform_factor_;
  double deform_time_start_;
  double deform_time_end_;

  // A few options for advance
  bool advance_analytic_(double dt);

  // misc setup information
  Teuchos::ParameterList prescribed_volumetric_deformation_plist_;

  // factory registration
  static RegisteredPKFactory<PrescribedVolumetricDeformation> reg_;
};

} // namespace
} // namespace

#endif
