/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Ugly hackjob to enable direct evaluation of the full model, on a single
  WRM/region.  This is bypassing much of the "niceness" of the framework, but
  seems necessary for solving a cell-wise correction equation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_PERMAFROST_MODEL_
#define AMANZI_PERMAFROST_MODEL_

#include "Teuchos_ParameterList.hpp"

#include "tensor.hh"
#include "Point.hh"

#include "ewc_model.hh"

namespace Amanzi {

namespace Flow { namespace FlowRelations { class WRMPermafrostModel; } }
namespace Flow { namespace FlowRelations { class PCIceWater;
                                           class PCLiqAtm;
                                           class CompressiblePorosityModel; } }
namespace Energy { namespace EnergyRelations { class IEM; class IEMWaterVapor;} }
namespace Relations { class EOS; class VaporPressureRelation; }


class PermafrostModel : public EWCModel {

 public:
  PermafrostModel() {}
  virtual void InitializeModel(const Teuchos::Ptr<State>& S);
  virtual void UpdateModel(const Teuchos::Ptr<State>& S);
  virtual int Evaluate(double T, double p, double base_poro, double& energy, double& wc);
  virtual int InverseEvaluate(double energy, double wc, double base_poro, double& T, double& p, bool verbose=false);
  virtual int InverseEvaluateEnergy(double energy, double p, double base_poro, double& T);
  virtual int EvaluateSaturations(double T, double p, double base_poro, double& s_gas, double& s_liq, double& s_ice);

 protected:
  bool IsSetUp_();

  int EvaluateEnergyAndWaterContent_(double T, double p, double base_poro,
          AmanziGeometry::Point& result);

  int EvaluateEnergyAndWaterContentAndJacobian_(double T, double p, double base_poro,
          AmanziGeometry::Point& result, WhetStone::Tensor& jac);

  int EvaluateEnergyAndWaterContentAndJacobian_FD_(double T, double p, double base_poro,
          AmanziGeometry::Point& result, WhetStone::Tensor& jac);

 protected:
  Teuchos::RCP<Flow::FlowRelations::WRMPermafrostModel> wrm_;
  Teuchos::RCP<Relations::EOS> liquid_eos_;
  Teuchos::RCP<Relations::EOS> gas_eos_;
  Teuchos::RCP<Relations::EOS> ice_eos_;
  Teuchos::RCP<Flow::FlowRelations::PCIceWater> pc_i_;
  Teuchos::RCP<Flow::FlowRelations::PCLiqAtm> pc_l_;
  Teuchos::RCP<Relations::VaporPressureRelation> vpr_;
  Teuchos::RCP<Energy::EnergyRelations::IEM> liquid_iem_;
  Teuchos::RCP<Energy::EnergyRelations::IEMWaterVapor> gas_iem_;
  Teuchos::RCP<Energy::EnergyRelations::IEM> ice_iem_;
  Teuchos::RCP<Energy::EnergyRelations::IEM> rock_iem_;
  Teuchos::RCP<Flow::FlowRelations::CompressiblePorosityModel> poro_model_;
  double rho_rock_;
  double p_atm_;

};



} //namespace


#endif
