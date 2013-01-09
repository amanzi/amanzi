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

#include "tensor.hpp"
#include "Point.hh"

namespace Amanzi {

namespace Flow { namespace FlowRelations { class WRM; } }
namespace Flow { namespace FlowRelations { class PCIceWater; } }
namespace Energy { namespace EnergyRelations { class IEM; class IEMWaterVapor;} }
namespace Relations { class EOS; class VaporPressureRelation; }


class PermafrostModel {

 public:
  PermafrostModel() :
      p_atm_(-1.e99),
      rho_rock_(-1.e99) {}

  void set_WRM(const Teuchos::RCP<Flow::FlowRelations::WRM>& wrm) { wrm_ = wrm; }
  void set_liquid_EOS(const Teuchos::RCP<Relations::EOS>& eos) { liquid_eos_ = eos; }
  void set_gas_EOS(const Teuchos::RCP<Relations::EOS>& eos) { gas_eos_ = eos; }
  void set_ice_EOS(const Teuchos::RCP<Relations::EOS>& eos) { ice_eos_ = eos; }
  void set_vapor_pressure_relation(const Teuchos::RCP<Relations::VaporPressureRelation>& vpr) { vpr_ = vpr; }
  void set_pc_ice_water(const Teuchos::RCP<Flow::FlowRelations::PCIceWater>& pc_i) { pc_i_ = pc_i; }
  void set_p_atm(double p_atm) { p_atm_ = p_atm; }
  void set_liquid_IEM(const Teuchos::RCP<Energy::EnergyRelations::IEM>& iem) { liquid_iem_ = iem; }
  void set_gas_IEM(const Teuchos::RCP<Energy::EnergyRelations::IEMWaterVapor>& iem) { gas_iem_ = iem; }
  void set_ice_IEM(const Teuchos::RCP<Energy::EnergyRelations::IEM>& iem) { ice_iem_ = iem; }
  void set_rock_IEM(const Teuchos::RCP<Energy::EnergyRelations::IEM>& iem) { rock_iem_ = iem; }
  void set_rock_density(double rho_rock) { rho_rock_ = rho_rock; }

  bool IsSetUp();

  // required methods from the base class
  virtual int EvaluateEnergyAndWaterContent(double T, double p, double poro,
          AmanziGeometry::Point& result);

  virtual int EvaluateEnergyAndWaterContentAndJacobian(double T, double p, double poro,
          AmanziGeometry::Point& result, WhetStone::Tensor& jac);

  virtual int InverseEvaluate(double energy, double wc, double poro, double& T, double& p);

 protected:

  virtual int EvaluateEnergyAndWaterContentAndJacobian_FD_(double T, double p, double poro,
          AmanziGeometry::Point& result, WhetStone::Tensor& jac);

  Teuchos::RCP<Flow::FlowRelations::WRM> wrm_;
  Teuchos::RCP<Relations::EOS> liquid_eos_;
  Teuchos::RCP<Relations::EOS> gas_eos_;
  Teuchos::RCP<Relations::EOS> ice_eos_;
  Teuchos::RCP<Flow::FlowRelations::PCIceWater> pc_i_;
  Teuchos::RCP<Relations::VaporPressureRelation> vpr_;
  Teuchos::RCP<Energy::EnergyRelations::IEM> liquid_iem_;
  Teuchos::RCP<Energy::EnergyRelations::IEMWaterVapor> gas_iem_;
  Teuchos::RCP<Energy::EnergyRelations::IEM> ice_iem_;
  Teuchos::RCP<Energy::EnergyRelations::IEM> rock_iem_;

  double rho_rock_;
  double p_atm_;

};



} //namespace


#endif
