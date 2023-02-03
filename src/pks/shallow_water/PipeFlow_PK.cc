/*
  Pipe Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Pipe flow model inherited from shallow water model. 
*/

#include "PipeFlow_PK.hh"
#include "NumericalFluxFactory.hh"

namespace Amanzi {
namespace ShallowWater {

//PipeFlow_PK::PipeFlow_PK(Teuchos::ParameterList& pk_tree,
//                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
//                  const Teuchos::RCP<State>& S,
//                  const Teuchos::RCP<TreeVector>& soln)
//                  : PK(pk_tree, glist, S, soln),
//                  ShallowWater_PK(pk_tree, glist, S, soln){
PipeFlow_PK::PipeFlow_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln)
                  : ShallowWater_PK(pk_tree, glist, S, soln), 
                  PK(pk_tree, glist, S, soln){
    

  hydrostatic_pressure_force_type_ = sw_list_->get<std::string>("hydrostatic pressure force type", "pipe flow");

  pipe_diameter_ = sw_list_->get<double>("pipe diameter", 1.0);

  Manning_coeff_ = sw_list_->get<double>("Manning coefficient", 0.005);

  celerity_ = sw_list_->get<double>("celerity", 100); // m/s

}

//--------------------------------------------------------------------
// Discretization of the friction source term
//--------------------------------------------------------------------
double PipeFlow_PK::NumericalSourceFriction(double htc, double Bc, double qx)
{

  double WettedAngle = 2.0 * acos(1.0 - 2.0 * (htc - Bc) / pipe_diameter_);
  double WettedPerimeter = 0.5 * pipe_diameter_ * WettedAngle;
  double num = - g_ * Manning_coeff_ * Manning_coeff_ * pow(WettedPerimeter, 4.0/3.0) * fabs(qx) * qx;
  double denom = pow( (htc - Bc), 7.0/3.0);
  double S1 = num / denom;

  return S1;
}

void PipeFlow_PK::Initialize(){

  //// gravity
  //g_ = norm(S_->Get<AmanziGeometry::Point>("gravity"));

  ShallowWater_PK::Initialize();

  // numerical flux
  Teuchos::ParameterList model_list;
  model_list.set<std::string>("numerical flux", sw_list_->get<std::string>("numerical flux", "central upwind"))
    .set<double>("gravity", g_);
  model_list.set<std::string>("numerical flux", sw_list_->get<std::string>("numerical flux", "central upwind"))
    .set<double>("pipe diameter", pipe_diameter_);
  model_list.set<std::string>("numerical flux", sw_list_->get<std::string>("numerical flux", "central upwind"))
    .set<std::string>("hydrostatic pressure force type", hydrostatic_pressure_force_type_);
  model_list.set<std::string>("numerical flux", sw_list_->get<std::string>("numerical flux", "central upwind"))
    .set<double>("celerity", celerity_);
  NumericalFluxFactory nf_factory;
  numerical_flux_ = nf_factory.Create(model_list);

}


}  // namespace ShallowWater
}  // namespace Amanzi

