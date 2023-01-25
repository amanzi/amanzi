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

PipeFlow_PK::PipeFlow_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln)
                  : ShallowWater_PK(pk_tree, glist, S, soln){

  hydrostatic_pressure_force_type_ = sw_list_->get<std::string>("hydrostatic pressure force type", "pipe flow");

  pipe_diameter_ = sw_list_->get<double>("pipe diameter", 1.0);

}

void PipeFlow_PK::Initialize(){

  // gravity
  g_ = norm(S_->Get<AmanziGeometry::Point>("gravity"));

  // numerical flux
  Teuchos::ParameterList model_list;
  model_list.set<std::string>("numerical flux", sw_list_->get<std::string>("numerical flux", "central upwind"))
    .set<double>("gravity", g_);
  model_list.set<std::string>("numerical flux", sw_list_->get<std::string>("numerical flux", "central upwind"))
    .set<double>("pipe diameter", pipe_diameter_);
  model_list.set<std::string>("numerical flux", sw_list_->get<std::string>("numerical flux", "central upwind"))
    .set<std::string>("hydrostatic pressure force type", hydrostatic_pressure_force_type_);
  NumericalFluxFactory nf_factory;
  numerical_flux_ = nf_factory.Create(model_list);

}


}  // namespace ShallowWater
}  // namespace Amanzi

