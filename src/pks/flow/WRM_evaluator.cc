/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The WRM Evaluator simply calls the WRM with the correct arguments.
*/

#include "FlowDefs.hh"
#include "WRM_BrooksCorey.hh"
#include "WRM_evaluator.hh"
#include "WRM_fake.hh"
#include "WRM_vanGenuchten.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist,
                           Teuchos::RCP<Teuchos::ParameterList> wrm_list,
                           Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    mesh_(mesh),
    SecondaryVariablesFieldEvaluator(plist)
{
  CreateWRM_(*wrm_list);
  CreateCell2Region_();

  InitializeFromPlist_();
}


WRMEvaluator::WRMEvaluator(const WRMEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    pressure_key_(other.pressure_key_),
    wrm_(other.wrm_),
    cell2region_(other.cell2region_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> WRMEvaluator::Clone() const {
  return Teuchos::rcp(new WRMEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void WRMEvaluator::InitializeFromPlist_() {
  // my key is for saturation
  my_keys_.push_back(plist_.get<std::string>("saturation key", "saturation_liquid"));

  // my dependency is pressure.
  pressure_key_ = plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(pressure_key_);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void WRMEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  Epetra_MultiVector& sat_c = *results[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pressure_key_)->ViewComponent("cell", false);

  int ncells = sat_c.MyLength();
  double patm = FLOW_PRESSURE_ATMOSPHERIC;

  for (int c = 0; c != ncells; ++c) {
    sat_c[0][c] = wrm_[(*cell2region_)[c]]->saturation(patm - pres_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void WRMEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const std::vector<Teuchos::Ptr<CompositeVector> > & results)
{
  Epetra_MultiVector& sat_c = *results[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pressure_key_)->ViewComponent("cell", false);

  int ncells = sat_c.MyLength();
  double patm = FLOW_PRESSURE_ATMOSPHERIC;

  for (int c = 0; c != ncells; ++c) {
    sat_c[0][c] = wrm_[(*cell2region_)[c]]->dSdPc(patm - pres_c[0][c]);
  }
}


/* ******************************************************************
* Create list of water retention models ordered by regions.
****************************************************************** */
void WRMEvaluator::CreateWRM_(Teuchos::ParameterList& plist)
{
  Errors::Message msg;

  // Find out how many WRM entries there are.
  int nblocks = 0;
  for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++) {
    if (plist.isSublist(plist.name(i))) nblocks++;
  }

  wrm_.resize(nblocks);

  int iblock = 0;
  for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++) {
    if (plist.isSublist(plist.name(i))) {
      Teuchos::ParameterList& wrm_list = plist.sublist(plist.name(i));

      std::string region;
      if (wrm_list.isParameter("region")) {
        region = wrm_list.get<std::string>("region");  // associated mesh block
      } else {
        msg << "Flow PK: WMR sublist \"" << plist.name(i).c_str() << "\" has no parameter \"region\".\n";
        Exceptions::amanzi_throw(msg);
      }

      if (wrm_list.get<std::string>("water retention model") == "van Genuchten") {
        double m = wrm_list.get<double>("van Genuchten m", FLOW_WRM_EXCEPTION);
        double alpha = wrm_list.get<double>("van Genuchten alpha", FLOW_WRM_EXCEPTION);
        double l = wrm_list.get<double>("van Genuchten l", FLOW_WRM_VANGENUCHTEN_L);
        double sr = wrm_list.get<double>("residual saturation", FLOW_WRM_EXCEPTION);
        double pc0 = wrm_list.get<double>("regularization interval", FLOW_WRM_REGULARIZATION_INTERVAL);
        std::string krel_function = wrm_list.get<std::string>("relative permeability model", "Mualem");

        wrm_[iblock] = Teuchos::rcp(new WRM_vanGenuchten(region, m, l, alpha, sr, krel_function, pc0));

      } else if (wrm_list.get<std::string>("water retention model") == "Brooks Corey") {
        double lambda = wrm_list.get<double>("Brooks Corey lambda", FLOW_WRM_EXCEPTION);
        double alpha = wrm_list.get<double>("Brooks Corey alpha", FLOW_WRM_EXCEPTION);
        double l = wrm_list.get<double>("Brooks Corey l", FLOW_WRM_BROOKS_COREY_L);
        double sr = wrm_list.get<double>("residual saturation", FLOW_WRM_EXCEPTION);
        double pc0 = wrm_list.get<double>("regularization interval", FLOW_WRM_REGULARIZATION_INTERVAL);
        std::string krel_function = wrm_list.get<std::string>("relative permeability model", "Mualem");

        wrm_[iblock] = Teuchos::rcp(new WRM_BrooksCorey(region, lambda, l, alpha, sr, krel_function, pc0));

      } else if (wrm_list.get<std::string>("water retention model") == "fake") {
        wrm_[iblock] = Teuchos::rcp(new WRM_fake(region));

      } else {
        msg << "Flow PK: unknown water retention model.";
        Exceptions::amanzi_throw(msg);
      }
      iblock++;
    }
  }
}


/* ******************************************************************
* Auxiliary map from cells to regions.  
****************************************************************** */
void WRMEvaluator::CreateCell2Region_()
{
  cell2region_ = Teuchos::rcp(new Epetra_IntVector(mesh_->cell_map(true)));

  Epetra_IntVector& map = *cell2region_;
  map.PutValue(-1);

  for (int mb = 0; mb < wrm_.size(); mb++) {
    std::string region = wrm_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) map[*i] = mb;
  }
  
  // pp_ = Teuchos::rcp(new ParallelCommunication(mesh_));
  // pp_->CopyMasterCell2GhostCell(map);

  // internal check
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    if (map[c] < 0) {
      Errors::Message msg;
      msg << "Flow PK: water retention models do not cover the whole domain.";
      Exceptions::amanzi_throw(msg);  
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi
