/*
This is the multiphase component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
         Quan Bui (mquanbui@math.umd.edu)
This file implements the main methods of the class for CapillaryPressure. 
Unlike the one in Flow, this is specialized for multiphase. The capillary pressure 
is computed from water saturation, instead of pressure (as in Richards).
*/

#include "errors.hh"

#include "MultiphaseDefs.hh"
#include "CapillaryPressure.hh"
#include "WaterRetentionModel.hh"
#include "WRM_Simple.hh"
#include "WRM_VanGenuchten.hh"
#include "WRM_BrooksCorey.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Initialize internal data.
****************************************************************** */
void CapillaryPressure::Init(Teuchos::ParameterList& plist)
{
  //atm_pressure = p0;
  ProcessParameterList_(plist);

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);

  Pc_ = Teuchos::rcp(new CompositeVector(cvs));
  dPc_dS_ = Teuchos::rcp(new CompositeVector(cvs));

  Pc_->PutScalarMasterAndGhosted(0.0);
  dPc_dS_->PutScalarMasterAndGhosted(0.0);

  pp_ = Teuchos::rcp(new ParallelCommunication(mesh_));
  map_c2mb_ = Teuchos::rcp(new Epetra_IntVector(mesh_->cell_map(true)));
  PopulateMapC2MB_();
}


/* ******************************************************************
* Defines relative permeability ONLY for cells.
****************************************************************** */
void CapillaryPressure::Compute(const CompositeVector& saturation_water)
{
  const Epetra_MultiVector& Sw_cell = *saturation_water.ViewComponent("cell");
  Epetra_MultiVector& Pc_cell = *Pc_->ViewComponent("cell");
  Epetra_MultiVector& dPc_dS_cell = *dPc_dS_->ViewComponent("cell");

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();

    std::vector<AmanziMesh::Entity_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double Sw = Sw_cell[0][*i];
      Pc_cell[0][*i] = WRM_[mb]->capillaryPressure(Sw);
      dPc_dS_cell[0][*i] = WRM_[mb]->dPc_dS(Sw);
    }
  }
}


/* ******************************************************************
* List to WRM models has to be provided.
****************************************************************** */
void CapillaryPressure::ProcessParameterList_(Teuchos::ParameterList& plist)
{
  Errors::Message msg;

  // create verbosity object
  Teuchos::ParameterList vlist;
  vo_ = new VerboseObject("Multiphase::CapillaryPressure", vlist); 

  int nblocks = 0;  // Find out how many WRM entries there are.
  for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++) {
    if (plist.isSublist(plist.name(i))) nblocks++;
  }

  WRM_.resize(nblocks);

  int iblock = 0;
  for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++) {
    if (plist.isSublist(plist.name(i))) {
      Teuchos::ParameterList& wrm_list = plist.sublist(plist.name(i));

      std::string region;
      if (wrm_list.isParameter("region")) {
        region = wrm_list.get<std::string>("region");  // associated mesh block
      } else {
        msg << "Multiphase PK: WMR sublist \"" << plist.name(i).c_str() << "\" has no parameter \"region\".\n";
        Exceptions::amanzi_throw(msg);
      }

      if (wrm_list.get<std::string>("water retention model") == "Simple") {
        double S_rw = wrm_list.get<double>("residual saturation wet", FLOW_WRM_EXCEPTION);
        double S_rn = wrm_list.get<double>("residual saturation nonwet", FLOW_WRM_EXCEPTION);
        double coef = wrm_list.get<double>("coefficient", FLOW_WRM_EXCEPTION);

        WRM_[iblock] = Teuchos::rcp(new WRM_Simple(region, S_rw, S_rn, coef));
      } else if (wrm_list.get<std::string>("water retention model") == "Van Genuchten") {
        double S_rw = wrm_list.get<double>("residual saturation wet", FLOW_WRM_EXCEPTION);
        double S_rn = wrm_list.get<double>("residual saturation nonwet", FLOW_WRM_EXCEPTION);
        double n = wrm_list.get<double>("van genuchten n", FLOW_WRM_EXCEPTION);
        double Pr = wrm_list.get<double>("van genuchten entry pressure", FLOW_WRM_EXCEPTION);
        //double alpha = wrm_list.get<double>("van genuchten alpha", FLOW_WRM_EXCEPTION);
        //double epsilon = wrm_list.get<double>("van genuchten epsilon", FLOW_WRM_EXCEPTION);
        //double gamma = wrm_list.get<double>("van genuchten gamma", FLOW_WRM_EXCEPTION);

        WRM_[iblock] = Teuchos::rcp(new WRM_VanGenuchten(region, S_rw, S_rn, n, Pr));
      } else if (wrm_list.get<std::string>("water retention model") == "Brooks Corey") {
        double S_rw = wrm_list.get<double>("residual saturation wet", FLOW_WRM_EXCEPTION);
        double S_rn = wrm_list.get<double>("residual saturation nonwet", FLOW_WRM_EXCEPTION);
        double pd = wrm_list.get<double>("entry pressure", FLOW_WRM_EXCEPTION);
        double lambda = wrm_list.get<double>("brooks corey lambda", FLOW_WRM_EXCEPTION);

        WRM_[iblock] = Teuchos::rcp(new WRM_BrooksCorey(region, S_rw, S_rn, pd, lambda));
      } else {
        msg << "Multiphase PK: unknown water retention model.";
        Exceptions::amanzi_throw(msg);
      }
      iblock++;
    }
  }
}


/* ******************************************************************
* Auxiliary map from cells to WRM models.                                               
****************************************************************** */
void CapillaryPressure::PopulateMapC2MB_()
{
  Epetra_IntVector& map = *map_c2mb_;
  map.PutValue(-1);

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      map[*i] = mb;
    }
  }
  
  pp_->CopyMasterCell2GhostCell(map);

  // internal check
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    if (map[c] < 0) {
      Errors::Message msg;
      msg << "Multiphase PK: water retention models do not cover the whole domain.\n";
      Exceptions::amanzi_throw(msg);  
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi

