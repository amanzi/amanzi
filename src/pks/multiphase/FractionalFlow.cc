/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Quan Bui (mquanbui@math.umd.edu)
  This file implements the FractionalFLow class, which is unique to
  IMPES formulation. It is used in the saturation equation.
*/

#include "errors.hh"

#include "MultiphaseDefs.hh"
#include "FractionalFlow.hh"
#include "WaterRetentionModel.hh"
#include "WRM_Simple.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Initialize internal data.
****************************************************************** */
void FractionalFlow::Init(Teuchos::ParameterList& plist)
{
  ProcessParameterList_(plist);
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);

  total_mobility_ = Teuchos::rcp(new CompositeVector(cvs));
  fractional_flow_w_ = Teuchos::rcp(new CompositeVector(cvs));
  dfw_dS_ = Teuchos::rcp(new CompositeVector(cvs));

  method_ = FLOW_RELATIVE_PERM_NONE;
  total_mobility_->PutScalarMasterAndGhosted(1.0);
  fractional_flow_w_->PutScalarMasterAndGhosted(0.0);

  pp_ = Teuchos::rcp(new ParallelCommunication(mesh_));
  map_c2mb_ = Teuchos::rcp(new Epetra_IntVector(mesh_->cell_map(true)));
  PopulateMapC2MB_();
}


/* ******************************************************************
* Defines fractional flow ONLY for cells.
****************************************************************** */
void FractionalFlow::Compute(const CompositeVector& Sw)
{
  const Epetra_MultiVector& Sw_cell = *Sw.ViewComponent("cell");
  Epetra_MultiVector& frac_flow_cell = *fractional_flow_w_->ViewComponent("cell");
  Epetra_MultiVector& total_mobility_cell = *total_mobility_->ViewComponent("cell");
  Epetra_MultiVector& dfw_dS_c = *dfw_dS_->ViewComponent("cell");

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();

    std::vector<AmanziMesh::Entity_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double lambda_w, lambda_nw, lambda_tot, dlamdaw_ds, dlamdanw_ds, df_dS;
      double Sw = Sw_cell[0][*i];
      lambda_w = WRM_[mb]->k_relative(Sw, "wetting")/mu_w_;
      lambda_nw = WRM_[mb]->k_relative(Sw, "non wetting")/mu_nw_;
      lambda_tot = lambda_w + lambda_nw;
      total_mobility_cell[0][*i] = lambda_tot;
      frac_flow_cell[0][*i] = lambda_w/lambda_tot;

      dlamdaw_ds = WRM_[mb]->dKdS(Sw, "wetting")/mu_w_;
      dlamdanw_ds = WRM_[mb]->dKdS(Sw, "non wetting")/mu_nw_;
      df_dS = (dlamdaw_ds * lambda_tot - lambda_w * (dlamdanw_ds + dlamdaw_ds)) / pow(lambda_tot, 2.0);
      dfw_dS_c[0][*i] = df_dS;
    }
  }

  //std::cout<<"Sw\n"<<Sw_cell<<"\n";
  //std::cout<<"frac flow "<<frac_flow_cell<<"\n";
}


/* ******************************************************************
* Auxiliary map from cells to WRM models.                                               
****************************************************************** */
void FractionalFlow::PopulateMapC2MB_()
{
  Epetra_IntVector& map = *map_c2mb_;
  map.PutValue(-1);

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) map[*i] = mb;
  }
  
  pp_->CopyMasterCell2GhostCell(map);

  // internal check
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    if (map[c] < 0) {
      Errors::Message msg;
      msg << "Multiphase PK: water retention models do not cover the whole domain.";
      Exceptions::amanzi_throw(msg);  
    }
  }
}

/* ******************************************************************
* Routine to process the necessary data to compute fractional flow
* from XML input file                   
****************************************************************** */
void FractionalFlow::ProcessParameterList_(Teuchos::ParameterList& plist)
{
  Errors::Message msg;

  // create verbosity object
  Teuchos::ParameterList vlist;
  vo_ = new VerboseObject("Flow::RelativePerm", vlist); 

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

      } else {
        msg << "Multiphase PK: unknown water retention model.";
        Exceptions::amanzi_throw(msg);
      }
      iblock++;
    }
  }

  // optional debug output
  //PlotWRMcurves(plist);
}

}  // namespace Flow
}  // namespace Amanzi

