/*
  MultiPhase

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Quan Bui (mquanbui@math.umd.edu)

  This file implements the main methods of the class for MPCoeff. 
  Unlike the one in Flow, this is specialized for multiphase. The relative 
  permeability is computed from water saturation, instead of pressure (as in Richards).
*/

#include "errors.hh"

#include "MultiphaseDefs.hh"
#include "MPCoeff.hh"
#include "WaterRetentionModel.hh"
#include "WRM_Simple.hh"
#include "WRM_VanGenuchten.hh"
#include "WRM_BrooksCorey.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Initialize internal data.
****************************************************************** */
void MPCoeff::Init(std::string phase_name, Teuchos::ParameterList& plist)
{
  //atm_pressure = p0;
  phase_ = phase_name;
  ProcessParameterList_(plist);

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);
  cvs.AddComponent("dirichlet_faces", AmanziMesh::BOUNDARY_FACE, 1);

  Krel_ = Teuchos::rcp(new CompositeVector(cvs));
  dKdS_ = Teuchos::rcp(new CompositeVector(cvs));
  rho_ = Teuchos::rcp(new CompositeVector(cvs));
  mpCoeff_ = Teuchos::rcp(new CompositeVector(cvs));
  rhoDerivKrel_ = Teuchos::rcp(new CompositeVector(cvs));
  dPc_ = Teuchos::rcp(new CompositeVector(cvs));

  //method_ = FLOW_RELATIVE_PERM_NONE;
  Krel_->PutScalarMasterAndGhosted(1.0);
  dKdS_->PutScalarMasterAndGhosted(0.0);
  rho_->PutScalarMasterAndGhosted(1.0);
  mpCoeff_->PutScalarMasterAndGhosted(1.0);
  rhoDerivKrel_->PutScalarMasterAndGhosted(1.0);
  dPc_->PutScalarMasterAndGhosted(1.0);

  Cg_ = 1.0;

  pp_ = Teuchos::rcp(new ParallelCommunication(mesh_));
  map_c2mb_ = Teuchos::rcp(new Epetra_IntVector(mesh_->cell_map(true)));
  PopulateMapC2MB_();
}


void MPCoeff::Init(std::string phase_name, Teuchos::ParameterList& plist, double Cg)
{
  Init(phase_name, plist);
  Cg_ = Cg;
}


/* ******************************************************************
* Defines relative permeability ONLY for cells.
****************************************************************** */
void MPCoeff::Compute(const CompositeVector& Sw)
{
  const Epetra_MultiVector& Sw_cell = *Sw.ViewComponent("cell");
  Epetra_MultiVector& Krel_cell = *Krel_->ViewComponent("cell");
  Epetra_MultiVector& dKdS_cell = *dKdS_->ViewComponent("cell");

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();

    std::vector<AmanziMesh::Entity_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double Sw = Sw_cell[0][*i];
      Krel_cell[0][*i] = WRM_[mb]->k_relative(Sw, phase_);
      dKdS_cell[0][*i] = WRM_[mb]->dKdS(Sw, phase_);
    }
  }

  // add boundary face component
  Krel_->ViewComponent("dirichlet_faces", true)->PutScalar(1.0);
  dKdS_->ViewComponent("dirichlet_faces", true)->PutScalar(0.0);
}


/* ******************************************************************
* Defines relative permeability ONLY for cells.
****************************************************************** */
void MPCoeff::Compute(const CompositeVector& primary_var, const CompositeVector& Sw)
{
  const Epetra_MultiVector& Sw_cell = *Sw.ViewComponent("cell");
  const Epetra_MultiVector& primary_var_cell = *primary_var.ViewComponent("cell");
  Epetra_MultiVector& Krel_cell = *Krel_->ViewComponent("cell");
  Epetra_MultiVector& dKdS_cell = *dKdS_->ViewComponent("cell");
  Epetra_MultiVector& rho_cell = *rho_->ViewComponent("cell");
  Epetra_MultiVector& mpCoeff_cell = *mpCoeff_->ViewComponent("cell");
  Epetra_MultiVector& rhoDerivKrel_cell = *rhoDerivKrel_->ViewComponent("cell");
  Epetra_MultiVector& dPc_dS_cell = *dPc_->ViewComponent("cell");

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();

    std::vector<AmanziMesh::Entity_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double Sw = Sw_cell[0][*i];
      double tmpRhoCell;
      if (phase_ == "non wetting")
        tmpRhoCell = Cg_ * (primary_var_cell[0][*i] + WRM_[mb]->capillaryPressure(Sw));
      else
        tmpRhoCell = primary_var_cell[0][*i];
      double krel = WRM_[mb]->k_relative(Sw, phase_);
      double dKdS = WRM_[mb]->dKdS(Sw, phase_);
      Krel_cell[0][*i] = krel;
      dKdS_cell[0][*i] = dKdS;
      rho_cell[0][*i] = tmpRhoCell;
      mpCoeff_cell[0][*i] = tmpRhoCell * krel;
      rhoDerivKrel_cell[0][*i] = tmpRhoCell * dKdS;
      dPc_dS_cell[0][*i] = WRM_[mb]->dPc_dS(Sw);
    }
  }

  // add boundary face component
  Krel_->ViewComponent("dirichlet_faces", true)->PutScalar(1.0);
  dKdS_->ViewComponent("dirichlet_faces", true)->PutScalar(0.0);
  mpCoeff_->ViewComponent("cell", true)->PutScalar(1.0);
  rhoDerivKrel_->ViewComponent("cell", true)->PutScalar(0.0);
  dPc_->ViewComponent("cell", true)->PutScalar(0.0);
}


/* ******************************************************************
* Auxiliary map from cells to WRM models.                                               
****************************************************************** */
void MPCoeff::PopulateMapC2MB_()
{
  Epetra_IntVector& map = *map_c2mb_;
  map.PutValue(-1);

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      map[*i] = mb;
    }
  }
  
  pp_->CopyMasterCell2GhostCell(map);

  // internal check
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    if (map[c] < 0) {
      Errors::Message msg;
      msg << "Multiphase PK: water retention models do not cover the whole domain.\n";
      Exceptions::amanzi_throw(msg);  
    }
  }
}

void MPCoeff::ProcessParameterList_(Teuchos::ParameterList& plist)
{
  Errors::Message msg;

  // create verbosity object
  Teuchos::ParameterList vlist;
  vo_ = new VerboseObject("Multiphase::Coeff", vlist); 

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

}  // namespace Flow
}  // namespace Amanzi

