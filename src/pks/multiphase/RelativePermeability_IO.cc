/*
This is the multiphase component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
         Quan Bui (mquanbui@math.umd.edu)
This file implements the I/O components of the class for RelativePermeability. 
Unlike the one in Flow, this is specialized for multiphase. The relative 
permeability is computed from water saturation, instead of pressure (as in Richards).
*/

#include "errors.hh"

#include "WaterRetentionModel.hh"
#include "WRM_Simple.hh"
#include "WRM_VanGenuchten.hh"
#include "WRM_BrooksCorey.hh"

#include "RelativePermeability.hh"
#include "MultiphaseDefs.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* List to WRM models has to be provided.
****************************************************************** */
void RelativePermeability::ProcessParameterList_(Teuchos::ParameterList& plist)
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

  // optional debug output
  PlotWRMcurves(plist);
}

/* ****************************************************************
* Process string for the relative permeability
**************************************************************** */
/*
void RelativePermeability::ProcessStringRelativePermeability(const std::string name)
{
  Errors::Message msg;
  if (name == "upwind with gravity") {
    method_ = Multiphase::FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  } else if (name == "cell centered") {
    method_ = Multiphase::FLOW_RELATIVE_PERM_CENTERED;
  } else if (name == "upwind with Darcy flux") {
    method_ = Multiphase::FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  } else if (name == "arithmetic mean") {
    method_ = Multiphase::FLOW_RELATIVE_PERM_ARITHMETIC_MEAN;
  } else if (name == "upwind amanzi") {
    method_ = Multiphase::FLOW_RELATIVE_PERM_AMANZI;
  } else {
    msg << "Multiphase PK: unknown relative permeability method has been specified.";
    Exceptions::amanzi_throw(msg);
  }
}
*/

/* ****************************************************************
* Plot water retention curves.
**************************************************************** */
/*
void RelativePermeability::PlotWRMcurves(Teuchos::ParameterList& plist)
{
  int MyPID = mesh_->cell_map(false).Comm().MyPID();
  if (MyPID == 0) {
    int mb(0); 
    for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++) {
      if (plist.isSublist(plist.name(i))) {
        Teuchos::ParameterList& wrm_list = plist.sublist(plist.name(i));

        if (wrm_list.isSublist("output")) {
          Teuchos::ParameterList& out_list = wrm_list.sublist("output");

          std::string fname = out_list.get<std::string>("file");
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "saving sat-krel-pc date in file \"" << fname << "\"..." << std::endl;
          std::ofstream ofile;
          ofile.open(fname.c_str());

          int ndata = out_list.get<int>("number of points", 100);
          ndata = std::max(ndata, 1);

          double sr = WRM_[mb]->residualSaturation();
          double ds = (1.0 - sr) / ndata;

          for (int i = 0; i < ndata; i++) {
            double sat = sr + ds * (i + 0.5);
            double pc = WRM_[mb]->capillaryPressure(sat);
            double krel = WRM_[mb]->k_relative(pc);
            ofile << sat << " " << krel << " " << pc << std::endl;
          }
          ofile << std::endl;
          ofile.close();
        }
        mb++; 
      }
    }
  }
}
*/

}  // namespace Flow
}  // namespace Amanzi

