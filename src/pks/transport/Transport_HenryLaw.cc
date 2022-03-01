/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

// Amanzi
#include "StringExt.hh"

#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Re-partition components between liquid and gas phases.
******************************************************************* */
void Transport_PK::PrepareAirWaterPartitioning_()
{
  henry_law_ = true;
  for (int i = 0; i < num_gaseous; i++) {
    int ig = num_aqueous + i;
    std::string name_l(component_names_[ig]);
    Amanzi::replace_all(name_l, "(g)", "(l)");

    int il = FindComponentNumber(name_l);
    air_water_map_.push_back(il);

    if (il < 0 || il >= num_aqueous) {
      henry_law_ = false;
      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "Gas component \"" << component_names_[ig] 
                   << "\" has no matching liquid component \"" << name_l << "\"\n";
      }
      break;
    }
  }

  if (henry_law_) {
    Teuchos::Array<double> empty;
    kH_ = tp_list_->sublist("molecular diffusion")
        .get<Teuchos::Array<double> >("air-water partitioning coefficient", empty).toVector();
  } else {
    air_water_map_.clear();
  }
}


/* *******************************************************************
* Re-partition components between liquid and gas phases.
******************************************************************* */
void Transport_PK::MakeAirWaterPartitioning_()
{
  auto& tcc_c = *tcc_tmp->ViewComponent("cell");
  const auto& sat_l = *ws;

  int ig, il;
  double sl, total;
  for (int i = 0; i < num_gaseous; ++i) {
    ig = num_aqueous + i;
    il = air_water_map_[i];

    for (int c = 0; c < ncells_owned; c++) {
      sl = sat_l[0][c];
      total = tcc_c[il][c] * sl + tcc_c[ig][c] * (1.0 - sl);
      tcc_c[ig][c] = total / (1.0 + (kH_[i] - 1.0) * sl);
      tcc_c[il][c] = tcc_c[ig][c] * kH_[i];
    } 
  }
}

}  // namespace Transport
}  // namespace Amanzi


