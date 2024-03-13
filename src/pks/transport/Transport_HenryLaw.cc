/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

// Amanzi
#include "StringExt.hh"

#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Re-partition components between liquid and gas phases.
******************************************************************* */
void
Transport_PK::PrepareAirWaterPartitioning_()
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
            .get<Teuchos::Array<double>>("air-water partitioning coefficient", empty)
            .toVector();
  } else {
    air_water_map_.clear();
  }
}


/* *******************************************************************
* Re-partition components between liquid and gas phases.
******************************************************************* */
void
Transport_PK::MakeAirWaterPartitioning_()
{
  auto& tcc_c = *tcc_tmp->ViewComponent("cell");
  const auto& sat_l =
    *S_->Get<CompositeVector>(saturation_liquid_key_, Tags::DEFAULT).ViewComponent("cell");

  int ig, il;
  double sl, sg, total;
  for (int i = 0; i < num_gaseous; ++i) {
    ig = num_aqueous + i;
    il = air_water_map_[i];

    for (int c = 0; c < ncells_owned; c++) {
      sl = sat_l[0][c];
      sg = 1.0 - sl;
      total = tcc_c[il][c] * sl + tcc_c[ig][c] * (1.0 - sl);
      if (sg > 0) {
        tcc_c[il][c] = total / (1.0 + sg / kH_[i]);
        tcc_c[ig][c] = (total - sl * tcc_c[il][c]) / sg;
      } else {
        tcc_c[il][c] = total;
        tcc_c[ig][c] = 0.0;
      }
    }
  }
}

} // namespace Transport
} // namespace Amanzi
