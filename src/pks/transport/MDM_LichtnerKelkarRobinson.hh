/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/
/*!

Anisotropic mechanical dispersion model by Lichtner, Kelkar, and Robinson, that
splits vertical from horizontal and longitudinal from transverse.

`"mechanical dispersion type`" = `"Lichtner-Kelkar-Robinson`"

.. _mdm-lichtner-kelkar-robinson-spec
.. admonition:: mdm-lichtner-kelkar-robinson-spec

   * `"alpha_lv`" ``[double]`` **0.** Vertical longitudinal
     dispersivity/dispersion coefficient.
   * `"alpha_lh`" ``[double]`` **0.** Horizontal longitudinal
     dispersivity/dispersion coefficient.
   * `"alpha_tv`" ``[double]`` **0.** Vertical transverse
     dispersivity/dispersion coefficient.
   * `"alpha_th`" ``[double]`` **0.** Horizontal transverse
     dispersivity/dispersion coefficient.

*/

#ifndef AMANZI_MDM_LICHTNER_KELKAR_ROBINSON_HH_
#define AMANZI_MDM_LICHTNER_KELKAR_ROBINSON_HH_

// TPLs
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Factory.hh"
#include "Point.hh"
#include "Tensor.hh"

// Transport
#include "MDM.hh"

namespace Amanzi {
namespace Transport {

class MDM_LichtnerKelkarRobinson : public MDM {
 public:
  explicit MDM_LichtnerKelkarRobinson(Teuchos::ParameterList& plist);
  ~MDM_LichtnerKelkarRobinson(){};

  // Required methods from the base class
  // -- dispersion tensor of rank 2
  WhetStone::Tensor mech_dispersion(double t,
                                    const AmanziGeometry::Point& xc,
                                    const AmanziGeometry::Point& u,
                                    int axi_symmetry,
                                    double wc,
                                    double phi) const;

  // -- the model is valid if at least one parameter is not zero.
  bool is_valid() const { return (alphaLH_ + alphaLV_ + alphaTH_ + alphaTV_ != 0.0); }

  // -- simple check of model applicability
  void set_dim(int dim)
  {
    AMANZI_ASSERT(dim == 3);
    dim_ = dim;
  }

 private:
  double alphaLH_, alphaLV_, alphaTH_, alphaTV_;

  static Utils::RegisteredFactory<MDM, MDM_LichtnerKelkarRobinson> reg_;
};

} // namespace Transport
} // namespace Amanzi

#endif
