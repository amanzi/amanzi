/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/
/*!

Anisotropic mechanical dispersion model by Burnett & Frind that splits the
transverse component of dispersion between the horizontal and vertical.

`"mechanical dispersion type`" = `"Burnett-Frind`"

.. _mdm-burnett-frind-spec
.. admonition:: mdm-burnett-frind-spec

   * `"alpha_l`" ``[double]`` Longitudinal dispersivity/dispersion coefficient.

   * `"alpha_tv`" ``[double]`` Vertical transverse dispersivity/dispersion
     coefficient.

   * `"alpha_th`" ``[double]`` Horizontal transverse dispersivity/dispersion
     coefficient.
*/

#ifndef AMANZI_MDM_BURNETT_FRIND_HH_
#define AMANZI_MDM_BURNETT_FRIND_HH_

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

class MDM_BurnettFrind : public MDM {
 public:
  explicit MDM_BurnettFrind(Teuchos::ParameterList& plist);
  ~MDM_BurnettFrind() {};

  // Required methods from the base class
  // -- dispersion tensor of rank 2
  WhetStone::Tensor mech_dispersion(double t,
                                    const AmanziGeometry::Point& xc,
                                    const AmanziGeometry::Point& u,
                                    int axi_symmetry,
                                    double wc,
                                    double phi) const;

  // -- the model is valid if at least one parameter is not zero.
  bool is_valid() const { return (alphaL_ + alphaTH_ + alphaTV_ != 0.0); }

  // -- check model applicability
  void set_dim(int dim)
  {
    AMANZI_ASSERT(dim == 3);
    dim_ = dim;
  }

 private:
  double alphaL_, alphaTH_, alphaTV_;

  static Utils::RegisteredFactory<MDM, MDM_BurnettFrind> reg_;
};

} // namespace Transport
} // namespace Amanzi

#endif
