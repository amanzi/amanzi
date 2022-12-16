/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

  Base class for mechanical dispersion models.
*/

#ifndef AMANZI_MECHANICAL_DISPERSION_MODEL_HH_
#define AMANZI_MECHANICAL_DISPERSION_MODEL_HH_

// Amanzi
#include "Point.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace Transport {

class MDM {
 public:
  virtual ~MDM(){};

  // returns dispersion tensor.
  virtual WhetStone::Tensor mech_dispersion(const AmanziGeometry::Point& u,
                                            int axi_symmetry,
                                            double wc,
                                            double phi) const = 0;

  // The model is valid if at least one parameter is not zero.
  virtual bool is_valid() const = 0;

  // This allows us to set space dimension which could be used for estimating
  // model applicability.
  virtual void set_dim(int dim) { dim_ = dim; }

 protected:
  int dim_;
};

} // namespace Transport
} // namespace Amanzi

#endif
