/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! FunctionStaticHead: f = p0 + rho*g*(z-z0)


/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!
:math:`f(z) = p0 + rho * g * (z0 - z)`

Note that dimension 0 is usually time.

* `"p0`" ``[double]`` Pressure at z0
* `"density`" ``[double]`` Density of water
* `"gravity`" ``[double]`` Gravity
* `"space dimension`" ``[int]`` Dimensionality, usually 3
* `"water table elevation`" ``[function-spec]`` Water table elevation function.


Example:
.. code-block:: xml

  <ParameterList name="function-static-head">
    <Parameter name="p0" type="double" value="101325.0"/>
    <Parameter name="density" type="double" value="1000.0"/>
    <Parameter name="gravity" type="double" value="9.8"/>
    <Parameter name="space dimension" type="int" value="3"/>
    <ParameterList name="water table elevation">
      <ParameterList name="function-constant">
        <Parameter name="value" type="double" value="1.0"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_STATIC_HEAD_FUNCTION_HH_
#define AMANZI_STATIC_HEAD_FUNCTION_HH_

namespace Amanzi {

#include <memory>

#include "Function.hh"

class FunctionStaticHead : public Function {
 public:
  FunctionStaticHead(double patm, double rho, double g, std::unique_ptr<Function> h, int dim)
      : patm_(patm), rho_g_(rho*g), h_(std::move(h)), dim_(dim) {}
  FunctionStaticHead(double patm, double rho, double g, const Function& h, int dim)
      : patm_(patm), rho_g_(rho*g), h_(h.Clone()), dim_(dim) {}
  FunctionStaticHead(const FunctionStaticHead& src)
      : patm_(src.patm_), rho_g_(src.rho_g_), h_(src.h_->Clone()), dim_(src.dim_) {}
  ~FunctionStaticHead() {}
  std::unique_ptr<Function> Clone() const {
    return std::make_unique<FunctionStaticHead>(*this);
  }
  // The array (t,x,y,z) is passed as *x, so that x[dim_] is z in 3D, y in 2D
  double operator()(const std::vector<double>& x) const { return patm_+rho_g_ * ((*h_)(x)-x[dim_]); }
 
 private:
  int dim_;
  double patm_, rho_g_;
  std::unique_ptr<Function> h_;
};

} // namespace Amanzi

#endif // AMANZI_STATIC_HEAD_FUNCTION_HH_
