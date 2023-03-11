/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

A two-scale porosity model (fracture + matrix) aka dual porosity
model. Current naming convention is that the fields used in the
single-porosity model correspond now to the fracture continuum.
Example: pressure = pressure in the fracture continuum;
pressure_msp = pressure in the matrix continuum.

* `"number of matrix nodes`" [int] defines number of matrix layers.
* `"matrix depth`" [double] is the characteristic length for matrix continuum.
* `"matrix volume fraction`" [double] defines relative volume of matrix continuum.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="multiscale models"> 
    <ParameterList name="_SOIL1">
      <Parameter name="regions" type="Array(string)" value="{_TOP HALF}"/>
      <Parameter name="multiscale model" type="string" value="dual porosity"/> 
      <ParameterList name="dual porosity parameters">
        <Paramater name="mass transfer coefficient" type="double" value="4.0e-5"/>
        <Paramater name="tolerance" type="double" value="1e-8"/>
      </ParameterList>

      <Parameter name="water retention model" type="string" value="van Genuchten"/>
      <Parameter name="van Genuchten alpha" type="double" value="0.000194"/>
      <Parameter name="van Genuchten m" type="double" value="0.28571"/>
      <Parameter name="van Genuchten l" type="double" value="0.5"/>
      <Parameter name="residual saturation" type="double" value="0.103"/>
      <Parameter name="relative permeability model" type="string" value="Mualem"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/

#ifndef MULTISCALE_FLOW_POROSITY_GDPM_HH_
#define MULTISCALE_FLOW_POROSITY_GDPM_HH_

#include <memory>

// TPLs
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "DenseVector.hh"
#include "Factory.hh"
#include "Mini_Diffusion1D.hh"
#include "SolverNewton.hh"

// Flow
#include "MultiscaleFlowPorosity.hh"
#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class MultiscaleFlowPorosity_GDPM : public MultiscaleFlowPorosity,
                                    public AmanziSolvers::SolverFnBase<WhetStone::DenseVector> {
 public:
  MultiscaleFlowPorosity_GDPM(Teuchos::ParameterList& plist);
  ~MultiscaleFlowPorosity_GDPM(){};

  // interface for porosity models
  // -- calculate field water storage assuming pressure equilibrium
  virtual double ComputeField(double phi, double n_l, double pcm) override;

  // -- local (cell-based) solver returns water storage and capilalry
  //   pressure in the matrix. max_itrs is input/output parameter
  virtual double WaterContentMatrix(double pcf0,
                                    WhetStone::DenseVector& pcm,
                                    double wcm0,
                                    double dt,
                                    double phi,
                                    double n_l,
                                    int& max_itrs) override;

  // -- number of matrix nodes
  virtual int NumberMatrixNodes() override { return matrix_nodes_; }

  // interface for nonlinear solvers
  // -- calculate nonlinear residual
  virtual void Residual(const Teuchos::RCP<WhetStone::DenseVector>& u,
                        const Teuchos::RCP<WhetStone::DenseVector>& f) override;

  // -- delegeating inversion to the diffusion operator
  virtual int ApplyPreconditioner(const Teuchos::RCP<const WhetStone::DenseVector>& u,
                                  const Teuchos::RCP<WhetStone::DenseVector>& hu) override
  {
    op_->ApplyInverse(*u, *hu);
    return 0;
  }

  // -- error estimate for convergence criteria
  virtual double ErrorNorm(const Teuchos::RCP<const WhetStone::DenseVector>& u,
                           const Teuchos::RCP<const WhetStone::DenseVector>& du) override;

  // -- calculate preconditioner
  virtual void UpdatePreconditioner(const Teuchos::RCP<const WhetStone::DenseVector>& u) override;

  // -- other required functions
  virtual void ChangedSolution() override{};

  // modifiers
  void set_op(std::shared_ptr<Operators::Mini_Diffusion1D> op) { op_ = op; }
  void set_bcl(double bcl) { bcl_ = bcl; }

 private:
  Teuchos::RCP<WRM> wrm_;

  int matrix_nodes_;
  double depth_, tau_, tol_, Ka_;

  double bcl_;
  std::shared_ptr<Operators::Mini_Diffusion1D> op_;

  static Utils::RegisteredFactory<MultiscaleFlowPorosity, MultiscaleFlowPorosity_GDPM> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
