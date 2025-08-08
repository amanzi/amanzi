/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Mimetic Finite Difference class of methods for diffusion on a set of
possibly non-matching meshes

Example:

.. code-block:: xml

  <ParameterList name="diffusion operator">
    <Parameter name="discretization primary" type="string" value="mfd: optimized for sparsity"/>
    <Parameter name="discretization secondary" type="string" value="mfd: optimized for sparsity"/>
    <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
    <Parameter name="preconditioner schema" type="Array(string)" value="{face, cell}"/>
    <Parameter name="gravity" type="bool" value="false"/>
    <Parameter name="meshes" type="Array(string)" value="{mesh1, mesh2, mesh3}"/>
    <ParameterList name="interfaces">
      <ParameterList name="interface 1">
        <Parameter name="meshes" type="Array(string)" value="{mesh1, mesh2}"/>
        <Parameter name="common surfacei regions" type="Array(string)" value="{surface12, surface21}"/>
      </ParameterList>
      <ParameterList name="interface 2">
        <Parameter name="meshes" type="Array(string)" value="{mesh2, mesh3}"/>
        <Parameter name="common surface regions" type="Array(string)" value="{surface23, surface32}"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_MULTI_MESH_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_MULTI_MESH_HH_

#include <map>
#include <utility>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "DenseMatrix.hh"
#include "exceptions.hh"
#include "Point.hh"
#include "State.hh"
#include "VerboseObject.hh"

#include "BCs.hh"
#include "OperatorDefs.hh"
#include "Operator_MultiMesh.hh"
#include "PDE_Diffusion.hh"
#include "TreeOperator.hh"
#include "TreeVector.hh"

namespace Amanzi {
namespace Operators {

class PDE_DiffusionMultiMesh {
 public:
  const int REFINE = 9;

 public:
  PDE_DiffusionMultiMesh(Teuchos::ParameterList& plist);

  void Init(const Teuchos::RCP<State>& S);

  // connect mesh domain with tensor coefficients
  void SetVariableTensorCoefficient(const std::string& name,
                                    const Teuchos::RCP<std::vector<WhetStone::Tensor>>& K)
  {
    K_[name] = K;
  }

  // connect mesh domain with boundary conditions
  void SetBCs(const std::string& name,
              const Teuchos::RCP<const BCs>& bc_trial,
              const Teuchos::RCP<const BCs>& bc_test)
  {
    int i = std::distance(names_.begin(), std::find(names_.begin(), names_.end(), name));
    pdes_[i]->SetBCs(bc_trial, bc_test);
  }

  // update linearized problem
  void UpdateMatrices(const Teuchos::RCP<const TreeVector>& flux,
                      const Teuchos::RCP<const TreeVector>& u);

  // modify matrix due to boundary conditions
  //    primary=true indicates that the operator updates both matrix and right-hand
  //      side using BC data. If primary=false, only matrix is changed.
  //    eliminate=true indicates that we eliminate essential BCs for a trial
  //      function, i.e. zeros go in the corresponding matrix columns and
  //      right-hand side is modified using BC values. This is the optional
  //      parameter that enforces symmetry for a symmetric tree operators.
  //    essential_eqn=true indicates that the operator places a positive number on
  //      the main matrix diagonal for the case of essential BCs. This is the
  //      implementation trick.
  void ApplyBCs(bool primary, bool eliminate, bool essential_eqn);

  // calculate fluxes after solving the problem
  void UpdateFlux(const Teuchos::Ptr<const TreeVector>& u, const Teuchos::Ptr<TreeVector>& flux);

  // access
  Teuchos::RCP<Operators::TreeOperator> get_matrix() { return matrix_; }

 protected:
  // utilities to define surface-surface intersection
  void meshToMeshMapParticles_(const AmanziMesh::Mesh& mesh1,
                               const std::string& rgn1,
                               const AmanziMesh::Mesh& mesh2,
                               const std::string& rgn2,
                               InterfaceData& data12);

  // Note: if no initial guess, set f2_guess to -1
  int findFace_(const AmanziGeometry::Point& xf1,
                const AmanziGeometry::Point& ray,
                const AmanziMesh::Mesh& mesh2,
                const std::string& rgn2,
                int f2_guess,
                int* stage);

  void ModifyMatrices_(int ib, const InterfaceData& data);

 protected:
  Teuchos::ParameterList plist_;
  std::map<std::string, Teuchos::RCP<const AmanziMesh::Mesh>> meshes_;
  int d_;

  double stability_;
  std::vector<std::vector<std::string>> interface_meshes_;
  std::vector<std::vector<int>> interface_block_;
  std::vector<InterfaceData> interface_data_;

  std::vector<std::string> names_;
  std::map<std::string, Teuchos::RCP<std::vector<WhetStone::Tensor>>> K_;
  std::vector<std::vector<WhetStone::DenseMatrix>> matrices_grad_;

  std::vector<Teuchos::RCP<PDE_Diffusion>> pdes_;
  std::vector<Teuchos::RCP<Operator_MultiMesh>> ops_;
  Teuchos::RCP<Operators::TreeOperator> matrix_;

  Teuchos::RCP<VerboseObject> vo_;
};

} // namespace Operators
} // namespace Amanzi

#endif
