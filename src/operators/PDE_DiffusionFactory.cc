/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Base factory for diffusion operators.
*/

#include "errors.hh"

#include "BCs.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFVwithGravity.hh"
#include "PDE_DiffusionFVonManifolds.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_DiffusionMFDwithGravity.hh"
#include "PDE_DiffusionNLFV.hh"
#include "PDE_DiffusionNLFVwithBndFaces.hh"
#include "PDE_DiffusionNLFVwithBndFacesGravity.hh"
#include "PDE_DiffusionNLFVwithGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Constructors
****************************************************************** */
PDE_DiffusionFactory::PDE_DiffusionFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : mesh_(mesh), const_k_(1.0), gravity_(false), const_b_(0.0){};


PDE_DiffusionFactory::PDE_DiffusionFactory(Teuchos::ParameterList& oplist,
                                           const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : oplist_(oplist), mesh_(mesh), const_k_(1.0), gravity_(false), const_b_(0.0)
{
  if (oplist.isParameter("diffusion coefficient")) {
    const_k_ = oplist.get<double>("diffusion coefficient");
  }
}


/* ******************************************************************
* Setup the problem
****************************************************************** */
void
PDE_DiffusionFactory::SetConstantTensorCoefficient(const WhetStone::Tensor& K)
{
  const_K_ = K;
  K_ = Teuchos::null;
}


void
PDE_DiffusionFactory::SetConstantScalarCoefficient(double k)
{
  const_k_ = k;
  k_ = Teuchos::null;
  dkdu_ = Teuchos::null;
}


void
PDE_DiffusionFactory::SetConstantGravitationalTerm(const AmanziGeometry::Point& g, double b)
{
  gravity_ = true;
  g_ = g;
  const_b_ = b;
  b_ = Teuchos::null;
  dbdu_ = Teuchos::null;
}


/* ******************************************************************
* Setup the problem
****************************************************************** */
Teuchos::RCP<PDE_Diffusion>
PDE_DiffusionFactory::Create()
{
  std::string name = oplist_.get<std::string>("discretization primary");
  bool fractured_matrix = oplist_.isParameter("fracture");

  if (oplist_.isSublist("gravity")) gravity_ = oplist_.get<bool>("gravity");

  if (gravity_ && norm(g_) == 0.0) {
    double tmp = oplist_.get<double>("gravity magnitude");
    g_[mesh_->space_dimension() - 1] = tmp;
  }

  Teuchos::RCP<PDE_Diffusion> op;

  // FV methods
  if (name == "fv: default" && !gravity_) {
    op = Teuchos::rcp(new PDE_DiffusionFV(oplist_, mesh_));
  } else if (name == "fv: default" && gravity_) {
    op = Teuchos::rcp(new PDE_DiffusionFVwithGravity(oplist_, mesh_, g_));

    // NLFV methods
  } else if (name == "nlfv: default" && !gravity_) {
    op = Teuchos::rcp(new PDE_DiffusionNLFV(oplist_, mesh_));
  } else if (name == "nlfv: default" && gravity_) {
    op = Teuchos::rcp(new PDE_DiffusionNLFVwithGravity(oplist_, mesh_));
  } else if (name == "nlfv: bnd_faces" && !gravity_) {
    op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFaces(oplist_, mesh_));
  } else if (name == "nlfv: bnd_faces" && gravity_) {
    op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFacesGravity(oplist_, mesh_));

    // MFD methods with non-uniform DOFs
  } else if (fractured_matrix) {
    auto op_tmp = Teuchos::rcp(new PDE_DiffusionFracturedMatrix(oplist_, mesh_, const_b_, g_));
    op_tmp->Init(oplist_);
    op = op_tmp;

    // MFD methods
  } else if (!gravity_) {
    auto op_tmp = Teuchos::rcp(new PDE_DiffusionMFD(oplist_, mesh_));
    op_tmp->Init(oplist_);
    op = op_tmp;

  } else {
    auto op_tmp = Teuchos::rcp(new PDE_DiffusionMFDwithGravity(oplist_, mesh_));
    op_tmp->Init(oplist_);
    op = op_tmp;
  }

  // setup problem coefficients
  if (K_ != Teuchos::null) {
    op->SetTensorCoefficient(K_);
  } else if (const_K_.rank() > 0) {
    op->SetConstantTensorCoefficient(const_K_);
  }

  if (k_ != Teuchos::null)
    op->SetScalarCoefficient(k_, dkdu_);
  else
    op->SetConstantScalarCoefficient(const_k_);

  if (gravity_) {
    auto op_tmp = Teuchos::rcp_dynamic_cast<PDE_DiffusionWithGravity>(op);
    op_tmp->SetGravity(g_);

    if (b_ != Teuchos::null)
      op_tmp->SetDensity(b_);
    else
      op_tmp->SetDensity(const_b_);
  }

  return op;
}


/* ******************************************************************
* Initialization of diffusion operator with optional gravity.
* This is the constructor used by Amanzi.
****************************************************************** */
Teuchos::RCP<PDE_Diffusion>
PDE_DiffusionFactory::Create(Teuchos::ParameterList& oplist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             const Teuchos::RCP<BCs>& bc,
                             double rho,
                             const AmanziGeometry::Point& g)
{
  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);
  bool fractured_matrix = oplist.isParameter("fracture");

  // FV methods
  if (name == "fv: default" && !flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "fv: default" && flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionFVwithGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);
    return op;

    // NLFV methods
  } else if (name == "nlfv: default" && !flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFV(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: default" && flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: bnd_faces" && !flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFaces(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: bnd_faces" && flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFacesGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);
    return op;

    // MFD methods with non-uniform DOFs
  } else if (fractured_matrix) {
    auto op = Teuchos::rcp(new PDE_DiffusionFracturedMatrix(oplist, mesh, rho, g));
    op->Init(oplist);
    op->SetBCs(bc, bc);
    return op;

    // MFD methods
  } else if (!flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionMFD(oplist, mesh));
    op->Init(oplist);
    op->SetBCs(bc, bc);
    return op;

  } else {
    auto op = Teuchos::rcp(new PDE_DiffusionMFDwithGravity(oplist, mesh, rho, g));
    op->Init(oplist);
    op->SetBCs(bc, bc);
    return op;
  }
}


/* ******************************************************************
* Initialization of diffusion operator with optional gravity.
* This is the factory used by Amanzi, though it makes life difficult
* for time-varying density.
****************************************************************** */
Teuchos::RCP<PDE_Diffusion>
PDE_DiffusionFactory::Create(Teuchos::ParameterList& oplist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             const Teuchos::RCP<BCs>& bc,
                             const Teuchos::RCP<const CompositeVector>& rho,
                             const AmanziGeometry::Point& g)
{
  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);

  // FV methods
  if (name == "fv: default" && !flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "fv: default" && flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionFVwithGravity(oplist, mesh, g));
    op->SetBCs(bc, bc);
    op->SetDensity(rho);
    return op;

    // NLFV methods
  } else if (name == "nlfv: default" && !flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFV(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: default" && flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: bnd_faces" && !flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFaces(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: bnd_faces" && flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFacesGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);
    return op;

    // MFD methods
  } else if (!flag) {
    auto op = Teuchos::rcp(new PDE_DiffusionMFD(oplist, mesh));
    op->Init(oplist);
    op->SetBCs(bc, bc);
    return op;

  } else {
    auto op = Teuchos::rcp(new PDE_DiffusionMFDwithGravity(oplist, mesh, g));
    op->Init(oplist);
    op->SetBCs(bc, bc);
    op->SetDensity(rho);
    return op;
  }
}


/* ******************************************************************
* Initialization of straight diffusion operator: method 1.
****************************************************************** */
Teuchos::RCP<PDE_Diffusion>
PDE_DiffusionFactory::Create(Teuchos::ParameterList& oplist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             const Teuchos::RCP<BCs>& bc)
{
  std::string name = oplist.get<std::string>("discretization primary");
  bool fractured_matrix = oplist.isParameter("fracture");

  // FV methods
  if (name == "fv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFV(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: bnd_faces") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFaces(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

    // MFD methods with non-uniform DOFs
  } else if (fractured_matrix) {
    AmanziGeometry::Point g(mesh->space_dimension()); // gravity should be turned-off in PList
    auto op = Teuchos::rcp(new PDE_DiffusionFracturedMatrix(oplist, mesh, 0.0, g));
    op->Init(oplist);
    op->SetBCs(bc, bc);
    return op;

    // MFD methods
  } else {
    auto op = Teuchos::rcp(new PDE_DiffusionMFD(oplist, mesh));
    op->Init(oplist);
    op->SetBCs(bc, bc);
    return op;
  }
}


/* ******************************************************************
* Initialization of straight diffusion operator: method 2.
****************************************************************** */
Teuchos::RCP<PDE_Diffusion>
PDE_DiffusionFactory::Create(Teuchos::ParameterList& oplist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  std::string name = oplist.get<std::string>("discretization primary");
  bool fractured_matrix = oplist.isParameter("fracture");

  if (name == "fv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionFV(oplist, mesh));
    return op;

  } else if (name == "nlfv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFV(oplist, mesh));
    return op;

  } else if (name == "nlfv: bnd_faces") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFaces(oplist, mesh));
    return op;

    // MFD methods with non-uniform DOFs
  } else if (fractured_matrix) {
    AmanziGeometry::Point g(mesh->space_dimension());
    auto op = Teuchos::rcp(new PDE_DiffusionFracturedMatrix(oplist, mesh, 0.0, g));
    op->Init(oplist);
    return op;

    // MFD methods
  } else {
    auto op = Teuchos::rcp(new PDE_DiffusionMFD(oplist, mesh));
    op->Init(oplist);
    return op;
  }
}


/* ******************************************************************
* Initialization of straight diffusion operator: method 3.
****************************************************************** */
Teuchos::RCP<PDE_Diffusion>
PDE_DiffusionFactory::Create(Teuchos::ParameterList& oplist,
                             const Teuchos::RCP<Operator>& global_op)
{
  std::string name = oplist.get<std::string>("discretization primary");

  if (name == "fv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionFV(oplist, global_op));
    return op;

    // NLFV methods
  } else if (name == "nlfv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFV(oplist, global_op));
    return op;

  } else if (name == "nlfv: bnd_faces") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFaces(oplist, global_op));
    return op;

    // MFD methods
  } else {
    auto op = Teuchos::rcp(new PDE_DiffusionMFD(oplist, global_op));
    op->Init(oplist);
    return op;
  }
}


/* ******************************************************************
* Initialization of diffusion operator with gravity: method 1.
*
* With gravity, assumed vector, temporally varying density.
* Used by ATS.
****************************************************************** */
Teuchos::RCP<PDE_DiffusionWithGravity>
PDE_DiffusionFactory::CreateWithGravity(Teuchos::ParameterList& oplist,
                                        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                        const Teuchos::RCP<BCs>& bc)
{
  std::string name = oplist.get<std::string>("discretization primary");

  if (name == "fv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionFVwithGravity(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

    // NLFV methods
  } else if (name == "nlfv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithGravity(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: bnd_faces") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFacesGravity(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;

    // MFD methods
  } else {
    auto op = Teuchos::rcp(new PDE_DiffusionMFDwithGravity(oplist, mesh));
    op->Init(oplist);
    op->SetBCs(bc, bc);
    return op;
  }
}


/* ******************************************************************
* Initialization of diffusion operator with gravity: method 2.
* With gravity, assumed vector, temporally varying density.
* Used by ATS.
****************************************************************** */
Teuchos::RCP<PDE_DiffusionWithGravity>
PDE_DiffusionFactory::CreateWithGravity(Teuchos::ParameterList& oplist,
                                        const Teuchos::RCP<Operator>& global_op,
                                        const Teuchos::RCP<BCs>& bc)
{
  std::string name = oplist.get<std::string>("discretization primary");
  bool manifolds = oplist_.isParameter("manifolds");

  if (name == "fv: default" && manifolds) {
    auto op = Teuchos::rcp(new PDE_DiffusionFVonManifolds(oplist_, global_op));
    return op;

  } else if (name == "fv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionFVwithGravity(oplist, global_op));
    op->SetBCs(bc, bc);
    return op;

    // NLFV
  } else if (name == "nlfv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithGravity(oplist, global_op));
    op->SetBCs(bc, bc);
    return op;

  } else if (name == "nlfv: bnd_faces") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFacesGravity(oplist, global_op));
    op->SetBCs(bc, bc);
    return op;

    // MFD methods
  } else {
    auto op = Teuchos::rcp(new PDE_DiffusionMFDwithGravity(oplist, global_op));
    op->Init(oplist);
    op->SetBCs(bc, bc);
    return op;
  }
}


/* ******************************************************************
* Initialization of diffusion operator with gravity: method 3.
* With gravity, assumed vector, temporally varying density.
* Used by ATS.
****************************************************************** */
Teuchos::RCP<PDE_DiffusionWithGravity>
PDE_DiffusionFactory::CreateWithGravity(Teuchos::ParameterList& oplist,
                                        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  std::string name = oplist.get<std::string>("discretization primary");

  if (name == "fv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionFVwithGravity(oplist, mesh));
    return op;

    // NLFV methods
  } else if (name == "nlfv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithGravity(oplist, mesh));
    return op;

  } else if (name == "nlfv: bnd_faces") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFacesGravity(oplist, mesh));
    return op;

    // MFD methods
  } else {
    auto op = Teuchos::rcp(new PDE_DiffusionMFDwithGravity(oplist, mesh));
    op->Init(oplist);
    return op;
  }
}


/* ******************************************************************
* Initialization of diffusion operator with gravity: method 4.
* With gravity, assumed vector, temporally varying density.
* Used by ATS.
****************************************************************** */
Teuchos::RCP<PDE_DiffusionWithGravity>
PDE_DiffusionFactory::CreateWithGravity(Teuchos::ParameterList& oplist,
                                        const Teuchos::RCP<Operator>& global_op)
{
  std::string name = oplist.get<std::string>("discretization primary");

  if (name == "fv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionFVwithGravity(oplist, global_op));
    return op;

    // NLFV methods
  } else if (name == "nlfv: default") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithGravity(oplist, global_op));
    return op;

  } else if (name == "nlfv: bnd_faces") {
    auto op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFacesGravity(oplist, global_op));
    return op;

    // MFD methods
  } else {
    auto op = Teuchos::rcp(new PDE_DiffusionMFDwithGravity(oplist, global_op));
    op->Init(oplist);
    return op;
  }
}

} // namespace Operators
} // namespace Amanzi
