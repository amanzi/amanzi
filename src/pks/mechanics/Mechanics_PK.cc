/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

*/

#include <vector>

#include "InverseFactory.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PorosityEvaluator.hh"
#include "StateArchive.hh"

#include "HydrostaticStressEvaluator.hh"
#include "MechanicsElasticity_PK.hh"

namespace Amanzi {
namespace Mechanics {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Define structure of this PK. We request physical fields and their
* evaluators. Selection of a few models is available and driven by
* model factories, evaluator factories, and parameters of the list
* "physical models and assumptions".
****************************************************************** */
void
Mechanics_PK::Setup()
{
  dt_ = 1e+98;
  mesh_ = S_->GetMesh();
  dim_ = mesh_->getSpaceDimension();

  hydrostatic_stress_key_ = Keys::getKey(domain_, "hydrostatic_stress");
  vol_strain_key_ = Keys::getKey(domain_, "volumetric_strain");

  young_modulus_key_ = Keys::getKey(domain_, "young_modulus");
  poisson_ratio_key_ = Keys::getKey(domain_, "poisson_ratio");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");
  undrained_split_coef_key_ = Keys::getKey(domain_, "undrained_split_coef");

  // constant fields
  S_->Require<AmanziGeometry::Point>("gravity", Tags::DEFAULT, "state");

  // primary fields
  // -- displacement
  const auto& list = ec_list_->sublist("operators").sublist("elasticity operator");
  auto schema = Operators::schemaFromPList(list, mesh_);

  if (!S_->HasRecord(displacement_key_)) {
    auto cvs = Operators::cvsFromSchema(schema, mesh_, true);
    *S_->Require<CV_t, CVS_t>(displacement_key_, Tags::DEFAULT, passwd_)
       .SetMesh(mesh_)
       ->SetGhosted(true) = cvs;

    eval_ = Teuchos::rcp_static_cast<EvaluatorPrimary<CV_t, CVS_t>>(
      S_->GetEvaluatorPtr(displacement_key_, Tags::DEFAULT));
  }

  if (!S_->HasRecord(hydrostatic_stress_key_)) {
    S_->Require<CV_t, CVS_t>(hydrostatic_stress_key_, Tags::DEFAULT, hydrostatic_stress_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  }
  {
    Teuchos::ParameterList elist(hydrostatic_stress_key_);
    elist.set<std::string>("tag", "");
    eval_hydro_stress_ = Teuchos::rcp(new HydrostaticStressEvaluator(elist));
    S_->SetEvaluator(hydrostatic_stress_key_, Tags::DEFAULT, eval_hydro_stress_);
  }

  if (!S_->HasRecord(vol_strain_key_)) {
    S_->Require<CV_t, CVS_t>(vol_strain_key_, Tags::DEFAULT, vol_strain_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1); // copy of states needs it
    S_->RequireEvaluator(vol_strain_key_, Tags::DEFAULT);
  }

  // -- rock properties
  S_->Require<CV_t, CVS_t>(young_modulus_key_, Tags::DEFAULT, passwd_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  S_->Require<CV_t, CVS_t>(poisson_ratio_key_, Tags::DEFAULT, passwd_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  if (!S_->HasRecord(particle_density_key_)) {
    S_->Require<CV_t, CVS_t>(particle_density_key_, Tags::DEFAULT, particle_density_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(particle_density_key_, Tags::DEFAULT);
  }

  // set units
  S_->GetRecordSetW(displacement_key_).set_units("m");
  S_->GetRecordSetW(poisson_ratio_key_).set_units("-");
  S_->GetRecordSetW(young_modulus_key_).set_units("Pa");
}


/* ******************************************************************
* Function goes through flow parameter list and initializes various
* objects including those created during the setup step.
****************************************************************** */
void
Mechanics_PK::Initialize()
{
  // Create verbosity object to print out initialiation statistics.
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "\nPK initialization started...\n";
  }

  // -- mesh dimensions
  ncells_owned_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  ncells_wghost_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  nfaces_owned_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nfaces_wghost_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  nnodes_owned_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  nnodes_wghost_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);
}


/* *******************************************************************
* Initialize boundary conditions
******************************************************************* */
void
Mechanics_PK::InitializeBCs()
{
  Teuchos::RCP<Teuchos::ParameterList> bc_list =
    Teuchos::rcp(new Teuchos::ParameterList(ec_list_->sublist("boundary conditions", true)));

  bcs_.clear();

  // -- displacement
  if (bc_list->isSublist("displacement")) {
    PK_DomainFunctionFactory<MechanicsBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("displacement");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        // nodal dofs
        auto bc =
          bc_factory.Create(spec, "no slip", AmanziMesh::NODE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("no slip");
        bc->set_type(WhetStone::DOF_Type::POINT);
        bc->set_kind(AmanziMesh::NODE);
        bcs_.push_back(bc);

        // bubble dofs
        auto bc2 =
          bc_factory.Create(spec, "no slip", AmanziMesh::FACE, Teuchos::null, Tags::DEFAULT, true);
        bc2->set_bc_name("no slip");
        bc2->set_type(WhetStone::DOF_Type::POINT);
        bc2->set_kind(AmanziMesh::FACE);
        bcs_.push_back(bc2);
      }
    }
  }

  if (bc_list->isSublist("kinematic")) {
    PK_DomainFunctionFactory<MechanicsBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("kinematic");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        // nodal dofs
        auto bc = bc_factory.Create(
          spec, "kinematic", AmanziMesh::NODE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("kinematic");
        bc->set_type(WhetStone::DOF_Type::NORMAL_COMPONENT);
        bc->set_kind(AmanziMesh::NODE);
        bcs_.push_back(bc);

        // bubble dofs
        auto bc2 = bc_factory.Create(
          spec, "kinematic", AmanziMesh::FACE, Teuchos::null, Tags::DEFAULT, true);
        bc2->set_bc_name("kinematic");
        bc2->set_type(WhetStone::DOF_Type::NORMAL_COMPONENT);
        bc2->set_kind(AmanziMesh::FACE);
        bcs_.push_back(bc2);
      }
    }
  }

  if (bc_list->isSublist("traction")) {
    PK_DomainFunctionFactory<MechanicsBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("traction");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        auto bc =
          bc_factory.Create(spec, "traction", AmanziMesh::FACE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("traction");
        bc->set_type(WhetStone::DOF_Type::POINT);
        bc->set_kind(AmanziMesh::FACE);
        bcs_.push_back(bc);
      }
    }
  }

  if (bc_list->isSublist("normal traction")) {
    PK_DomainFunctionFactory<MechanicsBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("normal traction");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        auto bc = bc_factory.Create(
          spec, "normal traction", AmanziMesh::FACE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("normal traction");
        bc->set_type(WhetStone::DOF_Type::NORMAL_COMPONENT);
        bc->set_kind(AmanziMesh::FACE);
        bcs_.push_back(bc);
      }
    }
  }

  // initialize boundary conditions (memory allocation)
  auto bc = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT));
  op_bcs_.push_back(bc);

  bc = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::SCALAR));
  op_bcs_.push_back(bc);

  bc = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::POINT));
  op_bcs_.push_back(bc);

  bc = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  op_bcs_.push_back(bc);
}


/* ******************************************************************
* We have four BCs on two locations (Node and Face) times two types
* (Point or Scalar).
* NOTE: The order of BCs is fixed and important: N+P, N+S, F+P, F+S
****************************************************************** */
void
Mechanics_PK::ComputeOperatorBCs()
{
  int d(mesh_->getSpaceDimension());
  dirichlet_bc_ = 0;

  // two nodal BCs are the first on the list
  for (int i = 0; i < op_bcs_.size(); ++i) {
    std::vector<int>& bc_model = op_bcs_[i]->bc_model();
    for (int n = 0; n < bc_model.size(); n++) {
      bc_model[n] = Operators::OPERATOR_BC_NONE;
    }
  }

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "no slip" && bcs_[i]->type() == WhetStone::DOF_Type::POINT &&
        bcs_[i]->kind() == AmanziMesh::Entity_kind::NODE) {
      std::vector<int>& bc_model = op_bcs_[0]->bc_model();
      std::vector<AmanziGeometry::Point>& bc_value = op_bcs_[0]->bc_value_point();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;
        bc_model[n] = Operators::OPERATOR_BC_DIRICHLET;
        for (int k = 0; k < d; ++k) {
          bc_value[n][k] = it->second[k];
        }
        dirichlet_bc_++;
      }
    }

    if (bcs_[i]->get_bc_name() == "no slip" && bcs_[i]->type() == WhetStone::DOF_Type::POINT &&
        bcs_[i]->kind() == AmanziMesh::Entity_kind::FACE) {
      std::vector<int>& bc_model = op_bcs_[3]->bc_model();
      std::vector<double>& bc_value = op_bcs_[3]->bc_value();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;

        double dot(0.0);
        const auto& normal = mesh_->getFaceNormal(n);
        for (int k = 0; k < d; ++k) {
          dot += it->second[k] * normal[k];
        }

        bc_model[n] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[n] = dot / mesh_->getFaceArea(n);
      }
    }

    if (bcs_[i]->get_bc_name() == "kinematic" &&
        bcs_[i]->type() == WhetStone::DOF_Type::NORMAL_COMPONENT &&
        bcs_[i]->kind() == AmanziMesh::Entity_kind::NODE) {
      std::vector<int>& bc_model = op_bcs_[1]->bc_model();
      std::vector<double>& bc_value = op_bcs_[1]->bc_value();

      if (bcs_[i]->plane_strain_direction() == "x") {
        for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
          int n = it->first;
          bc_model[n] |= Operators::OPERATOR_BC_PLANE_STRAIN_X;
          bc_value[n] = it->second[0];
        }
      } else if (bcs_[i]->plane_strain_direction() == "y") {
        for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
          int n = it->first;
          bc_model[n] |= Operators::OPERATOR_BC_PLANE_STRAIN_Y;
          bc_value[n] = it->second[0];
        }
      } else if (bcs_[i]->plane_strain_direction() == "z") {
        for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
          int n = it->first;
          bc_model[n] |= Operators::OPERATOR_BC_PLANE_STRAIN_Z;
          bc_value[n] = it->second[0];
        }
      } else {
        for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
          int n = it->first;
          bc_model[n] = Operators::OPERATOR_BC_KINEMATIC;
          bc_value[n] = it->second[0];
        }
      }
    }
  }

  // facial BCs make the second group on the list
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "traction") {
      std::vector<int>& bc_model = op_bcs_[2]->bc_model();
      std::vector<AmanziGeometry::Point>& bc_value = op_bcs_[2]->bc_value_point();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;
        bc_model[n] = Operators::OPERATOR_BC_NORMAL_STRESS;
        for (int k = 0; k < d; ++k) {
          bc_value[n][k] = it->second[k];
        }
      }
    }
  }
}


/* ******************************************************************
* Calculate sources and boundary conditions for operators.
****************************************************************** */
void
Mechanics_PK::UpdateSourceBoundaryData(double t_old, double t_new)
{
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(t_old, t_new);
    // bcs_[i]->ComputeSubmodel(mesh_);
  }

  // populate global arrays with data
  ComputeOperatorBCs();
}


/* ******************************************************************
* Add forcing term due to gravity
****************************************************************** */
void
Mechanics_PK::AddGravityTerm(CompositeVector& rhs)
{
  int d = mesh_->getSpaceDimension();
  double g = (S_->Get<AmanziGeometry::Point>("gravity"))[d - 1];

  const auto& rho_c = *S_->Get<CV_t>(particle_density_key_, Tags::DEFAULT).ViewComponent("cell");
  auto& rhs_v = *rhs.ViewComponent("node");

  rhs.PutScalarGhosted(0.0);

  for (int c = 0; c < ncells_owned_; ++c) {
    auto nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();

    double add = g * rho_c[0][c] * mesh_->getCellVolume(c) / nnodes;
    for (int n = 0; n < nnodes; ++n) rhs_v[d - 1][nodes[n]] += add;
  }

  rhs.GatherGhostedToMaster("node");
}


/* ******************************************************************
* Add pressure gradient in poroelasticity model
****************************************************************** */
void
Mechanics_PK::AddPressureGradient(CompositeVector& rhs)
{
  int d = mesh_->getSpaceDimension();
  const auto& p = S_->Get<CV_t>("pressure", Tags::DEFAULT);
  const auto& b_c = *S_->Get<CV_t>("biot_coefficient", Tags::DEFAULT).ViewComponent("cell");

  if (p.HasComponent("face")) {
    const auto& p_f = *p.ViewComponent("face", true);
    auto& rhs_v = *rhs.ViewComponent("node");

    p.ScatterMasterToGhosted("face");
    rhs.PutScalarGhosted(0.0);

    for (int c = 0; c < ncells_owned_; ++c) {
      // pressure gradient
      double vol = mesh_->getCellVolume(c);
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
      const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
      int nfaces = faces.size();

      AmanziGeometry::Point grad(d);
      for (int i = 0; i < nfaces; ++i) {
        int f = faces[i];
        const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
        const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
        grad += normal * (p_f[0][f] * dirs[i]);
      }
      grad /= vol;

      // pressure gradient
      auto nodes = mesh_->getCellNodes(c);
      int nnodes = nodes.size();

      for (int i = 0; i < d; ++i) {
        double add = b_c[0][c] * grad[i] * vol / nnodes;
        for (int n = 0; n < nnodes; ++n) rhs_v[i][nodes[n]] -= add;
      }
    }
  }

  rhs.GatherGhostedToMaster("node");
}


/* ******************************************************************
* Add temperature gradient in thermoelasticity model
****************************************************************** */
void
Mechanics_PK::AddTemperatureGradient(CompositeVector& rhs)
{
  int d = mesh_->getSpaceDimension();
  const auto& temp = S_->Get<CV_t>("temperature", Tags::DEFAULT);

  const auto& E = *S_->Get<CV_t>(young_modulus_key_, Tags::DEFAULT).ViewComponent("cell");
  const auto& nu = *S_->Get<CV_t>(poisson_ratio_key_, Tags::DEFAULT).ViewComponent("cell");

  auto eval = Teuchos::rcp_dynamic_cast<Evaluators::PorosityEvaluator>(
    S_->GetEvaluatorPtr("porosity", Tags::DEFAULT));

  if (temp.HasComponent("face")) {
    const auto& temp_f = *temp.ViewComponent("face", true);
    auto& rhs_v = *rhs.ViewComponent("node");

    temp.ScatterMasterToGhosted("face");
    rhs.PutScalarGhosted(0.0);

    for (int c = 0; c < ncells_owned_; ++c) {
      // compute temperature gradient
      double vol = mesh_->getCellVolume(c);
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
      const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
      int nfaces = faces.size();

      AmanziGeometry::Point grad(d);
      for (int i = 0; i < nfaces; ++i) {
        int f = faces[i];
        const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
        const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
        grad += normal * (temp_f[0][f] * dirs[i]);
      }
      grad /= vol;

      double a = eval->getThermalCoefficients(c).second;
      a *= (d == 3) ? E[0][c] / (1 - 2 * nu[0][c]) / 3 : E[0][c] / (1 - nu[0][c]) / 2;

      // temperature gradient
      auto nodes = mesh_->getCellNodes(c);
      int nnodes = nodes.size();

      for (int i = 0; i < d; ++i) {
        double add = a * grad[i] * vol / nnodes;
        for (int n = 0; n < nnodes; ++n) rhs_v[i][nodes[n]] -= add;
      }
    }
  }

  rhs.GatherGhostedToMaster("node");
}

} // namespace Mechanics
} // namespace Amanzi
