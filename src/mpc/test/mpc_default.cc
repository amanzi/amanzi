#include <iostream>
#include "stdlib.h"
#include "math.h"

#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "CycleDriver.hh"
#include "IO.hh"
#include "MeshFactory.hh"
#include "Mesh.hh"
#include "mpc_pks_registration.hh"
#include "lapack.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "PK_PhysicalBDF.hh"
#include "State.hh"

static int NEQN = 12;
static double AEQN[12][12] = { { 1.0, 0.4, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0 },
                               { 0.4, 1.0, 0.4, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0 },
                               { 0.3, 0.4, 1.0, 0.4, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0 },
                               { 0.2, 0.3, 0.4, 1.0, 0.4, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.0 },
                               { 0.1, 0.2, 0.3, 0.4, 1.0, 0.4, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1 },
                               { 0.1, 0.1, 0.2, 0.3, 0.4, 1.0, 0.4, 0.3, 0.2, 0.1, 0.1, 0.1 },
                               { 0.1, 0.1, 0.1, 0.2, 0.3, 0.4, 1.0, 0.4, 0.3, 0.2, 0.1, 0.1 },
                               { 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4, 1.0, 0.4, 0.3, 0.2, 0.1 },
                               { 0.0, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4, 1.0, 0.4, 0.3, 0.2 },
                               { 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4, 1.0, 0.4, 0.3 },
                               { 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4, 1.0, 0.4 },
                               { 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4, 1.0 } };
static double SOL[12] = { 14.23243562366807, 17.49990474999346, 20.22979170544159,
                          22.23675460623164, 23.37850928985512, 23.65632688437041,
                          23.65632688437041, 23.37850928985512, 22.23675460623161,
                          20.22979170544158, 17.49990474999348, 14.23243562366806 };

/*
static int NEQN = 4;
static double AEQN[4][4] = { { 1.0, 0.0, 0.0, 0.0 },
                             { 0.0, 1.0, 0.0, 0.0 },
                             { 0.0, 0.0, 1.0, 0.0 },
                             { 0.0, 0.0, 0.0, 1.0 } };
*/

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;


// get/set functionality
double&
FieldTimeStamp(const std::string& field, const Tag& tag, State& S)
{
  return S.GetW<double>("time", Tag(field + "@" + tag.get()), "time");
}


// interpolation functionality
void
InterpolateField(const std::string& field, State& S, double tint, CompositeVector& result)
{
  std::vector<double> times;
  std::vector<Tag> tags;

  for (auto& r : S.GetRecordSet(field)) {
    auto tag = r.first;
    double tnew = FieldTimeStamp(field, tag, S);

    bool flag(true);
    for (auto t : times) {
      if (fabs(tnew - t) < 1e-6) flag = false;
    }

    if (flag) {
      times.push_back(tnew);
      tags.push_back(tag);
    }
  }

  int ntimes = times.size();
  if (ntimes == 1) {
    result = S.Get<CompositeVector>(field, tags[0]);
  } else if (ntimes == 2) {
    const auto& v0 = S.Get<CV_t>(field, tags[0]);
    const auto& v1 = S.Get<CV_t>(field, tags[1]);

    double a = (tint - times[0]) / (times[1] - times[0]);
    result.Update(1.0 - a, v0, a, v1, 0.0);
  } else {
    AMANZI_ASSERT(false);
  }
}


class ImplicitPK : public Amanzi::PK_PhysicalBDF {
 public:
  ImplicitPK(Teuchos::ParameterList& pk_tree,
             const Teuchos::RCP<Teuchos::ParameterList>& glist,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& soln)
    : PK(pk_tree, glist, S, soln),
      PK_PhysicalBDF(pk_tree, glist, S, soln),
      dt_(0.0),
      cfl_(1.0),
      soln_(soln){};

  virtual void Setup() override;
  virtual void Initialize() override;

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  virtual double get_dt() override { return dt_ * cfl_; }
  virtual void set_dt(double dt) override { dt_ = dt; }

  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f) override;

  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu) override;
  virtual void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt) override{};

  virtual bool
  ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override
  {
    return false;
  }

  virtual double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;

  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override { return true; }

  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double dt,
                   Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override
  {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

  virtual void ChangedSolution() override{};

 private:
  double dt_, cfl_;
  int id0_, id1_, my_group_;
  std::string my_field_;
  Tag tag_next_, tag_current_;

  const Teuchos::RCP<TreeVector> soln_;
  static RegisteredPKFactory<ImplicitPK> reg_;
};


RegisteredPKFactory<ImplicitPK> ImplicitPK::reg_("implicit pk");


/* ****************************************************************
* Require fields for all groups.
**************************************************************** */
void
ImplicitPK::Setup()
{
  id0_ = plist_->get<int>("start id");
  id1_ = plist_->get<int>("end id");
  cfl_ = plist_->get<double>("cfl", 1.0);

  int n = id1_ - id0_ + 1;
  my_group_ = (int)id0_ / n;
  my_field_ = "field_group_" + std::to_string(my_group_);

  tag_current_ = Tags::CURRENT;
  tag_next_ = Tags::NEXT;

  int ngroups = NEQN / n;
  for (int i = 0; i < ngroups; ++i) {
    std::string field = "field_group_" + std::to_string(i);

    S_->Require<CV_t, CVS_t>(field, tag_next_, "state")
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, n);

    S_->Require<CV_t, CVS_t>(field, tag_current_, "state")
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, n);

    S_->Require<double>("time", Tag(field + "@" + tag_current_.get()), "time");
    S_->Require<double>("time", Tag(field + "@" + tag_next_.get()), "time");
  }
}


/* ****************************************************************
* Initialize owened field
**************************************************************** */
void
ImplicitPK::Initialize()
{
  dt_ = S_->Get<double>("dt", tag_next_);

  auto solution = S_->GetPtrW<CV_t>(my_field_, tag_next_, "state");
  soln_->SetData(solution);

  int n = id1_ - id0_ + 1;
  auto& uc0 = *S_->GetW<CV_t>(my_field_, tag_current_, "state").ViewComponent("cell");
  auto& uc1 = *S_->GetW<CV_t>(my_field_, tag_next_, "state").ViewComponent("cell");

  for (int c = 0; c < 4; ++c) {
    for (int i = 0; i < n; ++i) {
      uc0[i][c] = 1.0;
      uc1[i][c] = 1.0;
    }
  }

  S_->GetRecordW(my_field_, tag_next_, "state").set_initialized();

  S_->GetW<CV_t>(my_field_, tag_current_, "state") = S_->Get<CV_t>(my_field_, tag_next_);
  S_->GetRecordW(my_field_, tag_current_, "state").set_initialized();
}


/* ****************************************************************
* Advance my field using exact linear solver.
**************************************************************** */
bool
ImplicitPK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  S_->GetW<double>("time", Tag(my_field_ + "@" + tag_current_.get()), "time") = t_old;
  S_->GetW<double>("time", Tag(my_field_ + "@" + tag_next_.get()), "time") = t_new;

  double my_dt = t_new - t_old;
  // auto& uc = *soln_->Data()->ViewComponent("cell");
  auto& uc = *S_->GetW<CV_t>(my_field_, tag_next_, "state").ViewComponent("cell");

  int n(id1_ - id0_ + 1), nrhs(1), ierr;
  std::vector<int> ipiv(n);

  WhetStone::DenseMatrix A(n, n);
  WhetStone::DenseVector rhs(n);

  // interpolate all fields except my field at time t_new
  int ngroups = NEQN / n;
  std::vector<CompositeVector> result(ngroups, S_->Get<CV_t>(my_field_));

  for (int i = 0; i < ngroups; ++i) {
    std::string field = "field_group_" + std::to_string(i);
    if (field != my_field_) { InterpolateField(field, *S_, t_new, result[i]); }
  }

  // implicit solver (1/dt - A0) u0^{n+1} = u0^n/dt - A1 u1^{n + a}
  for (int c = 0; c < 4; ++c) {
    for (int i = 0; i < n; ++i) {
      rhs(i) = uc[i][c] / my_dt;

      for (int j = 0; j < NEQN; ++j) {
        int ngroup = (int)j / n;
        int jloc = j - ngroup * n;

        if (ngroup == my_group_) {
          A(i, jloc) = -AEQN[id0_ + i][j];
        } else {
          auto result_c = *result[ngroup].ViewComponent("cell");
          rhs(i) += AEQN[id0_ + i][j] * result_c[jloc][c];
        }
      }
      A(i, i) += 1.0 / my_dt;
    }

    WhetStone::DGESV_F77(&n, &nrhs, A.Values(), &n, &(ipiv[0]), rhs.Values(), &n, &ierr);

    for (int i = 0; i < n; ++i) { uc[i][c] = rhs(i); }
  }

  return false;
}


/* ****************************************************************
* All fields are used to calculate residual.
**************************************************************** */
void
ImplicitPK::FunctionalResidual(double t_old,
                               double t_new,
                               Teuchos::RCP<TreeVector> u_old,
                               Teuchos::RCP<TreeVector> u_new,
                               Teuchos::RCP<TreeVector> f)
{
  S_->GetW<double>("time", Tag(my_field_ + "@" + tag_current_.get()), "time") = t_old;
  S_->GetW<double>("time", Tag(my_field_ + "@" + tag_next_.get()), "time") = t_new;

  double my_dt = t_new - t_old;
  auto& u0c = *u_old->Data()->ViewComponent("cell");
  auto& u1c = *u_new->Data()->ViewComponent("cell");
  auto& fc = *f->Data()->ViewComponent("cell");

  int n(id1_ - id0_ + 1);

  // interpolate all fields except my field at time t_new
  int ngroups = NEQN / n;
  std::vector<CompositeVector> result(ngroups, S_->Get<CV_t>(my_field_));

  for (int i = 0; i < ngroups; ++i) {
    std::string field = "field_group_" + std::to_string(i);
    if (field != my_field_) { InterpolateField(field, *S_, t_new, result[i]); }
  }

  for (int c = 0; c < 4; ++c) {
    for (int i = 0; i < n; ++i) {
      fc[i][c] = (u1c[i][c] - u0c[i][c]) / my_dt;

      for (int j = 0; j < NEQN; ++j) {
        int ngroup = (int)j / n;
        int jloc = j - ngroup * n;

        if (ngroup == my_group_) {
          fc[i][c] -= AEQN[id0_ + i][j] * u1c[jloc][c];
        } else {
          auto result_c = *result[ngroup].ViewComponent("cell");
          fc[i][c] -= AEQN[id0_ + i][j] * result_c[jloc][c];
        }
      }
    }
  }
}


/* ****************************************************************
* Identity preconsitioner since the system is not stiff.
**************************************************************** */
void
ImplicitPK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  S_->GetW<CV_t>(my_field_, tag_current_, "state") = S_->Get<CV_t>(my_field_, tag_next_);
  S_->GetW<double>("time", Tag(my_field_ + "@" + tag_current_.get()), "time") = t_new;
}


/* ****************************************************************
* Identity preconsitioner since the system is not stiff.
**************************************************************** */
int
ImplicitPK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu)
{
  *pu = *u;
  pu->Scale(dt_);
  return 0;
}


double
ImplicitPK::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  return 0.0;
}


/* ****************************************************************
* The test.
**************************************************************** */
double
RunTest(int test, double dt)
{
  auto comm = Amanzi::getDefaultComm();

  // create mesh
  std::string xmlname = "test/mpc_default.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlname);

  Teuchos::ParameterList rlist = plist->sublist("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, rlist, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 2, 2);

  // create state
  Teuchos::ParameterList slist = plist->sublist("state");
  auto S = Teuchos::rcp(new Amanzi::State(slist));
  S->RegisterMesh("domain", mesh);

  S->Require<double>("time", Tags::CURRENT, "time");

  // allow CD to drive initialization
  plist->sublist("cycle driver").sublist("time periods").sublist("TP 0").sublist("PK tree") =
    plist->sublist("PK tree " + std::to_string(test));

  plist->sublist("cycle driver")
    .sublist("time periods")
    .sublist("TP 0")
    .set<double>("initial time step", dt);

  Amanzi::ObservationData obs_data;
  CycleDriver cd(plist, S, comm, obs_data);
  cd.Go();

  // checking that we created only two pks
  auto vo = Teuchos::rcp(new VerboseObject("MPC_DEFAULT", *plist));
  WriteStateStatistics(*S, *vo);

  // error
  auto& field = *S->Get<CV_t>("field_group_0", Tags::CURRENT).ViewComponent("cell");
  double err(0.0), ref(0.0);
  for (int i = 0; i < field.NumVectors(); ++i) {
    err += std::pow(field[i][0] - SOL[i], 2.0);
    ref += std::pow(SOL[i], 2.0);
  }
  err = std::pow(err / ref, 0.5);
  std::cout << "ERR=" << err << std::endl;
  CHECK_CLOSE(0.0, err, 0.01);

  return err;
}


TEST(MPC_DEFAULT_PK)
{
  double e0 = RunTest(1, 0.002);
  double e1 = RunTest(1, 0.001);
  CHECK_CLOSE(2.0, e0 / e1, 0.01);
}

TEST(MPC_DEFAULT_MPC_WEAK_2)
{
  double e0 = RunTest(2, 0.002);
  double e1 = RunTest(2, 0.001);
  CHECK_CLOSE(2.0, e0 / e1, 0.01);
}

TEST(MPC_DEFAULT_MPC_WEAK_2_FLIPPED)
{
  RunTest(3, 0.002);
}

TEST(MPC_DEFAULT_MPC_WEAK_3)
{
  double e0 = RunTest(4, 0.002);
  double e1 = RunTest(4, 0.001);
  CHECK_CLOSE(2.0, e0 / e1, 0.01);
}

TEST(MPC_DEFAULT_MPC_SUBCYCLED)
{
  double e0 = RunTest(5, 0.002);
  double e1 = RunTest(5, 0.001);
  CHECK_CLOSE(2.0, e0 / e1, 0.01);
}

TEST(MPC_DEFAULT_MPC_STRONG)
{
  RunTest(6, 0.002);
}

TEST(MPC_DEFAULT_MPC_TWO_LEVEL_STRONG)
{
  double e0 = RunTest(7, 0.002);
  double e1 = RunTest(7, 0.001);
  CHECK_CLOSE(2.0, e0 / e1, 0.01);
}

TEST(MPC_DEFAULT_MPC_TWO_LEVEL_SUBCYCLED)
{
  double e0 = RunTest(8, 0.002);
  double e1 = RunTest(8, 0.001);
  CHECK_CLOSE(2.0, e0 / e1, 0.01);
}
