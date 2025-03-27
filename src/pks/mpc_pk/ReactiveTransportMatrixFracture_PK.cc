/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Process kernel for coupling Transport and Chemistry PKs in the
  matrix and fracture network.
*/

#include "Chemistry_PK.hh"

#include "ReactiveTransportMatrixFracture_PK.hh"
#include "Transport_PK.hh"
#include "TransportMatrixFracture_PK.hh"
#include "TransportMatrixFractureImplicit_PK.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

// -----------------------------------------------------------------------------
// Standard constructor
// -----------------------------------------------------------------------------
ReactiveTransportMatrixFracture_PK::ReactiveTransportMatrixFracture_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : PK_MPCSubcycled(pk_tree, global_list, S, soln)
{
  coupled_chemistry_pk_ = Teuchos::rcp_dynamic_cast<ChemistryMatrixFracture_PK>(sub_pks_[0]);

  subcycling_ = my_list_->get<bool>("subcycle chemistry", true);
  AMANZI_ASSERT(master_ == 1);

  vo_ = Teuchos::rcp(new VerboseObject("CoupledRT_PK", *global_list));

  // tell each chemistry pk it is operator split, which means it will use TCC next
  // instead of current.  Also tell it to use Amanzi generic passwd
  for (const auto& chem_pk : *coupled_chemistry_pk_) {
    auto chem_pk_plist = Teuchos::sublist(Teuchos::sublist(global_list_, "PKs"), chem_pk->name());
    chem_pk_plist->set("operator split", true);
    chem_pk_plist->set("operator split tag", Tags::COPY.get());
    chem_pk_plist->set("primary variable password", "state");
  }
}


// -----------------------------------------------------------------------------
// Setup delegates work to base PK
// -----------------------------------------------------------------------------
void
ReactiveTransportMatrixFracture_PK::Setup()
{
  mesh_domain_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  // darcy fluxes use non-uniform distribution of DOFs
  // -- darcy flux for matrix
  if (!S_->HasRecord("volumetric_flow_rate")) {
    auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>("volumetric_flow_rate", Tags::DEFAULT, "transport")
      .SetMesh(mesh_domain_)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, mmap, gmap, 1);
  }

  // -- darcy flux for fracture
  if (!S_->HasRecord("fracture-volumetric_flow_rate")) {
    auto cvs = Operators::CreateManifoldCVS(mesh_fracture_);
    *S_->Require<CV_t, CVS_t>("fracture-volumetric_flow_rate", Tags::DEFAULT, "transport")
       .SetMesh(mesh_fracture_)
       ->SetGhosted(true) = *cvs;
  }

  tcc_matrix_key_ = "total_component_concentration";
  tcc_fracture_key_ = "fracture-total_component_concentration";

  // evaluators in fracture
  if (!S_->HasRecord("fracture-mass_density_liquid")) {
    Teuchos::ParameterList& elist = S_->FEList();
    Teuchos::ParameterList& ilist = S_->ICList();

    double rho = ilist.sublist("const_fluid_density").get<double>("value");
    elist.sublist("fracture-mass_density_liquid")
      .sublist("function")
      .sublist("All")
      .set<std::string>("region", "All")
      .set<std::string>("component", "cell")
      .sublist("function")
      .sublist("function-constant")
      .set<double>("value", rho);
    elist.sublist("fracture-mass_density_liquid")
      .set<std::string>("evaluator type", "independent variable");
  }

  // communicate chemistry engine to transport.
  auto ic = coupled_chemistry_pk_->begin();

  auto mpc_implicit = Teuchos::rcp_dynamic_cast<PK_MPC<PK_BDF>>(sub_pks_[1]);
  if (mpc_implicit.get() != NULL) {
    for (auto it = mpc_implicit->begin(); it != mpc_implicit->end(); ++it, ++ic) {
      auto it1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(*it);
      auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(*ic);
      it1->SetupChemistry(ic1);
    }
  } else {
    auto mpc_explicit = Teuchos::rcp_dynamic_cast<PK_MPC<PK>>(sub_pks_[1]);
    for (auto it = mpc_explicit->begin(); it != mpc_explicit->end(); ++it, ++ic) {
      auto it1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(*it);
      auto ic1 = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(*ic);
      it1->SetupChemistry(ic1);
    }
  }

  Amanzi::PK_MPCSubcycled::Setup();
}


// -----------------------------------------------------------------------------
// Calculate minimum of sub PKs timestep sizes
// -----------------------------------------------------------------------------
double
ReactiveTransportMatrixFracture_PK::get_dt()
{
  double dTchem = sub_pks_[0]->get_dt();
  double dTtran = sub_pks_[1]->get_dt();

  if (subcycling_) return dTtran;

  if (dTchem / dTtran > 0.99) dTchem *= 0.5;
  if (dTtran > dTchem) dTtran = dTchem;

  return dTchem;
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible
// -----------------------------------------------------------------------------
bool
ReactiveTransportMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = sub_pks_[1]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // advance the slave, subcycling if needed
  double dTchem = sub_pks_[0]->get_dt();
  double dt_next(dTchem), dt_done(0.0);

  try {
    S_->set_intermediate_time(t_old);
    bool done = false;

    while (!done) {
      // do not overstep
      if (t_old + dt_done + dt_next > t_new) { dt_next = t_new - t_old - dt_done; }
      S_->set_intermediate_time(t_old + dt_done + dt_next);

      fail = coupled_chemistry_pk_->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

      if (fail) {
        dt_next /= 2;
      } else {
        dt_done += dt_next;
        dt_next = sub_pks_[0]->get_dt();
      }

      // no state recovery is made, so the only option is to fail
      if (dt_next < min_dt_)
        Exceptions::amanzi_throw("Failure in coupled reactive transport PK: small timestep.");

      // check for subcycling condition
      done = std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1 * min_dt_;
    }

  } catch (const Errors::Message& chem_error) {
    fail = true;
  } catch (...) {
    fail = true;
  }

  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "\nStep failed. Recovered two tcc fields. Last dt=" << dt_next << "\n\n";
  }

  return fail;
};

} // namespace Amanzi
