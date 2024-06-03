/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "Reductions.hh"
#include "PK_Helpers.hh"
#include "PK_PhysicalBDF_Default.hh"

namespace Amanzi {

// constructor
PK_PhysicalBDF_Default::PK_PhysicalBDF_Default(const Comm_ptr_type& comm,
        Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& glist,
        const Teuchos::RCP<State>& S)
  : PK_Physical_Default<PK_BDF_Default>(comm, pk_tree, glist, S)
{}


void
PK_PhysicalBDF_Default::parseParameterList()
{
  PK_Physical_Default<PK_BDF_Default>::parseParameterList();

  // keys used here
  conserved_key_ = Keys::readKey(*plist_, domain_, "conserved quantity");
  cell_vol_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");

  // tolerances for ErrorNorm
  atol_ = plist_->get<double>("absolute error tolerance", 1.0);
  rtol_ = plist_->get<double>("relative error tolerance", 1.0);
  fluxtol_ = plist_->get<double>("flux error tolerance", 1.0);
}


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void
PK_PhysicalBDF_Default::setup()
{
  // call the meat of the base constructurs via Setup methods
  PK_Physical_Default<PK_BDF_Default>::setup();

  // convergence criteria is based on a conserved quantity
  PKHelpers::requireAtNext(conserved_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, true);

  // we also use a copy of the conserved quantity, as this is a better choice
  // in the error norm
  PKHelpers::requireAtCurrent(conserved_key_, tag_current_, *S_, name_, true);

  // cell volume used for ErrorNorm
  PKHelpers::requireAtNext(cell_vol_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, true);
};


int
PK_PhysicalBDF_Default::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                            Teuchos::RCP<TreeVector> Pu)
{
  // default preconditioner is the identify
  Pu->assign(*u);
  return 0;
}


// -----------------------------------------------------------------------------
// Default enorm that uses an abs and rel tolerance to monitor convergence.
// -----------------------------------------------------------------------------
double
PK_PhysicalBDF_Default::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<const TreeVector> res)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "ENorm (Infnorm) of: " << conserved_key_ << ": " << std::endl;

  // Abs tol based on old conserved quantity -- we know these have been vetted
  // at some level whereas the new quantity is some iterate, and may be
  // anything from negative to overflow.
  auto conserved = S_->Get<CompositeVector>(conserved_key_, tag_current_).viewComponent("cell", true);
  auto cv = S_->Get<CompositeVector>(cell_vol_key_, tag_next_).viewComponent("cell", true);

  Teuchos::RCP<const CompositeVector> dvec = res->getData();
  double h = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  double enorm_loc = 0;
  for (const auto& comp : *dvec) {
    Reductions::MaxLoc<double, GO> enorm_comp;
    const auto dvec_v = dvec->viewComponent(comp, false);

    if (comp == "cell") {
      // error done relative to extensive, conserved quantity
      int ncells = dvec_v.extent(0);
      Kokkos::parallel_reduce("PK_PhysicalBDF_Default::ErrorNorm", ncells,
              KOKKOS_LAMBDA(const int& c, Reductions::MaxLoc<double, GO>& lval) {
                double enorm_c = std::abs(h * dvec_v(c,0)) / (atol_ * cv(c,0)
                        + rtol_ * std::abs(conserved(c,0)));
                Reductions::MaxLoc<double,GO> myval(enorm_c, c);
                lval += myval;
              }, enorm_comp);

    } else if (comp == "face") {
      // error in flux -- relative to cell's extensive conserved quantity
      int nfaces = dvec_v.extent(0);
      Kokkos::parallel_reduce("PK_PhysicalBDF_Default::ErrorNorm", nfaces,
              KOKKOS_LAMBDA(const int& f, Reductions::MaxLoc<double, GO>& lval) {
                auto cells = mesh_->getFaceCells(f);
                double cv_min =
                  cells.size() == 1 ? cv(cells[0],0) : std::min(cv(cells[0],0), cv(cells[1],0));
                double conserved_min = cells.size() == 1 ?
                  conserved(cells[0],0) :
                  std::min(conserved(cells[0],0), conserved(cells[1],0));

                double enorm_f = fluxtol_ * h * std::abs(dvec_v(f,0)) /
                         (atol_ * cv_min + rtol_ * std::abs(conserved_min));
                lval += Reductions::MaxLoc<double, GO>(enorm_f, f);
              }, enorm_comp);

    } else {
      // pass...
    }

    // Write out Inf norms too.
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      double infnorm = dvec->getComponent(comp)->getVector(0)->normInf();

      // now do the global reduction in Teuchos, first converting to GID from LID
      enorm_comp.loc = dvec->getMap()->getComponentMap(comp)->getGlobalElement(enorm_comp.loc);
      auto global_enorm_comp = Reductions::reduceAllLoc<double>(*dvec->getMap()->getComm(), enorm_comp);
      *vo_->os() << "  ENorm (" << comp << ") = " << global_enorm_comp.val << "[" << global_enorm_comp.loc << "] ("
                 << infnorm << ")" << std::endl;
    }

    enorm_loc = std::max(enorm_loc, enorm_comp.val);
  }

  double enorm_global = 0;
  Teuchos::reduceAll(*dvec->getMap()->getComm(), Teuchos::REDUCE_MAX, 1, &enorm_loc, &enorm_global);
  return enorm_global;
};


void
PK_PhysicalBDF_Default::commitStep(double t_old, double t_new, const Tag& tag_next)
{
  PK_Physical_Default<PK_BDF_Default>::commitStep(t_old, t_new, tag_next);

  // copy over conserved quantity
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  S_->Get<CompositeVector>(conserved_key_, tag_next).scatterMasterToGhosted();
  PKHelpers::assign(conserved_key_, tag_current, tag_next, *S_);
}


} // namespace Amanzi
