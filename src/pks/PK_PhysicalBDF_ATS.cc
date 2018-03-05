/* -------------------------------------------------------------------------
  Author: Daniil Svyatsky

  Default base with a few methods implemented for ATS
------------------------------------------------------------------------- */

#include "PK_PhysicalBDF_ATS.hh"
#include "State.hh"
#include "TreeVector.hh"

namespace Amanzi {

PK_PhysicalBDF_ATS::PK_PhysicalBDF_ATS(Teuchos::ParameterList& FElist,
                                       const Teuchos::RCP<Teuchos::ParameterList>& plist,
                                       const Teuchos::RCP<State>& S,
                                       const Teuchos::RCP<TreeVector>& solution) :
    PK_PhysicalBDF(plist, FElist, solution)
{
  // domain -- default is the entire mesh, no prefix
  if (domain_.empty()) {
    domain_ = plist_->get<std::string>("domain name", std::string("domain"));
  }
  
  if (key_.empty()) {
    key_ = plist_->get<std::string>("primary variable");
  }

  // set up the primary variable solution, and its evaluator
  Teuchos::ParameterList& pv_sublist = FElist.sublist(key_);
  pv_sublist.set("evaluator name", key_);
  pv_sublist.set("field evaluator type", "primary variable");
}


PK_PhysicalBDF_ATS::PK_PhysicalBDF_ATS(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                                       Teuchos::ParameterList& FElist,
                                       const Teuchos::RCP<TreeVector>& solution):
  PK_PhysicalBDF(plist, FElist, solution)
{
  // domain -- default is the entire mesh, no prefix
  if (domain_.empty()) {
    domain_ = plist_->get<std::string>("domain name", std::string("domain"));
  }
  
  if (key_.empty()) {
    key_ = plist_->get<std::string>("primary variable");
  }

  // set up the primary variable solution, and its evaluator
  Teuchos::ParameterList& pv_sublist = FElist.sublist(key_);
  pv_sublist.set("evaluator name", key_);
  pv_sublist.set("field evaluator type", "primary variable");
}


void PK_PhysicalBDF_ATS::Setup(const Teuchos::Ptr<State>& S) {
  // THIS MAY BE CALLED MORE THAN ONCE!
  name_ = plist_->get<std::string>("PK name");

  // set up the VerboseObject
  vo_ = Teuchos::rcp(new VerboseObject(name_, *plist_));

  // get the mesh
  mesh_ = S->GetMesh(domain_);

  // set up the debugger
  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *plist_));

  // require primary variable evaluator
  S->RequireFieldEvaluator(key_);
  Teuchos::RCP<FieldEvaluator> fm = S->GetFieldEvaluator(key_);
  solution_evaluator_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  ASSERT(solution_evaluator_ != Teuchos::null);

  // initial timestep
  dt_ = plist_->get<double>("initial time step", 1.);

  // preconditioner assembly
  assemble_preconditioner_ = plist_->get<bool>("assemble preconditioner", true);

  if (!plist_->get<bool>("strongly coupled PK", false)) {
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");
    // -- check if continuation method
    // -- ETC Note this needs fixed if more than one continuation method used
    if (bdf_plist.isSublist("continuation parameters")) {
      S->RequireScalar("continuation_parameter", name_);
    }
  }

  // convergence criteria
  if (conserved_key_.empty()) {
    if (plist_->isParameter("conserved quantity suffix")) {
      Key conserved_default = getKey(domain_, plist_->get<std::string>("conserved quantity suffix"));
      conserved_key_ = plist_->get<std::string>("conserved quantity key", conserved_default);
    } else {
      conserved_key_ = plist_->get<std::string>("conserved quantity key");
    }
  }
  S->RequireField(conserved_key_)->SetMesh(mesh_)
      ->AddComponent("cell",AmanziMesh::CELL,true);
  S->RequireFieldEvaluator(conserved_key_);

  if (cell_vol_key_.empty()) {
    cell_vol_key_ = plist_->get<std::string>("cell volume key",
            getKey(domain_, "cell_volume"));
  }
  S->RequireField(cell_vol_key_)->SetMesh(mesh_)
      ->AddComponent("cell",AmanziMesh::CELL,true);
  S->RequireFieldEvaluator(cell_vol_key_);
  
  atol_ = plist_->get<double>("absolute error tolerance",1.0);
  rtol_ = plist_->get<double>("relative error tolerance",1.0);
  fluxtol_ = plist_->get<double>("flux error tolerance",1.0);
}


void PK_PhysicalBDF_ATS::Initialize(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<Field> field = S->GetField(key_, name_);

  if (!field->initialized()) {
    // initial conditions
    // -- Get the IC function plist.
    if (!plist_->isSublist("initial condition")) {
      std::stringstream messagestream;
      messagestream << name_ << " has no initial condition parameter list.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }

    // -- Calculate the IC.
    Teuchos::ParameterList ic_plist = plist_->sublist("initial condition");
    field->Initialize(ic_plist);

    // -- Update faces from cells if needed.
    // if (ic_plist.get<bool>("initialize faces from cells", false)) {
    //   DeriveFaceValuesFromCellValues_(field->GetFieldData().ptr());
    // }

    // communicate just to make sure values are initialized for valgrind's sake
    field->GetFieldData()->ScatterMasterToGhosted();
    solution_evaluator_->SetFieldAsChanged(S);
  }

  // -- Push the data into the solution.
  solution_->SetData(field->GetFieldData());

  // set up the timestepping algorithm
  if (!plist_->get<bool>("strongly coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");
    bdf_plist.set("initial time", S->time());
    time_stepper_ = Teuchos::rcp(new BDF1_TI<TreeVector,TreeVectorSpace>(*this, bdf_plist, solution_));

    // initialize continuation parameter if needed.
    if (bdf_plist.isSublist("continuation parameters")) {
      *S->GetScalarData("continuation_parameter", name_) = 1.;
      S->GetField("continuation_parameter", name_)->set_initialized();
    }

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->SetInitialState(S->time(), solution_, solution_dot);
  }
}


void PK_PhysicalBDF_ATS::ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u) {
  if (u->HasComponent("face")) {
    Epetra_MultiVector& u_f = *u->ViewComponent("face",false);
    unsigned int nfaces = u_f.MyLength();
    for (unsigned int f=0; f!=nfaces; ++f) {
      if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_f[0][f] = bc_values_[f];
      }
    }
  } else if (u->HasComponent("boundary_face")) {
    const Epetra_Map& vandalay_map = mesh_->exterior_face_map();
    const Epetra_Map& face_map = mesh_->face_map(false);

    Epetra_MultiVector& u_bf = *u->ViewComponent("boundary_face",false);
    unsigned int nfaces = u_bf.MyLength();
    for (unsigned int bf=0; bf!=nfaces; ++bf) {
      AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));
      if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_bf[0][bf] = bc_values_[f];
      }
    }
  }    
}


double PK_PhysicalBDF_ATS::BoundaryValue(
    const Teuchos::RCP<const Amanzi::CompositeVector>& solution, int face_id) {
  double value = 0.0;

  if (solution->HasComponent("face")){
    const Epetra_MultiVector& u = *solution -> ViewComponent("face",false);
    value = u[0][face_id];
  }
  else if  (solution->HasComponent("boundary_face")){
    const Epetra_MultiVector& u = *solution -> ViewComponent("boundary_face",false);
    const Epetra_Map& fb_map = mesh_->exterior_face_map();
    const Epetra_Map& f_map = mesh_->face_map(false);

    int face_gid = f_map.GID(face_id);
    int face_lbid = fb_map.LID(face_gid);

    value =  u[0][face_lbid];
  }
  else{
    Errors::Message msg("No component is defined for boundary faces\n");
    Exceptions::amanzi_throw(msg);
  }

  return value;
}


void PK_PhysicalBDF_ATS::ChangedSolution() {
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
}


double PK_PhysicalBDF_ATS::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                     Teuchos::RCP<const TreeVector> du) 
{

  S_next_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& conserved = *S_next_->GetFieldData(conserved_key_)
      ->ViewComponent("cell",true);
  const Epetra_MultiVector& cv = *S_next_->GetFieldData(cell_vol_key_)
      ->ViewComponent("cell",true);

  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "ENorm (Infnorm) of: " << conserved_key_ << ": " << std::endl;

  Teuchos::RCP<const CompositeVector> dvec = du->Data();
  double h = S_next_->time() - S_inter_->time();

  double enorm_val = 0.0;
  for (CompositeVector::name_iterator comp=dvec->begin();
       comp!=dvec->end(); ++comp) {
    double enorm_comp = 0.0;
    int enorm_loc = -1;
    const Epetra_MultiVector& dvec_v = *dvec->ViewComponent(*comp, false);

    if (*comp == std::string("cell")) {
      // error done relative to extensive, conserved quantity
      int ncells = dvec->size(*comp,false);
      for (unsigned int c=0; c!=ncells; ++c) {
        double enorm_c = std::abs(h * dvec_v[0][c])
            / (atol_*cv[0][c] + rtol_*std::abs(conserved[0][c]));

        if (enorm_c > enorm_comp) {
          enorm_comp = enorm_c;
          enorm_loc = c;
        }
      }
      
    } else if (*comp == std::string("face")) {
      // error in flux -- relative to cell's extensive conserved quantity
      int nfaces = dvec->size(*comp, false);

      for (unsigned int f=0; f!=nfaces; ++f) {
        AmanziMesh::Entity_ID_List cells;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &cells);
        double cv_min = cells.size() == 1 ? cv[0][cells[0]]
            : std::min(cv[0][cells[0]],cv[0][cells[1]]);
        double conserved_min = cells.size() == 1 ? conserved[0][cells[0]]
            : std::min(conserved[0][cells[0]],conserved[0][cells[1]]);
      
        double enorm_f = fluxtol_ * h * std::abs(dvec_v[0][f])
            / (atol_*cv_min + rtol_*std::abs(conserved_min));
        if (enorm_f > enorm_comp) {
          enorm_comp = enorm_f;
          enorm_loc = f;
        }
      }

    } else {
      double norm;
      dvec_v.Norm2(&norm);
      ASSERT(norm < 1.e-15);
    }

    // Write out Inf norms too.
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      double infnorm(0.);
      dvec_v.NormInf(&infnorm);

      ENorm_t err;
      ENorm_t l_err;
      l_err.value = enorm_comp;
      l_err.gid = dvec_v.Map().GID(enorm_loc);

      int ierr = MPI_Allreduce(&l_err, &err, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      ASSERT(!ierr);
      *vo_->os() << "  ENorm (" << *comp << ") = " << err.value << "[" << err.gid << "] (" << infnorm << ")" << std::endl;
    }

    enorm_val = std::max(enorm_val, enorm_comp);
  }

  double enorm_val_l = enorm_val;
  int ierr = MPI_Allreduce(&enorm_val_l, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  ASSERT(!ierr);
  return enorm_val;
};


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void PK_PhysicalBDF_ATS::DeriveFaceValuesFromCellValues_(const Teuchos::Ptr<CompositeVector>& cv) {

  cv->ScatterMasterToGhosted("cell");
  Teuchos::Ptr<const CompositeVector> cv_const(cv);
  const Epetra_MultiVector& cv_c = *cv_const->ViewComponent("cell",true);
  Epetra_MultiVector& cv_f = *cv->ViewComponent("face",false);

  int f_owned = cv_f.MyLength();
  for (int f=0; f!=f_owned; ++f) {
    AmanziMesh::Entity_ID_List cells;
    cv->Mesh()->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += cv_c[0][cells[n]];
    }
    cv_f[0][f] = face_value / ncells;
  }
};


bool PK_PhysicalBDF_ATS::AdvanceStep(double t_old, double t_new, bool reinit) {
    double dt = t_new - t_old;

    Teuchos::OSTab out = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "----------------------------------------------------------------" << std::endl
                 << "Advancing: t0 = " << S_inter_->time()
                 << " t1 = " << S_next_->time() << " h = " << dt << std::endl
                 << "----------------------------------------------------------------" << std::endl;

    State_to_Solution(S_next_, *solution_);

    // take a bdf timestep
    double dt_solver;
    bool fail;
    if (true) { // this is here simply to create a context for timer,
      // which stops the clock when it is destroyed at the
      // closing brace.
      fail = time_stepper_->TimeStep(dt, dt_solver, solution_);
    }

    if (!fail) {
      // commit the step as successful
      //    time_stepper_->CommitSolution(dt, solution_);
      //    commit_state(dt, S_next_);

      // update the timestep size
      if (dt_solver < dt_ && dt_solver >= dt) {
        // We took a smaller step than we recommended, and it worked fine (not
        // suprisingly).  Likely this was due to constraints from other PKs or
        // vis.  Do not reduce our recommendation.
      } else {
        dt_ = dt_solver;
      }
    } else {
      // take the decreased timestep size
      dt_ = dt_solver;
    }

    return fail;
  }


void PK_PhysicalBDF_ATS::CommitStep (double t_old, double t_new, const Teuchos::RCP<State>& S) {

    double dt = t_new - t_old;

    if (dt > 0. && time_stepper_ != Teuchos::null)
      time_stepper_->CommitSolution(dt, solution_);

  }

  // update the continuation parameter
  void PK_PhysicalBDF_ATS::UpdateContinuationParameter(double lambda) {
    *S_next_->GetScalarData("continuation_parameter", name_) = lambda;
    ChangedSolution();
  }
}

