/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MatrixMFD.hh"
#include "MatrixMFD_Factory.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

#include "implicit_snow_distribution_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

ImplicitSnowDistributionEvaluator::ImplicitSnowDistributionEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    assembled_(false) {
  my_key_ = plist_.get<std::string>("precipitation snow field", "precipitation_snow");

  // depenedencies
  mesh_name_ = plist_.get<std::string>("domain name", "surface");
  if (mesh_name_ == "domain") {
    cell_vol_key_ = "cell_volume";
  } else {
    cell_vol_key_ = mesh_name_+std::string("_cell_volume");
  }    
  dependencies_.insert(cell_vol_key_);

  elev_key_ = plist_.get<std::string>("elevation key", "elevation");
  dependencies_.insert(elev_key_);
  slope_key_ = plist_.get<std::string>("slope key", "slope_magnitude");
  dependencies_.insert(slope_key_);
  pd_key_ = plist_.get<std::string>("ponded depth key", "ponded_depth");
  dependencies_.insert(pd_key_);
  snow_height_key_ = plist_.get<std::string>("snow height key", "snow_depth");
  dependencies_.insert(snow_height_key_);

  // constant parameters
  kL_ = plist_.get<double>("snow distribution length");
  kdx_ = plist_.get<double>("characteristic horizontal grid size", 0.25);
  ktmax_ = plist_.get<double>("faux integration time", 86400.);
  kS_ = plist_.get<double>("characteristic slope", 1.);
  kCFL_ = plist_.get<double>("Courant number", 1.);
  kSWE_conv_ = plist_.get<double>("SWE-to-snow conversion ratio", 10.);
  tol_ = plist_.get<double>("solver tolerance", 1.e-2);
  atol_ = plist_.get<double>("absolute solver tolerance", 1.e-10);
  max_it_ = plist_.get<int>("max iterations", 10);

  FunctionFactory fac;
  precip_func_ = Teuchos::rcp(fac.Create(plist_.sublist("precipitation function")));
}

ImplicitSnowDistributionEvaluator::ImplicitSnowDistributionEvaluator(const ImplicitSnowDistributionEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    elev_key_(other.elev_key_),
    slope_key_(other.slope_key_),
    pd_key_(other.pd_key_),
    snow_height_key_(other.snow_height_key_),
    cell_vol_key_(other.cell_vol_key_),
    precip_func_(other.precip_func_),
    kL_(other.kL_),
    kdx_(other.kdx_),
    ktmax_(other.ktmax_),
    kS_(other.kS_),
    kCFL_(other.kCFL_),
    kSWE_conv_(other.kSWE_conv_),
    tol_(other.tol_),
    atol_(other.atol_),
    max_it_(other.max_it_),
    assembled_(other.assembled_),
    mesh_name_(other.mesh_name_),
    matrix_(other.matrix_),
    matrix_linsolve_(other.matrix_linsolve_) {}


void ImplicitSnowDistributionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (!assembled_) AssembleOperator_(S);

  double time = S->time();
  double Qe = (*precip_func_)(&time);
  double dt_sim = *S->GetScalarData("dt");

  // NOTE: snow precip comes in SWE, must convert it to snow depth!
  if (Qe * dt_sim * kSWE_conv_ > 0.) {
    // determine scaling of flow

    const double kV = kL_/ktmax_;
    double dt = kCFL_ * kdx_ / kV;

    const double kh0 = Qe * ktmax_ * kSWE_conv_;
    const double nm = std::pow(kh0, 1./3) * std::sqrt(kS_) / kV;

    int nsteps = std::ceil(ktmax_ / dt);
    dt = ktmax_ / nsteps;

    // scale report
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "Snow Distribution: taking " << nsteps << " steps of size " << dt
                 << " s for travel length " << kL_ << " m." << std::endl
                 << "  L     = " << kL_ << std::endl
                 << "  t_max = " << ktmax_ << std::endl
                 << "  V     = " << kV << std::endl
                 << "  Q     = " << kV*kh0 << std::endl
                 << "  -------" << std::endl
                 << "  h0    = " << kh0 << std::endl
                 << "  nm    = " << nm << std::endl
                 << "  V(man)= " << std::pow(kh0,1./3) * std::sqrt(kS_) / nm << std::endl
                 << "  -------" << std::endl;
    }

    // Gather mesh entities
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh(mesh_name_);
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    int nfaces_g = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

    // Gather null boundary conditions (no flux of snow into or out of domain)
    std::vector<Operators::MatrixBC> bc_markers(nfaces_g, Operators::MATRIX_BC_NULL);
    std::vector<double> bc_values(nfaces_g, 0.0);

    // Create temporary work space
    // -- not necessarily ghosted workspace
    Teuchos::RCP<CompositeVector> residual = Teuchos::rcp(new CompositeVector(*result));
    Teuchos::RCP<CompositeVector> result_prev = Teuchos::rcp(new CompositeVector(*result));
    Teuchos::RCP<CompositeVector> dresult = Teuchos::rcp(new CompositeVector(*result));
    
    // -- necessarily ghosted workspace
    CompositeVectorSpace ghosted_space(result->Map());
    ghosted_space.SetGhosted();
    Teuchos::RCP<CompositeVector> hz = Teuchos::rcp(new CompositeVector(ghosted_space));
    Teuchos::RCP<CompositeVector> Krel_c = Teuchos::rcp(new CompositeVector(ghosted_space));

    CompositeVectorSpace Krel_f_space;
    Krel_f_space.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);
    Teuchos::RCP<CompositeVector> Krel_uw = Teuchos::rcp(new CompositeVector(Krel_f_space));

    // Gather dependencies
    // NOTE: this is incorrect approximation... should be | sqrt( grad( z+h_pd+h_sd ) ) |
    const Epetra_MultiVector& slope = *S->GetFieldData(slope_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& cv = *S->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
    Teuchos::RCP<const CompositeVector> elev = S->GetFieldData(elev_key_);
    Teuchos::RCP<const CompositeVector> pd = S->GetFieldData(pd_key_);
    Teuchos::RCP<const CompositeVector> snow_height = S->GetFieldData(snow_height_key_);

    // initialize and begin timestep loop
    result->PutScalar(Qe * ktmax_ * kSWE_conv_);
    for (int istep=0; istep!=nsteps; ++istep) {
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Snow distribution inner timestep " << istep << " with size " << dt << std::endl
                   << "   Qe*t_total (min,max) = ";
        double min, max;
        result->ViewComponent("cell",false)->MinValue(&min);
        result->ViewComponent("cell",false)->MaxValue(&max);
        *vo_->os() << min << ", " << max << std::endl;
      }

      *result_prev = *result;
      double norm0 = 0.;
      bool done = false;
      int ncycle = 0;
      while(!done) {
        // update snow potential, z + h_pd + h_s + Qe*dt
        *hz = *result;
        hz->Update(1., *elev, 1., *pd, 1.);
        hz->Update(1., *snow_height, 1.);

        if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
          *vo_->os() << "  z + h_pd + h_s potential = " << std::endl;
          hz->Print(*vo_->os());
        }

        // update Krel
        {
          Epetra_MultiVector& Krel_c_vec = *Krel_c->ViewComponent("cell",false);
          const Epetra_MultiVector& result_c = *result->ViewComponent("cell",false);
          for (int c=0; c!=ncells; ++c) {
            Krel_c_vec[0][c] = std::pow(std::max(result_c[0][c],0.), 5./3)
              / (std::sqrt(std::max(slope[0][c], 1.e-6)) * nm);
          }
        }

        // communicate and upwind
        Krel_c->ScatterMasterToGhosted();
        hz->ScatterMasterToGhosted();
        {
          Epetra_MultiVector& Krel_uw_vec = *Krel_uw->ViewComponent("face",false);
          const Epetra_MultiVector& Krel_c_vec = *Krel_c->ViewComponent("cell",true);
          const Epetra_MultiVector& hz_c = *hz->ViewComponent("cell",true);
          for (int f=0; f!=nfaces; ++f) {
            AmanziMesh::Entity_ID_List cells;
            mesh->face_get_cells(f,AmanziMesh::USED,&cells);
            if (cells.size() == 1) {
            Krel_uw_vec[0][f] = Krel_c_vec[0][cells[0]];
            } else {
              if (hz_c[0][cells[0]] > hz_c[0][cells[1]]) {
                Krel_uw_vec[0][f] = Krel_c_vec[0][cells[0]];
              } else {
                Krel_uw_vec[0][f] = Krel_c_vec[0][cells[1]];
              }
            }
          }
        }

        // Re-assemble
        Krel_uw->ScatterMasterToGhosted();
        matrix_->CreateMFDstiffnessMatrices(Krel_uw.ptr());
        matrix_->CreateMFDrhsVectors();
        matrix_->ApplyBoundaryConditions(bc_markers, bc_values);

        // Apply the operator to get div flux
        matrix_->ComputeNegativeResidual(*hz, residual.ptr());
        if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
          *vo_->os() << "  DIFFUSION RESIDUAL = " << std::endl;
          residual->Print(*vo_->os());
        }

        // Accumulation
        {
          Epetra_MultiVector& residual_c = *residual->ViewComponent("cell",false);
          const Epetra_MultiVector& result_prev_c = *result_prev->ViewComponent("cell",false);
          const Epetra_MultiVector& result_c = *result->ViewComponent("cell",false);
          unsigned int ncells = residual_c.MyLength();
          for (unsigned int c=0; c!=ncells; ++c) {
            residual_c[0][c] += (result_c[0][c] - result_prev_c[0][c]) * cv[0][c] / dt;
          }
        }
        if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
          *vo_->os() << "  ACCUMULATION RESIDUAL = " << std::endl;
          residual->Print(*vo_->os());
        }

        double norm = 0.;
        residual->NormInf(&norm);
        if (ncycle == 0) norm0 = norm;
        ncycle++;
        if (vo_->os_OK(Teuchos::VERB_HIGH)) {
          *vo_->os() << "  inner iterate " << ncycle << " has error norm " << norm << std::endl;
        }
        if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
          residual->Print(*vo_->os());
        }

        if ((norm0 < atol_) || (norm / norm0 < tol_) || (ncycle > max_it_)) {
          done = true;
          continue;
        }

        // Apply the preconditioner
        matrix_->CreateMFDstiffnessMatrices(Krel_uw.ptr());
        matrix_->CreateMFDrhsVectors();
        matrix_->ApplyBoundaryConditions(bc_markers, bc_values);
        {
          std::vector<double>& Acc_cells = matrix_->Acc_cells();
          unsigned int ncells = Acc_cells.size();
          for (unsigned int c=0; c!=ncells; ++c) {
            Acc_cells[c] += cv[0][c] / dt;
          }
        }

        dresult->PutScalar(0.);
        matrix_linsolve_->ApplyInverse(*residual, *dresult);
        if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
          *vo_->os() << "  precon'd update = " << std::endl;
          dresult->Print(*vo_->os());
        }

        // Update
        result->Update(-1., *dresult, 1.);
        if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
          *vo_->os() << "  new snow depth = ";
          result->Print(*vo_->os());
        }
      }
    }

    // get Q back
    result->Scale(1./(ktmax_ * kSWE_conv_));
  } else {
    result->PutScalar(0.);
  }
}

// This is hopefully never called?
void ImplicitSnowDistributionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& results) {
  ASSERT(0);
}

void
ImplicitSnowDistributionEvaluator::AssembleOperator_(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh(mesh_name_);
  Teuchos::ParameterList& diffusion_plist = plist_.sublist("Diffusion");
  diffusion_plist.set("TPFA", true);
  diffusion_plist.set("scaled constraint equation", true);
  diffusion_plist.set<std::string>("MFD method", "two point flux approximation");
  diffusion_plist.set("TPFA use cells only", true);

  matrix_ = Operators::CreateMatrixMFD(diffusion_plist, mesh);
  matrix_->set_symmetric(true);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Teuchos::null);
  matrix_->InitPreconditioner();

  AmanziSolvers::LinearOperatorFactory<CompositeMatrix,CompositeVector,CompositeVectorSpace> fac;
  matrix_linsolve_ = fac.Create(plist_.sublist("Diffusion Solver"), matrix_);

  assembled_ = true;
}


} //namespace
} //namespace
} //namespace
