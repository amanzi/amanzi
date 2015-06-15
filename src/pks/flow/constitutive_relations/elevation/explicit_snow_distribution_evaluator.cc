/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MatrixMFD.hh"
#include "MatrixMFD_Factory.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

#include "explicit_snow_distribution_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

ExplicitSnowDistributionEvaluator::ExplicitSnowDistributionEvaluator(Teuchos::ParameterList& plist) :
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
  kCFL_ = plist_.get<double>("Courant number", 0.5);
  kSWE_conv_ = plist_.get<double>("SWE-to-snow conversion ratio", 10.);

  FunctionFactory fac;
  precip_func_ = Teuchos::rcp(fac.Create(plist_.sublist("precipitation function")));
}

ExplicitSnowDistributionEvaluator::ExplicitSnowDistributionEvaluator(const ExplicitSnowDistributionEvaluator& other) :
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
    assembled_(other.assembled_),
    mesh_name_(other.mesh_name_),
    matrix_(other.matrix_) {}


void ExplicitSnowDistributionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (!assembled_) AssembleOperator_(S);

  std::vector<double> time(1,S->time());
  double Qe = (*precip_func_)(time);
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
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
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

    // Gather null boundary conditions (no flux of snow into or out of domain)
    std::vector<Operators::MatrixBC> bc_markers(nfaces, Operators::OPERATOR_BC_NONE);
    std::vector<double> bc_values(nfaces, 0.0);

    // Create temporary work space
    Teuchos::RCP<CompositeVector> hz = Teuchos::rcp(new CompositeVector(*result));
    Teuchos::RCP<CompositeVector> divq = Teuchos::rcp(new CompositeVector(*result));
    Teuchos::RCP<CompositeVector> Krel_c = Teuchos::rcp(new CompositeVector(*result));
    CompositeVectorSpace Krel_f_space;
    Krel_f_space.SetMesh(mesh)->SetComponent("face",AmanziMesh::FACE,1);
    Teuchos::RCP<CompositeVector> Krel_uw = Teuchos::rcp(new CompositeVector(Krel_f_space));

    // Gather dependencies
    // NOTE: this is incorrect approximation... should be | sqrt( grad( z+h_pd ) ) |
    const Epetra_MultiVector& slope = *S->GetFieldData(slope_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& cv = *S->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
    Teuchos::RCP<const CompositeVector> elev = S->GetFieldData(elev_key_);
    Teuchos::RCP<const CompositeVector> pd = S->GetFieldData(pd_key_);
    Teuchos::RCP<const CompositeVector> snow_height = S->GetFieldData(snow_height_key_);

    // initialize and begin timestep loop
    result->PutScalar(Qe * ktmax_ * kSWE_conv_);
    for (int istep=0; istep!=nsteps; ++istep) {
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "Snow distribution inner timestep " << istep << " with size " << dt << std::endl
                   << "   Qe*t_total =" << std::endl;

      // update snow potential, z + h_pd + h_s + Qe*dt
      *hz = *result;
      hz->Update(1., *elev, 1., *pd, 1.);
      hz->Update(1., *snow_height, 1.);

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
      matrix_->CreateMFDstiffnessMatrices(Krel_uw.ptr());
      matrix_->CreateMFDrhsVectors();
      matrix_->ApplyBoundaryConditions(bc_markers, bc_values);

      // Apply the operator to get div flux
      matrix_->ComputeResidual(*hz, divq.ptr());

      // scale by 1/cell volume
      {
        Epetra_MultiVector& divq_c = *divq->ViewComponent("cell",false);
        for (unsigned int c=0; c!=ncells; ++c) {
          divq_c[0][c] /= cv[0][c];
        }
      }

      // update the timestep
      result->Update(dt, *divq, 1.);

      result->Print(std::cout);

      // Ensure non-negativity via clipping
      {
        Epetra_MultiVector& result_c = *result->ViewComponent("cell",false);
        for (int c=0; c!=ncells; ++c) {
          if (result_c[0][c] < 0. && vo_->os_OK(Teuchos::VERB_HIGH)) {
            *vo_->os() << "  Lost snow mass cell " << c << ": h = " << result_c[0][c] << std::endl
                       << "      slope = " << slope[0][c] << std::endl
                       << "      cellv = " << cv[0][c] << std::endl;
          }
          result_c[0][c] = std::max(0., result_c[0][c]);
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
void ExplicitSnowDistributionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& results) {
  ASSERT(0);
}

void
ExplicitSnowDistributionEvaluator::AssembleOperator_(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh(mesh_name_);
  matrix_ = Operators::CreateMatrixMFD(plist_.sublist("Diffusion"), mesh);

  matrix_->set_symmetric(true);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Teuchos::null);
  matrix_->InitPreconditioner();
}


} //namespace
} //namespace
} //namespace
