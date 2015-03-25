/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MatrixMFD.hh"
#include "MatrixMFD_Factory.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

#include "snow_distribution_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

SnowDistributionEvaluator::SnowDistributionEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    assembled_(false) {
  my_key_ = plist_.get<std::string>("precipitation snow field", "precipitation_snow");

  elev_key_ = plist_.get<std::string>("elevation key", "elevation");
  dependencies_.insert(elev_key_);
  slope_key_ = plist_.get<std::string>("slope key", "slope_magnitude");
  dependencies_.insert(slope_key_);
  pd_key_ = plist_.get<std::string>("ponded depth key", "ponded_depth");
  dependencies_.insert(pd_key_);
  snow_height_key_ = plist_.get<std::string>("snow height key", "snow_depth");
  dependencies_.insert(snow_height_key_);

  mesh_name_ = plist_.get<std::string>("domain name", "surface");
  manning_ = plist_.get<double>("manning coefficient", 1.);
  swe_conv_ = plist_.get<double>("SWE-to-snow conversion ratio", 10.);
  
  if (mesh_name_ == "domain") {
    cell_vol_key_ = "cell_volume";
  } else {
    cell_vol_key_ = mesh_name_+std::string("_cell_volume");
  }    
  dependencies_.insert(cell_vol_key_);
  
  FunctionFactory fac;
  precip_func_ = Teuchos::rcp(fac.Create(plist_.sublist("precipitation function")));
}

SnowDistributionEvaluator::SnowDistributionEvaluator(const SnowDistributionEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    elev_key_(other.elev_key_),
    slope_key_(other.slope_key_),
    pd_key_(other.pd_key_),
    snow_height_key_(other.snow_height_key_),
    cell_vol_key_(other.cell_vol_key_),
    precip_func_(other.precip_func_),
    manning_(other.manning_),
    swe_conv_(other.swe_conv_),
    assembled_(other.assembled_),
    mesh_name_(other.mesh_name_),
    matrix_(other.matrix_) {}


void SnowDistributionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  if (!assembled_) AssembleOperator_(S);

  result->PutScalar(0.);

  double dt = *S->GetScalarData("dt");
  std::vector<double> time(1,S->time());
  double Qe = (*precip_func_)(time);
  if (Qe * dt * swe_conv_ > 0.) {
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh(mesh_name_);
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    std::vector<Operators::MatrixBC> bc_markers(nfaces, Operators::OPERATOR_BC_NONE);
    std::vector<double> bc_values(nfaces, 0.0);
    
    Teuchos::RCP<CompositeVector> hz = Teuchos::rcp(new CompositeVector(*result));
    Teuchos::RCP<CompositeVector> residual = Teuchos::rcp(new CompositeVector(*result));
    Teuchos::RCP<CompositeVector> du = Teuchos::rcp(new CompositeVector(*result));

    Teuchos::RCP<CompositeVector> Krel_c = Teuchos::rcp(new CompositeVector(*result));
    CompositeVectorSpace Krel_f_space;
    Krel_f_space.SetMesh(mesh)->SetComponent("face",AmanziMesh::FACE,1);
    Teuchos::RCP<CompositeVector> Krel_uw = Teuchos::rcp(new CompositeVector(Krel_f_space));

    const Epetra_MultiVector& slope = *S->GetFieldData(slope_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& cv = *S->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
            
    int n_cycles = 10;
    int i_cycle = 0;
    double norm;
    result->PutScalar(Qe*dt*swe_conv_);
    for (int n=0; n!=n_cycles; ++n) {
      i_cycle = n;
      std::cout << "SNOW INNER ITERATE " << n << std::endl;

      // update snow potential, z + h_pd + h_s + Qe*dt
      *hz = *result;
      hz->Update(1., *S->GetFieldData(elev_key_),
                 1., *S->GetFieldData(pd_key_), 1.);
      hz->Update(1., *S->GetFieldData(snow_height_key_), 1.);

      // update Krel
      {
        Epetra_MultiVector& Krel_c_vec = *Krel_c->ViewComponent("cell",false);
        const Epetra_MultiVector& result_c = *result->ViewComponent("cell",false);
        for (int c=0; c!=ncells; ++c) {
          Krel_c_vec[0][c] = std::pow(std::max(result_c[0][c],0.), 5./3)
              / (std::sqrt(std::max(slope[0][c], 1.e-4)) * manning_);
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

      // Apply the operator to get a residual
      matrix_->ComputeNegativeResidual(*hz, residual.ptr());

      // Apply the accumulation term
      // Apply the source term
      Epetra_MultiVector& residual_c = *residual->ViewComponent("cell",false);
      const Epetra_MultiVector& result_c = *result->ViewComponent("cell",false);
      for (int c=0; c!=ncells; ++c) {
        residual_c[0][c] += cv[0][c] * result_c[0][c] > 0 ? result_c[0][c] : 0.;
        residual_c[0][c] -= cv[0][c] * dt * Qe * swe_conv_;
      }

      // smooth to find a correction
      residual->Norm2(&norm);
      std::cout << "Snow Distribution: itr=" << i_cycle << ", norm=" << norm << std::endl;
      std::cout << "Qs:" << std::endl;
      result->Print(std::cout);
      std::cout << "Residual:" << std::endl;
      residual->Print(std::cout);

      if (norm < 1.e-4) {
        break;
      }

      du->PutScalar(0.);
      // PC
      std::vector<double>& Acc_cells = matrix_->Acc_cells();
      for (int c=0; c!=ncells; ++c) {
        Acc_cells[c] += cv[0][c];
      }
      matrix_->ApplyBoundaryConditions(bc_markers, bc_values);
      matrix_->ApplyInverse(*residual, *du);

      // apply correction
      result->Update(-1., *du, 1.);

      // // check for negative results and backtrack
      // double factor = 1.;
      // double total_scaling = 1.;
      // double minval = 0;
      // result->ViewComponent("cell",false)->MinValue(&minval);
      // while (minval < 0.) {
      //   factor *= 0.5;
      //   result->Update(factor, *du, 1.);
      //   total_scaling -= factor;
      //   result->ViewComponent("cell",false)->MinValue(&minval);
      // }
      
      // log
      std::cout << "Correction:" << std::endl;
      // du->Scale(total_scaling);
      du->Print(std::cout);

    }
    std::cout << "Snow Distribution finished: itr=" << i_cycle << ", norm=" << norm << std::endl;
    
    // get Q back
    result->Scale(1./(dt * swe_conv_));
  }
}

// This is hopefully never called?
void SnowDistributionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& results) {
  ASSERT(0);
}

void
SnowDistributionEvaluator::AssembleOperator_(const Teuchos::Ptr<State>& S) {
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
