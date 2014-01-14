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
  pd_key_ = plist_.get<std::string>("ponded depth key", "ponded_depth");
  dependencies_.insert(pd_key_);
  snow_height_key_ = plist_.get<std::string>("snow height key", "snow_depth");
  dependencies_.insert(snow_height_key_);

  mesh_name_ = plist_.get<std::string>("mesh name", "surface");

  FunctionFactory fac;
  precip_func_ = Teuchos::rcp(fac.Create(plist_.sublist("precipitation function")));
}

SnowDistributionEvaluator::SnowDistributionEvaluator(const SnowDistributionEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    elev_key_(other.elev_key_),
    pd_key_(other.pd_key_),
    snow_height_key_(other.snow_height_key_),
    precip_func_(other.precip_func_),
    assembled_(other.assembled_),
    mesh_name_(other.mesh_name_),
    matrix_(other.matrix_) {}


void SnowDistributionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  if (!assembled_) AssembleOperator_(S);

  // result <-- Q_snow * dt
  double dt = *S->GetScalarData("dt");
  double time = S->time();
  double Qe = (*precip_func_)(&time);
  if (Qe * dt > 0.) {
    result->PutScalar(Qe * dt);
    CompositeVector res(*result);
    res = *result;

    // result += z + h_pd + h_snow
    std::cout << "elev = " << std::endl;
    S->GetFieldData(elev_key_)->Print(std::cout);
    res.Update(1., *S->GetFieldData(elev_key_),
               1., *S->GetFieldData(pd_key_), 1.);
    res.Update(1., *S->GetFieldData(snow_height_key_), 1.);

    // Apply the operator to get a residual
    CompositeVector du(*result);
    du.PutScalar(0.);
    matrix_->Apply(res, du);
    du.Print(std::cout);

    // smooth to find a correction
    res.PutScalar(0.);
    matrix_->ApplyInverse(du, res);
    res.Print(std::cout);

    // scale to ensure no negative precip
    double scaling = 0;
    const Epetra_MultiVector& du_c = *res.ViewComponent("cell",false);
    const Epetra_MultiVector& result_c = *result->ViewComponent("cell",false);
    for (int c=0; c!=du_c.MyLength(); ++c) {
      scaling = std::max(scaling, du_c[0][c]/result_c[0][c]);
    }

    double my_scaling(scaling);
    S->GetMesh(mesh_name_)->get_comm()->MaxAll(&my_scaling, &scaling, 1);
    if (scaling > 1) {
      du.Scale(1./scaling);
    }

    // apply correction
    std::cout << "result, pre-correction = " << std::endl;
    result->Print(std::cout);
    result->Update(-1., du, 1.);
    std::cout << "result, post-correction = " << std::endl;
    result->Print(std::cout);


    // get Q
    result->Scale(1./dt);
  } else {
    result->PutScalar(0.);
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
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  std::vector<Operators::MatrixBC> bc_markers(nfaces, Operators::MATRIX_BC_NULL);
  std::vector<double> bc_values(nfaces, 0.0);

  matrix_ = Operators::CreateMatrixMFD(plist_.sublist("smoothing operator"), mesh);

  matrix_->set_symmetric(true);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Teuchos::null);
  matrix_->InitPreconditioner();
  matrix_->CreateMFDstiffnessMatrices(Teuchos::null);
  matrix_->CreateMFDrhsVectors();
  std::vector<double>& Acc_cells = matrix_->Acc_cells();
  for (int c=0; c!=ncells; ++c) {
    Acc_cells[c] += 1.e-8;
  }

  matrix_->ApplyBoundaryConditions(bc_markers, bc_values);
  matrix_->AssembleGlobalMatrices();
  matrix_->ComputeSchurComplement(bc_markers, bc_values);
  matrix_->UpdatePreconditioner();
}


} //namespace
} //namespace
} //namespace
