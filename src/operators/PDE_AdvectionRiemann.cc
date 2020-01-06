/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// Amanzi
#include "Basis_Regularized.hh"
#include "BilinearFormFactory.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "Operator_Schema.hh"
#include "Op_Cell_Schema.hh"
#include "Op_Face_Schema.hh"
#include "PDE_AdvectionRiemann.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialize operator from parameter list.
****************************************************************** */
void PDE_AdvectionRiemann::InitAdvection_(Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList& schema_list = plist.sublist("schema");

  // parse discretization parameters
  auto base = global_schema_row_.StringToKind(schema_list.get<std::string>("base"));
  auto mfd = WhetStone::BilinearFormFactory::Create(schema_list, mesh_);

  matrix_ = plist.get<std::string>("matrix type");
  method_ = schema_list.get<std::string>("method");

  if (method_ == "dg modal") {
    dg_ = Teuchos::rcp_dynamic_cast<WhetStone::DG_Modal>(mfd);
  } else {
    Errors::Message msg;
    msg << "Advection operator: method \"" << method_ << "\" is invalid.";
    Exceptions::amanzi_throw(msg);
  }

  // -- fluxes
  flux_ = plist.get<std::string>("flux formula", "Rusanov");
  jump_on_test_ = plist.get<bool>("jump operator on test function", true);

  // constructor was given a mesh
  if (global_op_ == Teuchos::null) {
    local_schema_row_.Init(mfd, mesh_, base);
    global_schema_row_ = local_schema_row_;

    Teuchos::RCP<CompositeVectorSpace> cvs_row = Teuchos::rcp(new CompositeVectorSpace());
    cvs_row->SetMesh(mesh_)->SetGhosted(true);

    int num;
    AmanziMesh::Entity_kind kind;

    for (auto it = global_schema_row_.begin(); it != global_schema_row_.end(); ++it) {
      std::tie(kind, std::ignore, num) = *it;
      std::string name(local_schema_row_.KindToString(kind));
      cvs_row->AddComponent(name, kind, num);
    }

    // -- domain schema and cvs
    local_schema_col_.Init(mfd, mesh_, base);
    global_schema_col_ = local_schema_col_;

    Teuchos::RCP<CompositeVectorSpace> cvs_col = Teuchos::rcp(new CompositeVectorSpace());
    cvs_col->SetMesh(mesh_)->SetGhosted(true);

    for (auto it = global_schema_col_.begin(); it != global_schema_col_.end(); ++it) {
      std::tie(kind, std::ignore, num) = *it;
      std::string name(local_schema_col_.KindToString(kind));
      cvs_col->AddComponent(name, kind, num);
    }

    global_op_ = Teuchos::rcp(new Operator_Schema(cvs_row, cvs_col, plist, global_schema_row_, global_schema_col_));
    if (local_schema_col_.base() == AmanziMesh::CELL) {
      local_op_ = Teuchos::rcp(new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));
    } else if (local_schema_col_.base() == AmanziMesh::FACE) {
      local_op_ = Teuchos::rcp(new Op_Face_Schema(global_schema_row_, global_schema_col_, mesh_));
    }

  // constructor was given an Operator
  } else {
    global_schema_row_ = global_op_->schema_row();
    global_schema_col_ = global_op_->schema_col();

    mesh_ = global_op_->DomainMap().Mesh();
    local_schema_row_.Init(mfd, mesh_, base);
    local_schema_col_.Init(mfd, mesh_, base);

    if (local_schema_col_.base() == AmanziMesh::CELL) {
      local_op_ = Teuchos::rcp(new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));
    } else if (local_schema_col_.base() == AmanziMesh::FACE) {
      local_op_ = Teuchos::rcp(new Op_Face_Schema(global_schema_row_, global_schema_col_, mesh_));
    }
  }

  // register the advection Op
  global_op_->OpPushBack(local_op_);
}


/* ******************************************************************
* A simple first-order transport method of the form div (u C), where 
* u is the given velocity field and C is the advected field.
****************************************************************** */
void PDE_AdvectionRiemann::UpdateMatrices(
    const Teuchos::Ptr<const std::vector<WhetStone::Polynomial> >& u)
{
  double flux;
  WhetStone::DenseMatrix Aface;
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;

  if (matrix_ == "flux" && flux_ == "upwind") {
    for (int f = 0; f < nfaces_owned; ++f) {
      dg_->FluxMatrix(f, (*u)[f], Aface, true, jump_on_test_, &flux);
      matrix[f] = Aface;
    }
  } else if (matrix_ == "flux" && flux_ == "downwind") {
    for (int f = 0; f < nfaces_owned; ++f) {
      dg_->FluxMatrix(f, (*u)[f], Aface, false, jump_on_test_, &flux);
      matrix[f] = Aface;
    }
  } else if (matrix_ == "flux" && flux_ == "upwind at gauss points") {
    for (int f = 0; f < nfaces_owned; ++f) {
      dg_->FluxMatrixGaussPoints(f, (*u)[f], Aface, true, jump_on_test_);
      matrix[f] = Aface;
    }
  } else if (matrix_ == "flux" && flux_ == "downwind at gauss points") {
    for (int f = 0; f < nfaces_owned; ++f) {
      dg_->FluxMatrixGaussPoints(f, (*u)[f], Aface, false, jump_on_test_);
      matrix[f] = Aface;
    }
  } else if (matrix_ == "flux" && flux_ == "Rusanov") {
    // Polynomial Kc should be distributed here
    AmanziMesh::Entity_ID_List cells;
    for (int f = 0; f < nfaces_owned; ++f) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int c1 = cells[0];
      int c2 = (cells.size() == 2) ? cells[1] : c1;
      dg_->FluxMatrixRusanov(f, (*Kc_)[c1], (*Kc_)[c2], (*Kf_)[f], Aface);
      matrix[f] = Aface;
    }
  } else {
    Errors::Message msg;
    msg << "Unsupported matrix type=" << matrix_ << " and/or flux type=" << flux_ << "\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Flux matrices for the case of space-time polynomial velocity.
****************************************************************** */
void PDE_AdvectionRiemann::UpdateMatrices(double t)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;

  for (int f = 0; f < nfaces_owned; ++f) {
    int size = (*uc_)[f].size();

    // calculate dynamic surface flux 
    double flux_t(0.0), tmp(1.0);
    for (int i = 0; i < size; ++i) {
      flux_t += tmp * static_matrices_[f][i].uflux; 
      tmp *= t;
    }

    // sum up matrices with matching signs of point fluxes
    int nrows = static_matrices_[f][0].Uface.NumRows();
    int ncols = static_matrices_[f][0].Uface.NumRows();

    matrix[f].Reshape(nrows, ncols);
    matrix[f].PutScalar(0.0);

    tmp = 1.0;
    for (int i = 0; i < size; ++i) {
      double flux_i = static_matrices_[f][i].uflux;
      double sign = flux_i * flux_t;

      if (matrix_ == "flux" && flux_ == "upwind") {
        if (sign > 0.0) 
          matrix[f] += static_matrices_[f][i].Uface * tmp;
        else 
          matrix[f] += static_matrices_[f][i].Dface * tmp;
      }
      else if (matrix_ == "flux" && flux_ == "downwind") {
        if (sign > 0.0) 
          matrix[f] += static_matrices_[f][i].Dface * tmp;
        else 
          matrix[f] += static_matrices_[f][i].Uface * tmp;
      }
      else {
        Errors::Message msg;
        msg << "Unsupported matrix type=" << matrix_ << " and/or flux type=" << flux_ << "\n";
        Exceptions::amanzi_throw(msg);
      }

      tmp *= t;
    }
  }
}


/* *******************************************************************
* Apply boundary condition to the local matrices
******************************************************************* */
void PDE_AdvectionRiemann::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<std::vector<double> >& bc_value = bcs_trial_[0]->bc_value_vector();
  int nk = bc_value[0].size();

  Epetra_MultiVector& rhs_c = *global_op_->rhs()->ViewComponent("cell", true);

  AmanziMesh::Entity_ID_List cells;

  int dir, d = mesh_->space_dimension();
  std::vector<AmanziGeometry::Point> tau(d - 1);

  // create integration object for all mesh cells
  WhetStone::NumericalIntegration<AmanziMesh::Mesh> numi(mesh_);

  for (int f = 0; f != nfaces_owned; ++f) {
    if (bc_model[f] == OPERATOR_BC_DIRICHLET ||
        bc_model[f] == OPERATOR_BC_DIRICHLET_TYPE2) {
      // common section
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int c = cells[0];

      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c, &dir);

      // --set polynomial with Dirichlet data
      WhetStone::DenseVector coef(nk);
      for (int i = 0; i < nk; ++i) {
        coef(i) = bc_value[f][i];
      }

      WhetStone::Polynomial pf(d, dg_->order(), coef);
      pf.set_origin(xf);

      // -- convert boundary polynomial to regularized space polynomial
      pf.ChangeOrigin(mesh_->cell_centroid(c));

      // -- extract coefficients and update right-hand side 
      WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
      int nrows = Aface.NumRows();
      int ncols = Aface.NumCols();

      WhetStone::DenseVector v(nrows), av(ncols);
      v = pf.coefs();
      dg_->cell_basis(c).ChangeBasisNaturalToMy(v);

      Aface.Multiply(v, av, false);

      // now fork the work flow
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        for (int i = 0; i < ncols; ++i) {
          rhs_c[i][c] -= av(i);
        }
        local_op_->matrices_shadow[f] = Aface;
        Aface.PutScalar(0.0);
      } else {
        for (int i = 0; i < ncols; ++i) {
          rhs_c[i][c] += av(i);
        }
      }
    } 
    else if (bc_model[f] == OPERATOR_BC_REMOVE) {
      local_op_->matrices[f].PutScalar(0.0);
    }
    else if (bc_model[f] != OPERATOR_BC_NONE) {
      AMANZI_ASSERT(false);
    }
  } 
}


/* *******************************************************************
* Identify the advected flux of u
******************************************************************* */
void PDE_AdvectionRiemann::UpdateFlux(
    const Teuchos::Ptr<const CompositeVector>& h,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::RCP<BCs>& bc, const Teuchos::Ptr<CompositeVector>& flux)
{
  h->ScatterMasterToGhosted("cell");
  
  const Epetra_MultiVector& h_f = *h->ViewComponent("face", true);
  const Epetra_MultiVector& u_f = *u->ViewComponent("face", false);
  Epetra_MultiVector& flux_f = *flux->ViewComponent("face", false);

  flux->PutScalar(0.0);
  for (int f = 0; f < nfaces_owned; ++f) {
    flux_f[0][f] = u_f[0][f] * h_f[0][f];
  }  
}


/* *******************************************************************
* Space-time coefficients can be pre-processed.
******************************************************************* */
void PDE_AdvectionRiemann::CreateStaticMatrices_()
{
  AMANZI_ASSERT(matrix_ == "flux");

  static_matrices_.resize(nfaces_owned);

  for (int f = 0; f < nfaces_owned; ++f) {
    int size = (*uc_)[f].size();
    static_matrices_[f].clear();
 
    SurfaceFluxData data;
    for (int i = 0; i < size; ++i) {
      dg_->FluxMatrix(f, (*uc_)[f][i], data.Uface, true, jump_on_test_, &data.uflux);
      dg_->FluxMatrix(f, (*uc_)[f][i], data.Dface, false, jump_on_test_, &data.dflux);
      static_matrices_[f].push_back(data);
    }
  }

  static_matrices_initialized_ = true;
}

}  // namespace Operators
}  // namespace Amanzi
