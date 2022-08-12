/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "errors.hh"
#include "BilinearFormFactory.hh"

#include "Operator_Schema.hh"
#include "Op_Cell_Schema.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialize operator from parameter list.
****************************************************************** */
void PDE_Abstract::Init_(Teuchos::ParameterList& plist)
{
  Errors::Message msg;

  // parse parameters
  // -- discretization details
  matrix_ = plist.get<std::string>("matrix type");
  grad_on_test_ = plist.get<bool>("gradient operator on test function", true);

  // domain and rang of the operator
  bool symmetric(true);
  Teuchos::ParameterList range, domain;
  if (plist.isSublist("schema range") && plist.isSublist("schema domain")) {
    symmetric = false;
    range = plist.sublist("schema range");
    domain = plist.sublist("schema domain");
  }
  else if (plist.isSublist("schema")) {
    range = plist.sublist("schema");
    domain = range;
  }
  else {
    msg << "Schema mismatch for abstract operator.\n" 
        << "  Use \"schema\" for a square operator.\n"
        << "  Use \"schema range\" and \"schema domain\" for a general operator.\n";
    Exceptions::amanzi_throw(msg);
  }

  // compatibility of two schemas
  auto base = global_schema_row_.StringToKind(domain.get<std::string>("base"));
  auto tmp = global_schema_col_.StringToKind(domain.get<std::string>("base"));
  if (tmp != base) {
    msg << "Schema's base mismatch for abstract operator.\n";
    Exceptions::amanzi_throw(msg);
  }

  // discretization method:
  auto mfd_domain = WhetStone::BilinearFormFactory::Create(domain, mesh_);
  Teuchos::RCP<WhetStone::BilinearForm> mfd_range;
  if (!symmetric) 
    mfd_range = WhetStone::BilinearFormFactory::Create(range, mesh_);
  else
    mfd_range = mfd_domain;

  // At the moment, a bilinear form is based on one
  // element, so we need to specify its non-default location
  std::string location = plist.get<std::string>("factory", "schema domain");
  mfd_ = (location == "schema range") ? mfd_range : mfd_domain;

  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    // -- range schema and cvs
    local_schema_row_.Init(mfd_range, mesh_, base);
    global_schema_row_ = local_schema_row_;
    Teuchos::RCP<CompositeVectorSpace> cvs_row = Teuchos::rcp(new CompositeVectorSpace(cvsFromSchema(global_schema_row_, mesh_, true)));

    // -- domain schema and cvs
    local_schema_col_.Init(mfd_domain, mesh_, base);
    global_schema_col_ = local_schema_col_;
    Teuchos::RCP<CompositeVectorSpace> cvs_col = Teuchos::rcp(new CompositeVectorSpace(cvsFromSchema(global_schema_col_, mesh_, true)));

    global_op_ = Teuchos::rcp(new Operator_Schema(cvs_row, cvs_col, plist, global_schema_row_, global_schema_col_));

  } else {
    // constructor was given an Operator
    global_schema_row_ = global_op_->schema_row();
    global_schema_col_ = global_op_->schema_col();

    mesh_ = global_op_->DomainMap().Mesh();
    local_schema_row_.Init(mfd_range, mesh_, base);
    local_schema_col_.Init(mfd_domain, mesh_, base);
  }

  // register the advection Op
  local_op_ = Teuchos::rcp(new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));
  global_op_->OpPushBack(local_op_);

  // default values
  const auto coef = std::make_shared<CoefficientModel<WhetStone::Tensor> >(nullptr);
  interface_ = Teuchos::rcp(new InterfaceWhetStoneMFD<
      WhetStone::BilinearForm, CoefficientModel<WhetStone::Tensor> >(mfd_, coef));
}


/* ******************************************************************
* Populate containers of elemental matrices using MFD factory.
* NOTE: input parameters are not yet used.
****************************************************************** */
void PDE_Abstract::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                  const Teuchos::Ptr<const CompositeVector>& p)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;
  int d(mesh_->space_dimension());

  WhetStone::DenseMatrix Mcell, Acell, AcellT;
  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;

  if (matrix_ == "mass") {
    for (int c = 0; c < ncells_owned; ++c) {
      interface_->MassMatrix(c, Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "mass inverse") {
    for (int c = 0; c < ncells_owned; ++c) {
      interface_->MassMatrixInverse(c, Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "stiffness") {
    for (int c = 0; c < ncells_owned; ++c) {
      interface_->StiffnessMatrix(c, Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "divergence") {
    for (int c = 0; c < ncells_owned; ++c) {
      mfd_->DivergenceMatrix(c, Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "divergence transpose") {
    for (int c = 0; c < ncells_owned; ++c) {
      mfd_->DivergenceMatrix(c, Acell);
      AcellT.Transpose(Acell);
      matrix[c] = AcellT;
    }
  } else if (matrix_ == "advection" && coef_type_ == CoefType::CONSTANT) {
    if (u->HasComponent("cell")) {
      const Epetra_MultiVector& u_c = *u->ViewComponent("cell", false);

      for (int c = 0; c < ncells_owned; ++c) {
        AmanziGeometry::Point vc(d);
        for (int i = 0; i < d; ++i) vc[i] = u_c[i][c];

        mfd_->AdvectionMatrix(c, vc, Acell, grad_on_test_);
        matrix[c] = Acell;
      }
    } else if (u->HasComponent("node")) {
      u->ScatterMasterToGhosted();
      const Epetra_MultiVector& u_c = *u->ViewComponent("node", true);

      AmanziMesh::Entity_ID_List nodes;

      for (int c = 0; c < ncells_owned; ++c) {
        mesh_->cell_get_nodes(c, &nodes);

        AmanziGeometry::Point vn(d);
        std::vector<AmanziGeometry::Point> vec;

        for (int n = 0; n < nodes.size(); ++n) {
          int v = nodes[n];
          for (int i = 0; i < d; ++i) vn[i] = u_c[i][v];
          vec.push_back(vn);
        }
        
        mfd_->AdvectionMatrix(c, vec, Acell);
        matrix[c] = Acell;
      }
    } else {
      AMANZI_ASSERT(false);
    }
  } else if (matrix_ == "advection" && coef_type_ == CoefType::VECTOR_POLYNOMIAL) {
    for (int c = 0; c < ncells_owned; ++c) {
      mfd_->AdvectionMatrix(c, (*Kvec_poly_)[c], Acell, grad_on_test_);
      matrix[c] = Acell;
    }
  } else {
    Errors::Message msg;
    msg << "Unsupported matrix type = " << matrix_ << "\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Populate containers of elemental matrices using MFD factory.
****************************************************************** */
void PDE_Abstract::UpdateMatrices(double t) 
{
  // verify type of coefficient
  AMANZI_ASSERT(Kvec_stpoly_.get());

  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;

  if (matrix_ == "advection") {
    for (int c = 0; c < ncells_owned; ++c) {
      double tmp(t);
      int size = static_matrices_[c].size();

      matrix[c] = static_matrices_[c][0];
      for (int i = 1; i < size; ++i) {
        matrix[c] += tmp * static_matrices_[c][i]; 
        tmp *= t;
      }
    }
  } else {
    Errors::Message msg;
    msg << "Unsupported value of matrix=" << matrix_ << "\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* *******************************************************************
* Space-time coefficients can be pre-processed.
******************************************************************* */
void PDE_Abstract::CreateStaticMatrices_()
{
  AMANZI_ASSERT(matrix_ == "advection");

  int d(mesh_->space_dimension());
  WhetStone::DenseMatrix Acell;
  WhetStone::VectorPolynomial poly(d, d, 0);

  static_matrices_.resize(ncells_owned);

  for (int c = 0; c < ncells_owned; ++c) {
    int size = (*Kvec_stpoly_)[c][0].size();
    static_matrices_[c].clear();
 
    for (int i = 0; i < size; ++i) {
      for (int k = 0; k < d; ++k) poly[k] = (*Kvec_stpoly_)[c][k][i];
      mfd_->AdvectionMatrix(c, poly, Acell, grad_on_test_);
      static_matrices_[c].push_back(Acell);
    }
  }

  static_matrices_initialized_ = true;
}

}  // namespace Operators
}  // namespace Amanzi

