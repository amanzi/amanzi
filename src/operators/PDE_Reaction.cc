/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "BilinearFormFactory.hh"

#include "Op_Cell_Schema.hh"
#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "Operator_Schema.hh"
#include "PDE_Reaction.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialize operator from parameter list.
****************************************************************** */
void
PDE_Reaction::InitReaction_(Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList& schema_list = plist.sublist("schema");

  // parse discretization  parameters
  auto base = global_schema_row_.StringToKind(schema_list.get<std::string>("base"));
  mfd_ = WhetStone::BilinearFormFactory::Create(schema_list, mesh_);

  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    local_schema_row_.Init(mfd_, mesh_, base);
    global_schema_row_ = local_schema_row_;

    local_schema_col_ = local_schema_row_;
    global_schema_col_ = global_schema_row_;

    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);

    for (auto it = global_schema_row_.begin(); it != global_schema_row_.end(); ++it) {
      int num;
      AmanziMesh::Entity_kind kind;
      std::tie(kind, std::ignore, num) = *it;
      std::string name(local_schema_row_.KindToString(kind));

      const auto& maps = Amanzi::getMaps(*mesh_, kind);
      cvs->AddComponent(name, kind, maps.first, maps.second, num);
    }

    global_op_ =
      Teuchos::rcp(new Operator_Schema(cvs, cvs, plist, global_schema_row_, global_schema_col_));

  } else {
    // constructor was given an Operator
    global_schema_row_ = global_op_->schema_row();
    global_schema_col_ = global_op_->schema_col();

    mesh_ = global_op_->DomainMap().Mesh();
    local_schema_row_.Init(mfd_, mesh_, base);
    local_schema_col_ = local_schema_row_;
  }

  // register the advection Op
  local_op_ = Teuchos::rcp(new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));
  global_op_->OpPushBack(local_op_);
}


/* ******************************************************************
* Collection of local matrices.
* NOTE: Not all input parameters are used yet.
****************************************************************** */
void
PDE_Reaction::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                             const Teuchos::Ptr<const CompositeVector>& p)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;

  AmanziMesh::Entity_ID_List nodes;
  int d = mesh_->space_dimension();

  WhetStone::DenseMatrix Mcell;
  WhetStone::Tensor Kc(d, 1);

  for (int c = 0; c < ncells_owned; ++c) {
    if (Kpoly_.get()) {
      mfd_->MassMatrix(c, (*Kpoly_)[c], Mcell);
    } else {
      Kc(0, 0) = K_.get() ? (*K_)[0][c] : 1.0;
      mfd_->MassMatrix(c, Kc, Mcell);
    }

    matrix[c] = Mcell;
  }
}


/* ******************************************************************
* Populate containers of elemental matrices using MFD factory.
****************************************************************** */
void
PDE_Reaction::UpdateMatrices(double t)
{
  // verify type of coefficient
  AMANZI_ASSERT(Kpoly_st_.get());

  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;

  for (int c = 0; c < ncells_owned; ++c) {
    double tmp(t);
    int size = static_matrices_[c].size();

    matrix[c] = static_matrices_[c][0];
    for (int i = 1; i < size; ++i) {
      matrix[c] += tmp * static_matrices_[c][i];
      tmp *= t;
    }
  }
}


/* *******************************************************************
* Apply boundary condition to the local matrices
******************************************************************* */
void
PDE_Reaction::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  for (auto bc = bcs_trial_.begin(); bc != bcs_trial_.end(); ++bc) {}
}


/* *******************************************************************
* Space-time coefficients can be pre-processed.
******************************************************************* */
void
PDE_Reaction::CreateStaticMatrices_()
{
  WhetStone::DenseMatrix Mcell;

  static_matrices_.resize(ncells_owned);

  for (int c = 0; c < ncells_owned; ++c) {
    int size = (*Kpoly_st_)[c].size();
    static_matrices_[c].clear();

    for (int i = 0; i < size; ++i) {
      mfd_->MassMatrix(c, (*Kpoly_st_)[c][i], Mcell);
      static_matrices_[c].push_back(Mcell);
    }
  }

  static_matrices_initialized_ = true;
}

} // namespace Operators
} // namespace Amanzi
