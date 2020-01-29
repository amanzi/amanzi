/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#include <sstream>

// Amanzi
#include "dbc.hh"
#include "DenseVector.hh"
#include "MatrixFE.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"

// Operators
#include "Op.hh"
#include "Op_Cell_Cell.hh"
//#include "Op_Cell_Edge.hh"
//#include "Op_Cell_FaceCell.hh"
//#include "Op_Cell_Face.hh"
//#include "Op_Cell_Node.hh"
//#include "Op_Cell_Schema.hh"
//#include "Op_Diagonal.hh"
//#include "Op_Edge_Edge.hh"
#include "Op_Face_Cell.hh"
#include "Op_Face_CellBndFace.hh"
//#include "Op_Face_Schema.hh"
//#include "Op_Node_Node.hh"
//#include "Op_Node_Schema.hh"
#include "Op_SurfaceCell_SurfaceCell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "OperatorUtils.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Deprecated constructor: still supported for compatability
****************************************************************** */
Operator::Operator(const Teuchos::RCP<const CompositeSpace>& cvs,
                   Teuchos::ParameterList& plist,
                   int schema) :
    Operator(cvs, cvs, plist, Schema(schema), Schema(schema)) {}


/* ******************************************************************
* New default constructor for general (rectangular) operator.
* Code of two constructors can be optimized.
****************************************************************** */
Operator::Operator(const Teuchos::RCP<const CompositeSpace>& cvs_row,
                   const Teuchos::RCP<const CompositeSpace>& cvs_col,
                   Teuchos::ParameterList& plist,
                   const Schema& schema_row,
                   const Schema& schema_col) :
    cvs_row_(cvs_row),
    cvs_col_(cvs_col),
    schema_row_(schema_row),
    schema_col_(schema_col),
    shift_(0.0)
{
  mesh_ = cvs_col_->Mesh();
  rhs_ = Teuchos::rcp(new CompositeVector(cvs_row_));

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  if (mesh_->valid_edges()) {
    nedges_owned = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::OWNED);
    nedges_wghost = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::ALL);
  } else {
    nedges_owned = 0;
    nedges_wghost = 0;
  }

  Teuchos::ParameterList vo_list = plist.sublist("Verbose Object");
  vo_ = Teuchos::rcp(new VerboseObject("Operators", vo_list));

  shift_ = plist.get<double>("diagonal shift", 0.0);

  apply_calls_ = 0; 
}



/* ******************************************************************
* Create structure of a global square matrix.
****************************************************************** */
void Operator::SymbolicAssembleMatrix()
{
  // NOTE: not yet implemented, as we have no preconditioner, and the user
  // shouldn't call this anyway (only the PC should!)
  //
  // What still needs to be implemented is GraphFE and MatrixFE (or ditch if no
  // longer necessary) and the SymbolicAssemble sub-calls.
  AMANZI_ASSERT(false);
  // // Create the supermap given a space (set of possible schemas) and a
  // // specific schema (assumed/checked to be consistent with the space).
  //smap_ = createSuperMap(cvs_col_.ptr());

  // // create the graph
  //int row_size = MaxRowSize(*mesh_, schema(), 1);
  //Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->getMap(),
  //    smap_->getGhostedMap(), smap_->getGhostedMap(), row_size));

  // // fill the graph
  //SymbolicAssembleMatrix(*smap_, *graph, 0, 0);

  // // Completing and optimizing the graphs
  //int ierr = graph->FillComplete(smap_->getMap(), smap_->getMap());
  //AMANZI_ASSERT(!ierr);

  // // create global matrix
  //Amat_ = Teuchos::rcp(new MatrixFE(graph));
  //A_ = Amat_->Matrix();
}


/* ******************************************************************
* Create structure of a global matrix.
****************************************************************** */
// void Operator::SymbolicAssembleMatrix(const SuperMap& map, GraphFE& graph,
//                                       int my_block_row, int my_block_col) const
// {
//   // first of double dispatch via Visitor pattern
//   for (auto& it : *this) {
//     it->SymbolicAssembleMatrixOp(this, map, graph, my_block_row, my_block_col);
//   }
// }


// /* ******************************************************************
// * Default visit method for symbolic assemble: Coupling
// ****************************************************************** */
// void Operator::SymbolicAssembleMatrixOp(const Op_Diagonal& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const
// {
//   AMANZI_ASSERT(cvs_col_->HasComponent(op.col_compname()));
//   AMANZI_ASSERT(cvs_row_->HasComponent(op.row_compname()));

//   const std::vector<int>& row_gids = map.GhostIndices(my_block_row, op.row_compname(), 0);
//   const std::vector<int>& col_gids = map.GhostIndices(my_block_col, op.col_compname(), 0);

//   const auto& col_lids = op.col_inds();
//   const auto& row_lids = op.row_inds();

//   std::vector<int> lid_r, lid_c;

//   int ierr(0);
//   for (int n = 0; n != col_lids.size(); ++n) {
//     int ndofs = col_lids[n].size();

//     lid_r.clear();
//     lid_c.clear();

//     for (int i = 0; i != ndofs; ++i) {
//       lid_r.push_back(row_gids[row_lids[n][i]]);
//       lid_c.push_back(col_gids[col_lids[n][i]]);
//     }
//     ierr |= graph.InsertMyIndices(ndofs, lid_r.data(), ndofs, lid_c.data());
//   }
//   AMANZI_ASSERT(!ierr);
// }


/* ******************************************************************
* Populate matrix entries.
****************************************************************** */
void Operator::AssembleMatrix()
{
  // NOTE: not yet implemented, as we have no preconditioner, and the user
  // shouldn't call this anyway (only the PC should!)
  AMANZI_ASSERT(false);

  // if (Amat_ == Teuchos::null) {
  //   Errors::Message msg("Symbolic assembling was not performed.");
  //   Exceptions::amanzi_throw(msg);
  // }

  // Amat_->Zero();
  // AssembleMatrix(*smap_, *Amat_, 0, 0);
  // Amat_->FillComplete();

  // if (shift_ != 0.0) {
  //   Amat_->DiagonalShift(shift_);
  // }
}


/* ******************************************************************
* Populates matrix entries.
****************************************************************** */
// void Operator::AssembleMatrix(const SuperMap& map, MatrixFE& matrix,
//                               int my_block_row, int my_block_col) const
// {
//   // first of double dispatch via Visitor pattern
//   for (auto& it : *this) {
//     it->AssembleMatrixOp(this, map, matrix, my_block_row, my_block_col);
//   }
// }


/* ******************************************************************
* Default visit methods for assemble: Coupling
****************************************************************** */
// void Operator::AssembleMatrixOp(const Op_Diagonal& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const
// {
//   AMANZI_ASSERT(cvs_col_->HasComponent(op.col_compname()));
//   AMANZI_ASSERT(cvs_row_->HasComponent(op.row_compname()));

//   const std::vector<int>& row_gids = map.GhostIndices(my_block_row, op.row_compname(), 0);
//   const std::vector<int>& col_gids = map.GhostIndices(my_block_col, op.col_compname(), 0);

//   const auto& col_lids = op.col_inds();
//   const auto& row_lids = op.row_inds();

//   std::vector<int> lid_r, lid_c;

//   int ierr(0);
//   for (int n = 0; n != col_lids.size(); ++n) {
//     int ndofs = col_lids[n].size();

//     lid_r.clear();
//     lid_c.clear();

//     for (int i = 0; i != ndofs; ++i) {
//       lid_r.push_back(row_gids[row_lids[n][i]]);
//       lid_c.push_back(col_gids[col_lids[n][i]]);
//     }

//     ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[n]);
//   }
//   AMANZI_ASSERT(!ierr);
// }


/* ******************************************************************
* Linear algebra operations with matrices: r = f - A * u.
****************************************************************** */
int Operator::ComputeResidual(const CompositeVector& u, CompositeVector& r, bool zero)
{
  int ierr;
  if (zero) {
    ierr = apply(u, r);
  } else {
    ierr = apply(u, r, -1.0);
  }
  r.update(1.0, *rhs_, -1.0);
  return ierr;
}


/* ******************************************************************
* Linear algebra operations with matrices: r = A * u - f.
****************************************************************** */
int Operator::ComputeNegativeResidual(const CompositeVector& u, CompositeVector& r, bool zero)
{
  int ierr;
  if (zero) {
    ierr = apply(u, r);
  } else {
    ierr = apply(u, r, 1.0);
  }    
  
  r.update(-1.0, *rhs_, 1.0);

  return ierr;
}


/* ******************************************************************
* Parallel matvec product Y = A * X.
******************************************************************* */
int Operator::apply(const CompositeVector& X, CompositeVector& Y, double scalar) const
{
  X.ScatterMasterToGhosted();

  // initialize ghost elements
  if (scalar == 0.0) {
    Y.putScalarMasterAndGhosted(0.0);
  } else if (scalar == 1.0) {
    Y.putScalarGhosted(0.0);
  } else {
    Y.scale(scalar);
    Y.putScalarGhosted(0.0);
  }

  // apply local matrices
  apply_calls_++;
  for (auto& it : *this) it->ApplyMatrixFreeOp(this, X, Y);

  // gather off-process contributions
  Y.GatherGhostedToMaster();
  return 0;
}


/* ******************************************************************
* Parallel matvec product Y = A * X.
* This method is mainly for debugging! Matrix-free apply could better.
******************************************************************* */
int Operator::applyAssembled(const CompositeVector& X, CompositeVector& Y, double scalar) const
{
  X.ScatterMasterToGhosted();

  // initialize ghost elements
  if (scalar == 0.0) {
    Y.putScalarMasterAndGhosted(0.0);
  } else if (scalar == 1.0) {
    Y.putScalarGhosted(0.0);
  } else {
    Y.scale(scalar);
    Y.putScalarGhosted(0.0);
  }

  Vector_type Xcopy(A_->getRowMap());
  Vector_type Ycopy(A_->getRowMap());

  int ierr = copyToSuperVector(*smap_, X, Xcopy);
  A_->apply(Xcopy, Ycopy);
  ierr |= addFromSuperVector(*smap_, Ycopy, Y);

  if (ierr) {
    Errors::Message msg;
    msg << "Operators: ApplyAssemble failed.\n";
    Exceptions::amanzi_throw(msg);
  }

  apply_calls_++;

  return ierr;
}


/* ******************************************************************
* Parallel matvec product Y = A * X.
******************************************************************* */
int Operator::applyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  if (preconditioner_.get() == nullptr) {
    Errors::Message msg("Operator did not initialize a preconditioner.\n");
    Exceptions::amanzi_throw(msg);
  }
  int ierr = preconditioner_->applyInverse(X, Y);
  if (ierr) {
    Errors::Message msg("Operator: ApplyInverse failed.\n");
    Exceptions::amanzi_throw(msg);
  }
  return ierr;
}


/* ******************************************************************
* Default visit method for apply 
****************************************************************** */
// int Operator::ApplyMatrixFreeOp(const Op_Diagonal& op,
//                                 const CompositeVector& X, CompositeVector& Y) const
// {
//   AMANZI_ASSERT(cvs_col_->HasComponent(op.col_compname()));
//   AMANZI_ASSERT(cvs_row_->HasComponent(op.row_compname()));

//   const Epetra_MultiVector& Xf = *X.ViewComponent(op.col_compname(), true);
//   Epetra_MultiVector& Yf = *Y.ViewComponent(op.row_compname(), true);
  
//   const auto& col_lids = op.col_inds();
//   const auto& row_lids = op.row_inds();

//   int ierr(0);
//   for (int n = 0; n != col_lids.size(); ++n) {
//     int ndofs = col_lids[n].size();
//     AMANZI_ASSERT(ndofs == 1);

//     Yf[0][row_lids[n][0]] += op.matrices[n](0, 0) * Xf[0][col_lids[n][0]];
//   }
//   return 0;
// }


/* ******************************************************************
*                       DEPRECATED
* Initialization of the preconditioner. Note that boundary conditions
* may be used in re-implementation of this virtual function.
****************************************************************** */
void Operator::InitPreconditioner(const std::string& prec_name,
                                  const Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory<Operator,CompositeVector> factory;
  preconditioner_ = factory.Create(prec_name, plist);
  UpdatePreconditioner();
}


/* ******************************************************************
*                       DEPRECATED
* Initialization of the preconditioner. Note that boundary conditions
* may be used in re-implementation of this virtual function.
****************************************************************** */
void Operator::InitPreconditioner(Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory<Operator,CompositeVector> factory;
  preconditioner_ = factory.Create(plist);
  UpdatePreconditioner();
}

/* ******************************************************************
* Init owned local operators.
****************************************************************** */
void Operator::Zero()
{
  rhs_->putScalarMasterAndGhosted(0.0);
  int nops = ops_.size();
  for (int i = 0; i < nops; ++i) {
    if (! (ops_properties_[i] & OPERATOR_PROPERTY_DATA_READ_ONLY))
       ops_[i]->Zero();
  }
}

/* ******************************************************************
* Copy the diagonal.
****************************************************************** */
void Operator::getLocalDiagCopy(CompositeVector& X) const
{
  X.putScalarMasterAndGhosted(0.);
  for (auto& op : *this) {
    op->GetLocalDiagCopy(X);
  }
  X.GatherGhostedToMaster();
}

/* ******************************************************************
* Two-stage initialization of preconditioner, part 1.
* Create the preconditioner and set options. Symbolic assemble of 
* operator's matrix must have been called.
****************************************************************** */
void Operator::InitializePreconditioner(Teuchos::ParameterList& plist)
{
  if (smap_.get() == nullptr) {
    smap_ = createSuperMap(getDomainMap().ptr());
  }

  // provide block ids for block strategies.
  //
  // NOTE: not implemented in SuperMap (but ignored for now since we can't use HYPRE anyway)
  // if (plist.isParameter("preconditioner type") &&
  //     plist.get<std::string>("preconditioner type") == "boomer amg" &&
  //     plist.isSublist("boomer amg parameters")) {

  //   // NOTE: Hypre frees this
  //   auto block_ids = smap_->BlockIndices();

  //   plist.sublist("boomer amg parameters").set("number of unique block indices", block_ids.first);

  //   // Note, this passes a raw pointer through a ParameterList.  I was surprised
  //   // this worked too, but ParameterList is a boost::any at heart... --etc
  //   plist.sublist("boomer amg parameters").set("block indices", block_ids.second);
  // }
  AmanziPreconditioners::PreconditionerFactory<Operator,CompositeVector> factory;
  preconditioner_ = factory.Create(plist);
}


/* ******************************************************************
* Two-stage initialization of preconditioner, part 2.
* Set the preconditioner structure. Operator's matrix must have been
* assembled.
****************************************************************** */
void Operator::UpdatePreconditioner()
{
  if (preconditioner_.get() == NULL) {
    Errors::Message msg("Operator has no matrix or preconditioner for update.\n");
    Exceptions::amanzi_throw(msg);
  }

  // pass the preconditioner a non-owning RCP of this
  preconditioner_->Update(Teuchos::rcpFromRef(*this));
}



/* ******************************************************************
* Update the RHS with this vector.
* Note that derived classes may reimplement this with a volume term.
****************************************************************** */
void Operator::UpdateRHS(const CompositeVector& source, bool volume_included) {
  if (!volume_included) {
    Errors::Message msg("Operator::UpdateRHS expected volume_included == True in default implementation.");
    Exceptions::amanzi_throw(msg);
  }
  for (auto& it : *rhs_) {
    if (source.HasComponent(it)) {
      rhs_->GetComponent(it, false)->update(1.0, *source.GetComponent(it, false), 1.0);
    }
  }
}


/* ******************************************************************
* Rescale the local matrices via dispatch.
****************************************************************** */
void Operator::Rescale(double scaling)
{
  for (auto& it : *this) it->Rescale(scaling);
}


/* ******************************************************************
* Rescale the local matrices via dispatch.
****************************************************************** */
void Operator::Rescale(const CompositeVector& scaling)
{
  scaling.ScatterMasterToGhosted();
  for (auto& it : *this) it->Rescale(scaling);
}


/* ******************************************************************
* Rescale the local matrices for particular operator.
****************************************************************** */
void Operator::Rescale(const CompositeVector& scaling, int iops)
{
  AMANZI_ASSERT(iops < ops_.size());
  scaling.ScatterMasterToGhosted();
  ops_[iops]->Rescale(scaling);
}


/* ******************************************************************
* Check points allows us to revert boundary conditions, source terms,
* and accumulation terms. They are useful for operators with constant
* coefficients and varying boundary conditions, e.g. for modeling
* saturated flows.
****************************************************************** */
void Operator::CreateCheckPoint()
{
  rhs_checkpoint_ = Teuchos::rcp(new CompositeVector(*rhs_));
}


void Operator::RestoreCheckPoint()
{
  // The routine should be called after checkpoint is created.
  AMANZI_ASSERT(rhs_checkpoint_ != Teuchos::null);

  // restore accumulation and source terms
  *rhs_ = *rhs_checkpoint_;

  // restore local matrices without boundary conditions
  for (auto& it : *this) {
    it->CopyShadowToMaster();
  }
}


/* ******************************************************************
* New implementation of check-point algorithm.
****************************************************************** */
int Operator::CopyShadowToMaster(int iops) 
{
  int nops = ops_.size();
  AMANZI_ASSERT(iops < nops);
  ops_[iops]->CopyShadowToMaster();

  return 0;
} 


/* ******************************************************************
* Find a matrix block matching the given pattern.
****************************************************************** */
Operator::const_op_iterator
Operator::FindMatrixOp(int schema_dofs, int matching_rule, bool action) const
{
  for (const_op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->Matches(schema_dofs, matching_rule)) {
      return it;
    }
  }

  if (action) {
    Errors::Message msg;
    msg << "Operators: Matching rule " << matching_rule << " not found.\n";
    Exceptions::amanzi_throw(msg);
  }

  return end();
}


/* ******************************************************************
* Find a matrix block matching the given pattern.
****************************************************************** */
Operator::op_iterator
Operator::FindMatrixOp(int schema_dofs, int matching_rule, bool action)
{
  for (op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->Matches(schema_dofs, matching_rule)) {
      return it;
    }
  }

  if (action) {
    Errors::Message msg;
    msg << "Operators: Matching rule " << matching_rule << " not found.\n";
    Exceptions::amanzi_throw(msg);
  }

  return end();
}


/* ******************************************************************
* Push back.
****************************************************************** */
void Operator::OpPushBack(const Teuchos::RCP<Op>& block, int properties) {
  ops_.push_back(block);
  ops_properties_.push_back(properties);
}


/* ******************************************************************
* Add more operators to the existing list. The added operators have
* no special properties. 
****************************************************************** */
void Operator::OpExtend(op_iterator begin, op_iterator end)
{
  int nops = ops_.size();
  int nnew = nops + std::distance(begin, end);

  ops_.reserve(nnew);
  ops_.insert(ops_.end(), begin, end);
  ops_properties_.resize(nnew, 0);  
}


/* ******************************************************************
* Generic error message.
****************************************************************** */
int Operator::SchemaMismatch_(const std::string& schema1, const std::string& schema2) const
{
  std::stringstream err;
  err << "Schemas mismatch: " << schema1 << " != " << schema2;
  Errors::Message message(err.str());
  Exceptions::amanzi_throw(message);
  return 1;
}


/* ******************************************************************
* Populates matrix entries.
****************************************************************** */
std::string Operator::PrintDiagnostics() const
{
  std::stringstream msg;
  for (const_op_iterator it = begin(); it != end(); ++it) {
    msg << "<" << (*it)->schema_string << "> ";
  }
  return msg.str();
}


/* ******************************************************************
* Visit methods for Apply: Cell.
****************************************************************** */
// int Operator::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
//                                 const CompositeVector& X, CompositeVector& Y) const {
//   return SchemaMismatch_(op.schema_string, schema_string_);
// }


// int Operator::ApplyMatrixFreeOp(const Op_Cell_Face& op,
//                                 const CompositeVector& X, CompositeVector& Y) const {
//   return SchemaMismatch_(op.schema_string, schema_string_);
// }


// int Operator::ApplyMatrixFreeOp(const Op_Cell_Node& op,
//                                 const CompositeVector& X, CompositeVector& Y) const {
//   return SchemaMismatch_(op.schema_string, schema_string_);
// }


// int Operator::ApplyMatrixFreeOp(const Op_Cell_Edge& op,
//                                 const CompositeVector& X, CompositeVector& Y) const {
//   return SchemaMismatch_(op.schema_string, schema_string_);
// }


int Operator::ApplyMatrixFreeOp(const Op_Cell_Cell& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


// int Operator::ApplyMatrixFreeOp(const Op_Cell_Schema& op,
//                                 const CompositeVector& X, CompositeVector& Y) const {
//   return SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for Apply: Face
****************************************************************** */
int Operator::ApplyMatrixFreeOp(const Op_Face_Cell& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}

int Operator::ApplyMatrixFreeOp(const Op_Face_CellBndFace& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


// int Operator::ApplyMatrixFreeOp(const Op_Face_Schema& op,
//                                 const CompositeVector& X, CompositeVector& Y) const {
//   return SchemaMismatch_(op.schema_string, schema_string_);
// }


// /* ******************************************************************
// * Visit methods for Apply: Edges
// ****************************************************************** */
// int Operator::ApplyMatrixFreeOp(const Op_Edge_Edge& op,
//                                 const CompositeVector& X, CompositeVector& Y) const {
//   return SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for Apply: Node
****************************************************************** */
// int Operator::ApplyMatrixFreeOp(const Op_Node_Node& op,
//                                 const CompositeVector& X, CompositeVector& Y) const {
//   return SchemaMismatch_(op.schema_string, schema_string_);
// }


// int Operator::ApplyMatrixFreeOp(const Op_Node_Schema& op,
//                                 const CompositeVector& X, CompositeVector& Y) const {
//   return SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for Apply: SurfaceCell
****************************************************************** */
int Operator::ApplyMatrixFreeOp(const Op_SurfaceCell_SurfaceCell& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for Apply: SurfaceFace
****************************************************************** */
int Operator::ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Cell.
****************************************************************** */
// void Operator::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


// void Operator::SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


// void Operator::SymbolicAssembleMatrixOp(const Op_Cell_Node& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


// void Operator::SymbolicAssembleMatrixOp(const Op_Cell_Edge& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


void Operator::SymbolicAssembleMatrixOp(const Op_Cell_Cell& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


// void Operator::SymbolicAssembleMatrixOp(const Op_Cell_Schema& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for symbolic assemble: Face.
****************************************************************** */
void Operator::SymbolicAssembleMatrixOp(const Op_Face_Cell& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}

void Operator::SymbolicAssembleMatrixOp(const Op_Face_CellBndFace& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


// void Operator::SymbolicAssembleMatrixOp(const Op_Face_Schema& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for symbolic assemble: Edge.
****************************************************************** */
// void Operator::SymbolicAssembleMatrixOp(const Op_Edge_Edge& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for symbolic assemble: Node.
****************************************************************** */
// void Operator::SymbolicAssembleMatrixOp(const Op_Node_Node& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


// void Operator::SymbolicAssembleMatrixOp(const Op_Node_Schema& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for symbolic assemble: SurfaceCell
****************************************************************** */
void Operator::SymbolicAssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for symbolic assemble: SurfaceFace.
****************************************************************** */
void Operator::SymbolicAssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for assemble: Cell.
****************************************************************** */
// void Operator::AssembleMatrixOp(const Op_Cell_FaceCell& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


// void Operator::AssembleMatrixOp(const Op_Cell_Face& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


// void Operator::AssembleMatrixOp(const Op_Cell_Node& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


// void Operator::AssembleMatrixOp(const Op_Cell_Edge& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


void Operator::AssembleMatrixOp(const Op_Cell_Cell& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


// void Operator::AssembleMatrixOp(const Op_Cell_Schema& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for assemble: Face.
****************************************************************** */
void Operator::AssembleMatrixOp(const Op_Face_Cell& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}

void Operator::AssembleMatrixOp(const Op_Face_CellBndFace& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


// void Operator::AssembleMatrixOp(const Op_Face_Schema& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for assemble: Edge.
****************************************************************** */
// void Operator::AssembleMatrixOp(const Op_Edge_Edge& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for assemble: Node.
****************************************************************** */
// void Operator::AssembleMatrixOp(const Op_Node_Node& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }

// void Operator::AssembleMatrixOp(const Op_Node_Schema& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const {
//   SchemaMismatch_(op.schema_string, schema_string_);
// }


/* ******************************************************************
* Visit methods for assemble: Surface Cell
****************************************************************** */
void Operator::AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const
{
  std::stringstream err;
  err << "Assemble matrix: invalid schema combination -- " << op.schema_string
      << " cannot be used with a matrix on " << schema_string_;
  Errors::Message message(err.str());
  Exceptions::amanzi_throw(message);
}


/* ******************************************************************
* Visit methods for assemble: Surface Face
****************************************************************** */
void Operator::AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const
{
  std::stringstream err;
  err << "Assemble matrix: invalid schema combination -- " << op.schema_string
      << " cannot be used with a matrix on " << schema_string_;
  Errors::Message message(err.str());
  Exceptions::amanzi_throw(message);
}


/* ******************************************************************
* Local assemble routines (for new schema)
// ****************************************************************** */
// void Operator::ExtractVectorCellOp(int c, const Schema& schema,
//                                    WhetStone::DenseVector& v, const CompositeVector& X) const {
//   Errors::Message msg("Extracton fo local cell-based vector is missing for this operator");
//   Exceptions::amanzi_throw(msg);
// }


void Operator::AssembleVectorCellOp(int c, const Schema& schema,
                                    const WhetStone::DenseVector& v, CompositeVector& X) const {
  Errors::Message msg("Assembly fo local cell-based vector is missing for this operator");
  Exceptions::amanzi_throw(msg);
}


// void Operator::ExtractVectorFaceOp(int c, const Schema& schema,
//                                    WhetStone::DenseVector& v, const CompositeVector& X) const {
//   Errors::Message msg("Extracton fo local cell-based vector is missing for this operator");
//   Exceptions::amanzi_throw(msg);
// }


// void Operator::AssembleVectorFaceOp(int c, const Schema& schema,
//                                     const WhetStone::DenseVector& v, CompositeVector& X) const {
//   Errors::Message msg("Assembly fo local cell-based vector is missing for this operator");
//   Exceptions::amanzi_throw(msg);
// }


// void Operator::ExtractVectorNodeOp(int n, const Schema& schema,
//                                    WhetStone::DenseVector& v, const CompositeVector& X) const {
//   Errors::Message msg("Extracton fo local node-based vector is missing for this operator");
//   Exceptions::amanzi_throw(msg);
// }


// void Operator::AssembleVectorNodeOp(int n, const Schema& schema,
//                                     const WhetStone::DenseVector& v, CompositeVector& X) const {
//   Errors::Message msg("Assembly fo local node-based vector is missing for this operator");
//   Exceptions::amanzi_throw(msg);
// }

}  // namespace Operators
}  // namespace Amanzi

