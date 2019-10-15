/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <sstream>

// TPLs
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Export.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

// Amanzi
#include "dbc.hh"
#include "DenseVector.hh"
#include "MatrixFE.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"

// Operators
#include "Op.hh"
#include "Op_Cell_Cell.hh"
#include "Op_Cell_Edge.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Cell_Face.hh"
#include "Op_Cell_Node.hh"
#include "Op_Cell_Schema.hh"
#include "Op_Diagonal.hh"
#include "Op_Edge_Edge.hh"
#include "Op_Face_Cell.hh"
#include "Op_Face_CellBndFace.hh"
#include "Op_Face_Schema.hh"
#include "Op_Node_Node.hh"
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
Operator::Operator(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                   Teuchos::ParameterList& plist, int schema)
  : cvs_row_(cvs),
    cvs_col_(cvs),
    schema_row_(schema),
    schema_col_(schema),
    shift_(0.0)
{
  mesh_ = cvs_col_->Mesh();
  rhs_ = Teuchos::rcp(new CompositeVector(*cvs_row_, true));

  ncells_owned =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  nfaces_owned =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nnodes_owned =
    mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);

  ncells_wghost =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_wghost =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  nnodes_wghost =
    mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  if (mesh_->valid_edges()) {
    nedges_owned =
      mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::OWNED);
    nedges_wghost =
      mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::ALL);
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
 * New default constructor for general (rectangular) operator.
 * Code of two constructors can be optimized.
 ****************************************************************** */
Operator::Operator(const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
                   const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
                   Teuchos::ParameterList& plist, const Schema& schema_row,
                   const Schema& schema_col)
  : cvs_row_(cvs_row),
    cvs_col_(cvs_col),
    schema_row_(schema_row),
    schema_col_(schema_col),
    shift_(0.0)
{
  mesh_ = cvs_col_->Mesh();
  rhs_ = Teuchos::rcp(new CompositeVector(*cvs_row_, true));

  ncells_owned =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  nfaces_owned =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nnodes_owned =
    mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);

  ncells_wghost =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_wghost =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  nnodes_wghost =
    mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  if (mesh_->valid_edges()) {
    nedges_owned =
      mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::OWNED);
    nedges_wghost =
      mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::ALL);
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
 * Init owned local operators.
 ****************************************************************** */
void
Operator::Init()
{
  rhs_->putScalarMasterAndGhosted(0.0);
  int nops = ops_.size();
  for (int i = 0; i < nops; ++i) {
    if (!(ops_properties_[i] & OPERATOR_PROPERTY_DATA_READ_ONLY))
      ops_[i]->Init();
  }
}


/* ******************************************************************
 * Create structure of a global square matrix.
 ****************************************************************** */
void
Operator::SymbolicAssembleMatrix()
{
  // Create the supermap given a space (set of possible schemas) and a
  // specific schema (assumed/checked to be consistent with the space).
  smap_ = CreateSuperMap(*cvs_col_, schema(), 1);

  // create the graph
  int row_size = MaxRowSize(*mesh_, schema(), 1);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(
    smap_->getMap(), smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  // fill the graph
  SymbolicAssembleMatrix(*smap_, *graph, 0, 0);

  // Completing and optimizing the graphs
  int ierr = graph->FillComplete(smap_->getMap(), smap_->getMap());
  AMANZI_ASSERT(!ierr);

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


/* ******************************************************************
 * Create structure of a global matrix.
 ****************************************************************** */
void
Operator::SymbolicAssembleMatrix(const SuperMap& map, GraphFE& graph,
                                 int my_block_row, int my_block_col) const
{
  // first of double dispatch via Visitor pattern
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->SymbolicAssembleMatrixOp(
      this, map, graph, my_block_row, my_block_col);
  }
}


/* ******************************************************************
 * Populate matrix entries.
 ****************************************************************** */
void
Operator::AssembleMatrix()
{
  if (Amat_ == Teuchos::null) {
    Errors::Message msg("Symbolic assembling was not performed.");
    Exceptions::amanzi_throw(msg);
  }

  Amat_->Zero();
  AssembleMatrix(*smap_, *Amat_, 0, 0);
  Amat_->FillComplete();

  if (shift_ != 0.0) { Amat_->DiagonalShift(shift_); }

  // std::stringstream filename_s2;
  // filename_s2 << "assembled_matrix" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *Amat_
  // ->Matrix());
  // exit(0);
}


/* ******************************************************************
 * Populates matrix entries.
 ****************************************************************** */
void
Operator::AssembleMatrix(const SuperMap& map, MatrixFE& matrix,
                         int my_block_row, int my_block_col) const
{
  // first of double dispatch via Visitor pattern
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->AssembleMatrixOp(this, map, matrix, my_block_row, my_block_col);
  }
}


/* ******************************************************************
 * Linear algebra operations with matrices: r = f - A * u.
 ****************************************************************** */
int
Operator::ComputeResidual(const CompositeVector& u, CompositeVector& r,
                          bool zero)
{
  int ierr;
  if (zero) {
    ierr = Apply(u, r);
  } else {
    ierr = Apply(u, r, -1.0);
  }
  r.update(1.0, *rhs_, -1.0);
  return ierr;
}


/* ******************************************************************
 * Linear algebra operations with matrices: r = A * u - f.
 ****************************************************************** */
int
Operator::ComputeNegativeResidual(const CompositeVector& u, CompositeVector& r,
                                  bool zero)
{
  int ierr;
  if (zero) {
    ierr = Apply(u, r);
  } else {
    ierr = Apply(u, r, 1.0);
  }

  r.update(-1.0, *rhs_, 1.0);

  return ierr;
}


/* ******************************************************************
 * Parallel matvec product Y = A * X.
 ******************************************************************* */
int
Operator::Apply(const CompositeVector& X, CompositeVector& Y,
                double scalar) const
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

  apply_calls_++;

  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it)
    (*it)->ApplyMatrixFreeOp(this, X, Y);

  return 0;
}


/* ******************************************************************
 * Parallel matvec product Y = A^T * X.
 ******************************************************************* */
int
Operator::ApplyTranspose(const CompositeVector& X, CompositeVector& Y,
                         double scalar) const
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

  apply_calls_++;

  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->ApplyTransposeMatrixFreeOp(this, X, Y);
  }

  return 0;
}


/* ******************************************************************
 * Parallel matvec product Y = A * X.
 * This method is mainly for debugging! Matrix-free apply could better.
 ******************************************************************* */
int
Operator::ApplyAssembled(const CompositeVector& X, CompositeVector& Y,
                         double scalar) const
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

  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());

  int ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy);
  ierr |= A_->apply(Xcopy, Ycopy);
  ierr |= AddSuperVectorToCompositeVector(*smap_, Ycopy, Y);

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
int
Operator::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  int ierr(1);

  Epetra_Vector Xcopy(*smap_->getMap());
  Epetra_Vector Ycopy(*smap_->getMap());

  ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy);
  ierr |= preconditioner_->applyInverse(Xcopy, Ycopy);
  ierr |= CopySuperVectorToCompositeVector(*smap_, Ycopy, Y);

  if (ierr) {
    Errors::Message msg("Operator: ApplyInverse failed.\n");
    Exceptions::amanzi_throw(msg);
  }

  return ierr;
}


/* ******************************************************************
 *                       DEPRECATED
 * Initialization of the preconditioner. Note that boundary conditions
 * may be used in re-implementation of this virtual function.
 ****************************************************************** */
void
Operator::InitPreconditioner(const std::string& prec_name,
                             const Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory factory;

  preconditioner_ = factory.Create(prec_name, plist);
  UpdatePreconditioner();
}


/* ******************************************************************
 *                       DEPRECATED
 * Initialization of the preconditioner. Note that boundary conditions
 * may be used in re-implementation of this virtual function.
 ****************************************************************** */
void
Operator::InitPreconditioner(Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(plist);
  preconditioner_->update(A_);
}


/* ******************************************************************
 * Two-stage initialization of preconditioner, part 1.
 * Create the preconditioner and set options. Symbolic assemble of
 * operator's matrix must have been called.
 ****************************************************************** */
void
Operator::InitializePreconditioner(Teuchos::ParameterList& plist)
{
  if (smap_.get() == NULL) {
    if (plist.isParameter("preconditioner type") &&
        plist.get<std::string>("preconditioner type") == "identity") {
      smap_ = CreateSuperMap(*cvs_col_, schema(), 1);
    } else {
      Errors::Message msg("Operator has no super map to be initialized.\n");
      Exceptions::amanzi_throw(msg);
    }
  }

  // provide block ids for block strategies.
  if (plist.isParameter("preconditioner type") &&
      plist.get<std::string>("preconditioner type") == "boomer amg" &&
      plist.isSublist("boomer amg parameters")) {
    // NOTE: Hypre frees this
    auto block_ids = smap_->BlockIndices();

    plist.sublist("boomer amg parameters")
      .set("number of unique block indices", block_ids.first);

    // Note, this passes a raw pointer through a ParameterList.  I was surprised
    // this worked too, but ParameterList is a boost::any at heart... --etc
    plist.sublist("boomer amg parameters")
      .set("block indices", block_ids.second);
  }

  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(plist);
}


/* ******************************************************************
 * Two-stage initialization of preconditioner, part 2.
 * Set the preconditioner structure. Operator's matrix must have been
 * assembled.
 ****************************************************************** */
void
Operator::UpdatePreconditioner()
{
  if (preconditioner_.get() == NULL) {
    Errors::Message msg(
      "Operator has no matrix or preconditioner for update.\n");
    Exceptions::amanzi_throw(msg);
  }
  preconditioner_->update(A_);
}


/* ******************************************************************
 * Update the RHS with this vector.
 * Note that derived classes may reimplement this with a volume term.
 ****************************************************************** */
void
Operator::UpdateRHS(const CompositeVector& source, bool volume_included)
{
  for (auto it = rhs_->begin(); it != rhs_->end(); ++it) {
    if (source.HasComponent(*it)) {
      rhs_->ViewComponent(*it, false)
        ->update(1.0, *source.ViewComponent(*it, false), 1.0);
    }
  }
}


/* ******************************************************************
 * Rescale the local matrices via dispatch.
 ****************************************************************** */
void
Operator::Rescale(double scaling)
{
  for (op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->Rescale(scaling);
  }
}


/* ******************************************************************
 * Rescale the local matrices via dispatch.
 ****************************************************************** */
void
Operator::Rescale(const CompositeVector& scaling)
{
  scaling.ScatterMasterToGhosted();
  for (op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->Rescale(scaling);
  }
}


/* ******************************************************************
 * Rescale the local matrices for particular operator.
 ****************************************************************** */
void
Operator::Rescale(const CompositeVector& scaling, int iops)
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
void
Operator::CreateCheckPoint()
{
  rhs_checkpoint_ = Teuchos::rcp(new CompositeVector(*rhs_));
}


void
Operator::RestoreCheckPoint()
{
  // The routine should be called after checkpoint is created.
  AMANZI_ASSERT(rhs_checkpoint_ != Teuchos::null);

  // restore accumulation and source terms
  *rhs_ = *rhs_checkpoint_;

  // restore local matrices without boundary conditions
  for (op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->RestoreCheckPoint();
  }
}


/* ******************************************************************
 * New implementation of check-point algorithm.
 ****************************************************************** */
int
Operator::CopyShadowToMaster(int iops)
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
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    if ((*it)->Matches(schema_dofs, matching_rule)) { return it; }
  }

  if (action) {
    Errors::Message msg;
    msg << "Operators: Matching rule " << matching_rule << " not found.\n";
    Exceptions::amanzi_throw(msg);
  }

  return OpEnd();
}


/* ******************************************************************
 * Find a matrix block matching the given pattern.
 ****************************************************************** */
Operator::op_iterator
Operator::FindMatrixOp(int schema_dofs, int matching_rule, bool action)
{
  for (op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    if ((*it)->Matches(schema_dofs, matching_rule)) { return it; }
  }

  if (action) {
    Errors::Message msg;
    msg << "Operators: Matching rule " << matching_rule << " not found.\n";
    Exceptions::amanzi_throw(msg);
  }

  return OpEnd();
}


/* ******************************************************************
 * Push back.
 ****************************************************************** */
void
Operator::OpPushBack(const Teuchos::RCP<Op>& block, int properties)
{
  ops_.push_back(block);
  ops_properties_.push_back(properties);
}


/* ******************************************************************
 * Add more operators to the existing list. The added operators have
 * no special properties.
 ****************************************************************** */
void
Operator::OpExtend(op_iterator begin, op_iterator end)
{
  int nops = ops_.size();
  int nnew = nops + std::distance(begin, end);

  ops_.reserve(nnew);
  ops_.insert(ops_.end(), begin, end);
  ops_properties_.resize(nnew, 0);
}


/* ******************************************************************
 * Copies to/from SuperVector for use by Amesos.
 ****************************************************************** */
void
Operator::CopyVectorToSuperVector(const CompositeVector& cv,
                                  Epetra_Vector& sv) const
{
  CopyCompositeVectorToSuperVector(*smap_, cv, sv);
}

void
Operator::CopySuperVectorToVector(const Epetra_Vector& sv,
                                  CompositeVector& cv) const
{
  CopySuperVectorToCompositeVector(*smap_, sv, cv);
}


/* ******************************************************************
 * Generic error message.
 ****************************************************************** */
int
Operator::SchemaMismatch_(const std::string& schema1,
                          const std::string& schema2) const
{
  std::stringstream err;
  err << "Schemas mismatch " << schema1 << " != " << schema2;
  Errors::Message message(err.str());
  Exceptions::amanzi_throw(message);
  return 1;
}


/* ******************************************************************
 * Populates matrix entries.
 ****************************************************************** */
std::string
Operator::PrintDiagnostics() const
{
  std::stringstream msg;
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    msg << "<" << (*it)->schema_string << "> ";
  }
  return msg.str();
}


/* ******************************************************************
 * Visit methods for Apply: Cell.
 ****************************************************************** */
int
Operator::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                            const CompositeVector& X, CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int
Operator::ApplyMatrixFreeOp(const Op_Cell_Face& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int
Operator::ApplyMatrixFreeOp(const Op_Cell_Node& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int
Operator::ApplyMatrixFreeOp(const Op_Cell_Edge& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int
Operator::ApplyMatrixFreeOp(const Op_Cell_Cell& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int
Operator::ApplyMatrixFreeOp(const Op_Cell_Schema& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for Apply: Face
 ****************************************************************** */
int
Operator::ApplyMatrixFreeOp(const Op_Face_Cell& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}

int
Operator::ApplyMatrixFreeOp(const Op_Face_CellBndFace& op,
                            const CompositeVector& X, CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int
Operator::ApplyMatrixFreeOp(const Op_Face_Schema& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for Apply: Edges
 ****************************************************************** */
int
Operator::ApplyMatrixFreeOp(const Op_Edge_Edge& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for Apply: Node
 ****************************************************************** */
int
Operator::ApplyMatrixFreeOp(const Op_Node_Node& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for Apply: SurfaceCell
 ****************************************************************** */
int
Operator::ApplyMatrixFreeOp(const Op_SurfaceCell_SurfaceCell& op,
                            const CompositeVector& X, CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for Apply: SurfaceFace
 ****************************************************************** */
int
Operator::ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
                            const CompositeVector& X, CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for Apply: Coupling
 ****************************************************************** */
int
Operator::ApplyMatrixFreeOp(const Op_Diagonal& op, const CompositeVector& X,
                            CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for ApplyTranspose: Cell
 ****************************************************************** */
int
Operator::ApplyTransposeMatrixFreeOp(const Op_Cell_Schema& op,
                                     const CompositeVector& X,
                                     CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int
Operator::ApplyTransposeMatrixFreeOp(const Op_Face_Schema& op,
                                     const CompositeVector& X,
                                     CompositeVector& Y) const
{
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: Cell.
 ****************************************************************** */
void
Operator::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                   const SuperMap& map, GraphFE& graph,
                                   int my_block_row, int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::SymbolicAssembleMatrixOp(const Op_Cell_Face& op, const SuperMap& map,
                                   GraphFE& graph, int my_block_row,
                                   int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::SymbolicAssembleMatrixOp(const Op_Cell_Node& op, const SuperMap& map,
                                   GraphFE& graph, int my_block_row,
                                   int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::SymbolicAssembleMatrixOp(const Op_Cell_Edge& op, const SuperMap& map,
                                   GraphFE& graph, int my_block_row,
                                   int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::SymbolicAssembleMatrixOp(const Op_Cell_Cell& op, const SuperMap& map,
                                   GraphFE& graph, int my_block_row,
                                   int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::SymbolicAssembleMatrixOp(const Op_Cell_Schema& op,
                                   const SuperMap& map, GraphFE& graph,
                                   int my_block_row, int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: Face.
 ****************************************************************** */
void
Operator::SymbolicAssembleMatrixOp(const Op_Face_Cell& op, const SuperMap& map,
                                   GraphFE& graph, int my_block_row,
                                   int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}

void
Operator::SymbolicAssembleMatrixOp(const Op_Face_CellBndFace& op,
                                   const SuperMap& map, GraphFE& graph,
                                   int my_block_row, int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::SymbolicAssembleMatrixOp(const Op_Face_Schema& op,
                                   const SuperMap& map, GraphFE& graph,
                                   int my_block_row, int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: Edge.
 ****************************************************************** */
void
Operator::SymbolicAssembleMatrixOp(const Op_Edge_Edge& op, const SuperMap& map,
                                   GraphFE& graph, int my_block_row,
                                   int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: Node.
 ****************************************************************** */
void
Operator::SymbolicAssembleMatrixOp(const Op_Node_Node& op, const SuperMap& map,
                                   GraphFE& graph, int my_block_row,
                                   int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: SurfaceCell
 ****************************************************************** */
void
Operator::SymbolicAssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                   const SuperMap& map, GraphFE& graph,
                                   int my_block_row, int my_block_col) const
{
  std::stringstream err;
  err << "Symbolic assemble: invalid schema combination -- " << op.schema_string
      << " cannot be used with a matrix on " << schema_string_;
  Errors::Message message(err.str());
  Exceptions::amanzi_throw(message);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: SurfaceFace.
 ****************************************************************** */
void
Operator::SymbolicAssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                   const SuperMap& map, GraphFE& graph,
                                   int my_block_row, int my_block_col) const
{
  std::stringstream err;
  err << "Symbolic assemble: invalid schema combination -- " << op.schema_string
      << " cannot be used with a matrix on " << schema_string_;
  Errors::Message message(err.str());
  Exceptions::amanzi_throw(message);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: Coupling
 ****************************************************************** */
void
Operator::SymbolicAssembleMatrixOp(const Op_Diagonal& op, const SuperMap& map,
                                   GraphFE& graph, int my_block_row,
                                   int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for assemble: Cell.
 ****************************************************************** */
void
Operator::AssembleMatrixOp(const Op_Cell_FaceCell& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::AssembleMatrixOp(const Op_Cell_Face& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::AssembleMatrixOp(const Op_Cell_Node& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::AssembleMatrixOp(const Op_Cell_Edge& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::AssembleMatrixOp(const Op_Cell_Cell& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::AssembleMatrixOp(const Op_Cell_Schema& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for assemble: Face.
 ****************************************************************** */
void
Operator::AssembleMatrixOp(const Op_Face_Cell& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}

void
Operator::AssembleMatrixOp(const Op_Face_CellBndFace& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


void
Operator::AssembleMatrixOp(const Op_Face_Schema& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for assemble: Edge.
 ****************************************************************** */
void
Operator::AssembleMatrixOp(const Op_Edge_Edge& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for assemble: Node.
 ****************************************************************** */
void
Operator::AssembleMatrixOp(const Op_Node_Node& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Visit methods for assemble: Surface Cell
 ****************************************************************** */
void
Operator::AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                           const SuperMap& map, MatrixFE& mat, int my_block_row,
                           int my_block_col) const
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
void
Operator::AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                           const SuperMap& map, MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  std::stringstream err;
  err << "Assemble matrix: invalid schema combination -- " << op.schema_string
      << " cannot be used with a matrix on " << schema_string_;
  Errors::Message message(err.str());
  Exceptions::amanzi_throw(message);
}


/* ******************************************************************
 * Visit methods for assemble: Coupling
 ****************************************************************** */
void
Operator::AssembleMatrixOp(const Op_Diagonal& op, const SuperMap& map,
                           MatrixFE& mat, int my_block_row,
                           int my_block_col) const
{
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
 * Local assemble routines (for new schema)
 ****************************************************************** */
void
Operator::ExtractVectorCellOp(int c, const Schema& schema,
                              WhetStone::DenseVector& v,
                              const CompositeVector& X) const
{
  Errors::Message msg(
    "Extracton fo local cell-based vector is missing for this operator");
  Exceptions::amanzi_throw(msg);
}


void
Operator::AssembleVectorCellOp(int c, const Schema& schema,
                               const WhetStone::DenseVector& v,
                               CompositeVector& X) const
{
  Errors::Message msg(
    "Extracton fo local cell-based vector is missing for this operator");
  Exceptions::amanzi_throw(msg);
}


void
Operator::ExtractVectorFaceOp(int c, const Schema& schema,
                              WhetStone::DenseVector& v,
                              const CompositeVector& X) const
{
  Errors::Message msg(
    "Extracton fo local cell-based vector is missing for this operator");
  Exceptions::amanzi_throw(msg);
}


void
Operator::AssembleVectorFaceOp(int c, const Schema& schema,
                               const WhetStone::DenseVector& v,
                               CompositeVector& X) const
{
  Errors::Message msg(
    "Extracton fo local cell-based vector is missing for this operator");
  Exceptions::amanzi_throw(msg);
}

} // namespace Operators
} // namespace Amanzi
