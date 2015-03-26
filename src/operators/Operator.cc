/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
  Ethan Coon (ecoon@lanl.gov)
*/

#include "EpetraExt_RowMatrixOut.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_FECrsGraph.h"

#include "DenseVector.hh"
#include "PreconditionerFactory.hh"

#include "SuperMap.hh"
#include "MatrixFE.hh"

#include "Op.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Cell_Face.hh"
#include "Op_Cell_Cell.hh"
#include "Op_Cell_Node.hh"
#include "Op_Cell_Edge.hh"
#include "Op_Face_Cell.hh"
#include "Op_Node_Node.hh"

#include "OperatorDefs.hh"
#include "OperatorUtils.hh"

#include "Operator.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Default constructor.
 ****************************************************************** */
Operator::Operator(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                   Teuchos::ParameterList& plist,
                   int schema) :
    cvs_(cvs),
    schema_(schema),
    data_validity_(true),
    symbolic_assembled_(false),
    assembled_(false)
{
  mesh_ = cvs_->Mesh();
  rhs_ = Teuchos::rcp(new CompositeVector(*cvs_, true));

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  Teuchos::ParameterList vo_list = plist.sublist("Verbose Object");
  vo_ = Teuchos::rcp(new VerboseObject("Operators", vo_list));
}


/* ******************************************************************
* Zeros everything out
****************************************************************** */
void Operator::Init()
{
  rhs_->PutScalarMasterAndGhosted(0.0);
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->Init();
  }
}


/* ******************************************************************
* Create structure of a global matrix.
****************************************************************** */
void Operator::SymbolicAssembleMatrix()
{
  // Create the supermap given a space (set of possible schemas) and a
  // specific schema (assumed/checked to be consistent with the sapce).
  smap_ = createSuperMap(*cvs_, schema(), 1);

  // create the graph
  int row_size = MaxRowSize(*mesh_, schema(), 1);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->Map(),
          smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  // fill the graph
  SymbolicAssembleMatrix(*smap_, *graph, 0, 0);

  // Completing and optimizing the graphs
  int ierr = graph->FillComplete(smap_->Map(), smap_->Map());
  ASSERT(!ierr);

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


/* ******************************************************************
* Create structure of a global matrix.
****************************************************************** */
void Operator::SymbolicAssembleMatrix(const SuperMap& map, GraphFE& graph,
                                      int my_block_row, int my_block_col) const
{
  // first of double dispatch via Visitor pattern
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->SymbolicAssembleMatrixOp(this, map, graph, my_block_row, my_block_col);
  }
}


/* ******************************************************************
* Populate matrix entries.
****************************************************************** */
void Operator::AssembleMatrix()
{
  if (Amat_ == Teuchos::null) {
    Errors::Message msg("Symbolic assembling was not performed.");
    Exceptions::amanzi_throw(msg);
  }

  Amat_->Zero();
  AssembleMatrix(*smap_, *Amat_, 0, 0);
  Amat_->FillComplete();
}


/* ******************************************************************
* Populates matrix entries.
****************************************************************** */
void Operator::AssembleMatrix(const SuperMap& map, MatrixFE& matrix,
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
int Operator::ComputeResidual(const CompositeVector& u, CompositeVector& r, bool zero)
{
  int ierr;
  if (zero) {
    ierr = Apply(u, r);
  } else {
    ierr = Apply(u, r, -1.0);
  }
  r.Update(1.0, *rhs_, -1.0);
  return ierr;
}


/* ******************************************************************
* Linear algebra operations with matrices: r = A * u - f.
****************************************************************** */
int Operator::ComputeNegativeResidual(const CompositeVector& u, CompositeVector& r, bool zero)
{
  int ierr;
  if (zero) {
    ierr = Apply(u, r);
  } else {
    ierr = Apply(u, r, 1.0);
  }    
  r.Update(-1.0, *rhs_, 1.0);
  return ierr;
}


/* ******************************************************************
* Parallel matvec product Y = A * X.
******************************************************************* */
int
Operator::Apply(const CompositeVector& X, CompositeVector& Y, double scalar) const
{
  X.ScatterMasterToGhosted();

  // initialize ghost elements
  if (scalar == 0.0) {
    Y.PutScalarMasterAndGhosted(0.0);
  } else if (scalar == 1.0) {
    Y.PutScalarGhosted(0.0);
  } else {
    Y.Scale(scalar);
    Y.PutScalarGhosted(0.0);
  }

  int ierr(0);
  // THIS NEEDS TESTING: which should be the preferred execution pathway?
  // Apply via assembled matrix or via local matrices (assuming the assembled
  // matrix is available).

  // if (assembled_) {
  //   Epetra_Vector Xcopy(A_->RowMap());
  //   Epetra_Vector Ycopy(A_->RowMap());
  //   int ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy, 0);
  //   ierr |= A_->Apply(Xcopy, Ycopy);
  //   ierr |= AddSuperVectorToCompositeVector(*smap_, Ycopy, Y, 0);
  //   ASSERT(!ierr);
  // } else {
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->ApplyMatrixFreeOp(this, X, Y);
  }
  // }

  return ierr;
}


/* ******************************************************************
* Parallel matvec product Y = A * X.
******************************************************************* */
int Operator::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  // Y = X;
  // return 0;
  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());
  int ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy, 0);

  // dump the schur complement
  // std::stringstream filename_s2;
  // filename_s2 << "schur_PC_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *A_);

  ierr |= preconditioner_->ApplyInverse(Xcopy, Ycopy);
  ierr |= CopySuperVectorToCompositeVector(*smap_, Ycopy, Y, 0);
  ASSERT(!ierr);

  return ierr;
}


/* ******************************************************************
* Initialization of the preconditioner. Note that boundary conditions
* may be used in re-implementation of this virtual function.
****************************************************************** */
void Operator::InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(prec_name, plist);
  preconditioner_->Update(A_);
}


/* ******************************************************************
* Initialization of the preconditioner. Note that boundary conditions
* may be used in re-implementation of this virtual function.
****************************************************************** */
void Operator::InitPreconditioner(Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(plist);
  preconditioner_->Update(A_);
}


/* ******************************************************************
* Update the RHS with this vector.
* Note that derived classes may reimplement this with a volume term.
****************************************************************** */
void Operator::UpdateRHS(const CompositeVector& source, bool volume_included) {
  for (CompositeVector::name_iterator it = rhs_->begin();
       it != rhs_->end(); ++it) {
    if (source.HasComponent(*it)) {
      rhs_->ViewComponent(*it, false)->Update(1.0, *source.ViewComponent(*it, false), 1.0);
    }
  }
}


/* ******************************************************************
* Rescale the local matrices.
****************************************************************** */
void
Operator::Rescale(const CompositeVector& scaling) {
  // Dispatch Rescaling to the Ops.
  scaling.ScatterMasterToGhosted();
  for (op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->Rescale(scaling);
  }
}


// /* ******************************************************************
//  * Enforce the BCs on local matrices
//  ****************************************************************** */
// void
// Operator::ApplyBCs(const Teuchos::RCP<BCs>& bc) {
//   bc_ = bc;

//   // Dispatch BCs to the Ops.  Note that since we allow multiple Ops, the BC
//   // can be applied in one of several Ops, but all others must, for example,
//   // zero their appropriate row.  The flag here allows Operators to only set
//   // an identity row in one Op while zeroing all others.
//   bool bc_sucessfully_applied = false;
//   for (op_iterator it = OpBegin(); it != OpEnd(); ++it) {
//     bc_sucessfully_applied |=
//         (*it)->ApplyBC(*bc, rhs_.ptr(), bc_sucessfully_applied);
//   }
//   // This check is currently broken by adv + diffusion when diffusion is MFD
//   // -- need to reconsider whether this is useful or can be reformulated to be
//   // useful.  
//   // if (!bc_sucessfully_applied) {
//   //   Errors::Message msg("Operators: ApplyBC not sucessful for this bc and operator schema combinations.");
//   //   Exceptions::amanzi_throw(msg);
//   // }    
// }


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
  ASSERT(rhs_checkpoint_ != Teuchos::null);

  // restore accumulation and source terms
  *rhs_ = *rhs_checkpoint_;

  // restore local matrices without boundary conditions
  for (op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    (*it)->RestoreCheckPoint();
  }
}


/* ******************************************************************
* Find a matrix block matching the given pattern.
****************************************************************** */
Operator::const_op_iterator
Operator::FindMatrixOp(int schema_dofs, int matching_rule, bool action) const
{
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    if ((*it)->Matches(schema_dofs, matching_rule)) {
      return it;
    }
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
    if ((*it)->Matches(schema_dofs, matching_rule)) {
      return it;
    }
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
void Operator::OpPushBack(const Teuchos::RCP<Op>& block) {
  ops_.push_back(block);
}


/* ******************************************************************
* Extension.
****************************************************************** */
void Operator::OpExtend(op_iterator begin, op_iterator end)
{
  ops_.reserve(ops_.size() + std::distance(begin, end));
  ops_.insert(ops_.end(), begin, end);
}


/* ******************************************************************
* Generic error message.
****************************************************************** */
int Operator::SchemaMismatch_(const std::string& schema1, const std::string& schema2) const
{
  std::stringstream err;
  err << "Invalid schema combination -- " << schema1
      << " cannot be used with a matrix on " << schema2;
  Errors::Message message(err.str());
  Exceptions::amanzi_throw(message);
  return 1;
}


/* ******************************************************************
* Visit methods for Apply: Cell.
****************************************************************** */
int Operator::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int Operator::ApplyMatrixFreeOp(const Op_Cell_Face& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int Operator::ApplyMatrixFreeOp(const Op_Cell_Node& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int Operator::ApplyMatrixFreeOp(const Op_Cell_Edge& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int Operator::ApplyMatrixFreeOp(const Op_Cell_Cell& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for Apply: Face
****************************************************************** */
int Operator::ApplyMatrixFreeOp(const Op_Face_Cell& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for Apply: Node
****************************************************************** */
int Operator::ApplyMatrixFreeOp(const Op_Node_Node& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for RHS: Cell
****************************************************************** */
void Operator::AssembleRHSOp(const Op_Cell_FaceCell& op, CompositeVector& rhs) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::AssembleRHSOp(const Op_Cell_Face& op, CompositeVector& rhs) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::AssembleRHSOp(const Op_Cell_Node& op, CompositeVector& rhs) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::AssembleRHSOp(const Op_Cell_Cell& op, CompositeVector& rhs) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for RHS: Face
****************************************************************** */
void Operator::AssembleRHSOp(const Op_Face_Cell& op, CompositeVector& rhs) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for RHS: Node
****************************************************************** */
void Operator::AssembleRHSOp(const Op_Node_Node& op, CompositeVector& rhs) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Cell.
****************************************************************** */
void Operator::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::SymbolicAssembleMatrixOp(const Op_Cell_Node& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::SymbolicAssembleMatrixOp(const Op_Cell_Edge& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::SymbolicAssembleMatrixOp(const Op_Cell_Cell& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Face.
****************************************************************** */
void Operator::SymbolicAssembleMatrixOp(const Op_Face_Cell& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Node.
****************************************************************** */
void Operator::SymbolicAssembleMatrixOp(const Op_Node_Node& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for assemble: Cell.
****************************************************************** */
void Operator::AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::AssembleMatrixOp(const Op_Cell_Face& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::AssembleMatrixOp(const Op_Cell_Node& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::AssembleMatrixOp(const Op_Cell_Edge& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


void Operator::AssembleMatrixOp(const Op_Cell_Cell& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for assemble: Face.
****************************************************************** */
void Operator::AssembleMatrixOp(const Op_Face_Cell& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}


/* ******************************************************************
* Visit methods for assemble: Node.
****************************************************************** */
void Operator::AssembleMatrixOp(const Op_Node_Node& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const {
  SchemaMismatch_(op.schema_string, schema_string_);
}

}  // namespace Operators
}  // namespace Amanzi

