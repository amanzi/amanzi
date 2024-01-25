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
#include <fstream>

#include "MatrixMarket_Tpetra.hpp"

// Amanzi
#include "dbc.hh"
#include "DenseVector.hh"
#include "MatrixFE.hh"
#include "InverseFactory.hh"
#include "SuperMap.hh"

// Operators
#include "Op.hh"
#include "Op_Cell_Cell.hh"
//#include "Op_Cell_Edge.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Cell_Face.hh"
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
#include "GraphFE.hh"

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
    shift_(0.0),
    shift_min_(0.0),
    plist_(plist),
    num_colors_(0),
    coloring_(Teuchos::null),
    inverse_pars_set_(false),
    initialize_complete_(false),
    compute_complete_(false),
    assembly_complete_(false)
{
  mesh_ = cvs_col_->getMesh();
  rhs_ = Teuchos::rcp(new CompositeVector(cvs_row_));

  ncells_owned = mesh_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  nfaces_owned = mesh_->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
  nnodes_owned = mesh_->getNumEntities(AmanziMesh::NODE, AmanziMesh::Parallel_kind::OWNED);

  ncells_wghost = mesh_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::ALL);
  nfaces_wghost = mesh_->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);
  nnodes_wghost = mesh_->getNumEntities(AmanziMesh::NODE, AmanziMesh::Parallel_kind::ALL);

  if (mesh_->hasEdges()) {
    nedges_owned = mesh_->getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::OWNED);
    nedges_wghost = mesh_->getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::ALL);
  } else {
    nedges_owned = 0;
    nedges_wghost = 0;
  }

  Teuchos::ParameterList vo_list = plist.sublist("Verbose Object");
  vo_ = Teuchos::rcp(new VerboseObject("Operators", vo_list));

  shift_ = plist.get<double>("diagonal shift", 0.0);
  shift_min_ = plist.get<double>("diagonal shift minimum", 0.0);
  AMANZI_ASSERT(shift_min_ <= 0.0); // not yet implemented

  apply_calls_ = 0;

  if (plist_.isSublist("inverse")) {
    auto& inv_list = plist_.sublist("inverse");
    AmanziSolvers::setMakeOneIterationCriteria(inv_list);
    set_inverse_parameters(inv_list);
  }
}



/* ******************************************************************
* Create structure of a global square matrix.
****************************************************************** */
void Operator::SymbolicAssembleMatrix()
{
  // Create the supermap given a space (set of possible schemas) and a
  // specific schema (assumed/checked to be consistent with the space).
  if (!smap_.get()) smap_ = createSuperMap(cvs_col_.ptr());

  // create the graph
  int row_size = MaxRowSize(*mesh_, schema(), 1);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->getMap(),
     smap_->getGhostedMap(), smap_->getGhostedMap(), row_size));

  // fill the graph
  SymbolicAssembleMatrix(*smap_, *graph, 0, 0);

  // Completing and optimizing the graphs
  graph->fillComplete(smap_->getMap(), smap_->getMap());

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  Amat_->fillComplete(); // some inverse methods seem to need this fillComplete, even with no data.
  A_ = Amat_->getMatrix();
  assembly_complete_ = false;
}


/* ******************************************************************
* Create structure of a global matrix.
****************************************************************** */
void Operator::SymbolicAssembleMatrix(const SuperMap& map, GraphFE& graph,
                                      int my_block_row, int my_block_col) const
{
  // first of double dispatch via Visitor pattern
  for (auto& it : *this) {
    it->SymbolicAssembleMatrixOp(this, map, graph, my_block_row, my_block_col);
  }
}


// /* ******************************************************************
// * Default visit method for symbolic assemble: Coupling
// ****************************************************************** */
// void Operator::SymbolicAssembleMatrixOp(const Op_Diagonal& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row, int my_block_col) const
// {
//   AMANZI_ASSERT(cvs_col_->hasComponent(op.col_compname()));
//   AMANZI_ASSERT(cvs_row_->hasComponent(op.row_compname()));

//   const std::vector<int>& row_gids = map.viewGhostIndices(my_block_row, op.row_compname(), 0);
//   const std::vector<int>& col_gids = map.viewGhostIndices(my_block_col, op.col_compname(), 0);

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
//     ierr |= graph.insertLocalIndices(ndofs, lid_r.data(), ndofs, lid_c.data());
//   }
//   AMANZI_ASSERT(!ierr);
// }


/* ******************************************************************
* Populate matrix entries.
****************************************************************** */
void Operator::AssembleMatrix()
{
  if (Amat_ == Teuchos::null) SymbolicAssembleMatrix();

  // note, this is called prior to AssembleMatrix() because Schur complements
  // override this because ApplyAssembled is not valid for them.
  assembly_complete_ = true;

  Amat_->resumeFill();
  Amat_->zero();
  AssembleMatrix(*smap_, *Amat_, 0, 0);
  Amat_->fillComplete();

  if (shift_ != 0.0) {
    Amat_->diagonalShift(shift_);
  // } else if (shift_min_ != 0.0) {
  //   Amat_->diagonalShiftMin(shift_min_);
  }

  compute_complete_ = false;

  // WriteMatrix("assembled_matrix");
  // throw("assembed matrix written");
}

/* ******************************************************************
* Populates matrix entries.
****************************************************************** */
void Operator::AssembleMatrix(const SuperMap& map, MatrixFE& matrix,
                              int my_block_row, int my_block_col) const
{
  // first of double dispatch via Visitor pattern
  for (auto& it : *this) {
    it->AssembleMatrixOp(this, map, matrix, my_block_row, my_block_col);
  }
}


/* ******************************************************************
* Default visit methods for assemble: Coupling
****************************************************************** */
// void Operator::AssembleMatrixOp(const Op_Diagonal& op,
//                                 const SuperMap& map, MatrixFE& mat,
//                                 int my_block_row, int my_block_col) const
// {
//   AMANZI_ASSERT(cvs_col_->hasComponent(op.col_compname()));
//   AMANZI_ASSERT(cvs_row_->hasComponent(op.row_compname()));

//   const std::vector<int>& row_gids = map.viewGhostIndices(my_block_row, op.row_compname(), 0);
//   const std::vector<int>& col_gids = map.viewGhostIndices(my_block_col, op.col_compname(), 0);

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

//     ierr |= mat.sumIntoLocalValues(lid_r.data(), lid_c.data(), op.matrices[n]);
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
  if (assembly_complete_) {
    if (shift_ != 0.0) {
      // shift_ != 0 implies a diagonal shift is being applied.  That diagonal
      // shift is done at assemble-time.  This is used to ensure that the
      // matrix is non-singular, which often happens with zero rows.  When it
      // is used to make the operator non-singular, the intent is to allow it
      // to be invertible, but not to change the meaning of the forward
      // operator.  Therefore we have to unshift the assembled matrix prior to
      // doing the forward apply.
      Amat_->diagonalShift(-shift_);
      int ierr = applyAssembled(X, Y, scalar);
      Amat_->diagonalShift(shift_);
      return ierr;
    } else if (shift_min_ != 0.0) {
      // A shift_min_ is lossy -- we cannot "unshift" it because we cannot know
      // what the original diagonal value was.  Therefore must use the
      // unassembled forward apply (or re-assemble without the shift_min).
      return applyUnassembled(X, Y, scalar);
    } else {
      return applyAssembled(X, Y, scalar);
    }
  } else {
    return applyUnassembled(X, Y, scalar);
  }
}

int
Operator::applyUnassembled(const CompositeVector& X, CompositeVector& Y, double scalar) const
{
  X.scatterMasterToGhosted();

  if (scalar == 0.0) {
    Y.putScalarMasterAndGhosted(0.0);
  } else if (scalar == 1.0) {
    Y.putScalarGhosted(0.0);
  } else {
    Y.scale(scalar);
    Y.putScalarGhosted(0.0);
  }

  apply_calls_++;
  for (auto& it : *this) it->ApplyMatrixFreeOp(this, X, Y);
  Y.gatherGhostedToMaster();
  return 0;
}


/* ******************************************************************
* Parallel matvec product Y = A * X.
* This method is mainly for debugging! Matrix-free apply could better.
******************************************************************* */
int Operator::applyAssembled(const CompositeVector& X, CompositeVector& Y, double scalar) const
{
  X.scatterMasterToGhosted();

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
  if (!compute_complete_) const_cast<Operator*>(this)->computeInverse();
  if (preconditioner_.get() == nullptr) {
    Errors::Message msg("Operator did not initialize a preconditioner.\n");
    Exceptions::amanzi_throw(msg);
  }
  return preconditioner_->applyInverse(X, Y);
}


/* ******************************************************************
* Default visit method for apply 
****************************************************************** */
// int Operator::ApplyMatrixFreeOp(const Op_Diagonal& op,
//                                 const CompositeVector& X, CompositeVector& Y) const
// {
//   AMANZI_ASSERT(cvs_col_->hasComponent(op.col_compname()));
//   AMANZI_ASSERT(cvs_row_->hasComponent(op.row_compname()));

//   const Epetra_MultiVector& Xf = *X.viewComponent(op.col_compname(), true);
//   Epetra_MultiVector& Yf = *Y.viewComponent(op.row_compname(), true);
  
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
  compute_complete_ = false;
  assembly_complete_ = false;
}

/* ******************************************************************
* Copy the diagonal.
****************************************************************** */
void Operator::getLocalDiagCopy(CompositeVector& X) const
{
  X.putScalarMasterAndGhosted(0.);
  for (auto& op : *this) {
    op->SumLocalDiag(X);
  }
  X.gatherGhostedToMaster();
}


void Operator::set_inverse_parameters(const std::string& prec_name,
        const Teuchos::ParameterList& plist) {
  Teuchos::ParameterList inner_plist(plist.sublist(prec_name));
  set_inverse_parameters(inner_plist);
}


void Operator::set_inverse_parameters(const std::string& prec_name,
        const Teuchos::ParameterList& prec_list,
        const std::string& iter_name,
        const Teuchos::ParameterList& iter_list,
        bool make_one_iteration) {
  auto inv_plist = AmanziSolvers::mergePreconditionerSolverLists(
      prec_name, prec_list, iter_name, iter_list, make_one_iteration);
  inv_plist.setName(iter_name);
  set_inverse_parameters(inv_plist);
}


/* ******************************************************************
* Two-stage initialization of preconditioner, part 1.
* Create the preconditioner and set options. Symbolic assemble of
* operator's matrix must have been called.
****************************************************************** */
void Operator::set_inverse_parameters(Teuchos::ParameterList& plist)
{
  // delay pc construction until we know we have structure and can create the
  // coloring.
  inv_plist_ = plist;
  inverse_pars_set_ = true;
  initialize_complete_ = false;
  compute_complete_ = false;
}


/* ******************************************************************
* Two-stage initialization of preconditioner, part 1.
* Create the preconditioner and set options. Symbolic assemble of 
* operator's matrix must have been called.
****************************************************************** */
void Operator::initializeInverse()
{
  if (!inverse_pars_set_) {
    Errors::Message msg("No inverse parameter list.  Provide a sublist \"inverse\" or ensure set_inverse_parameters() is called.");
    msg << " In: " << typeid(*this).name() << "\n";
    Exceptions::amanzi_throw(msg);
  }

  smap_ = createSuperMap(getDomainMap().ptr());

  // provide block ids for block strategies.
  if (inv_plist_.isParameter("preconditioning method") &&
      inv_plist_.get<std::string>("preconditioning method") == "boomer amg" &&
      inv_plist_.isSublist("boomer amg parameters") &&
      inv_plist_.sublist("boomer amg parameters").get<bool>("use block indices", false)) {
    assert(false && "Not implemented");
    //if (coloring_ == Teuchos::null || num_colors_ == 0) {
    //  auto block_ids = smap_->BlockIndices();
    //  set_coloring(block_ids.first, block_ids.second);
    //}
    inv_plist_.sublist("boomer amg parameters").set("number of unique block indices", num_colors_);
    inv_plist_.sublist("boomer amg parameters").set("block indices", coloring_);
  }
  preconditioner_ = AmanziSolvers::createInverse(inv_plist_, Teuchos::rcpFromRef(*this));
  preconditioner_->initializeInverse(); // NOTE: calls this->SymbolicAssembleMatrix()
  initialize_complete_ = true;
  compute_complete_ = false;
}

void Operator::computeInverse()
{
  if (!initialize_complete_) initializeInverse();
  // assembly must be possible now
  AMANZI_ASSERT(preconditioner_.get());
  preconditioner_->computeInverse(); // NOTE: calls this->AssembleMatrix()
  compute_complete_ = true;
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
    if (source.hasComponent(it)) {
      rhs_->getComponent(it, false)->update(1.0, *source.getComponent(it, false), 1.0);
    }
  }
}


/* ******************************************************************
* Rescale the local matrices via dispatch.
****************************************************************** */
void Operator::Rescale(double scaling)
{
  for (auto& it : *this) it->Rescale(scaling);
  assembly_complete_ = false;
  compute_complete_ = false;
}


/* ******************************************************************
* Rescale the local matrices via dispatch.
****************************************************************** */
void Operator::Rescale(const CompositeVector& scaling)
{
  scaling.scatterMasterToGhosted();
  for (auto& it : *this) it->Rescale(scaling);
  assembly_complete_ = false;
  compute_complete_ = false;
}


/* ******************************************************************
* Rescale the local matrices for particular operator.
****************************************************************** */
void Operator::Rescale(const CompositeVector& scaling, int iops)
{
  AMANZI_ASSERT(iops < ops_.size());
  scaling.scatterMasterToGhosted();
  ops_[iops]->Rescale(scaling);
  assembly_complete_ = false;
  compute_complete_ = false;
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
int Operator::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


int Operator::ApplyMatrixFreeOp(const Op_Cell_Face& op,
                                const CompositeVector& X, CompositeVector& Y) const {
  return SchemaMismatch_(op.schema_string, schema_string_);
}


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
                                    const WhetStone::DenseVector<>& v, CompositeVector& X) const {
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


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Operator> Operator::clone() const {
  Errors::Message msg("Cloning of a derived Operator class is missing");
  Exceptions::amanzi_throw(msg);
  return Teuchos::null;
}

void Operator::WriteMatrix(const std::string& fname_base) const
{
  // write to disk
  std::string filename = fname_base + ".dat";
  Tpetra::MatrixMarket::Writer<Matrix_type> writer;
  writer.writeSparseFile(filename.c_str(), A_, "A", "A", true);

  const auto smap = smap_->getMap();
  for (int r=0; r!=mesh_->getComm()->getSize(); ++r) {
    if (r == mesh_->getComm()->getRank()) {
      std::ofstream os(fname_base+"_map.dat", std::ofstream::app);
      if (r == 0)
        os << "# Rank EType LID GID SMAP_LID SMAP_GID x y z" << std::endl;

      std::vector<std::string> compnames = {"cell", "face", "boundary_face", "edge", "node"};
      for (const auto& comp : compnames) {
        auto mesh_map = mesh_->getMap(AmanziMesh::createEntityKind(comp), false);
        int lcv_block = 0;
        bool not_done_block = smap_->hasComponent(lcv_block, comp, 0);
        while(not_done_block) {
          int lcv_dof = 0;
          while(smap_->hasComponent(lcv_block, comp, lcv_dof)) {
            auto ids = smap_->viewIndices<Amanzi::MirrorHost>(lcv_block, comp, lcv_dof);
            for (int i=0; i!=ids.extent(0); ++i) {
              os << r << " "
                 << comp << " "
                 << i << " "
                 << mesh_map->getGlobalElement(i) << " "
                 << ids(i) << " "
                 << smap->getGlobalElement(ids(i)) << " "
                 << mesh_->getCentroid(AmanziMesh::createEntityKind(comp), i) << std::endl;
            }
            lcv_dof++;
          }

          lcv_block++;
          not_done_block = smap_->hasComponent(lcv_block, comp, lcv_dof);
        }
      }
      os.close();
    }
    mesh_->getComm()->barrier();
  }
}


}  // namespace Operators
}  // namespace Amanzi

