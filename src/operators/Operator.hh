/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

//! Operator represents a linear map, and typically encapsulates a discretization.
/*!

Operators are discrete forms of linearized PDEs operators.
They form a layer between physical process kernels and solvers
and include accumulation, diffusion, advection, elasticity, reaction,
and source operators.
The residual associated with an operator :math:`L_h` helps to
understand the employed sign convention:

.. math::
  r = f - L_h u.

Operator represents a map from linear space X to linear space Y.  Typically,
this map is a linear map, and encapsulates much of the discretization involved
in moving from continuous to discrete equations. The spaces X and Y are described
by CompositeVectors (CV). A few maps X->Y are supported.

An operator provides an interface for applying both the forward and inverse
linear map (assuming the map is invertible).

Typically the Operator is never seen by the user; instead the user provides
input information for helper classes based on the continuous mathematical
operator and the desired discretization.  These helpers build the needed
``Operator``, which may include information from multiple helpers (i.e. in the
case of Jacobian Operators for a PDE).

However, one option may be provided by the user, which is related to dealing
with nearly singular operators:

.. admonition:: operators-spec

  * `"diagonal shift`" ``[double]`` **0.0** Adds a scalar shift to the diagonal
    of the ``Operator``, which can be useful if the ``Operator`` is singular or
    near-singular.

A PK decides how to bundle operators in a collection of operators.
For example, an advection-diffusion problem may benefit from using
a single operator that combines two operators representing diffusion and advection process.
Collection of operators must be used for implicit solvers and for building preconditioners.
In such a case, the collections acts as a single operator.

Operators use a few tools that are generic in nature and can be used independently by PKs.
The list includes reconstruction and limiting algorithms.


Schema
------

The operators use notion of schema to describe operator's abstract structure.
Old operators use a simple schema which is simply the list of geometric objects where
scalar degrees of freedom are defined.
New operators use a list to define location, type, and number of degrees of freedom.
In addition, the base of local stencil is either *face* or *cell*.
A rectangular operator needs two schemas do describe its domain (called `"schema domain`")
and its range (called `"schema range`").
A square operator may use either two identical schema lists or a single list called `"schema`".

.. code-block:: xml

  <ParameterList name="pks operator name">  <!-- parent list-->
  <ParameterList name="schema domain">
    <Parameter name="base" type="string" value="cell"/>
    <Parameter name="location" type="Array(string)" value="{node, face}"/>
    <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
    <Parameter name="number" type="Array(int)" value="{2, 1}"/>
  </ParameterList>
  <ParameterList name="schema domain">
    <Parameter name="base" type="string" value="cell"/>
    <Parameter name="location" type="Array(string)" value="{node, face}"/>
    <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
    <Parameter name="number" type="Array(int)" value="{2, 1}"/>
  </ParameterList>
  </ParameterList>

This example describes a square operator with two degrees of freedom per mesh node and one
degree of freedom per mesh face.
The face-based degree of freedom is the normal component of a vector field.
Such set of degrees of freedom is used in the Bernardi-Raugel element for discretizing
Stokes equations.
Parameter `"base`" indicates that local matrices are associated with mesh cells.

*/


/*
Developer and implementation notes:

An ``Operator`` consists of, potentially, both a single, global, assembled
matrix AND an un-ordered, additive collection of lower-rank (or equal) local
operators, here called ``Op``s. During its construction, an operator can grow
by assimilating more ``Op``s. An ``Operator`` knows how to Apply and Assemble
all of its local ``Op``s.

Typically the forward operator is applied using only local ``Op``s -- the
inverse operator typically requires assembling a matrix, which may represent
the entire operator or may be a Schur complement.

In all ``Operator``s, ``Op``s, and matrices, a key concept is the schema.
The retired schema includes, at the least, an enum representing the Degrees
of Fredom (DoFs) associated with the Operator/matrix's domain (and implied
equivalent range, X=Y).  A schema may also, in the case of ``Op``s, include
information on the base entity on which the local matrix lives.
The new schema specifies dofs for the operator domain and range. It also
includes location of the geometric base.

For instance, the global matrix schema may be CELL, indicating that the domain
and range is defined as the space of cells on the mesh provided at
construction.  For example, a local Op might be defined with a scheme of FACE
base entities and CELL DoFs, indicating that there is one local operator on
each face, and it consists of small dense matrices with a domain and range of
all cells connected to that face (two for interior faces, one for boundary
faces).

Notes for application code developers:

* We note that typically the Op object is not directly constructed; instead
  helper classes based on the discretization and continuous operator are used.
  For instance, OperatorDiffusion is a class that constructs local Ops that
  can be assembled into the global Operator to form a diffusion matrix.

* For simple operations (i.e. diffusion), Operators are not constructed
  directly.  Instead, a helper class that contains methods for creating and
  populating the Ops within the Operator is created. This helper class can
  create the appropriate Operator itself. More complex operations, i.e.
  advection-diffusion, can be generated by creating an Operator that is the
  union of all dof requirements, and then passing this Operator into the
  helper's constructor. When this is done, the helper simply checks to make
  sure the Operator contains the necessary dofs and adds local Ops to the
  Operator's list of Ops.

Note on implementation for discretization/framework developers:

* Ops work via a visitor pattern.  Assembly (resp Apply, BCs, and
  SymbolicAssembly) are implemented by the (base class) Operator calling a
  dispatch to the (base virtual class) Op, which then dispatches back to the
  (derived class) Operator so that type information of both the Operator
  (i.e. global matrix info) and the Op (i.e. local matrix info) are known.

  Adding new Ops is a fairly straightforward process of providing
  implementations of how to use the mesh to implement the forward action and
  how to assumble the local Op into the global Operator matrix.
*/

#ifndef AMANZI_OPERATOR_HH_
#define AMANZI_OPERATOR_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "Mesh.hh"
#include "CompositeVectorSpace.hh"
#include "CompositeVector.hh"
#include "DenseVector.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"
#include "Matrix.hh"

namespace Amanzi {

class CompositeVector;

namespace Operators {

class BCs;
class SuperMap;
class MatrixFE;
class GraphFE;
class Op;
class Op_Cell_FaceCell;
class Op_Cell_Face;
class Op_Cell_Cell;
class Op_Cell_Node;
class Op_Cell_Edge;
class Op_Cell_Schema;
class Op_Diagonal;
class Op_Face_Face;
class Op_Face_Cell;
class Op_Face_CellBndFace;
class Op_Face_Schema;
class Op_Edge_Edge;
class Op_Node_Node;
class Op_Node_Schema;
class Op_SurfaceCell_SurfaceCell;
class Op_SurfaceFace_SurfaceCell;
class Op_MeshInjection;


class Operator : public Matrix<CompositeVector, CompositeVectorSpace> {
 public:
  using Vector_t = CompositeVector;
  using VectorSpace_t = CompositeVector::VectorSpace_t;

  // constructors
  // At the moment CVS is the domain and range of the operator
  Operator() { apply_calls_ = 0; }

  // deprecated but yet supported
  Operator(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
           Teuchos::ParameterList& plist,
           int schema);

  // general operator (domain may differ from range)
  Operator(const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
           const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
           Teuchos::ParameterList& plist,
           const Schema& schema_row,
           const Schema& schema_col);

  // bijective operator (domain = range)
  Operator(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
           Teuchos::ParameterList& plist,
           const Schema& schema)
    : Operator(cvs, cvs, plist, schema, schema) {};

  virtual ~Operator() = default;

  // a specialized copy constructor is used to extend the operator, e.g.
  // by adding more operators to it.
  virtual Teuchos::RCP<Operator> Clone() const;

  void Init();

  // main members
  // -- virtual methods potentially altered by the schema, Schur complements
  virtual int Apply(const CompositeVector& X, CompositeVector& Y) const override
  {
    return Apply(X, Y, 0.0);
  }
  // -- icomputes Y = A * X + scalar * Y
  virtual int Apply(const CompositeVector& X, CompositeVector& Y, double scalar) const;
  virtual int ApplyAssembled(const CompositeVector& X,
                             CompositeVector& Y,
                             double scalar = 0.0) const;
  virtual int ApplyUnassembled(const CompositeVector& X,
                               CompositeVector& Y,
                               double scalar = 0.0) const;
  virtual int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const override;

  // versions that make it easier to deal with Amanzi input spec format
  void set_inverse_parameters(const std::string& prec_name, const Teuchos::ParameterList& plist);
  void set_inverse_parameters(const std::string& prec_name,
                              const Teuchos::ParameterList& prec_list,
                              const std::string& iter_name,
                              const Teuchos::ParameterList& iter_list,
                              bool make_one_iteration = true);

  // -- preferred methods -- set_parameters, initialize, compute
  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override;
  virtual void InitializeInverse() override;
  virtual void ComputeInverse() override;

  // symbolic assembly:
  // -- wrapper
  virtual void SymbolicAssembleMatrix();
  // -- first dispatch
  virtual void SymbolicAssembleMatrix(const SuperMap& map,
                                      GraphFE& graph,
                                      int my_block_row,
                                      int my_block_col) const;

  // actual assembly:
  // -- wrapper
  virtual void AssembleMatrix();
  // -- first dispatch
  virtual void AssembleMatrix(const SuperMap& map,
                              MatrixFE& matrix,
                              int my_block_row,
                              int my_block_col) const;

  // modifiers
  // -- add a vector to operator's rhs vector
  virtual void UpdateRHS(const CompositeVector& source, bool volume_included = true);
  // -- rescale elemental matrices
  virtual void Rescale(double scaling);
  virtual void Rescale(const CompositeVector& scaling);
  virtual void Rescale(const CompositeVector& scaling, int iops);

  // -- default functionality
  const CompositeVectorSpace& DomainMap() const override { return *get_domain_map(); }
  const CompositeVectorSpace& RangeMap() const override { return *get_range_map(); }
  const Teuchos::RCP<const CompositeVectorSpace>& get_domain_map() const { return cvs_col_; }
  const Teuchos::RCP<const CompositeVectorSpace>& get_range_map() const { return cvs_row_; }
  const Teuchos::RCP<const CompositeVectorSpace>& get_col_map() const { return cvs_col_; }
  const Teuchos::RCP<const CompositeVectorSpace>& get_row_map() const { return cvs_row_; }

  int ComputeResidual(const CompositeVector& u, CompositeVector& r, bool zero = true);
  int ComputeNegativeResidual(const CompositeVector& u, CompositeVector& r, bool zero = true);

  void CreateCheckPoint();
  void RestoreCheckPoint();

  // supporting members
  int CopyShadowToMaster(int iops);

  // access
  virtual std::string name() const override
  {
    return std::string("Operator (") + schema_string_ + ")";
  }
  int schema() const { return schema_col_.OldSchema(); }
  const Schema& schema_col() const { return schema_col_; }
  const Schema& schema_row() const { return schema_row_; }

  void set_schema_string(const std::string& schema_string) { schema_string_ = schema_string; }
  const std::string& get_schema_string() const { return schema_string_; }

  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() const { return mesh_; }
  Teuchos::RCP<SuperMap> get_supermap() const { return smap_; }

  Teuchos::RCP<Epetra_CrsMatrix> A() { return A_; }
  Teuchos::RCP<const Epetra_CrsMatrix> A() const { return A_; }
  Teuchos::RCP<CompositeVector> rhs() { return rhs_; }
  Teuchos::RCP<const CompositeVector> rhs() const { return rhs_; }
  void set_rhs(const Teuchos::RCP<CompositeVector>& rhs) { rhs_ = rhs; }
  void set_coloring(int num_colors, const Teuchos::RCP<std::vector<int>>& coloring)
  {
    num_colors_ = num_colors;
    coloring_ = coloring;
  }

  int apply_calls() { return apply_calls_; }

  // block access
  typedef std::vector<Teuchos::RCP<Op>>::const_iterator const_op_iterator;
  std::size_t size() const { return ops_.size(); }
  const_op_iterator begin() const { return ops_.begin(); }
  const_op_iterator end() const { return ops_.end(); }
  const_op_iterator FindMatrixOp(int schema_dofs, int matching_rule, bool action) const;

  typedef std::vector<Teuchos::RCP<Op>>::iterator op_iterator;
  op_iterator begin() { return ops_.begin(); }
  op_iterator end() { return ops_.end(); }
  op_iterator FindMatrixOp(int schema_dofs, int matching_rule, bool action);

  // block mutate
  void OpPushBack(const Teuchos::RCP<Op>& op) { ops_.push_back(op); }
  void OpExtend(op_iterator begin, op_iterator end);
  void OpReplace(const Teuchos::RCP<Op>& op, int index) { ops_[index] = op; }
  void OpErase(op_iterator begin, op_iterator end) { ops_.erase(begin, end); }

  // quality control
  void Verify() const;

 public:
  // visit methods for Apply
  virtual int ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Cell_Face& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Cell_Node& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Cell_Edge& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Cell_Cell& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Cell_Schema& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;

  virtual int ApplyMatrixFreeOp(const Op_Face_Cell& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Face_Face& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Face_CellBndFace& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Face_Schema& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;

  virtual int ApplyMatrixFreeOp(const Op_Edge_Edge& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;

  virtual int ApplyMatrixFreeOp(const Op_Node_Node& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Node_Schema& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;

  virtual int ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_SurfaceCell_SurfaceCell& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_MeshInjection& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;

  virtual int ApplyMatrixFreeOp(const Op_Diagonal& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;

  // visit methods for symbolic assemble
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Node& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Edge& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Cell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Schema& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;

  virtual void SymbolicAssembleMatrixOp(const Op_Face_Cell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Face_CellBndFace& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Face_Schema& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Face_Face& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;

  virtual void SymbolicAssembleMatrixOp(const Op_Edge_Edge& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;

  virtual void SymbolicAssembleMatrixOp(const Op_Node_Node& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Node_Schema& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;

  virtual void SymbolicAssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_MeshInjection& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;

  virtual void SymbolicAssembleMatrixOp(const Op_Diagonal& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;

  // visit methods for assemble
  virtual void AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Cell_Face& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Cell_Node& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Cell_Edge& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Cell_Cell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Cell_Schema& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;

  virtual void AssembleMatrixOp(const Op_Face_Cell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Face_CellBndFace& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Face_Schema& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Face_Face& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;

  virtual void AssembleMatrixOp(const Op_Edge_Edge& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;

  virtual void AssembleMatrixOp(const Op_Node_Node& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Node_Schema& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;

  virtual void AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_MeshInjection& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;

  virtual void AssembleMatrixOp(const Op_Diagonal& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;

  // deep copy for building interfaces to TPLs, mainly to solvers
  // -- composite vectors
  void CopyVectorToSuperVector(const CompositeVector& cv, Epetra_Vector& sv) const;
  void CopySuperVectorToVector(const Epetra_Vector& sv, CompositeVector& cv) const;

  // diagnostics
  std::string PrintDiagnostics() const;
  double residual() const override
  {
    AMANZI_ASSERT(preconditioner_.get());
    return preconditioner_->residual();
  }
  int num_itrs() const override
  {
    AMANZI_ASSERT(preconditioner_.get());
    return preconditioner_->num_itrs();
  }
  int returned_code() const override
  {
    AMANZI_ASSERT(preconditioner_.get());
    return preconditioner_->returned_code();
  }
  std::string returned_code_string() const override
  {
    AMANZI_ASSERT(preconditioner_.get());
    return preconditioner_->returned_code_string();
  }

 protected:
  int SchemaMismatch_(const std::string& schema1, const std::string& schema2) const;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const CompositeVectorSpace> cvs_row_;
  Teuchos::RCP<const CompositeVectorSpace> cvs_col_;
  Teuchos::ParameterList plist_;

  mutable std::vector<Teuchos::RCP<Op>> ops_;
  Teuchos::RCP<CompositeVector> rhs_, rhs_checkpoint_;

  int ncells_owned, nfaces_owned, nnodes_owned, nedges_owned;
  int ncells_wghost, nfaces_wghost, nnodes_wghost, nedges_wghost;

  int num_colors_;
  Teuchos::RCP<std::vector<int>> coloring_;
  Teuchos::ParameterList inv_plist_;
  bool inverse_pars_set_;
  bool initialize_complete_;
  bool compute_complete_;
  mutable bool assembly_complete_;

  Teuchos::RCP<Epetra_CrsMatrix> A_;
  Teuchos::RCP<Matrix<CompositeVector>> preconditioner_;

  Teuchos::RCP<MatrixFE> Amat_;
  Teuchos::RCP<SuperMap> smap_;

  Teuchos::RCP<VerboseObject> vo_;

  int schema_old_;
  Schema schema_row_, schema_col_;
  std::string schema_string_;
  double shift_, shift_min_;

  mutable int apply_calls_;

 private:
  // Operator(const Operator& op);
  Operator& operator=(const Operator& op);
};

} // namespace Operators
} // namespace Amanzi


#endif
