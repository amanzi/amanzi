/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_FECrsGraph.h"

#include "DenseVector.hh"
#include "PreconditionerFactory.hh"

#include "SuperMap.hh"
#include "MatrixFE.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "OperatorUtils.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Default constructor.
****************************************************************** */
Operator::Operator(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) 
    : cvs_(cvs), data_validity_(true)
{
  mesh_ = cvs_->Mesh();
  rhs_ = Teuchos::rcp(new CompositeVector(*cvs_, true));
  diagonal_ = Teuchos::rcp(new CompositeVector(*cvs_, true));

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  Teuchos::ParameterList plist;
  vo_ = Teuchos::rcp(new VerboseObject("Operators", plist));

  nonstandard_symbolic_ = 0;
}


/* ******************************************************************
* Second constructor.
****************************************************************** */
Operator::Operator(const Operator& op) 
    : data_validity_(true), 
      mesh_(op.mesh_),
      cvs_(op.cvs_), 
      blocks_(op.blocks_), 
      blocks_schema_(op.blocks_schema_), 
      blocks_shadow_(op.blocks_shadow_), 
      diagonal_(op.diagonal_),
      rhs_(op.rhs_),
      bc_(op.bc_),
      A_(op.A_),
      Amat_(op.Amat_),
      smap_(op.smap_),
      preconditioner_(op.preconditioner_)
{
  op.data_validity_ = false;
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  Teuchos::ParameterList plist;
  vo_ = Teuchos::rcp(new VerboseObject("Operators", plist));

  nonstandard_symbolic_ = 0;
}


/* ******************************************************************
* Copy data structures from another operator.
****************************************************************** */
Operator&
Operator::operator=(const Operator& op)
{
  if (this != &op) {
    data_validity_ = true; 
    op.data_validity_ = false;
    mesh_ = op.mesh_;
    cvs_ = op.cvs_; 
    blocks_ = op.blocks_;
    blocks_schema_ = op.blocks_schema_;
    blocks_shadow_ = op.blocks_shadow_; 
    diagonal_ = op.diagonal_;
    rhs_ = op.rhs_;
    bc_ = op.bc_;
    A_ = op.A_;
    Amat_ = op.Amat_;
    smap_ = op.smap_;
        
    preconditioner_ = op.preconditioner_;

    ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

    ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

    nonstandard_symbolic_ = op.nonstandard_symbolic_;
  }
  return *this;
}


/* ******************************************************************
* Initialization of the operator.                                           
****************************************************************** */
void Operator::Init()
{
  diagonal_->PutScalarMasterAndGhosted(0.0);
  rhs_->PutScalarMasterAndGhosted(0.0);

  WhetStone::DenseMatrix null_matrix;

  int n = blocks_.size();
  for (int i = 0; i < n; i++) { 
    std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[i];
    std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[i];
    int m = matrix.size();
    for (int k = 0; k < m; k++) {
      matrix[k] = 0.0;
      matrix_shadow[k] = null_matrix;
    }
  }
}


/* ******************************************************************
* Create a global matrix.
* The non-standard option has default value 0. This violation of the
* standards is due to the experimental nature of this option.
****************************************************************** */
void Operator::SymbolicAssembleMatrix(int schema, int nonstandard)
{
  // Create the supermap given a space (set of possible schemas) and a
  // specific schema (assumed/checked to be consistent with the sapce).
  smap_ = createSuperMap(*cvs_, schema, 1);

  // create the graph
  int row_size = MaxRowSize(*mesh_, schema, 1);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->Map(),
      smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  // fill the graph
  SymbolicAssembleMatrix(schema, nonstandard, *smap_, *graph, 0, 0);
  
  // Completing and optimizing the graphs
  int ierr = graph->FillComplete(smap_->Map(), smap_->Map());
  ASSERT(!ierr);

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


/* ******************************************************************
* Populate the sparsity structure of a global matrix
****************************************************************** */
void Operator::SymbolicAssembleMatrix(int schema, int nonstandard,
                                      const SuperMap& map, GraphFE& graph,
                                      int my_block_row, int my_block_col) const
{
  // populate matrix graph using blocks that fit the schema
  AmanziMesh::Entity_ID_List cells, faces, nodes;
  int lid_r[OPERATOR_MAX_NODES];
  int lid_c[OPERATOR_MAX_NODES];

  int nblocks = blocks_.size();
  for (int nb = 0; nb != nblocks; ++nb) {
    int subschema = blocks_schema_[nb] & schema;

    // Non-standard combinations of schemas 
    if ((nonstandard == 1 || nonstandard_symbolic_ == 1) && 
        (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_CELL) && 
        subschema == OPERATOR_SCHEMA_DOFS_CELL) {

      // ELEMENT: face, DOF: cell
      const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
      const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);

      int ierr(0);
      for (int f=0; f!=nfaces_owned; ++f) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

        int ncells = cells.size();
        for (int n=0; n!=ncells; ++n) {
          lid_r[n] = cell_row_inds[cells[n]];
          lid_c[n] = cell_col_inds[cells[n]];
        }

        ierr |= graph.InsertMyIndices(ncells, lid_r, ncells, lid_c);
      }
      ASSERT(!ierr);

    // Typical representatives of cell-based methods are MFD and FEM.
    } else if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_CELL) {

      if (subschema == OPERATOR_SCHEMA_DOFS_FACE) {
        // ELEMENT: cell, DOF: face
        const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
        const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

        int ierr(0);
        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_faces(c, &faces);
          int nfaces = faces.size();

          for (int n=0; n!=nfaces; ++n) {
            lid_r[n] = face_row_inds[faces[n]];
            lid_c[n] = face_col_inds[faces[n]];
          }
          ierr |= graph.InsertMyIndices(nfaces, lid_r, nfaces, lid_c);
        }
        ASSERT(!ierr);

      } else if (subschema == OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
        // ELEMENT: cell, DOFS: cell and face
        const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
        const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);
        const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
        const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);

        int ierr(0);
        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_faces(c, &faces);
          int nfaces = faces.size();

          for (int n=0; n!=nfaces; ++n) {
            lid_r[n] = face_row_inds[faces[n]];
            lid_c[n] = face_col_inds[faces[n]];
          }
          lid_r[nfaces] = cell_row_inds[c];
          lid_c[nfaces] = cell_col_inds[c];
          ierr |= graph.InsertMyIndices(nfaces+1, lid_r, nfaces+1, lid_c);
        }
        ASSERT(!ierr);

      } else if (subschema == OPERATOR_SCHEMA_DOFS_NODE) {
        // ELEMENT: cell, DOFS: node
        const std::vector<int>& node_row_inds = map.GhostIndices("node", my_block_row);
        const std::vector<int>& node_col_inds = map.GhostIndices("node", my_block_col);

        int ierr(0);
        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_nodes(c, &nodes);
          int nnodes = nodes.size();

          for (int n=0; n!=nnodes; ++n) {
            lid_r[n] = node_row_inds[nodes[n]];
            lid_c[n] = node_col_inds[nodes[n]];
          }
          ierr |= graph.InsertMyIndices(nnodes, lid_r, nnodes, lid_c);
        }
        ASSERT(!ierr);          

      } else {
        ASSERT(false);
      }

    // Typical representative of face-based methods is FV.
    } else if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_FACE) {
      // ELEMENT: face, DOF: cell
      const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
      const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);

      int ierr(0);
      for (int f=0; f!=nfaces_owned; ++f) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

        int ncells = cells.size();
        for (int n=0; n!=ncells; ++n) {
          lid_r[n] = cell_row_inds[cells[n]];
          lid_c[n] = cell_col_inds[cells[n]];
        }

        ierr |= graph.InsertMyIndices(ncells, lid_r, ncells, lid_c);
      }
      ASSERT(!ierr);          
      

    // Typical representative of node-based methods is MPFA.
    } else if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_NODE) {
      // not implemented yet
      ASSERT(false);
    }
  }
}


/* ******************************************************************
* Assemble elemental face-based matrices into four global matrices. 
****************************************************************** */
void Operator::AssembleMatrix(int schema)
{
  if (Amat_ == Teuchos::null) {
    Errors::Message msg;
    msg << "Symbolic assembling was not performed.";
    Exceptions::amanzi_throw(msg);
  }

  Amat_->Zero();
  AssembleMatrix(schema, *smap_, *Amat_, 0, 0);
  Amat_->FillComplete();
}

  

/* ******************************************************************
* Assemble elemental face-based matrices into four global matrices. 
****************************************************************** */
void Operator::AssembleMatrix(
    int schema, const SuperMap& map,
    MatrixFE& mat, int my_block_row, int my_block_col) const
{
  AmanziMesh::Entity_ID_List cells, faces, nodes;
  int lid_r[OPERATOR_MAX_NODES + 1];
  int lid_c[OPERATOR_MAX_NODES + 1];
  double values[OPERATOR_MAX_NODES + 1];

  int nblocks = blocks_.size();
  for (int nb = 0; nb != nblocks; ++nb) {
    int subschema = blocks_schema_[nb] & schema;
    std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[nb];

    if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_CELL) {
      if (subschema == OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
        // ELEMENT: cell, DOFS: face and cell
        const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
        const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);
        const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
        const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);        

        int ierr(0);
        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_faces(c, &faces);

          int nfaces = faces.size();
          for (int n=0; n!=nfaces; ++n) {
            lid_r[n] = face_row_inds[faces[n]];
            lid_c[n] = face_col_inds[faces[n]];
          }
          lid_r[nfaces] = cell_row_inds[c];
          lid_c[nfaces] = cell_col_inds[c];

          ierr |= mat.SumIntoMyValues(lid_r, lid_c, matrix[c]);
        }
        ASSERT(!ierr);
        
      } else if (subschema == OPERATOR_SCHEMA_DOFS_FACE) {
        // ELEMENT: cell, DOFS: face
        const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
        const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

        int ierr(0);
        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_faces(c, &faces);

          int nfaces = faces.size();
          for (int n=0; n!=nfaces; ++n) {
            lid_r[n] = face_row_inds[faces[n]];
            lid_c[n] = face_col_inds[faces[n]];
          }

          ierr |= mat.SumIntoMyValues(lid_r, lid_c, matrix[c]);
        }
        ASSERT(!ierr);

      } else if (subschema == OPERATOR_SCHEMA_DOFS_NODE) {
        // ELEMENT: cell, DOFS: node
        const std::vector<int>& node_row_inds = map.GhostIndices("node", my_block_row);
        const std::vector<int>& node_col_inds = map.GhostIndices("node", my_block_col);

        int ierr(0);
        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_nodes(c, &nodes);

          int nnodes = nodes.size();
          for (int n=0; n!=nnodes; ++n) {
            lid_r[n] = node_row_inds[nodes[n]];
            lid_c[n] = node_col_inds[nodes[n]];
          }

          ierr |= mat.SumIntoMyValues(lid_r, lid_c, matrix[c]);
        }
        ASSERT(!ierr);
      }

    } else if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_FACE) {
      // ELEMENT: face, DOFS: cell
      if (subschema == OPERATOR_SCHEMA_DOFS_CELL) {
        const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
        const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);

        int ierr(0);
        for (int f = 0; f != nfaces_owned; ++f) {
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

          int ncells = cells.size();
          for (int n=0; n!=ncells; ++n) {
            lid_r[n] = cell_row_inds[cells[n]];
            lid_c[n] = cell_col_inds[cells[n]];
          }

          ierr |= mat.SumIntoMyValues(lid_r, lid_c, matrix[f]);
        } 
        ASSERT(!ierr);
      }
    }
  }

  // Add diagonal
  if (diagonal_->HasComponent("face")) {
    const std::vector<int>& face_row_inds = map.Indices("face", my_block_row);
    const std::vector<int>& face_col_inds = map.Indices("face", my_block_col);
    Epetra_MultiVector& diag = *diagonal_->ViewComponent("face");
    int ierr(0);
    for (int f = 0; f != nfaces_owned; ++f)
      ierr |= mat.SumIntoMyValues(face_row_inds[f], 1, &diag[0][f], &face_col_inds[f]);
    ASSERT(!ierr);
  }

  if (diagonal_->HasComponent("cell")) {
    const std::vector<int>& cell_row_inds = map.Indices("cell", my_block_row);
    const std::vector<int>& cell_col_inds = map.Indices("cell", my_block_col);
    Epetra_MultiVector& diag = *diagonal_->ViewComponent("cell");
    int ierr(0);
    for (int c = 0; c != ncells_owned; ++c)
      ierr |= mat.SumIntoMyValues(cell_row_inds[c], 1, &diag[0][c], &cell_col_inds[c]);
    ASSERT(!ierr);
  }

  if (diagonal_->HasComponent("node")) {
    const std::vector<int>& node_row_inds = map.Indices("node", my_block_row);
    const std::vector<int>& node_col_inds = map.Indices("node", my_block_col);
    Epetra_MultiVector& diag = *diagonal_->ViewComponent("node");
    int ierr(0);
    for (int v = 0; v != nnodes_owned; ++v) 
      ierr |= mat.SumIntoMyValues(node_row_inds[v], 1, &diag[0][v], &node_col_inds[v]);
    ASSERT(!ierr);
  }
}


/* ******************************************************************
* Applies boundary conditions to matrix_blocks and update the
* right-hand side and the diagonal block.                                           
* NOTE. It will take to implement a few other PKs to realize the level
*       of abstraction of the default implementation below.
****************************************************************** */
void Operator::ApplyBCs()
{
  Errors::Message msg;
  bool applied_bc(false);

  // clean ghosted values
  for (CompositeVector::name_iterator name = rhs_->begin(); name != rhs_->end(); ++name) {
    Epetra_MultiVector& rhs_g = *rhs_->ViewComponent(*name, true);
    int n0 = rhs_->ViewComponent(*name, false)->MyLength();
    int n1 = rhs_g.MyLength();
    for (int i = n0; i < n1; i++) rhs_g[0][i] = 0.0;

    Epetra_MultiVector& diag_g = *diagonal_->ViewComponent(*name, true);
    for (int i = n0; i < n1; i++) diag_g[0][i] = 0.0;
  }

  int nblocks = blocks_.size();
  for (int nb = 0; nb < nblocks; nb++) {
    int schema = blocks_schema_[nb];

    if (schema & OPERATOR_SCHEMA_BASE_CELL) {
      applied_bc  = ApplyBC_Cell_Mixed_(nb);
      applied_bc |= ApplyBC_Cell_Nodal_(nb);
    } else if (schema & OPERATOR_SCHEMA_BASE_FACE) {
      applied_bc = ApplyBC_Face_(nb);
    }
  }

  if (!applied_bc) {
    msg << "Operator boundary conditions: The only supported schema is CELL.";
    Exceptions::amanzi_throw(msg);
  }

  // Account for the ghosted values
  diagonal_->GatherGhostedToMaster(Add);
  rhs_->GatherGhostedToMaster(Add);
}


/* ******************************************************************
* Add additional set of boundary conditions.
******************************************************************* */
void Operator::AddBCs(Teuchos::RCP<BCs> bc)
{
  bc_.push_back(bc);
}


/* ******************************************************************
* Process a cell-based schema.
******************************************************************* */
bool Operator::ApplyBC_Cell_Mixed_(int nb)
{
  Teuchos::RCP<BCs> bc_f = GetBCofType(OPERATOR_BC_TYPE_FACE);
  if (bc_f == Teuchos::null) return false;

  bool applied_bc(false);

  int schema = blocks_schema_[nb];
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[nb];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[nb];

  AmanziMesh::Entity_ID_List faces;

  const std::vector<int>& bc_model = bc_f->bc_model();
  const std::vector<double>& bc_value = bc_f->bc_value();
  const std::vector<double>& bc_mixed = bc_f->bc_mixed();

  if ((schema & OPERATOR_SCHEMA_DOFS_FACE) && (schema & OPERATOR_SCHEMA_DOFS_CELL)) {
    applied_bc = true;
    Epetra_MultiVector& rhs_face = *rhs_->ViewComponent("face", true);
    Epetra_MultiVector& rhs_cell = *rhs_->ViewComponent("cell");
    Epetra_MultiVector& diag = *diagonal_->ViewComponent("face");

    for (int c=0; c!=ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseMatrix& Acell = matrix[c];

      bool flag(true);
      for (int n=0; n!=nfaces; ++n) {
        int f = faces[n];
        double value = bc_value[f];

        if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
          if (flag) {  // make a copy of elemental matrix
            matrix_shadow[c] = Acell;
            flag = false;
          }
          for (int m = 0; m < nfaces; m++) {
            rhs_face[0][faces[m]] -= Acell(m, n) * value;
            Acell(n, m) = Acell(m, n) = 0.0;
          }
          rhs_face[0][f] = value;
          diag[0][f] = 1.0;

          rhs_cell[0][c] -= Acell(nfaces, n) * value;
          Acell(nfaces, n) = 0.0;
          Acell(n, nfaces) = 0.0;
        } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
          rhs_face[0][f] -= value * mesh_->face_area(f);
        } else if (bc_model[f] == OPERATOR_BC_MIXED) {
          if (flag) {  // make a copy of elemental matrix
            matrix_shadow[c] = Acell;
            flag = false;
          }
          double area = mesh_->face_area(f);
          rhs_face[0][f] -= value * area;
          Acell(n, n) += bc_mixed[f] * area;
        }
      }
    }
  }

  return applied_bc;
}


/* ******************************************************************
* Process a cell-based schema.
******************************************************************* */
bool Operator::ApplyBC_Cell_Nodal_(int nb)
{
  bool applied_bc(false);

  int schema = blocks_schema_[nb];
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[nb];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[nb];

  AmanziMesh::Entity_ID_List faces, nodes;

  Teuchos::RCP<BCs> bc_v = GetBCofType(OPERATOR_BC_TYPE_NODE);
  Teuchos::RCP<BCs> bc_f = GetBCofType(OPERATOR_BC_TYPE_FACE);

  if (schema & OPERATOR_SCHEMA_DOFS_NODE) {
    applied_bc = true;
    Epetra_MultiVector& rhs_node = *rhs_->ViewComponent("node", true);
    Epetra_MultiVector& diag = *diagonal_->ViewComponent("node", true);

    for (int c=0; c!=ncells_owned; ++c) {
      bool flag(true);
      WhetStone::DenseMatrix& Acell = matrix[c];

      if (bc_f != Teuchos::null) {
        const std::vector<int>& bc_model = bc_f->bc_model();
        const std::vector<double>& bc_value = bc_f->bc_value();
        const std::vector<double>& bc_mixed = bc_f->bc_mixed();

        mesh_->cell_get_faces(c, &faces);
        int nfaces = faces.size();

        for (int n=0; n!=nfaces; ++n) {
          int f = faces[n];

          if (bc_model[f] == OPERATOR_BC_NEUMANN) {
            double value = bc_value[f];
            double area = mesh_->face_area(f);

            mesh_->face_get_nodes(f, &nodes);
            int nnodes = nodes.size();

            for (int m = 0; m < nnodes; m++) {
              int v = nodes[m];
              rhs_node[0][v] -= value * area / nnodes;
            }
          } else if (bc_model[f] == OPERATOR_BC_MIXED) {
            if (flag) {  // make a copy of cell-based matrix
              matrix_shadow[c] = Acell;
              flag = false;
            }
            double value = bc_value[f];
            double area = mesh_->face_area(f);

            mesh_->face_get_nodes(f, &nodes);
            int nnodes = nodes.size();

            for (int m = 0; m < nnodes; m++) {
              int v = nodes[m];
              rhs_node[0][v] -= value * area / nnodes;
              Acell(n, n) += bc_mixed[f] * area / nnodes;
            }
          }
        }
      }

      if (bc_v != Teuchos::null) {
        const std::vector<int>& bc_model = bc_v->bc_model();
        const std::vector<double>& bc_value = bc_v->bc_value();

        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int n=0; n!=nnodes; ++n) {
          int v = nodes[n];
          double value = bc_value[v];

          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
            if (flag) {  // make a copy of cell-based matrix
              matrix_shadow[c] = Acell;
              flag = false;
            }
            for (int m = 0; m < nnodes; m++) {
              rhs_node[0][nodes[m]] -= Acell(m, n) * value;
              Acell(n, m) = Acell(m, n) = 0.0;
            }
            rhs_node[0][v] = value;
            diag[0][v] = 1.0;
          }
        }
      }
    }
  } 

  return applied_bc;
}


/* ******************************************************************
* Process a face-based schema.
******************************************************************* */
bool Operator::ApplyBC_Face_(int nb)
{
  bool applied_bc(false);

  int schema = blocks_schema_[nb];
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[nb];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[nb];

  AmanziMesh::Entity_ID_List cells, nodes;

  Teuchos::RCP<BCs> bc_f = GetBCofType(OPERATOR_BC_TYPE_FACE);
  const std::vector<int>& bc_model = bc_f->bc_model();
  const std::vector<double>& bc_value = bc_f->bc_value();
  const std::vector<double>& bc_mixed = bc_f->bc_mixed();

  if (schema & OPERATOR_SCHEMA_DOFS_CELL) {
    applied_bc = true;
    Epetra_MultiVector& rhs_cell = *rhs_->ViewComponent("cell");

    for (int f = 0; f != nfaces_owned; ++f) {
      WhetStone::DenseMatrix& Aface = matrix[f];

      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        rhs_cell[0][cells[0]] += bc_value[f] * Aface(0, 0);
      }
      else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        matrix_shadow[f] = Aface;

        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        rhs_cell[0][cells[0]] -= bc_value[f] * mesh_->face_area(f);
        Aface *= 0.0;
      }
      // solve system of two equations in three unknowns
      else if (bc_model[f] == OPERATOR_BC_MIXED) {
        matrix_shadow[f] = Aface;

        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        double area = mesh_->face_area(f);
        double factor = area / (1.0 + bc_mixed[f] * area / Aface(0, 0));
        rhs_cell[0][cells[0]] -= bc_value[f] * factor;
        Aface(0, 0) = bc_mixed[f] * factor;
      }
    }
  }

  return applied_bc;
}


/* ******************************************************************
* Parallel matvec product A * X.                                              
******************************************************************* */
int Operator::Apply(const CompositeVector& X, CompositeVector& Y, double scalar) const
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
    

  // Multiply by the diagonal block.
  Y.Multiply(1.0, *diagonal_, X, 0.0);

  // Multiply by the remaining matrix blocks.
  AmanziMesh::Entity_ID_List faces, cells, nodes;

  int nblocks = blocks_.size();
  for (int nb = 0; nb < nblocks; nb++) {
    std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[nb];
    int schema = blocks_schema_[nb];

    if (schema & OPERATOR_SCHEMA_BASE_CELL) {
      if (schema & OPERATOR_SCHEMA_DOFS_FACE && schema & OPERATOR_SCHEMA_DOFS_CELL) {
        const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
        const Epetra_MultiVector& Xc = *X.ViewComponent("cell");

        Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
        Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_faces(c, &faces);
          int nfaces = faces.size();

          WhetStone::DenseVector v(nfaces + 1), av(nfaces + 1);
          for (int n=0; n!=nfaces; ++n) {
            v(n) = Xf[0][faces[n]];
          }
          v(nfaces) = Xc[0][c];

          WhetStone::DenseMatrix& Acell = matrix[c];
          Acell.Multiply(v, av, false);

          for (int n=0; n!=nfaces; ++n) {
            Yf[0][faces[n]] += av(n);
          }
          Yc[0][c] += av(nfaces);
        } 
      } else if (schema & OPERATOR_SCHEMA_DOFS_NODE) {
        const Epetra_MultiVector& Xv = *X.ViewComponent("node", true);
        Epetra_MultiVector& Yv = *Y.ViewComponent("node", true);

        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_nodes(c, &nodes);
          int nnodes = nodes.size();

          WhetStone::DenseVector v(nnodes), av(nnodes);
          for (int n=0; n!=nnodes; ++n) {
            v(n) = Xv[0][nodes[n]];
          }

          WhetStone::DenseMatrix& Acell = matrix[c];
          Acell.Multiply(v, av, false);

          for (int n=0; n!=nnodes; ++n) {
            Yv[0][nodes[n]] += av(n);
          }
        } 
      }
    // new schema
    } else if (schema & OPERATOR_SCHEMA_BASE_FACE) {
      if (schema & OPERATOR_SCHEMA_DOFS_CELL) {
        const Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);
        Epetra_MultiVector& Yc = *Y.ViewComponent("cell", true);

        for (int f=0; f!=nfaces_owned; ++f) {
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
          int ncells = cells.size();

          WhetStone::DenseVector v(ncells), av(ncells);
          for (int n=0; n!=ncells; ++n) {
            v(n) = Xc[0][cells[n]];
          }

          WhetStone::DenseMatrix& Aface = matrix[f];
          Aface.Multiply(v, av, false);

          for (int n=0; n!=ncells; ++n) {
            Yc[0][cells[n]] += av(n);
          }
        }
      } 
    }
  }
  Y.GatherGhostedToMaster(Add);

  return 0;
}


/* ******************************************************************
* Linear algebra operations with matrices: r = f - A * u
****************************************************************** */
int Operator::ComputeResidual(const CompositeVector& u, CompositeVector& r)
{
  int ierr = Apply(u, r);
  r.Update(1.0, *rhs_, -1.0);
  return ierr;
}


/* ******************************************************************
* Linear algebra operations with matrices: r = A * u - f                                                 
****************************************************************** */
int Operator::ComputeNegativeResidual(const CompositeVector& u, CompositeVector& r)
{
  int ierr = Apply(u, r);
  r.Update(-1.0, *rhs_, 1.0);
  return ierr;
}


/* ******************************************************************
* Parallel matvec product A * X.                                              
******************************************************************* */
int Operator::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  // Y = X;
  // return 0;
  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());

  int ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy, 0);
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
* Adds time derivative (ss * u - s0 * u0) / dT.
****************************************************************** */
void Operator::AddAccumulationTerm(
    const CompositeVector& u0, const CompositeVector& s0, 
    const CompositeVector& ss, double dT, const std::string& name)
{
  AmanziMesh::Entity_ID_List nodes;

  CompositeVector entity_volume(ss);

  if (name == "cell" && ss.HasComponent("cell")) {
    Epetra_MultiVector& volume = *entity_volume.ViewComponent(name); 

    for (int c=0; c!=ncells_owned; ++c) {
      volume[0][c] = mesh_->cell_volume(c); 
    }
  } else if (name == "face" && ss.HasComponent("face")) {
    // Missing code.
    ASSERT(false);
  } else if (name == "node" && ss.HasComponent("node")) {
    Epetra_MultiVector& volume = *entity_volume.ViewComponent(name, true); 
    volume.PutScalar(0.0);

    for (int c=0; c!=ncells_owned; ++c) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int i = 0; i < nnodes; i++) {
        volume[0][nodes[i]] += mesh_->cell_volume(c) / nnodes; 
      }
    }

    entity_volume.GatherGhostedToMaster(name);
  } else {
    ASSERT(false);
  }

  const Epetra_MultiVector& u0c = *u0.ViewComponent(name);
  const Epetra_MultiVector& s0c = *s0.ViewComponent(name);
  const Epetra_MultiVector& ssc = *ss.ViewComponent(name);

  Epetra_MultiVector& volume = *entity_volume.ViewComponent(name); 
  Epetra_MultiVector& diag = *diagonal_->ViewComponent(name);
  Epetra_MultiVector& rhs = *rhs_->ViewComponent(name);

  int n = u0c.MyLength();
  for (int i = 0; i < n; i++) {
    double factor = volume[0][i] / dT;
    diag[0][i] += factor * ssc[0][i];
    rhs[0][i] += factor * s0c[0][i] * u0c[0][i];
  }
}


/* ******************************************************************
* Adds time derivative ss * (u - u0) / dT.
****************************************************************** */
void Operator::AddAccumulationTerm(
    const CompositeVector& u0, const CompositeVector& ss, 
    double dT, const std::string& name)
{
  AmanziMesh::Entity_ID_List nodes;

  CompositeVector entity_volume(ss);

  if (name == "cell" && ss.HasComponent("cell")) {
    Epetra_MultiVector& volume = *entity_volume.ViewComponent(name); 

    for (int c=0; c!=ncells_owned; ++c) {
      volume[0][c] = mesh_->cell_volume(c); 
    }
  } else if (name == "face" && ss.HasComponent("face")) {
    // Missing code.
    ASSERT(false);
  } else if (name == "node" && ss.HasComponent("node")) {
    Epetra_MultiVector& volume = *entity_volume.ViewComponent(name, true); 
    volume.PutScalar(0.0);

    for (int c=0; c!=ncells_owned; ++c) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int i = 0; i < nnodes; i++) {
        volume[0][nodes[i]] += mesh_->cell_volume(c) / nnodes; 
      }
    }

    entity_volume.GatherGhostedToMaster(name);
  } else {
    ASSERT(false);
  }

  const Epetra_MultiVector& u0c = *u0.ViewComponent(name);
  const Epetra_MultiVector& ssc = *ss.ViewComponent(name);

  Epetra_MultiVector& volume = *entity_volume.ViewComponent(name); 
  Epetra_MultiVector& diag = *diagonal_->ViewComponent(name);
  Epetra_MultiVector& rhs = *rhs_->ViewComponent(name);

  int n = u0c.MyLength();
  for (int i = 0; i < n; i++) {
    double factor = volume[0][i] * ssc[0][i] / dT;
    diag[0][i] += factor;
    rhs[0][i] += factor * u0c[0][i];
  }
}


/* ******************************************************************
* Adds time derivative ss * (u - u0).
****************************************************************** */
void Operator::AddAccumulationTerm(
    const CompositeVector& u0, const CompositeVector& ss,
    const std::string& name)
{
  if (!ss.HasComponent(name)) ASSERT(false);

  const Epetra_MultiVector& u0c = *u0.ViewComponent(name);
  const Epetra_MultiVector& ssc = *ss.ViewComponent(name);

  Epetra_MultiVector& diag = *diagonal_->ViewComponent(name);
  Epetra_MultiVector& rhs = *rhs_->ViewComponent(name);

  int n = u0c.MyLength();
  for (int i = 0; i < n; i++) {
    diag[0][i] += ssc[0][i];
    rhs[0][i] += ssc[0][i] * u0c[0][i];
  }
}


/* ******************************************************************
* Check points allows us to revert boundary conditions, source terms, 
* and accumulation terms. They are useful for operators with constant
* coefficients and varying boundary conditions, e.g. for modeling
* saturated flows.
****************************************************************** */
void Operator::CreateCheckPoint()
{
  diagonal_checkpoint_ = Teuchos::rcp(new CompositeVector(*diagonal_));
  rhs_checkpoint_ = Teuchos::rcp(new CompositeVector(*rhs_));
}

void Operator::RestoreCheckPoint()
{
  // The routine should be called after checkpoint is created.
  ASSERT(diagonal_checkpoint_ != Teuchos::null);
  ASSERT(rhs_checkpoint_ != Teuchos::null);

  // restore accumulation and source terms
  *diagonal_ = *diagonal_checkpoint_;
  *rhs_ = *rhs_checkpoint_;

  // restore boundary conditions
  int n = blocks_.size();
  for (int i = 0; i < n; i++) { 
    std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[i];
    std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[i];
    int m = matrix.size();
    for (int k = 0; k < m; k++) {
      if (matrix_shadow[k].NumRows() != 0) {
        matrix[k] = matrix_shadow[k];
      }
    }
  }
}


/* ******************************************************************
* Find a matrix block matching the given pattern.
****************************************************************** */
int Operator::FindMatrixBlock(int schema_dofs, int matching_rule, bool action) const
{
  int nblocks = blocks_.size();
  for (int nb = 0; nb < nblocks; nb++) {
    int schema = blocks_schema_[nb];
    if (matching_rule == OPERATOR_SCHEMA_RULE_EXACT) {
      if ((schema & schema_dofs) == schema_dofs) return nb;
    } else if (matching_rule == OPERATOR_SCHEMA_RULE_SUBSET) {
      if (schema & schema_dofs) return nb;
    }
  }

  if (action) {
    Errors::Message msg;
    msg << "Operators: Matching rule " << matching_rule << " not found.\n";
    Exceptions::amanzi_throw(msg);
  }

  return -1;
}


/* ******************************************************************
* Find BC matching the integer tag type. 
****************************************************************** */
const Teuchos::RCP<BCs> Operator::GetBCofType(int type) const
{
  int nbc = bc_.size();
  for (int i = 0; i < nbc; i++) {
    if (bc_[i]->type() == type) return bc_[i];
  }
  return Teuchos::null;
}


/* ******************************************************************
* Extension of Mesh API. 
****************************************************************** */
/*
int Operator::FindFacePositionInCell_(int f, int c)
{
  AmanziMesh::Entity_ID_List faces;
  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  for (int n = 0; n < nfaces; ++n) {
    if (f == faces[n]) return n;
  }
  return -1;
}
*/
}  // namespace Operators
}  // namespace Amanzi



