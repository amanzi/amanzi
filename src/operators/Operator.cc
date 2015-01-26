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
    : cvs_(cvs), data_validity_(true), my_dof_index_(0)
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
      my_dof_index_(op.my_dof_index_),
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
    my_dof_index_ = op.my_dof_index_;
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
  std::vector<std::string> compnames;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  std::vector<Teuchos::RCP<const Epetra_Map> > ghost_maps;
  if (schema & OPERATOR_SCHEMA_DOFS_FACE) {
    ASSERT(cvs_->HasComponent("face"));
    compnames.push_back("face");
    dofnums.push_back(1);
    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cvs_->Mesh(), AmanziMesh::FACE);
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }

  if (schema & OPERATOR_SCHEMA_DOFS_CELL) {
    ASSERT(cvs_->HasComponent("cell"));
    compnames.push_back("cell");
    dofnums.push_back(1);
    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cvs_->Mesh(), AmanziMesh::CELL);
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }
  if (schema & OPERATOR_SCHEMA_DOFS_NODE) {
    ASSERT(cvs_->HasComponent("node"));
    compnames.push_back("node");
    dofnums.push_back(1);
    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cvs_->Mesh(), AmanziMesh::NODE);
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }

  smap_ = Teuchos::rcp(new SuperMap(cvs_->Comm(),
          compnames, dofnums, maps, ghost_maps));

  // estimate size of the matrix graph
  int row_size(0), dim = mesh_->space_dimension();
  if (schema & OPERATOR_SCHEMA_DOFS_FACE) {
    int i = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;

    for (int c = 0; c < ncells_owned; c++) {
      int nfaces = mesh_->cell_get_num_faces(c);
      i = std::max(i, nfaces);
    }
    row_size += 2 * i;
  }
  if (schema & OPERATOR_SCHEMA_DOFS_CELL) {
    int i = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;
    row_size += i + 1;
  }
  if (schema & OPERATOR_SCHEMA_DOFS_NODE) {
    int i = (dim == 2) ? OPERATOR_QUAD_NODES : OPERATOR_HEX_NODES;
    row_size += 8 * i;
  }
    
  // create the graph
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->Map(),
          smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  // populate matrix graph using blocks that fit the schema
  AmanziMesh::Entity_ID_List cells, faces, nodes;
  int lid[OPERATOR_MAX_NODES];

  int nblocks = blocks_.size();
  for (int nb=0; nb!=nblocks; ++nb) {
    int subschema = blocks_schema_[nb] & schema;

    // Non-standard combinations of schemas 
    if ((nonstandard == 1 || nonstandard_symbolic_ == 1) && 
        (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_CELL) && 
        subschema == OPERATOR_SCHEMA_DOFS_CELL) {

      // ELEMENT: face, DOF: cell
      const std::vector<int>& cell_inds = smap_->GhostIndices("cell", my_dof_index_);

      for (int f=0; f!=nfaces_owned; ++f) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

        int ncells = cells.size();
        for (int n=0; n!=ncells; ++n) {
          lid[n] = cell_inds[cells[n]];
        }

        graph->InsertMyIndices(ncells, lid, ncells, lid);
      }

    // Typical representatives of cell-based methods are MFD and FEM.
    } else if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_CELL) {

      if (subschema == OPERATOR_SCHEMA_DOFS_FACE) {
        // ELEMENT: cell, DOF: face
        const std::vector<int>& face_inds = smap_->GhostIndices("face", my_dof_index_);
        
        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_faces(c, &faces);
          int nfaces = faces.size();

          for (int n=0; n!=nfaces; ++n) {
            lid[n] = face_inds[faces[n]];
          }
          graph->InsertMyIndices(nfaces, lid, nfaces, lid);
        }

      } else if (subschema == OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
        // ELEMENT: cell, DOFS: cell and face
        const std::vector<int>& face_inds = smap_->GhostIndices("face", my_dof_index_);
        const std::vector<int>& cell_inds = smap_->GhostIndices("cell", my_dof_index_);

        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_faces(c, &faces);
          int nfaces = faces.size();

          for (int n=0; n!=nfaces; ++n) {
            lid[n] = face_inds[faces[n]];
          }
          lid[nfaces] = cell_inds[c];
          graph->InsertMyIndices(nfaces+1, lid, nfaces+1, lid);
        }

      } else if (subschema == OPERATOR_SCHEMA_DOFS_NODE) {
        // ELEMENT: cell, DOFS: node
        const std::vector<int>& node_inds = smap_->GhostIndices("node", my_dof_index_);

        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_nodes(c, &nodes);
          int nnodes = nodes.size();

          for (int n=0; n!=nnodes; ++n) {
            lid[n] = node_inds[nodes[n]];
          }
          int ierr = graph->InsertMyIndices(nnodes, lid, nnodes, lid);
          ASSERT(!ierr);          
        }

      } else {
        ASSERT(false);
      }

    // Typical representative of face-based methods is FV.
    } else if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_FACE) {

      // ELEMENT: face, DOF: cell
      const std::vector<int>& cell_inds = smap_->GhostIndices("cell", my_dof_index_);

      for (int f=0; f!=nfaces_owned; ++f) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

        int ncells = cells.size();
        for (int n=0; n!=ncells; ++n) {
          lid[n] = cell_inds[cells[n]];
        }

        graph->InsertMyIndices(ncells, lid, ncells, lid);
      }

    // Typical representative of node-based methods is MPFA.
    } else if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_NODE) {
      // not implemented yet
      ASSERT(false);
    }
  }

  // Completing and optimizing the graphs
  graph->FillComplete(smap_->Map(), smap_->Map());

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


/* ******************************************************************
* Assemble elemental face-based matrices into four global matrices. 
****************************************************************** */
void Operator::AssembleMatrix(int schema)
{
  // populate matrix
  Amat_->Zero();

  AmanziMesh::Entity_ID_List cells, faces, nodes;
  int lid[OPERATOR_MAX_NODES + 1];
  double values[OPERATOR_MAX_NODES + 1];

  int nblocks = blocks_.size();
  for (int nb=0; nb!=nblocks; ++nb) {
    int subschema = blocks_schema_[nb] & schema;
    std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[nb];

    if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_CELL) {
      // Element-based assembly on cells.

      if (subschema == OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {

        // ELEMENT: cell, DOFS: face and cell
        const std::vector<int>& face_inds = smap_->GhostIndices("face", my_dof_index_);
        const std::vector<int>& cell_inds = smap_->GhostIndices("cell", my_dof_index_);
        
        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_faces(c, &faces);

          int nfaces = faces.size();
          for (int n=0; n!=nfaces; ++n) {
            lid[n] = face_inds[faces[n]];
          }
          lid[nfaces] = cell_inds[c];

          Amat_->SumIntoMyValues(lid, matrix[c]);
        }
        
      } else if (subschema == OPERATOR_SCHEMA_DOFS_FACE) {
        // ELEMENT: cell, DOFS: face
        const std::vector<int>& face_inds = smap_->GhostIndices("face", my_dof_index_);

        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_faces(c, &faces);

          int nfaces = faces.size();
          for (int n=0; n!=nfaces; ++n) {
            lid[n] = face_inds[faces[n]];
          }

          Amat_->SumIntoMyValues(lid, matrix[c]);
        }
        
      } else if (subschema == OPERATOR_SCHEMA_DOFS_NODE) {
        // ELEMENT: cell, DOFS: node
        const std::vector<int>& node_inds = smap_->GhostIndices("node", my_dof_index_);

        for (int c=0; c!=ncells_owned; ++c) {
          mesh_->cell_get_nodes(c, &nodes);

          int nnodes = nodes.size();
          for (int n=0; n!=nnodes; ++n) {
            lid[n] = node_inds[nodes[n]];
          }


          int ierr = Amat_->SumIntoMyValues(lid, matrix[c]);
          ASSERT(!ierr);
        }
      }

    } else if (blocks_schema_[nb] & OPERATOR_SCHEMA_BASE_FACE) {
      // element loop over faces

      if (subschema == OPERATOR_SCHEMA_DOFS_CELL) {
        // ELEMENT: face, DOF: cell
        const std::vector<int>& cell_inds = smap_->GhostIndices("cell", my_dof_index_);

        for (int f=0; f!=nfaces_owned; ++f) {
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

          int ncells = cells.size();
          for (int n=0; n!=ncells; ++n) {
            lid[n] = cell_inds[cells[n]];
          }

          Amat_->SumIntoMyValues(lid, matrix[f]);
        }
      }
    }
  }

  Amat_->FillComplete();

  // Add diagonal (a hack)
  Epetra_Vector tmp(A_->RowMap());
  A_->ExtractDiagonalCopy(tmp);

  if (diagonal_->HasComponent("face")) {
    const std::vector<int>& face_inds = smap_->Indices("face", my_dof_index_);
    Epetra_MultiVector& diag = *diagonal_->ViewComponent("face");
    for (int f=0; f!=nfaces_owned; ++f) tmp[face_inds[f]] += diag[0][f];
  }

  if (diagonal_->HasComponent("cell")) {
    const std::vector<int>& cell_inds = smap_->Indices("cell", my_dof_index_);
    Epetra_MultiVector& diag = *diagonal_->ViewComponent("cell");
    for (int c=0; c!=ncells_owned; ++c) tmp[cell_inds[c]] += diag[0][c];
  }

  if (diagonal_->HasComponent("node")) {
    const std::vector<int>& node_inds = smap_->Indices("node", my_dof_index_);
    Epetra_MultiVector& diag = *diagonal_->ViewComponent("node");
    for (int v=0; v!=nnodes_owned; ++v) tmp[node_inds[v]] += diag[0][v];
  }

  A_->ReplaceDiagonalValues(tmp);
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

    for (int f=0; f!=nfaces_owned; ++f) {
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
int Operator::Apply(const CompositeVector& X, CompositeVector& Y) const
{
  // initialize ghost elements 
  X.ScatterMasterToGhosted();
  Y.PutScalarMasterAndGhosted(0.0);

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
  int ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy, my_dof_index_);
  ierr |= preconditioner_->ApplyInverse(Xcopy, Ycopy);
  ierr |= CopySuperVectorToCompositeVector(*smap_, Ycopy, Y, my_dof_index_);
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



