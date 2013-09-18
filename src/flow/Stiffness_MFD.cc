/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <vector>

#include "Epetra_FECrsGraph.h"

#include "mfd3d_diffusion.hh"

#include "Flow_PK.hh"
#include "FlowDefs.hh"
#include "Stiffness_MFD.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Constructor                                      
****************************************************************** */
Stiffness_MFD::Stiffness_MFD(Teuchos::RCP<Flow_State> FS, const Epetra_Map& map)
    : FS_(FS), map_(map)
{ 
  mesh_ = FS_->mesh();
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.                                           
****************************************************************** */
void Stiffness_MFD::CreateMFDstiffnessMatrices(std::vector<WhetStone::Tensor>& K, double factor)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziMesh::Entity_ID_List nodes;

  mfd.ModifyStabilityScalingFactor(factor);

  Avv_cells_.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    WhetStone::DenseMatrix Avv(nnodes, nnodes);
    int ok = mfd.StiffnessMatrix(c, K[c], Avv);

    Avv_cells_.push_back(Avv);

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("Stiffness_MFD: unexpected failure of LAPACK in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Simply allocates memory.                                           
****************************************************************** */
void Stiffness_MFD::CreateMFDrhsVectors()
{
  Fv_cells_.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List nodes;

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    Epetra_SerialDenseVector Fv(nnodes);  // Entries are initialized to 0.
    Fv_cells_.push_back(Fv);
  }
}


/* ******************************************************************
* Applies boundary conditions to elemental stiffness matrices and
* creates elemental rigth-hand-sides.                                           
****************************************************************** */
void Stiffness_MFD::ApplyBoundaryConditions(
    std::vector<int>& bc_model, std::vector<double>& bc_values)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List nodes;

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    WhetStone::DenseMatrix& Bvv = Avv_cells_[c];  // B means elemental.
    Epetra_SerialDenseVector& Fv = Fv_cells_[c];

    for (int n = 0; n < nnodes; n++) {
      int v = nodes[n];
      double value = bc_values[v];

      if (bc_model[v] == FLOW_BC_FACE_PRESSURE) {
        for (int m = 0; m < nnodes; m++) {
          Fv[m] -= Bvv(m, n) * value;
          Bvv(n, m) = Bvv(m, n) = 0.0;
        }
        Bvv(n, n) = 1.0;
        Fv[n] = value;
      }
    }
  }
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Stiffness_MFD::SymbolicAssembleGlobalMatrices()
{
  const Epetra_Map& vmap_wghost = mesh_->node_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? FLOW_QUAD_NODES : FLOW_HEX_NODES;
  Epetra_FECrsGraph vv_graph(Copy, map_, 8*avg_entries_row);

  AmanziMesh::Entity_ID_List nodes;
  int nodes_LID[FLOW_MAX_NODES];  // Contigious memory is required.
  int nodes_GID[FLOW_MAX_NODES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int n = 0; n < nnodes; n++) {
      nodes_LID[n] = nodes[n];
      nodes_GID[n] = vmap_wghost.GID(nodes_LID[n]);
    }
    vv_graph.InsertGlobalIndices(nnodes, nodes_GID, nnodes, nodes_GID);
  }
  vv_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  Avv_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, vv_graph));
  Avv_->GlobalAssemble();

  rhs_ = Teuchos::rcp(new Epetra_Vector(map_));
}


/* ******************************************************************
* Assemble elemental mass matrices into four global matrices. 
* We need an auxiliary GHOST-based vector to assemble the RHS.
****************************************************************** */
void Stiffness_MFD::AssembleGlobalMatrices()
{
  Avv_->PutScalar(0.0);

  const Epetra_Map& vmap_wghost = mesh_->node_map(true);
  AmanziMesh::Entity_ID_List nodes;

  int nodes_LID[FLOW_MAX_NODES];
  int nodes_GID[FLOW_MAX_NODES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int n = 0; n < nnodes; n++) {
      nodes_LID[n] = nodes[n];
      nodes_GID[n] = vmap_wghost.GID(nodes_LID[n]);
    }
    Avv_->SumIntoGlobalValues(nnodes, nodes_GID, Avv_cells_[c].Values());
  }
  Avv_->GlobalAssemble();

  // We repeat some of the loops for code clarity.
  Epetra_Vector rhs_wghost(vmap_wghost);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int n = 0; n < nnodes; n++) {
      int v = nodes[n];
      rhs_wghost[v] += Fv_cells_[c][n];
    }
  }
  FS_->CombineGhostNode2MasterNode(rhs_wghost, Add);

  int nnodes = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  for (int v = 0; v < nnodes; v++) (*rhs_)[v] = rhs_wghost[v];
}


/* ******************************************************************
* Initialization of the preconditioner                                                 
****************************************************************** */
void Stiffness_MFD::InitPreconditioner(int method, Teuchos::ParameterList& prec_list)
{
  method_ = method;

  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    ML_list = prec_list;
    MLprec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Avv_, ML_list, false));
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
#ifdef HAVE_HYPRE
    hypre_ncycles = prec_list.get<int>("cycle applications", 5);  // Boomer AMG parameters
    hypre_nsmooth = prec_list.get<int>("smoother sweeps", 3);
    hypre_tol = prec_list.get<double>("tolerance", 0.0);
    hypre_strong_threshold = prec_list.get<double>("strong threshold", 0.0);
    hypre_verbosity = prec_list.get<int>("verbosity", 0);
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    ifp_plist_ = prec_list;
  }
}


/* ******************************************************************
* Rebuild the preconditioner.                                                 
****************************************************************** */
void Stiffness_MFD::UpdatePreconditioner()
{
  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    if (MLprec->IsPreconditionerComputed()) MLprec->DestroyPreconditioner();
    MLprec->SetParameterList(ML_list);
    MLprec->ComputePreconditioner();
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
#ifdef HAVE_HYPRE
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*Avv_));

    Teuchos::RCP<FunctionParameter> functs[8];
    functs[0] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 0));
    functs[1] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, hypre_verbosity)); 
    functs[2] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth));
    functs[3] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles));
    functs[4] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6)); 
    functs[5] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold)); 
    functs[6] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, hypre_tol)); 
    functs[7] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCycleType, 1));  

    Teuchos::ParameterList hypre_list("Preconditioner List");
    // hypre_list.set("Solver", PCG);
    hypre_list.set("Preconditioner", BoomerAMG);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 8);
    hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs); 

    IfpHypre_Sff_->SetParameters(hypre_list);
    IfpHypre_Sff_->Initialize();
    IfpHypre_Sff_->Compute();
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = ifp_plist_.get<int>("overlap", 0);
    ifp_plist_.set<std::string>("schwarz: combine mode", "Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*Avv_, ovl));
    ifp_prec_->SetParameters(ifp_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();
  }
}


/* ******************************************************************
* Parallel matvec product A * X.                                              
****************************************************************** */
int Stiffness_MFD::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  (*Avv_).Multiply(false, X, Y);
  return 0;
}


/* ******************************************************************
* The OWNED cell-based and face-based d.o.f. are packed together into 
* the X and Y Epetra vectors, with the cell-based in the first part.
****************************************************************** */
int Stiffness_MFD::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    MLprec->ApplyInverse(X, Y);
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) { 
#ifdef HAVE_HYPRE
    IfpHypre_Sff_->ApplyInverse(X, Y);
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    ifp_prec_->ApplyInverse(X, Y);
  }
  return 0;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

