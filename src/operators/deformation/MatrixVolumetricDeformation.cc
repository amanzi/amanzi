/*
  License:
  Authors: Ethan Coon (ecoon@lanl.gov) (ATS version)

  Deformation optimization matrix
*/

#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "composite_vector_factory.hh"
#include "MatrixVolumetricDeformation.hh"

#define MESH_TYPE 0 // 0 = HEXES, 1 = TRIANGULAR PRISMS

namespace Amanzi {
namespace Operators {


MatrixVolumetricDeformation::MatrixVolumetricDeformation(
          Teuchos::ParameterList& plist,
          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
          const Teuchos::RCP<std::vector<std::string> >& bottom_region_list) :
    plist_(plist),
    mesh_(mesh),
    bottom_region_list_(bottom_region_list) {
  InitializeFromOptions_();
  Assemble_();
  UpdateInverse_();
};


MatrixVolumetricDeformation::MatrixVolumetricDeformation(
          const MatrixVolumetricDeformation& other) :
    plist_(other.plist_),
    mesh_(other.mesh_),
    bottom_region_list_(other.bottom_region_list_) {
  InitializeFromOptions_();
  Assemble_();
  UpdateInverse_();
};


// Apply matrix, b <-- Ax
void MatrixVolumetricDeformation::Apply(const CompositeVector& x,
        const Teuchos::Ptr<CompositeVector>& b) {
  ASSERT(0);
}


// Apply the inverse, x <-- A^-1 b
void MatrixVolumetricDeformation::ApplyInverse(const CompositeVector& b,
        const Teuchos::Ptr<CompositeVector>& x) {
  int ierr(0);

  // ensure we have a solver
  if (prec_method_ == PREC_METHOD_NULL) {
    Errors::Message msg("MatrixVolumetricDeformation::ApplyInverse requires a specified method");
    Exceptions::amanzi_throw(msg);
  }

  // Equations of the form: dVdz x = b, solved via:
  //   dVdz^T * dVdz * x = dVdz^T * b, where operator_ = dVdz^T * dVdz
  // -- form dVdz^T * b
  Teuchos::RCP<CompositeVector> rhs = domain()->CreateVector(false);
  rhs->CreateData();
  dVdz_->SetUseTranspose(true);
  dVdz_->Apply(*b.ViewComponent("cell",false),
              *rhs->ViewComponent("node",false));

  // Must also apply the boundary condition, dz_bottom = 0
  // -- Fix the bottom nodes, they may not move
  {
    Epetra_MultiVector& rhs_n = *rhs->ViewComponent("node",false);

    for (std::vector<std::string>::const_iterator region_name=
             bottom_region_list_->begin();
         region_name!=bottom_region_list_->end(); ++region_name) {
      AmanziMesh::Entity_ID_List region_nodes;
      mesh_->get_set_entities(*region_name, AmanziMesh::NODE,
              AmanziMesh::OWNED, &region_nodes);

      for (AmanziMesh::Entity_ID_List::const_iterator n=region_nodes.begin();
           n!=region_nodes.end(); ++n) {
        rhs_n[0][*n] = 0;
      }
    }
  }

  // Solve the system x = operator_^-1 * rhs
  if (prec_method_ == TRILINOS_ML) {
    ierr |= ml_prec_->ApplyInverse(*rhs->ViewComponent("node",false),
            *x->ViewComponent("node", false));
  } else if (prec_method_ == TRILINOS_ILU) {
    ierr |= ilu_prec_->ApplyInverse(*rhs->ViewComponent("node",false),
            *x->ViewComponent("node", false));
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ierr |= ifp_prec_->ApplyInverse(*rhs->ViewComponent("node",false),
            *x->ViewComponent("node", false));
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG || prec_method_ == HYPRE_EUCLID) {
    ierr |= IfpHypre_->ApplyInverse(*rhs->ViewComponent("node",false),
            *x->ViewComponent("node", false));
#endif
  } else {
    ASSERT(0);
  }
  ASSERT(!ierr);
}

void MatrixVolumetricDeformation::InitializeFromOptions_() {
  // parameters for optimization
  smoothing_ = plist_.get<double>("smoothing coefficient");
  diagonal_shift_ = plist_.get<double>("diagonal shift", 1.e-6);

  // method for inversion
  prec_method_ = PREC_METHOD_NULL;
  if (plist_.isParameter("preconditioner")) {
    std::string precmethodstring = plist_.get<string>("preconditioner");
    if (precmethodstring == "ML") {
      prec_method_ = TRILINOS_ML;
    } else if (precmethodstring == "ILU" ) {
      prec_method_ = TRILINOS_ILU;
    } else if (precmethodstring == "Block ILU" ) {
      prec_method_ = TRILINOS_BLOCK_ILU;
#ifdef HAVE_HYPRE
    } else if (precmethodstring == "HYPRE AMG") {
      prec_method_ = HYPRE_AMG;
    } else if (precmethodstring == "HYPRE Euclid") {
      prec_method_ = HYPRE_EUCLID;
    } else if (precmethodstring == "HYPRE ParaSails") {
      prec_method_ = HYPRE_EUCLID;
#endif
    } else {
#ifdef HAVE_HYPRE
      Errors::Message msg("Matrix_MFD: The specified preconditioner "+precmethodstring+" is not supported, we only support ML, ILU, HYPRE AMG, HYPRE Euclid, and HYPRE ParaSails");
#else
      Errors::Message msg("Matrix_MFD: The specified preconditioner "+precmethodstring+" is not supported, we only support ML, and ILU");
#endif
      Exceptions::amanzi_throw(msg);
    }
  }
};



void MatrixVolumetricDeformation::Assemble_() {
  int ierr = 0;

  // HARD CODED CRUFT!
#if MESH_TYPE
  int nnz = 6;
  int indices[6];
  double values[6];
#else
  int nnz = 8;
  int indices[8];
  double values[8];
#endif


  double eps = 1.e-4;

  // Domain and Range: matrix inverse solves the problem, given dV, what is
  // dnode_z that results in that dV.  Therefore, A * dnode_z = dV
  const Epetra_Map& cell_map = mesh_->cell_epetra_map(false);
  const Epetra_Map& node_map = mesh_->node_epetra_map(false);
  const Epetra_Map& node_map_wghost = mesh_->node_epetra_map(true);

  range_ = Teuchos::rcp(new CompositeVectorFactory());
  range_->SetMesh(mesh_)->SetComponent("cell",AmanziMesh::CELL,1);
  domain_ = Teuchos::rcp(new CompositeVectorFactory());
  domain_->SetMesh(mesh_)->SetComponent("node",AmanziMesh::NODE,1);

  unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  unsigned int nnodes = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);


  // == Create dVdz ==
  // -- create the matrix
  dVdz_ =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy,cell_map,nnz,true));

  // -- Assemble
  for (unsigned int c=0; c!=ncells; ++c) {
    for (int n=0; n!=nnz; ++n) values[n] = 0.;

    AmanziMesh::Entity_ID_List nodes;
    mesh_->cell_get_nodes(c, &nodes);
    ASSERT(nodes.size() == nnz);

    // determine the upward/downward faces
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

    int my_up_n = -1;
    int my_down_n = -1;
    for (int n=0; n!=faces.size(); ++n) {
      if (mesh_->face_normal(faces[n],false,c)[2] > eps) {
        ASSERT(my_up_n < 0);
        my_up_n = n;
      } else if (mesh_->face_normal(faces[n],false,c)[2] < -eps) {
        ASSERT(my_down_n < 0);
        my_down_n = n;
      }
    }
    ASSERT(my_up_n >= 0);
    ASSERT(my_down_n >= 0);

    // Determine the perpendicular area (remember, face_normal() is weighted
    // by face area already).
    double perp_area = mesh_->face_normal(faces[my_up_n])[2];

    // loop over nodes in the top face, setting the Jacobian
    AmanziMesh::Entity_ID_List nodes_up;
    mesh_->face_get_nodes(faces[my_up_n], &nodes_up);
    for (int n=0; n!=nodes_up.size(); ++n) {
      values[n] = perp_area / (nnz/2.);
      indices[n] = node_map_wghost.GID(nodes_up[n]);
    }

    // loop over nodes in the bottom face, setting the Jacobian
    AmanziMesh::Entity_ID_List nodes_down;
    mesh_->face_get_nodes(faces[my_down_n], &nodes_down);
    int nnz_half = nnz/2;
    for (int n=0; n!=nodes_down.size(); ++n) {
      values[n+nnz_half] = -perp_area / (nnz/2.);
      indices[n+nnz_half] = node_map_wghost.GID(nodes_down[n]);
    }

    // ensure all are nonzero
    for (int n=0; n!=nnz; ++n)
      ASSERT(std::abs(values[n]) > 0.);

    // Assemble
    int ierr = dVdz_->InsertGlobalValues(cell_map.GID(c), nnz,
            values, indices);
    ASSERT(!ierr);
  }

  ierr = dVdz_->FillComplete(node_map, cell_map);
  ASSERT(!ierr);

  // dump dvdz
  EpetraExt::RowMatrixToMatlabFile("dvdz.txt", *dVdz_);

  // == Form the normal equations, dVdz^T * dVdz ==

  // calculate nnz
  int *nnz_op = new int[nnodes];
  int max_nnode_neighbors = 0;
  for (unsigned int no=0; no!=nnodes; ++no) {
    AmanziMesh::Entity_ID_List cells;
    mesh_->node_get_cells(no, AmanziMesh::USED,&cells);
    // THIS IS VERY SPECIFIC TO GEOMETRY!
#if MESH_TYPE
    int nnode_neighbors = (cells.size()/2 + 1) * 3;
#else
    int nnode_neighbors = 9*3;
#endif
    max_nnode_neighbors = std::max(nnode_neighbors, max_nnode_neighbors);
    nnz_op[no] = nnode_neighbors;
  }

  operator_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,node_map,nnz_op));
  delete[] nnz_op;

  // form explicitly dVdz^T * dVdz
  EpetraExt::MatrixMatrix::Multiply(*dVdz_,true, *dVdz_,false, *operator_);

  // dump Operator intermediate step
  EpetraExt::RowMatrixToMatlabFile("op_normal_eq.txt", *operator_);

  // == Add in optimization components ==
  // -- Add a small shift to the diagonal
  Epetra_Vector diag(node_map);
  ierr = operator_->ExtractDiagonalCopy(diag);  ASSERT(!ierr);
  for (int n=0; n!=nnodes; ++n) diag[n] += diagonal_shift_;
  ierr = operator_->ReplaceDiagonalValues(diag);  ASSERT(!ierr);

  // dump Operator intermediate step
  EpetraExt::RowMatrixToMatlabFile("op_diag.txt", *operator_);

  // -- Add in a diffusive term
  int *node_indices = new int[max_nnode_neighbors];
  double *node_values = new double[max_nnode_neighbors];

  for (int n=0; n!=nnodes; ++n) {
    int nneighbors;
    int n_gid = node_map.GID(n);
    // get the row
    ierr = operator_->ExtractGlobalRowCopy(n_gid,
            max_nnode_neighbors,nneighbors,
            node_values, node_indices);  ASSERT(!ierr);

    // apply the stencil, - nneighbors * smoothing_ on diag, smoothing_ on
    // off-diagonals
    for (int i=0; i!=nneighbors; ++i)
      node_values[i] += node_indices[i] == n_gid ?
          nneighbors * smoothing_ : -smoothing_;

    // replace the row
    ierr = operator_->ReplaceGlobalValues(n_gid, nneighbors,
            node_values, node_indices);
    ASSERT(!ierr);
  }


  // -- Fix the bottom nodes, they may not move
  for (std::vector<std::string>::const_iterator region_name=
           bottom_region_list_->begin();
       region_name!=bottom_region_list_->end(); ++region_name) {
    AmanziMesh::Entity_ID_List region_nodes;
    mesh_->get_set_entities(*region_name, AmanziMesh::NODE,
                            AmanziMesh::OWNED, &region_nodes);

    for (AmanziMesh::Entity_ID_List::const_iterator n=region_nodes.begin();
         n!=region_nodes.end(); ++n) {
      // extract the row
      int n_gid = node_map.GID(*n);
      int nneighbors;
      ierr = operator_->ExtractGlobalRowCopy(n_gid, max_nnode_neighbors,
              nneighbors, node_values, node_indices); ASSERT(!ierr);

      // zero the row, 1 on diagonal
      for (int i=0; i!=nneighbors; ++i)
        node_values[i] = node_indices[i] == n_gid ? 1. : 0.;

      // replace the row
      ierr = operator_->ReplaceGlobalValues(n_gid, nneighbors,
              node_values, node_indices); ASSERT(!ierr);
    }
  }

  ierr = operator_->FillComplete();  ASSERT(!ierr);

  // dump Operator intermediate step
  EpetraExt::RowMatrixToMatlabFile("op.txt", *operator_);

  // clean up
  delete[] node_indices;
  delete[] node_values;
}



void MatrixVolumetricDeformation::UpdateInverse_() {
  // Set up the solver
  if (prec_method_ == TRILINOS_ML) {
    if (ml_prec_->IsPreconditionerComputed()) ml_prec_->DestroyPreconditioner();
    ml_prec_->SetParameterList(ml_plist_);
    ml_prec_->ComputePreconditioner();
  } else if (prec_method_ == TRILINOS_ILU) {
    ilu_prec_ = Teuchos::rcp(new Ifpack_ILU(&*operator_));
    ilu_prec_->SetParameters(ilu_plist_);
    ilu_prec_->Initialize();
    ilu_prec_->Compute();
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = ifp_plist_.get<int>("overlap",0);
    ifp_plist_.set<std::string>("schwarz: combine mode","Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*operator_, ovl));
    ifp_prec_->SetParameters(ifp_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG) {
    IfpHypre_ = Teuchos::rcp(new Ifpack_Hypre(&*operator_));
    Teuchos::RCP<FunctionParameter> functs[8];
    functs[0] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 0));
    functs[1] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, 0));
    functs[2] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth_));
    functs[3] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles_));
    functs[4] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6));
    functs[5] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold_));
    functs[6] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, hypre_tol_));
    functs[7] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCycleType, 1));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", BoomerAMG);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 8);
    hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs);

    IfpHypre_->SetParameters(hypre_list);
    IfpHypre_->Initialize();
    IfpHypre_->Compute();
  } else if (prec_method_ == HYPRE_EUCLID) {
    IfpHypre_ = Teuchos::rcp(new Ifpack_Hypre(&*operator_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", Euclid);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);

    IfpHypre_->SetParameters(hypre_list);
    IfpHypre_->Initialize();
    IfpHypre_->Compute();
  } else if (prec_method_ == HYPRE_PARASAILS) {
    IfpHypre_ = Teuchos::rcp(new Ifpack_Hypre(&*operator_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", ParaSails);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);

    IfpHypre_->SetParameters(hypre_list);
    IfpHypre_->Initialize();
    IfpHypre_->Compute();
#endif
  }
};


}
}
