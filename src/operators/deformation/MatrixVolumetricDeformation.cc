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

#define MESH_TYPE 1 // 0 = HEXES, 1 = TRIANGULAR PRISMS

namespace Amanzi {
namespace Operators {


MatrixVolumetricDeformation::MatrixVolumetricDeformation(
          Teuchos::ParameterList& plist,
          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
          const Teuchos::RCP<const AmanziMesh::Entity_ID_List>& fixed_nodes) :
    plist_(plist),
    mesh_(mesh),
    fixed_nodes_(fixed_nodes) {
  InitializeFromOptions_();
  Assemble_();
  UpdateInverse_();
};


MatrixVolumetricDeformation::MatrixVolumetricDeformation(
          const MatrixVolumetricDeformation& other) :
    plist_(other.plist_),
    mesh_(other.mesh_),
    fixed_nodes_(other.fixed_nodes_) {
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
    for (AmanziMesh::Entity_ID_List::const_iterator n=fixed_nodes_->begin();
         n!=fixed_nodes_->end(); ++n) {
      rhs_n[0][*n] = 0;
    }
  }

  ierr |= IfpHypre_->ApplyInverse(*rhs->ViewComponent("node",false),
				  *x->ViewComponent("node", false));
  ASSERT(!ierr);

  // write a measure of error
  Epetra_MultiVector error(*b.ViewComponent("cell",false));
  dVdz_->SetUseTranspose(false);
  dVdz_->Apply(*x->ViewComponent("node",false), error);
  error.Update(-1., *b.ViewComponent("cell",false), 1.);
  double err(0.);
  error.NormInf(&err);
  std::cout << "ERROR IN dV = " << err << std::endl;

}

void MatrixVolumetricDeformation::InitializeFromOptions_() {
  // parameters for optimization
  smoothing_ = plist_.get<double>("smoothing coefficient");
  diagonal_shift_ = plist_.get<double>("diagonal shift", 1.e-6);

  if ( plist_.isSublist("HYPRE AMG Parameters") ) {
    Teuchos::ParameterList & hypre_list = plist_.sublist("HYPRE AMG Parameters");

    hypre_ncycles_ = hypre_list.get<int>("number of cycles",2);
    hypre_nsmooth_ = hypre_list.get<int>("number of smoothing iterations",2);
    hypre_tol_     = hypre_list.get<double>("tolerance", 1e-12);
    hypre_strong_threshold_ = hypre_list.get<double>("strong threshold", 0.5);
    hypre_verbose_ = hypre_list.get<int>("verbosity level",0);
    hypre_coarsen_type_ = hypre_list.get<int>("coarsen type",0);
    hypre_relax_type_ = hypre_list.get<int>("relax type",6);
    hypre_cycle_type_ = hypre_list.get<int>("cycle type",1);

  } else {
    // set reasonable defaults for HYPRE AMG
    hypre_ncycles_ = 100;
    hypre_nsmooth_ = 2;
    hypre_tol_ = 1e-12;
    hypre_strong_threshold_ = 0.5;
    hypre_verbose_ = 0;
    hypre_coarsen_type_ = 0;
    hypre_relax_type_ = 6;
    hypre_cycle_type_ = 1;

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
    double perp_area_up = mesh_->face_normal(faces[my_up_n])[2];
    double perp_area_down = mesh_->face_normal(faces[my_down_n])[2];
    ASSERT(std::abs(perp_area_up - perp_area_down) < 1.e-6);

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
    // determine the set of node neighbors
    std::set<int> neighbors;
    AmanziMesh::Entity_ID_List faces;
    mesh_->node_get_faces(n, AmanziMesh::USED, &faces);
    for (AmanziMesh::Entity_ID_List::const_iterator f=faces.begin();
         f!=faces.end(); ++f) {
      AmanziMesh::Entity_ID_List neighbor_nodes;
      mesh_->face_get_nodes(*f, &neighbor_nodes);
      neighbors.insert(neighbor_nodes.begin(), neighbor_nodes.end());
    }

    // remove my node
    neighbors.erase(n);

    // apply the stencil
    int nneighbors = neighbors.size() + 1;
    AmanziGeometry::Point center(3);
    mesh_->node_get_coordinates(n,&center);
    node_indices[0] = node_map.GID(n);
    node_values[0] = 0.;

    int lcv = 1;
    for (std::set<int>::const_iterator hn=neighbors.begin();
         hn!=neighbors.end(); ++hn) {
      node_indices[lcv] = node_map_wghost.GID(*hn);

      AmanziGeometry::Point coord(3);
      mesh_->node_get_coordinates(*hn,&coord);
      double smoothing = smoothing_ / AmanziGeometry::norm(center - coord);
      node_values[0] += smoothing;
      node_values[lcv] = - smoothing;
      lcv++;
    }

    // add into the row
    ierr = operator_->SumIntoGlobalValues(node_map.GID(n), nneighbors,
            node_values, node_indices);
    ASSERT(!ierr);
  }

  ierr = operator_->FillComplete();  ASSERT(!ierr);


  // -- Fix the bottom nodes, they may not move
  for (AmanziMesh::Entity_ID_List::const_iterator n=fixed_nodes_->begin();
       n!=fixed_nodes_->end(); ++n) {

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

  ierr = operator_->FillComplete();  ASSERT(!ierr);

  // dump Operator intermediate step
  EpetraExt::RowMatrixToMatlabFile("op.txt", *operator_);

  // clean up
  delete[] node_indices;
  delete[] node_values;
}



void MatrixVolumetricDeformation::UpdateInverse_() {

  IfpHypre_ = Teuchos::rcp(new Ifpack_Hypre(&*operator_));
  Teuchos::RCP<FunctionParameter> functs[8];
  functs[0] = Teuchos::rcp(new FunctionParameter(Solver, &HYPRE_BoomerAMGSetCoarsenType, hypre_coarsen_type_));
  functs[1] = Teuchos::rcp(new FunctionParameter(Solver, &HYPRE_BoomerAMGSetPrintLevel, hypre_verbose_));
  functs[2] = Teuchos::rcp(new FunctionParameter(Solver, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth_));
  functs[3] = Teuchos::rcp(new FunctionParameter(Solver, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles_));
  functs[4] = Teuchos::rcp(new FunctionParameter(Solver, &HYPRE_BoomerAMGSetRelaxType, hypre_relax_type_));
  functs[5] = Teuchos::rcp(new FunctionParameter(Solver, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold_));
  functs[6] = Teuchos::rcp(new FunctionParameter(Solver, &HYPRE_BoomerAMGSetTol, hypre_tol_));
  functs[7] = Teuchos::rcp(new FunctionParameter(Solver, &HYPRE_BoomerAMGSetCycleType, hypre_cycle_type_));

  Teuchos::ParameterList hypre_list;
  hypre_list.set("Solver", BoomerAMG);
  hypre_list.set("SolveOrPrecondition", Solver);
  hypre_list.set("SetPreconditioner", false);
  hypre_list.set("NumFunctions", 8);
  hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs);

  IfpHypre_->SetParameters(hypre_list);
  IfpHypre_->Initialize();
  IfpHypre_->Compute();

};


}
}
