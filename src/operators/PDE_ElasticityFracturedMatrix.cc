/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

// TPLs
#include "Epetra_BlockMap.h"

// Amanzi
#include "CompositeVector.hh"
#include "MeshFactory.hh"
#include "Op_Cell_Schema.hh"
#include "Operator_Schema.hh"
#include "PDE_ElasticityFracturedMatrix.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void
PDE_ElasticityFracturedMatrix::Init(Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList& schema_list = plist.sublist("schema");
  mfd_ = WhetStone::BilinearFormFactory::Create(schema_list, mesh_);

  Schema my_schema;
  base_ = my_schema.StringToKind(schema_list.get<std::string>("base"));
  my_schema.Init(mfd_, mesh_, base_);

  local_schema_col_ = my_schema;
  local_schema_row_ = my_schema;
  global_schema_col_ = my_schema;
  global_schema_row_ = my_schema;

  // extract mesh in fractures
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(*mesh_->getGeometricModel().get()));
  AmanziMesh::MeshFactory meshfactory(mesh_->getComm(), gm);
  meshfactory.set_preference(AmanziMesh::Preference({ AmanziMesh::Framework::MSTK }));

  std::vector<std::string> names = plist.sublist("schema").get<Teuchos::Array<std::string>>("fracture").toVector();
  fracture_ = meshfactory.create(mesh_, names, AmanziMesh::Entity_kind::FACE);

  // create global and local operators
  cvs_ = CreateFracturedMatrixCVS(mesh_, fracture_, my_schema.get_items());
  global_op_ = Teuchos::rcp(new Operator_Schema(cvs_, plist, my_schema));

  local_op_ = Teuchos::rcp(new Op_Cell_Schema(my_schema, my_schema, mesh_));
  global_op_->OpPushBack(local_op_);

  // create reverse node map
  int nnodes = fracture_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);

  node_to_node_.resize(nnodes_owned, -1);
  for (int n = 0; n < nnodes; ++n) {
    int v = fracture_->getEntityParent(AmanziMesh::Entity_kind::NODE, n);
    node_to_node_[v] = n;
  }

  // other parameters
  C_ = Teuchos::null;
  C_default_.Init(1, 1);
  C_default_(0, 0) = 1.0; // to run the code without a tensor
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void
PDE_ElasticityFracturedMatrix::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                              const Teuchos::Ptr<const CompositeVector>& u)
{
  PDE_Elasticity::UpdateMatrices(flux, u);

  const auto& nmap = *cvs_->Map("node", true);
  const auto& fmap = *cvs_->Map("face", true);

  for (int c = 0; c < ncells_owned; c++) {
    const auto& Acell = local_op_->matrices[c];
    int nrows = Acell.NumRows();
    std::vector<int> map(nrows);

    const auto& nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();

    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    // create maps from nodes to DOFs
    int np(0);
    for (int n = 0; n < nnodes; ++n) {
      int v = nodes[n];
      int ndofs = nmap.ElementSize(v);

      int shift(0);
      if (ndofs > 1) shift = UniqueIndexNodeToCells(*mesh_, *fracture_, v, node_to_node_[v], c);

      for (int k = 0; k < 3; ++k) {
        map[3 * n + k] = np + ndofs * k + shift;
      }
      np += 3 * ndofs;
    }

    // create maps from faces to DOFs
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      int ndofs = fmap.ElementSize(f);

      int shift(0);
      if (ndofs == 2) shift = UniqueIndexFaceToCells(*mesh_, f, c);

      map[3 * nnodes + n] = np + shift;
      np += ndofs;
    }

    // resize element matrices and right-hand side
    if (np != nrows) {
      WhetStone::DenseMatrix Anew(np, np);
      Anew.PutScalar(0.0);

      for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < nrows; ++j) { Anew(map[i], map[j]) = Acell(i, j); }
      }

      local_op_->matrices[c] = Anew;
    }
  }
}


/* ******************************************************************
* Apply boundary conditions to the local matrices. We always zero-out
* matrix rows for essential test BCs. As to trial BCs, there are
* options: eliminate them or not. Finally we may add the essential BC
* the the system of equations as the trivial equations.
*
* Supported BCs are Dirichlet BCs.
****************************************************************** */
void
PDE_ElasticityFracturedMatrix::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  const std::vector<int>& bcv_model_trial = bcs_trial_[0]->bc_model();
  const std::vector<int>& bcv_model_test = bcs_test_[0]->bc_model();
  const std::vector<AmanziGeometry::Point>& bcv_value = bcs_trial_[0]->bc_value_point();

  const std::vector<int>& bcf_model_trial = bcs_trial_[1]->bc_model();
  const std::vector<int>& bcf_model_test = bcs_test_[1]->bc_model();
  const std::vector<double>& bcf_value = bcs_trial_[1]->bc_value();

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);
  Epetra_MultiVector& rhs_node = *global_op_->rhs()->ViewComponent("node", true);

  const auto& nmap = *cvs_->Map("node", true);
  const auto& fmap = *cvs_->Map("face", true);

  for (int c = 0; c != ncells_owned; ++c) {
    const auto& nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();

    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
    int nrows = Acell.NumRows();

    // un-roll multiple DOFs into a linear array
    std::vector<int> lid(nrows), ids(nrows, 0), comp(nrows, 0), bctrial(nrows), bctest(nrows);
    std::vector<double> bcval(nrows);

    int np(0), np0;
    for (int n = 0; n != nnodes; ++n) {
      int v = nodes[n];
      int first = nmap.FirstPointInElement(v);
      int ndofs = nmap.ElementSize(v);

      for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < ndofs; ++k) {
          lid[np] = first + k;
          ids[np] = v;
          comp[np] = i;
          bctrial[np] = bcv_model_trial[v];
          bctest[np] = bcv_model_test[v];
          bcval[np] = bcv_value[v][i];
          np++;
        }
      }
    }
    np0 = np;

    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int first = fmap.FirstPointInElement(f);
      int ndofs = fmap.ElementSize(f);

      for (int k = 0; k < ndofs; ++k) {
        lid[np] = first + k;
        ids[np] = f;
        bctrial[np] = bcf_model_trial[f];
        bctest[np] = bcf_model_test[f];
        bcval[np] = bcf_value[f];
        np++;
      }
    }

    // essential conditions for test functions
    for (int n = 0; n != np; ++n) {
      if (bctest[n] == OPERATOR_BC_DIRICHLET) {
        if (flag) { // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        for (int m = 0; m < nrows; m++) Acell(n, m) = 0.0;
      }
      else if (bctest[n] != OPERATOR_BC_NONE) {
        AMANZI_ASSERT(true);
      }
    }

    // conditions for trial functions
    for (int n = 0; n != np; ++n) {
      double value = bcval[n];
      if (bctrial[n] == OPERATOR_BC_DIRICHLET) {
        // make a copy of elemental matrix for post-processing
        if (flag) {
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }

        if (eliminate) {
          for (int m = 0; m < np; m++) {
            if (m < np0) 
              rhs_node[comp[m]][lid[m]] -= Acell(m, n) * value;
            else 
              rhs_face[0][lid[m]] -= Acell(m, n) * value;
            Acell(m, n) = 0.0;
          }
        }

        if (essential_eqn) {
          if (n < np0) {
            rhs_node[comp[n]][lid[n]] = value;
            auto cells = mesh_->getNodeCells(ids[n]);
            Acell(n, n) = 1.0 / cells.size();
          } else {
            rhs_face[0][lid[n]] = value;
            Acell(n, n) = 1.0;
          }
        }
      }
    }
  }

  global_op_->rhs()->GatherGhostedToMaster(Add);
}


/* ******************************************************************
* WARNING: Since diffusive flux may be discontinuous (e.g. for
* Richards equation), we derive it only once (using flag) and in
* exactly the same manner as other routines.
* **************************************************************** */
void
PDE_ElasticityFracturedMatrix::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                          const Teuchos::Ptr<CompositeVector>& flux)
{
  // Initialize intensity in ghost faces.
  flux->PutScalarMasterAndGhosted(0.0);
  u->ScatterMasterToGhosted("face");

  const Epetra_MultiVector& u_cell = *u->ViewComponent("cell");
  const Epetra_MultiVector& u_face = *u->ViewComponent("face", true);
  Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);

  const auto& fmap = *cvs_->Map("face", true);

  int ndofs_owned = flux->ViewComponent("face")->MyLength();
  int ndofs_wghost = flux_data.MyLength();

  std::vector<int> hits(ndofs_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
    int nfaces = faces.size();

    // un-roll multiple DOFs in a linear array
    int nrows = 2 * nfaces; // pessimistic estimate
    std::vector<int> lid(nrows), mydir(nrows), nohit(nrows, 0), map(nfaces);
    WhetStone::DenseVector v(nrows), av(nrows);

    int np(0);
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int first = fmap.FirstPointInElement(f);
      int ndofs = fmap.ElementSize(f);

      int shift(0);
      if (ndofs == 2) shift = UniqueIndexFaceToCells(*mesh_, f, c);
      map[n] = np + shift;

      for (int k = 0; k < ndofs; ++k) {
        lid[np] = first + k;
        mydir[np] = dirs[n];
        if (ndofs > 1) nohit[np] = 1;

        v(np) = u_face[0][first + k];
        np++;
      }
    }
    v(np) = u_cell[0][c];

    v.Reshape(np + 1);
    av.Reshape(np + 1);

    if (local_op_->matrices_shadow[c].NumRows() == 0) {
      local_op_->matrices[c].Multiply(v, av, false);
    } else {
      local_op_->matrices_shadow[c].Multiply(v, av, false);
    }

    // points of the master/slave interface require special logic
    for (int n = 0; n < np; n++) {
      int g = lid[n];
      if (g < ndofs_owned || nohit[n] == 1) {
        flux_data[0][g] -= av(n) * mydir[n];
        hits[g]++;
      }
      if (nohit[n] == 1) hits[g] = 1;
    }
  }

  for (int g = 0; g != ndofs_owned; ++g) { flux_data[0][g] /= hits[g]; }

  flux->GatherGhostedToMaster();
}

} // namespace Operators
} // namespace Amanzi
