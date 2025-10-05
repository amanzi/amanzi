/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Operators

*/

#ifndef AMANZI_OPERATOR_VIRTUAL_HH_
#define AMANZI_OPERATOR_VIRTUAL_HH_

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "KDTree.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi
#include "MeshDefs.hh"
#include "MeshVirtual.hh"
#include "Op_Face_Cell.hh"
#include "Operator_Cell.hh"
#include "PDE_Diffusion.hh"
#include "Tensor.hh"

namespace Amanzi {

/* *****************************************************************
* Supporting class: geometry solver
***************************************************************** */
class PDE_Geometry : public virtual Operators::PDE_Diffusion {
 public:
  PDE_Geometry(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::ParameterList& plist);

  virtual void
  SetTensorCoefficient(const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K) override {};

  virtual void 
  SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                       const Teuchos::RCP<const CompositeVector>& dkdp) override {};

  virtual void
  UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                 const Teuchos::Ptr<const CompositeVector>& u) override;

  virtual void
  UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                 const Teuchos::Ptr<const CompositeVector>& u,
                                 double scalar_factor = 1.0) override {};

  virtual void
  UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                 const Teuchos::Ptr<const CompositeVector>& u,
                                 const Teuchos::Ptr<const CompositeVector>& factor) override {};

  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override {};

  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;
  virtual void ModifyMatrices(const CompositeVector& u) override;
  virtual void ScaleMassMatrices(double s) override {};

  virtual double ComputeTransmissibility(int f) const override { return 0.0; }
  virtual double ComputeGravityFlux(int f) const override { return 0.0; }

  virtual void ScaleMatricesColumns(const CompositeVector& s) override {};

 private:
  Teuchos::ParameterList plist_;
};


/* ******************************************************************
* Initialization
****************************************************************** */
PDE_Geometry::PDE_Geometry(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                           Teuchos::ParameterList& plist)
  : PDE_Diffusion(mesh),
    plist_(plist)
{
  // Define stencil for the FV-type method
  local_op_schema_ = Operators::OPERATOR_SCHEMA_BASE_FACE | Operators::OPERATOR_SCHEMA_DOFS_CELL;
  global_op_schema_ = Operators::OPERATOR_SCHEMA_DOFS_CELL;

  // global operator has cell DoFs
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh_)->SetGhosted(true);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);

  global_op_ = Teuchos::rcp(new Operators::Operator_Cell(cvs, plist, global_op_schema_));

  // local operator is face based
  std::string name = "Diffusion: FACE_CELL";
  local_op_ = Teuchos::rcp(new Operators::Op_Face_Cell(name, mesh));
  global_op_->OpPushBack(local_op_);
}


/* ******************************************************************
* 
****************************************************************** */
void
PDE_Geometry::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                             const Teuchos::Ptr<const CompositeVector>& u)
{
  int row = plist_.get<int>("vector component");
  int pos = plist_.get<int>("polynomial component");

  // updating matrix blocks
  for (int f = 0; f != nfaces_owned; ++f) {
    const auto& xf = mesh_->getFaceCentroid(f);
    const auto& cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface = 0.0;

    // Note: Neumann condition is imposed on all boundary faces and 
    //       the related matrices will be removed from the system
    double tij = (pos == 0) ? 1.0 : xf[pos - 1];
    if (ncells == 2) tij /= norm(mesh_->getCellCentroid(cells[0]) - mesh_->getCellCentroid(cells[1]));

    for (int i = 0; i != ncells; ++i) {
      Aface(i, i) = tij;
      for (int j = i + 1; j != ncells; ++j) {
        Aface(i, j) = -tij;
        Aface(j, i) = -tij;
      }
    }
    local_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Remove the first equations from the system
****************************************************************** */
void
PDE_Geometry::ModifyMatrices(const CompositeVector& u) 
{
  auto& rhs_c = *global_op_->rhs()->ViewComponent("cell", true);

  for (int f = 0; f != nfaces_owned; ++f) {
    const auto& cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    auto& Aface = local_op_->matrices[f];

    for (int i = 0; i != ncells; ++i) {
      int c = cells[i];
      if (mesh_->getEntityGID(AmanziMesh::Entity_kind::CELL, c) == 0) {
        for (int j = 0; j != ncells; ++j) {
          Aface(i, j) = 0.0;
          Aface(j, i) = 0.0;
        }
        Aface(i, i) = 1.0;
        rhs_c[0][c] = 0.0;
      }
    }
  }
}


/* ******************************************************************
* Special implementation of boundary conditions.
****************************************************************** */
void
PDE_Geometry::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  AMANZI_ASSERT(bcs_trial_.size() > 0);
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  auto& rhs_c = *global_op_->rhs()->ViewComponent("cell", true);

  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] != Operators::OPERATOR_BC_NONE) {
      const auto& cells = mesh_->getFaceCells(f);
      int c = cells[0];

      if (bc_model[f] == Operators::OPERATOR_BC_NEUMANN) {
        local_op_->matrices_shadow[f] = local_op_->matrices[f];
        local_op_->matrices[f](0, 0) = 0.0;
        if (primary) rhs_c[0][c] -= bc_value[f] * mesh_->getFaceArea(f);
      } else {
        AMANZI_ASSERT(false);
      }
    }
  }
}


/* *****************************************************************
* Supporting function: generator of random numbers
***************************************************************** */
double random(double range) {
  return 2 * range * (double(rand()) / RAND_MAX - 0.5);
}


/* *****************************************************************
* Supporting function: initialization of a cloud of points
***************************************************************** */
void
generatePointCloud(int nx, int ny,
                   double Lx, double Ly,
                   const std::vector<double>& origin,
                   AmanziMesh::Point_List& face_centroids_bnd, 
                   AmanziMesh::Point_List& face_normals_bnd, 
                   AmanziMesh::Point_List& cell_centroids,
                   std::vector<double>& cell_volumes)
{
  int d = origin.size();
  double hx(Lx / nx), hy(Ly / ny);
  double hmin = std::min(hx, hy);
  AmanziGeometry::Point xp(d);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      double x = hx * (i + 0.5);
      double y = hy * (j + 0.5);

      // xp[0] = x;
      // xp[1] = y;
      xp[0] = x + random(hmin / 12);
      xp[1] = y + random(hmin / 12);
      cell_centroids.push_back(xp);
      cell_volumes.push_back(Lx * Ly / nx / ny * (1.0 + random(0.0 / 12)));
    }
  }

  for (int n = 0; n < nx; ++n) {
    face_centroids_bnd.push_back(AmanziGeometry::Point(hx * n + hx / 2, 0.0));
    face_normals_bnd.push_back(AmanziGeometry::Point(0.0, -hx));

    face_centroids_bnd.push_back(AmanziGeometry::Point(Lx - hx * n - hx / 2, Ly));
    face_normals_bnd.push_back(AmanziGeometry::Point(0.0, hx));
  }
  for (int n = 0; n < ny; ++n) {
    face_centroids_bnd.push_back(AmanziGeometry::Point(Lx, hy * n + hy / 2));
    face_normals_bnd.push_back(AmanziGeometry::Point(hy, 0.0));

    face_centroids_bnd.push_back(AmanziGeometry::Point(0.0, Ly - hy * n - hy / 2));
    face_normals_bnd.push_back(AmanziGeometry::Point(-hy, 0.0));
  }

  // shift coordinates
  for (int n = 0; n < cell_centroids.size(); ++n) {
    for (int k = 0; k < d; ++k) cell_centroids[n][k] += origin[k];
  }
  for (int n = 0; n < face_centroids_bnd.size(); ++n) {
    for (int k = 0; k < d; ++k) face_centroids_bnd[n][k] += origin[k];
  }
}


/* *****************************************************************
* Supporting function: connectivity graph
***************************************************************** */
std::vector<AmanziMesh::Entity_ID_List>
createConnectivityGraph(std::vector<AmanziGeometry::Point>& face_centroids_bnd,
                        std::vector<AmanziGeometry::Point>& cell_centroids)
{
  double angleMin = 0.8660;  // 30 degress

  int d = cell_centroids[0].dim();
  int ncells = cell_centroids.size();
  int nfaces_bnd = face_centroids_bnd.size();

  std::vector<AmanziMesh::Entity_ID_List> graph_int;
  std::vector<AmanziMesh::Entity_ID_List> graph_bnd(nfaces_bnd);

  AmanziMesh::KDTree tree;
  tree.Init(&cell_centroids);

  // boundary <-> internal pairs
  // we take the first point inside
  std::vector<int> count(ncells, 0);

  for (int n = 0; n < nfaces_bnd; ++n) {
    auto [idx, dist2] = tree.SearchNearest(face_centroids_bnd[n], 20);
    int nresults = idx.size();

    for (int i = 0; i < nresults; ++i) {
      int m = idx[i];
      graph_bnd[n].push_back(m);
      count[m]++;
      break;
    }
  }

  // internal <-> internal pairs
  int nmin_glb(1000000), nmax_glb(0), nfilter_max(0), nfilter_avg(0);

  for (int n = 0; n < ncells; ++n) {
    std::vector<std::pair<size_t, double>> matches, filter0, filter1;

    nanoflann::SearchParams params;
    double query_pt[2] = { cell_centroids[n][0], cell_centroids[n][1] };

    int nresults, mresults(0), nmin(4 - count[n]), nmax(6 - count[n]);
    int num = nmin - 1;
    while (mresults < nmin || mresults > nmax) {
      num++;
      AMANZI_ASSERT(num < 100);

      auto [idx, dist2] = tree.SearchNearest(cell_centroids[n], num);
      int nresults = idx.size();

      // filter results
      // -- take only internal points
      filter0.clear();
      for (int k = 0; k < nresults; ++k) {
        int m = idx[k];
        if (m != n) {
          filter0.push_back(std::make_pair(m, std::sqrt(dist2[k])));
        }
      }
      mresults = filter0.size();

      // -- remove far-away points in small sectors (30 degrees)
      std::vector<int> flag(mresults); 
      for (int i = 0; i < mresults; ++i) {
        double di = filter0[i].second;
        AmanziGeometry::Point xi = cell_centroids[filter0[i].first] - cell_centroids[n];

        for (int j = i + 1; j < mresults; ++j) {
          double dj = filter0[j].second;
          AmanziGeometry::Point xj = cell_centroids[filter0[j].first] - cell_centroids[n];
          double angle = (xi * xj) / di / dj;
          if (angle > angleMin) {
            if (di > dj) flag[i] = 1;
            else flag[j] = 1;
          }
        }
      }

      filter1.clear();
      for (int i = 0; i < mresults; ++i) {
        if (flag[i] == 0) filter1.push_back(filter0[i]);
      }
      mresults = filter1.size();
      
      if (mresults >=nmin && mresults <= nmax) {
        for (auto it : filter1) {
          int m = it.first;
          if (n < m) {
            graph_int.push_back({ n, m });
          }
        }
        nmin_glb = std::min(nmin_glb, mresults + count[n]);
        nmax_glb = std::max(nmax_glb, mresults + count[n]);
        nfilter_max = std::max(nfilter_max, nresults - mresults);
        nfilter_avg += nresults - mresults;
      }
    }
  }

  // copy boundary connection
  for (auto it : graph_bnd) graph_int.push_back(it);
  int nfaces = graph_int.size();

  std::cout << "Number of cells: " << ncells << std::endl;
  std::cout << "Number of internal faces: " << nfaces - nfaces_bnd << std::endl;
  std::cout << "Number of boundary faces: " << nfaces_bnd << std::endl;
  std::cout << "  max links filtered: " << nfilter_max << std::endl;
  std::cout << "  avg links filtered: " << double(nfilter_avg) / ncells << std::endl;
  std::cout << "Number of neighbors: " << nmin_glb << " " << nmax_glb << std::endl;

  return graph_int;
}


/* *****************************************************************
* Supporting function: compute moments
***************************************************************** */
Teuchos::RCP<CompositeVector>
computeGeometryMoments(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                       int nfaces_int, int nfaces,
                       Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  int d = mesh->getSpaceDimension();

  auto preconditioner_list = Teuchos::sublist(plist, "preconditioners", true);
  auto solver_list = Teuchos::sublist(plist, "solvers", true);

  Amanzi::CompositeVectorSpace cvs;
  cvs.SetMesh(mesh);
  cvs.SetGhosted(true);
  cvs.SetComponent("face", AmanziMesh::FACE, d);
  auto moments = Teuchos::rcp(new CompositeVector(cvs, true));
  auto& moments_f = *moments->ViewComponent("face");

  // operator
  Teuchos::ParameterList pde_list;
  for (int i = 0; i < d; ++i) {
    pde_list.set<int>("vector component", i);
    pde_list.set<int>("polynomial component", 0);
    auto pde = Teuchos::rcp(new PDE_Geometry(mesh, pde_list));
    auto op = pde->global_operator();

    auto bc = Teuchos::rcp(new Operators::BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
    auto& bc_model = bc->bc_model();
    auto& bc_value = bc->bc_value();

    for (int f = nfaces_int; f < nfaces; ++f) {
      const auto& normal = mesh->getFaceNormal(f);
      double area = mesh->getFaceArea(f);

      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = normal[i] / area;
    }

    pde->SetBCs(bc, bc);
    pde->UpdateMatrices(Teuchos::null, Teuchos::null);
    pde->ApplyBCs(true, true, true);
    pde->ModifyMatrices(*op->rhs());

    op->SymbolicAssembleMatrix();
    op->AssembleMatrix();

    op->set_inverse_parameters("Hypre AMG", *preconditioner_list, "PCG", *solver_list);

    auto rhs = op->rhs();
    auto sol(*rhs);
    op->ApplyInverse(*rhs, sol);

    const auto& sol_c = *sol.ViewComponent("cell");

    for (int f = 0; f < nfaces; ++f) {
      const auto& cells = mesh->getFaceCells(f);
      int ncells = cells.size();

      if (ncells == 2) {
        int c1 = cells[0];
        int c2 = cells[1];
        moments_f[i][f] = -(sol_c[0][c2] - sol_c[0][c1]) / norm(mesh->getCellCentroid(c2) - mesh->getCellCentroid(c1));
      } else {
        const auto& normal = mesh->getFaceNormal(f);
        moments_f[i][f] = normal[i]; 
      } 
    }
  }

  return moments;
}

} // namespace Amanzi

#endif
