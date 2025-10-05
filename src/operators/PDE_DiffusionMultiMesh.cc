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

#include "nanoflann.hpp"

#include "KDTree.hh"
#include "MeshDefs.hh"
#include "Point.hh"

#include "Op_Cell_FaceCell.hh"
#include "Operator_MultiMesh.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_DiffusionMultiMesh.hh"
#include "TreeOperator.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Constructor
****************************************************************** */
PDE_DiffusionMultiMesh::PDE_DiffusionMultiMesh(Teuchos::ParameterList& plist)
  : plist_(plist), d_(0), stability_(1.0)
{
  vo_ = Teuchos::rcp(new VerboseObject("PDE", plist));
}


/* ******************************************************************
* Initialization of the PDE operator
****************************************************************** */
void
PDE_DiffusionMultiMesh::Init(const Teuchos::RCP<State>& S)
{
  // meshes
  names_ = plist_.get<Teuchos::Array<std::string>>("meshes").toVector();
  int nmeshes = names_.size();

  for (int i = 0; i < nmeshes; ++i) {
    meshes_[names_[i]] = S->GetMesh(names_[i]);
  }

  stability_ = plist_.get<double>("stability constant");
  method_ = plist_.get<std::string>("intersection method");
  nparticles_ = plist_.get<int>("number of particles", 10);
  assemble_stability_matrix_ = plist_.get<bool>("assemble stability matrix", false);
  d_ = meshes_[names_[0]]->getSpaceDimension();

  // parse interface data
  Teuchos::ParameterList interfaces = plist_.sublist("interfaces");
  std::vector<std::string> names, rgns;

  for (auto& it : interfaces) {
    auto sublist = interfaces.sublist(it.first);
    names = sublist.get<Teuchos::Array<std::string>>("meshes").toVector();
    AMANZI_ASSERT(names.size() == 2);

    AMANZI_ASSERT(std::find(names_.begin(), names_.end(), names[0]) != names_.end());
    AMANZI_ASSERT(std::find(names_.begin(), names_.end(), names[1]) != names_.end());

    auto mesh1 = meshes_[names[0]];
    auto mesh2 = meshes_[names[1]];
    std::vector<Teuchos::RCP<const AmanziMesh::Mesh>> submeshes = { mesh1, mesh2 };

    rgns = sublist.get<Teuchos::Array<std::string>>("common surface regions").toVector();
    AMANZI_ASSERT(rgns.size() == 2);

    double kappa = sublist.get<double>("contact conductance", 0.0);

    InterfaceData data;
    for (int k = 0; k < 2; ++k) {
      if (method_ == "particles") {
        meshToMeshMapParticles_(*submeshes[k], rgns[k], *submeshes[1 - k], rgns[1 - k], data);
      } else {
        meshToMeshMapConvexHull_(*submeshes[k], rgns[k], *submeshes[1 - k], rgns[1 - k], data);
      }

      interface_weights_.push_back(data);
      interface_meshes_.push_back({ names[k], names[1 - k] });
      interface_conductance_.push_back(kappa);
    }
  }
  int ninterfaces = interface_weights_.size();

  // create tree vector space
  auto tvs = Teuchos::rcp(new TreeVectorSpace());
  for (auto& name : names_) {
    Operators::PDE_DiffusionFactory opfactory(plist_, meshes_[name]);
    opfactory.SetVariableTensorCoefficient(K_[name]);

    auto op = opfactory.Create();
    pdes_.push_back(op);

    auto global_op = op->global_operator();
    auto tmp = Teuchos::rcp(new TreeVectorSpace(global_op->get_domain_map()));
    tvs->PushBack(tmp);
  }

  matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  // create "diagonal" operators
  boundary_block_.resize(nmeshes);
  boundary_data_.resize(nmeshes);
  boundary_conductance_.resize(nmeshes);

  for (int i0 = 0; i0 < nmeshes; ++i0) {
    int nfaces_wghost = meshes_[names_[i0]]->getNumEntities(AmanziMesh::Entity_kind::FACE,
                                                            AmanziMesh::Parallel_kind::ALL);
    boundary_block_[i0].resize(nfaces_wghost, -1);
    boundary_conductance_[i0].resize(nfaces_wghost);

    for (int k = 0; k < ninterfaces; ++k) {
      if (interface_meshes_[k][0] == names_[i0]) {
        int i1 = std::distance(names_.begin(),
                               std::find(names_.begin(), names_.end(), interface_meshes_[k][1]));
        for (auto it : interface_weights_[k]) {
          boundary_block_[i0][it.first] = i1;
          boundary_conductance_[i0][it.first] = interface_conductance_[k];
        }
        boundary_data_[i0].insert(interface_weights_[k].begin(), interface_weights_[k].end());
      }
    }

    auto op = Teuchos::rcp(new Operator_MultiMesh(plist_,
                                                  pdes_[i0]->global_operator(),
                                                  pdes_[i0]->local_op(),
                                                  boundary_block_[i0],
                                                  boundary_data_[i0]));
    matrix_->set_operator_block(i0, i0, op);
  }

  // interface matrices without stability term
  matrices_grad_.resize(nmeshes);
  for (int n = 0; n < nmeshes; ++n) {
    int ncells_owned = meshes_[names_[n]]->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                          AmanziMesh::Parallel_kind::OWNED);
    matrices_grad_[n].resize(ncells_owned);
  }
}


/* *****************************************************************
* Update matrices for PDEs and pealties
* **************************************************************** */
void
PDE_DiffusionMultiMesh::UpdateMatrices(const Teuchos::RCP<const TreeVector>& flux,
                                       const Teuchos::RCP<const TreeVector>& u)
{
  for (int i = 0; i < meshes_.size(); ++i) {
    auto op = pdes_[i]->global_operator();
    op->Init();

    pdes_[i]->UpdateMatrices(Teuchos::null, Teuchos::null);
    ModifyMatrices_(i, boundary_data_[i]);
  }
}


/* *****************************************************************
* Apply boundary conditions
* **************************************************************** */
void
PDE_DiffusionMultiMesh::ModifyMatrices_(int ib, const InterfaceData& data)
{
  auto op = *matrix_->get_operator_block(ib, ib)->begin();

  double beta = (d_ - 2.0) / 2;
  double alpha = std::sqrt(stability_);

  for (int c = 0; c < op->matrices.size(); ++c) {
    const auto& faces = meshes_[names_[ib]]->getCellFaces(c);
    int nfaces = faces.size();

    // count new matrix size
    int ncol(nfaces + 1);
    for (int n = 0; n < nfaces; ++n) {
      auto it = data.find(faces[n]);
      if (it != data.end()) ncol += it->second.size();
    }

    // populate a flux coupling matrix
    if (ncol > nfaces + 1) {
      WhetStone::DenseMatrix C(nfaces + 1, ncol), CT(ncol, nfaces + 1);
      WhetStone::DenseMatrix L(nfaces + 1, ncol), LT(ncol, nfaces + 1);
      C.PutScalar(0.0);
      L.PutScalar(0.0);

      ncol = nfaces + 1;
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        auto it = data.find(f);
        if (it == data.end()) {
          C(n, n) = 1.0;
        } else {
          double scale = boundary_conductance_[ib][f];
          if (scale > 0.0) {
            scale *= meshes_[names_[ib]]->getFaceArea(f);

            L(n, n) = -scale;
            for (auto coef : it->second) {
              L(n, ncol) = coef.second * scale;
              ncol++;
            }
          } else {
            scale = alpha * std::pow(meshes_[names_[ib]]->getFaceArea(f), beta);

            C(n, n) = 0.5;
            L(n, n) = -scale;
            for (auto coef : it->second) {
              C(n, ncol) = coef.second / 2;
              L(n, ncol) = coef.second * scale;
              ncol++;
            }
          }
        }
      }
      C(nfaces, nfaces) = 1.0;

      CT.Transpose(C);
      LT.Transpose(L);

      matrices_grad_[ib][c] = CT * (op->matrices)[c] * C;
      (op->matrices)[c] = (assemble_stability_matrix_) ? LT * L : matrices_grad_[ib][c] + LT * L;
    }
    else if (assemble_stability_matrix_) {
      (op->matrices)[c].PutScalar(0.0);
    }
  }
}


/* *****************************************************************
* Apply boundary conditions
* **************************************************************** */
void
PDE_DiffusionMultiMesh::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  if (assemble_stability_matrix_) return;

  for (int ib = 0; ib < matrix_->get_row_size(); ++ib) {
    auto op = *matrix_->get_operator_block(ib, ib)->begin();

    const auto& data = boundary_data_[ib];
    const auto& boundary_block = boundary_block_[ib];

    const std::vector<int>& bc_model_trial = (pdes_[ib]->get_bcs_trial())[0]->bc_model();
    const std::vector<int>& bc_model_test = (pdes_[ib]->get_bcs_test())[0]->bc_model();
    const std::vector<double>& bc_value = (pdes_[ib]->get_bcs_trial())[0]->bc_value();

    auto rhs = matrix_->get_operator_block(ib, ib)->rhs();
    auto& rhs_face = *rhs->ViewComponent("face", true);
    auto& rhs_cell = *rhs->ViewComponent("cell");
    rhs->PutScalarGhosted(0.0);

    for (int c = 0; c < op->matrices.size(); ++c) {
      const auto& faces = meshes_[names_[ib]]->getCellFaces(c);
      int nfaces = faces.size();

      bool flag(true);
      WhetStone::DenseMatrix& Acell = op->matrices[c];
      int ncols = Acell.NumCols();

      // essential conditions for test functions
      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        if (bc_model_test[f] == OPERATOR_BC_DIRICHLET) {
          if (flag) { // make a copy of elemental matrix
            op->matrices_shadow[c] = Acell;
            flag = false;
          }
          for (int m = 0; m < ncols; m++) Acell(n, m) = 0.0;
        }
      }

      // conditions for trial functions
      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        double value = bc_value[f];

        if (bc_model_trial[f] == OPERATOR_BC_DIRICHLET) {
          if (flag) { // make a copy of elemental matrix
            op->matrices_shadow[c] = Acell;
            flag = false;
          }

          if (eliminate) {
            int ncol(nfaces + 1);
            for (int m = 0; m < nfaces; m++) {
              int f2 = faces[m];
              rhs_face[0][f2] -= Acell(m, n) * value;
              Acell(m, n) = 0.0;

              auto it = data.find(f2);
              if (it != data.end()) {
                for (auto coef : it->second) {
                  int jb = boundary_block[f2];
                  AMANZI_ASSERT(jb >= 0);
                  auto& tmp_face =
                    *matrix_->get_operator_block(jb, jb)->rhs()->ViewComponent("face", true);

                  tmp_face[0][coef.first] -= Acell(ncol, n) * value;
                  Acell(ncol, n) = 0.0;
                  ncol++;
                }
              }
            }

            rhs_cell[0][c] -= Acell(nfaces, n) * value;
            Acell(nfaces, n) = 0.0;
          }

          if (essential_eqn) {
            rhs_face[0][f] = value;
            Acell(n, n) = 1.0;
          }
        }
      }
    }
  }
}


/* ******************************************************************
* Flux for each mesh face
* **************************************************************** */
void
PDE_DiffusionMultiMesh::UpdateFlux(const Teuchos::Ptr<const TreeVector>& u,
                                   const Teuchos::Ptr<TreeVector>& flux)
{
  for (int ib = 0; ib < matrix_->get_row_size(); ++ib) {
    auto op = *matrix_->get_operator_block(ib, ib)->begin();

    const auto& data = boundary_data_[ib];
    const auto& boundary_block = boundary_block_[ib];

    const auto& u_cell = *u->SubVector(ib)->Data()->ViewComponent("cell");
    const auto& u_face = *u->SubVector(ib)->Data()->ViewComponent("face");
    auto& flux_face = *flux->SubVector(ib)->Data()->ViewComponent("face");

    int nfaces_owned = flux_face.MyLength();
    std::vector<int> hits(nfaces_owned, 0);

    for (int c = 0; c < op->matrices.size(); ++c) {
      int nrows = op->matrices[c].NumRows();

      const auto& [faces, dirs] = meshes_[names_[ib]]->getCellFacesAndDirections(c);
      int nfaces = faces.size();
      int ncol(nfaces + 1);

      WhetStone::DenseVector v(nrows), av(nrows), scale(nfaces);
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        v(n) = u_face[0][f];
        scale(n) = 1.0;

        auto it = data.find(f);
        if (it != data.end()) {
          int ib_other = boundary_block[f];
          const auto& u_other = *u->SubVector(ib_other)->Data()->ViewComponent("face");
          for (auto coef : it->second) {
            v(ncol++) = u_other[0][coef.first];
          }
          scale(n) = 2.0;
        }
      }
      v(nfaces) = u_cell[0][c];

      if (ncol > nfaces + 1) {
        AMANZI_ASSERT(matrices_grad_[ib][c].NumCols() == v.NumRows());
        matrices_grad_[ib][c].Multiply(v, av, false);
      } else if (op->matrices_shadow[c].NumRows() == 0) {
        op->matrices[c].Multiply(v, av, false);
      } else {
        op->matrices_shadow[c].Multiply(v, av, false);
      }

      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        flux_face[0][f] -= av(n) * dirs[n] * scale(n);
        hits[f]++;
      }
    }

    for (int f = 0; f != nfaces_owned; ++f) {
      flux_face[0][f] /= hits[f];
    }
  }
}


/* *****************************************************************
* Approximate projection of mesh1 onto mesh2 using particles
* **************************************************************** */
void
PDE_DiffusionMultiMesh::meshToMeshMapParticles_(const AmanziMesh::Mesh& mesh1,
                                                const std::string& rgn1,
                                                const AmanziMesh::Mesh& mesh2,
                                                const std::string& rgn2,
                                                InterfaceData& data12)
{
  auto block1 =
    mesh1.getSetEntities(rgn1, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nblock1 = block1.size();

  data12.clear();

  int f2(-1), dir, stage, count(0);
  double s, r, t;

  // count particle hits from mesh1 to mesh2
  for (int n = 0; n < nblock1; ++n) {
    int f1 = block1[n];
    const auto& coords1 = mesh1.getFaceCoordinates(f1);

    int c1 = getFaceOnBoundaryInternalCell(mesh1, f1);
    const AmanziGeometry::Point& ray = mesh1.getFaceNormal(f1, c1, &dir);

    if (d_ == 2) {
      for (int k = 0; k < nparticles_; ++k) {
        s = double(k + 0.5) / nparticles_;
        auto xa1 = coords1[0] * (1.0 - s) + coords1[1] * s;

        f2 = findFace_(xa1, ray, mesh2, rgn2, f2, &stage);
        if (stage == 1) count++;

        data12[f1][f2] += mesh1.getFaceArea(f1) / nparticles_;
      }
    } else if (d_ == 3) {
      int nnodes = coords1.size();
      AmanziGeometry::Point xa1(d_);

      if (nnodes == 4) { // quadrilateral
        for (int k = 0; k < nparticles_; ++k) {
          s = double(k + 0.5) / nparticles_;

          for (int l = 0; l < nparticles_; ++l) {
            r = double(l + 0.5) / nparticles_;
            xa1 = coords1[0] * (1.0 - s) * (1.0 - r) + coords1[1] * s * (1.0 - r) +
                  coords1[2] * s * r + coords1[3] * (1.0 - s) * r;

            f2 = findFace_(xa1, ray, mesh2, rgn2, f2, &stage);
            if (stage == 1) count++;
            data12[f1][f2] += mesh1.getFaceArea(f1) / nparticles_ / nparticles_;
          }
        }
      } else if (nnodes == 3) { // triangle
        int m(0), ntri((nparticles_ + 1) * nparticles_ / 2);
        for (int k = 0; k < nparticles_; ++k) {
          s = double(k + 0.333) / nparticles_;

          for (int l = 0; l < nparticles_; ++l) {
            r = double(l + 0.333) / nparticles_;
            t = 1.0 - s - r;
            if (t > 0.0) {
              xa1 = coords1[0] * s + coords1[1] * r + coords1[2] * t;

              f2 = findFace_(xa1, ray, mesh2, rgn2, f2, &stage);
              if (stage == 1) count++;
              data12[f1][f2] += mesh1.getFaceArea(f1) / ntri;
              m++;
            }
          }
        }
        AMANZI_ASSERT(m == ntri);
      } else { 
        AMANZI_ASSERT(false);
      }
    }
  }

  // -- statistics
  double sum1(0.0), sum2(0.0);
  for (int f1 : block1) sum1 += mesh1.getFaceArea(f1);

  for (int n = 0; n < nblock1; ++n) {
    int f1 = block1[n];
    for (auto& data : data12[f1]) sum2 += data.second;
  }

  // compute interpolation coefficients for a scalar quantity
  for (int n = 0; n < nblock1; ++n) {
    int f1 = block1[n];
    double area = mesh1.getFaceArea(f1);
    for (auto& data : data12[f1]) data.second /= area;
  }

  // -- verify
  double sum, tol(1e-10);
  for (int n = 0; n < nblock1; ++n) {
    int f1 = block1[n];
    sum = 0.0;
    for (auto& data : data12[f1]) {
      sum += data.second;
      AMANZI_ASSERT(data.second <= 1.0 + tol);
    }
    AMANZI_ASSERT(sum <= 1.0 + tol);
  }

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    int count2 = nblock1 * std::pow(nparticles_, d_ - 1);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "interface names: \"" << rgn1 << "\" or \"" << rgn2 << "\", areas: " << sum1
               << " " << sum2 << "\n  " << count << " particles in the 2nd stage out of " << count2
               << " (" << 100.0 * count / count2 << "%)" << std::endl;
  }
}


/* *****************************************************************
* Convex interpolation on mesh2 for face cetroids on mesh1
* **************************************************************** */
void
PDE_DiffusionMultiMesh::meshToMeshMapConvexHull_(const AmanziMesh::Mesh& mesh1,
                                                 const std::string& rgn1,
                                                 const AmanziMesh::Mesh& mesh2,
                                                 const std::string& rgn2,
                                                 InterfaceData& data12)
{
  const auto& block1 = mesh1.getSetEntities(rgn1, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nblock1 = block1.size();

  const auto& block2 = mesh2.getSetEntities(rgn2, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nblock2 = block2.size();

  // construct a kd-tree index:
  std::vector<Amanzi::AmanziGeometry::Point> points;
  for (int f2 : block2) points.push_back(mesh2.getFaceCentroid(f2));

  AmanziMesh::KDTree tree;
  tree.Init(&points);

  // find closest face centroid on mesh2
  int f2a, f2b, f2c, count(0);
  bool flag, found;
  double d2a, d2b;

  data12.clear();

  for (int n1 = 0; n1 < nblock1; ++n1) {
    int f1 = block1[n1];
    const auto& xp1 = mesh1.getFaceCentroid(f1);

    int num(10);
    auto [idx, dist2] = tree.SearchNearest(xp1, num);
    num = idx.size();

    if (d_ == 2) {
      f2a = block2[idx[0]];
      d2a = std::sqrt(dist2[0]);
      const auto& xp2 = mesh2.getFaceCentroid(f2a);
      AmanziGeometry::Point v2a = xp2 - xp1;
      
      found = false;
      for (int n = 1; n < num; ++n) {
        f2b = block2[idx[n]];
        const auto& xp2 = mesh2.getFaceCentroid(f2b);
        AmanziGeometry::Point v2b = xp2 - xp1;

        if (v2a * v2b < 0.0) {
          d2b = std::sqrt(dist2[n]);
          found = true;
          break;
        }
      }

      if (!found) {
        if (dist2[0] < 1e-24) {
          data12[f1][f2a] = 1.0;
        } else {
          f2b = block2[idx[1]];
          d2b = std::sqrt(dist2[1]);

          data12[f1][f2a] = d2b / (d2b - d2a);
          data12[f1][f2b] = 1.0 - d2b / (d2b - d2a);
          count++;
        }
      } else {
        data12[f1][f2a] = d2b / (d2a + d2b);
        data12[f1][f2b] = d2a / (d2a + d2b);
      }

    } else if (d_ == 3) {
      std::array<double, 3> lambdas;

      f2a = block2[idx[0]];
      found = false;

      for (int i1 = 1; i1 < num; ++i1) {
        for (int i2 = i1 + 1; i2 < num; ++i2) {
          flag = pointInTriangle_(xp1, points[idx[0]], points[idx[i1]], points[idx[i2]], lambdas);
          if (flag) {
            f2b = block2[idx[i1]];    
            f2c = block2[idx[i2]];    
            data12[f1][f2a] = lambdas[0];
            data12[f1][f2b] = lambdas[1];
            data12[f1][f2c] = lambdas[2];

            found = true;
            break;
          }
        }
        if (found) break;
      }

      // boundary effect: take the closest point
      if (!found) {
        data12[f1][f2a] = 1.0;
        count++;
      }
    }
  }

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "interface names: \"" << rgn1 << "\" or \"" << rgn2 
               << "\n " << count << " points outside a simplex " 
               << " (" << 100.0 * count / nblock1 << "%)" << std::endl;
  }
}


/* *****************************************************************
* Given a point in mesh1, it returns face in mesh2
* **************************************************************** */
int
PDE_DiffusionMultiMesh::findFace_(const AmanziGeometry::Point& xf1,
                                  const AmanziGeometry::Point& ray,
                                  const AmanziMesh::Mesh& mesh2,
                                  const std::string& rgn2,
                                  int f2_guess,
                                  int* stage)
{
  const auto& block =
    mesh2.getSetEntities(rgn2, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nblock = block.size();

  // check of ray hits interior of a face
  *stage = 0;
  int n0 = (f2_guess >= 0) ? -1 : 0; // check if initial guess is known

  for (int n = n0; n < nblock; ++n) {
    int f2 = (n == -1) ? f2_guess : block[n];
    const auto& xf2 = mesh2.getFaceCentroid(f2);
    const auto& normal = mesh2.getFaceNormal(f2);

    double angle = ray * normal;
    if (std::fabs(angle) < 1e-12) AMANZI_ASSERT(false);
    double s = ((xf1 - xf2) * normal) / angle;
    auto xf1_proj = xf1 - s * ray;

    const auto& coords2 = mesh2.getFaceCoordinates(f2);
    if (d_ == 2 && (coords2[0] - xf1_proj) * (coords2[1] - xf1_proj) <= 0.0) return f2;
    if (d_ == 3 && point_in_polygon(xf1_proj, coords2) ) return f2;
  }

  // check for closest face
  int f2min(-1);
  double dist_min(1e+50);

  *stage = 1;
  for (int n = 0; n < nblock; ++n) {
    int f2 = block[n];
    const auto& xf2 = mesh2.getFaceCentroid(f2);
    const auto& normal = mesh2.getFaceNormal(f2);

    double angle = ray * normal;
    double s = ((xf1 - xf2) * normal) / angle;
    auto xf1_proj = xf1 - s * ray;

    const auto& coords2 = mesh2.getFaceCoordinates(f2);

    s = norm(coords2[0] - xf1_proj);
    for (int i = 1; i < coords2.size(); ++i) {
      s = std::min(s, norm(coords2[i] - xf1_proj));
    }

    if (s < dist_min) {
      dist_min = s;
      f2min = f2;
    }
  }

  AMANZI_ASSERT(f2min >= 0);
  return f2min;
}


/* *****************************************************************
* Check if point is inside 3D triangle 
* **************************************************************** */
bool
PDE_DiffusionMultiMesh::pointInTriangle_(const AmanziGeometry::Point& testpnt,
                                         const AmanziGeometry::Point& xa,
                                         const AmanziGeometry::Point& xb,
                                         const AmanziGeometry::Point& xc,
                                         std::array<double, 3>& lambdas,
                                         double tol)
{
  int i0(0), i1(1);
  if (d_ == 3) {
    double xmax = std::max(std::fabs(xb[0] - xa[0]), std::fabs(xc[0] - xa[0]));
    double ymax = std::max(std::fabs(xb[1] - xa[1]), std::fabs(xc[1] - xa[1]));
    double zmax = std::max(std::fabs(xb[2] - xa[2]), std::fabs(xc[2] - xa[2]));

    // projection on one of the three planes: XY, XZ, and YZ
    if (zmax <= std::min(xmax, ymax)) {
      i0 = 0;
      i1 = 1;
    } else if (ymax <= std::min(xmax, zmax)) {
      i0 = 0;
      i1 = 2;
    } else if (xmax <= std::min(ymax, zmax)) {
      i0 = 1;
      i1 = 2;
    }
  }

  double a[2], b[2], c[2];
  a[0] = xa[i0] - xc[i0];
  a[1] = xa[i1] - xc[i1];

  b[0] = xb[i0] - xc[i0];
  b[1] = xb[i1] - xc[i1];

  c[0] = testpnt[i0] - xc[i0];
  c[1] = testpnt[i1] - xc[i1];

  double det = a[0] * b[1] - a[1] * b[0];
  double norma = a[0] * a[0] + a[1] * a[1];
  double normb = b[0] * b[0] + b[1] * b[1];
  if (std::abs(det) < tol * std::sqrt(norma * normb)) return false;

  lambdas[0] = (b[1] * c[0] - b[0] * c[1]) / det;
  lambdas[1] = (a[0] * c[1] - a[1] * c[0]) / det;
  lambdas[2] = 1.0 - lambdas[0] - lambdas[1];

  if (lambdas[0] >= -tol && lambdas[1] >= -tol && lambdas[2] >= -tol) {
    for (int i = 0; i < 3; ++i) lambdas[i] = std::fabs(lambdas[i]);
    return true;
  }
  return false;
}

} // namespace Operators
} // namespace Amanzi
