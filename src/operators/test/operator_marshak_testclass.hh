#ifndef OPERATOR_MARSHAK_TESTCLASS_HH_
#define OPERATOR_MARSHAK_TESTCLASS_HH_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

namespace Amanzi{

class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh, double T0)
    : mesh_(mesh),
      TemperatureFloor(T0),
      TemperatureSource(100.0) { 
    auto cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::CELL, 1)
        ->AddComponent("face", AmanziMesh::FACE, 1)
        ->AddComponent("dirichlet_faces", AmanziMesh::BOUNDARY_FACE, 1);

    values_ = cvs->Create();
    derivatives_ = cvs->Create();
  }
  ~HeatConduction() {};

  // main members
  void UpdateValues(const CompositeVector& u,
                    Operators::BCs& bcs) {
    const auto uc = u.ViewComponent<Amanzi::MirrorHost>("cell", true); 
    auto values_c = values_->ViewComponent<Amanzi::MirrorHost>("cell", true);
    auto derivs_c = derivatives_->ViewComponent<Amanzi::MirrorHost>("cell", true);
    auto bc_value = bcs.bc_value();
    auto bc_model = bcs.bc_model();

    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    for (int c = 0; c < ncells; c++) {
      values_c(c,0) = std::pow(uc(c,0), 3.0);
      derivs_c(c,0) = 3.0 * std::pow(uc(c,0), 2.0);
    }

    // add boundary face component
    auto vbf = values_->ViewComponent<Amanzi::MirrorHost>("dirichlet_faces", true);
    auto dbf = derivatives_->ViewComponent<Amanzi::MirrorHost>("dirichlet_faces", true);
    auto ext_face_map = mesh_->exterior_face_map(true);
    auto face_map = mesh_->face_map(true);
    for (int f=0; f!=face_map->getLocalNumElements(); ++f) {
      if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
        AmanziMesh::Entity_ID_View cells;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
        AMANZI_ASSERT(cells.size() == 1);
        int bf = ext_face_map->getLocalElement(face_map->getGlobalElement(f));
        vbf(bf,0) = std::pow(bc_value[f], 3.0);
        dbf(bf,0) = 3.0 * std::pow(bc_value[f], 2.0);
      }
    }
    

  }

  double Conduction(int c, double T) const { return T * T * T; }

  Teuchos::RCP<CompositeVector> values() { return values_; }
  Teuchos::RCP<CompositeVector> derivatives() { return derivatives_; }
 
  double exact(double t, const Amanzi::AmanziGeometry::Point& p) {
    double x = p[0], c = 0.4;
    double xi = c * t - x;
    return (xi > 0.0) ? std::pow(3 * c * (c * t - x), 1.0 / 3) : TemperatureFloor;
  }

 public:
  double TemperatureSource;
  double TemperatureFloor;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<CompositeVector> values_, derivatives_;
};

}  // namespace Amanzi

#endif
