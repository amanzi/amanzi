/*
  WhetStone

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

// TPLs
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "Mesh.hh"
#include "MeshFactory.hh"

// Amanzi::WhetStone
#include "DenseMatrix.hh"
#include "DG_Modal.hh"
#include "MeshMaps_VEM.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Lagrange.hh"
#include "MFD3D_LagrangeSerendipity.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"


/* ****************************************************************
* Test of determinant of transformation
**************************************************************** */
TEST(DG_MAP_DETERMINANT_CELL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Convergence of determinant." << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // create two meshes
  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<Mesh> mesh0 = meshfactory("test/one_pentagon.exo");
  Teuchos::RCP<Mesh> mesh1 = meshfactory("test/one_pentagon.exo");

  // deform the second mesh
  int dim(2), cell(0), nnodes(5), nfaces(5);
  AmanziGeometry::Point xv(dim);
  Entity_ID_List nodeids, faces;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes; ++v) {
    mesh1->node_get_coordinates(v, &xv);

    xv[0] = (xv[0] + xv[0] * xv[0] + xv[0] * xv[0] * xv[0]) / 3;
    xv[1] = (xv[1] + xv[1] * xv[1]) / 2;

    nodeids.push_back(v);
    new_positions.push_back(xv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);

  // collect geometric data
  Teuchos::ParameterList plist;
  plist.set<std::string>("method", "unknown")
       .set<int>("method order", 1)
       .set<std::string>("projector", "H1");
  auto maps = std::make_shared<MeshMaps_VEM>(mesh0, mesh1, plist);

  std::vector<WhetStone::VectorPolynomial> vf(nfaces);
  for (int f = 0; f < nfaces; ++f) {
    maps->VelocityFace(f, vf[f]);
  }

  // cell-baced velocities and Jacobian matrices
  // test piecewise linear deformation (part II)
  VectorPolynomial det;
  VectorPolynomial uc;
  MatrixPolynomial J;

  auto moments = std::make_shared<WhetStone::DenseVector>();
  auto numi = std::make_shared<NumericalIntegration>(mesh0);
  std::vector<const char*> list = {"HarmonicCRk", "L2HarmonicPk", "HarmonicPk", "SerendipityPk"};
  
  for (auto name : list) {
    double fac(0.5), volume = mesh1->cell_volume(cell);
    for (int k = 1; k < 6; ++k) {
      if (std::strcmp(name, "HarmonicCRk") == 0) {
        MFD3D_CrouzeixRaviart mfd(mesh0);
        mfd.set_order(k);
        mfd.H1CellHarmonic(cell, vf, moments, uc);
      } else if (std::strcmp(name, "L2HarmonicPk") == 0) {
        if (k > 2) continue;
        MFD3D_Lagrange mfd(mesh0);
        mfd.set_order(k);
        mfd.L2CellHarmonic(cell, vf, moments, uc);
      } else if (std::strcmp(name, "SerendipityPk") == 0) {
        if (k > 3) continue;
        MFD3D_LagrangeSerendipity mfd(mesh0);
        mfd.set_order(k);
        mfd.L2Cell(cell, vf, moments, uc);
      } else if (std::strcmp(name, "HarmonicPk") == 0) {
        if (k > 3) continue;
        MFD3D_Lagrange mfd(mesh0);
        mfd.set_order(k);
        mfd.H1CellHarmonic(cell, vf, moments, uc);
      }
      maps->Jacobian(uc, J);
      maps->Determinant(1.0, J, det);

      double tmp = numi->IntegratePolynomialCell(cell, det[0]);
      double err = tmp - volume;
      fac /= (k + 1);
      CHECK(fabs(err) < fac);

      uc[0].ChangeOrigin(mesh0->cell_centroid(cell));
      uc[1].ChangeOrigin(mesh0->cell_centroid(cell));
      printf("k=%d  %14s  vol=%8.6g  err=%12.8f  |poly|=%9.6g %9.6g\n",
          k, name, tmp, err, uc[0].NormMax(), uc[1].NormMax());
    }
  }

  delete comm;
}

