/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "UnitTest++.h"

// Operators
#include "OperatorDefs.hh"
#include "Mini_Diffusion1D.hh"


/* *****************************************************************
* This test diffusion solver one dimension: u(x) = x^3, K = 2.
* **************************************************************** */
void MiniDiffusion1D_Constant(double bcl, int type_l, double bcr, int type_r) {
  using namespace Amanzi;
  using namespace Amanzi::Operators;

  std::cout << "\nTest: 1D elliptic solver: constant coefficient" << std::endl;

  double pl2_err[2], ph1_err[2];
  for (int loop = 0; loop < 2; ++loop) {
    int ncells = (loop + 1) * 30;
    double length(1.0);
    auto mesh = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(ncells + 1));
    // make a non-uniform mesh
    double h = length / ncells;
    for (int i = 0; i < ncells + 1; ++i) (*mesh)(i) = h * i;
    for (int i = 1; i < ncells; ++i) (*mesh)(i) += h * std::sin(3 * h * i) / 4;

    // initialize diffusion operator with constant coefficient
    Mini_Diffusion1D op;
    op.Init(mesh);

    if (loop == 0) {
      double K(2.0);
      op.Setup(K);
    } else {
      auto K = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(ncells));
      for (int i = 0; i < ncells; ++i) (*K)(i) = 2.0;
      op.Setup(K, NULL, NULL);
    }

    op.UpdateMatrices();

    // create right-hand side
    WhetStone::DenseVector& rhs = op.rhs();
    for (int i = 0; i < ncells; ++i) { 
      double xc = op.mesh_cell_centroid(i);
      double hc = op.mesh_cell_volume(i);
      rhs(i) = -12.0 * xc * hc;
    }

    // apply boundary condition
    op.ApplyBCs(bcl, type_l, bcr, type_r);

    // solve the problem
    WhetStone::DenseVector sol(rhs);
    op.ApplyInverse(rhs, sol);  

    // compute error
    double hc, xc, err, pnorm(1.0), hnorm(1.0);
    pl2_err[loop] = 0.0; 
    ph1_err[loop] = 0.0;

    for (int i = 0; i < ncells; ++i) {
      hc = op.mesh_cell_volume(i);
      xc = op.mesh_cell_centroid(i);
      err = xc * xc * xc - sol(i);

      pl2_err[loop] += err * err * hc;
      pnorm += xc * xc * xc * hc;
    }

    pl2_err[loop] = std::pow(pl2_err[loop] / pnorm, 0.5);
    ph1_err[loop] = std::pow(ph1_err[loop] / hnorm, 0.5);
    printf("BCs:%2d%2d  L2(p)=%9.6f H1(p)=%9.6f\n", type_l, type_r, pl2_err[loop], ph1_err[loop]);

    CHECK(pl2_err[loop] < 1e-3 / (loop + 1) && ph1_err[loop] < 1e-4 / (loop + 1));
  }
  CHECK(pl2_err[0] / pl2_err[1] > 3.7);
}


TEST(OPERATOR_MINI_DIFFUSION_CONSTANT) {
  int dir = Amanzi::Operators::OPERATOR_BC_DIRICHLET;
  int neu = Amanzi::Operators::OPERATOR_BC_NEUMANN;
  MiniDiffusion1D_Constant(0.0, dir, 1.0, dir);
  MiniDiffusion1D_Constant(0.0, dir, -6.0, neu);
  MiniDiffusion1D_Constant(0.0, neu, 1.0, dir);
}


/* *****************************************************************
* This test diffusion solver one dimension: u(x) = x^2, K(x) = x+1
* **************************************************************** */
void MiniDiffusion1D_Variable(double bcl, int type_l, double bcr, int type_r) {
  using namespace Amanzi;
  using namespace Amanzi::Operators;

  std::cout << "\nTest: 1D elliptic solver: variable coefficient" << std::endl;

  double pl2_err[2], ph1_err[2];
  for (int loop = 0; loop < 2; ++loop) {
    int ncells = (loop + 1) * 30;
    double length(1.0);
    auto mesh = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(ncells + 1));
    // make a non-uniform mesh
    double h = length / ncells;
    for (int i = 0; i < ncells + 1; ++i) (*mesh)(i) = h * i;
    for (int i = 1; i < ncells; ++i) (*mesh)(i) += h * std::sin(3 * h * i) / 4;

    // initialize diffusion operator with constant coefficient
    Mini_Diffusion1D op;
    op.Init(mesh);

    auto K = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(ncells));
    for (int i = 0; i < ncells; ++i) {
      double xc = op.mesh_cell_centroid(i);
      (*K)(i) = xc + 1.0;
    }
    op.Setup(K, NULL, NULL);

    op.UpdateMatrices();

    // create right-hand side
    WhetStone::DenseVector& rhs = op.rhs();
    for (int i = 0; i < ncells; ++i) { 
      double xc = op.mesh_cell_centroid(i);
      double hc = op.mesh_cell_volume(i);
      rhs(i) = -(4 * xc + 2.0) * hc;
    }

    // apply boundary condition
    op.ApplyBCs(bcl, type_l, bcr, type_r);

    // solve the problem
    WhetStone::DenseVector sol(rhs);
    op.ApplyInverse(rhs, sol);  

    // compute error
    double hc, xc, err, pnorm(1.0), hnorm(1.0);
    pl2_err[loop] = 0.0; 
    ph1_err[loop] = 0.0;

    for (int i = 0; i < ncells; ++i) {
      hc = op.mesh_cell_volume(i);
      xc = op.mesh_cell_centroid(i);
      err = xc * xc - sol(i);

      pl2_err[loop] += err * err * hc;
      pnorm += xc * xc * hc;
    }

    pl2_err[loop] = std::pow(pl2_err[loop] / pnorm, 0.5);
    ph1_err[loop] = std::pow(ph1_err[loop] / hnorm, 0.5);
    printf("BCs:%2d%2d  L2(p)=%9.6f H1(p)=%9.6f\n", type_l, type_r, pl2_err[loop], ph1_err[loop]);

    CHECK(pl2_err[loop] < 1e-3 / (loop + 1) && ph1_err[loop] < 1e-4 / (loop + 1));
  }
  CHECK(pl2_err[0] / pl2_err[1] > 3.7);
}


TEST(OPERATOR_MINI_DIFFUSION_VARIABLE) {
  int dir = Amanzi::Operators::OPERATOR_BC_DIRICHLET;
  int neu = Amanzi::Operators::OPERATOR_BC_NEUMANN;
  MiniDiffusion1D_Variable(0.0, dir, 1.0, dir);
  MiniDiffusion1D_Variable(0.0, dir, -4.0, neu);
  MiniDiffusion1D_Variable(0.0, neu, 1.0, dir);
}

