/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Implicit time discretization is unconditionally stable.
*/

#include <string>

#include "DenseVector.hh"
#include "Mini_Diffusion1D.hh"
#include "OperatorDefs.hh"

#include "MultiscaleTransportPorosity_GDPM.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Simple constructor
****************************************************************** */
MultiscaleTransportPorosity_GDPM::MultiscaleTransportPorosity_GDPM(
    Teuchos::ParameterList& plist)
{
  auto& sublist = plist.sublist("generalized dual porosity parameters");
  matrix_nodes_ = sublist.get<int>("number of matrix nodes");

  // depth is defined for each matrix block as A_m / V_m, so in general,
  // it depends on geometry 
  depth_ = sublist.get<double>("matrix depth");
  tau_ = sublist.get<double>("matrix tortuosity");
  mol_diff_ = plist.get<Teuchos::Array<double> >("molecular diffusion").toVector();

  // make uniform mesh inside matrix
  auto mesh = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_ + 1));
  double h = depth_ / matrix_nodes_;
  for (int i = 0; i < matrix_nodes_ + 1; ++i) (*mesh)(i) = h * i;

  // initialize diffusion operators for each species
  int ncomp = mol_diff_.size();
  op_diff_.resize(ncomp);
  for (int i = 0; i < ncomp; ++i) {
    op_diff_[i].Init(mesh);
    op_diff_[i].Setup(mol_diff_[i]);
    op_diff_[i].UpdateMatrices();
  }
}


/* ******************************************************************
* It should be called only once; otherwise, create an evaluator.
****************************************************************** */
double MultiscaleTransportPorosity_GDPM::ComputeSoluteFlux(
    double flux_liquid, double& tcc_f, WhetStone::DenseVector& tcc_m, int icomp,
    double dt, double wcf0, double wcf1, double wcm0, double wcm1, double phi)
{
  // make a copy of static mini-operator
  Operators::Mini_Diffusion1D op(op_diff_[icomp]);

  double scale = phi * tau_;
  op.ScaleMatrix(scale);
  op.AddAccumulationTerm(wcm0, wcm1, dt, tcc_m);

  // get Schur complement due to fracture equation. This is one of
  // a few possible implementations.
  double al, ad, ar, al_mod, ad_mod, tcc_f_mod, beta;
  op.GetMatrixRow(0, &al, &ad, &ar); 

  beta = al * dt / depth_;
  al_mod = al * wcf1 / (wcf1 - beta);
  ad_mod = ad + al - al_mod;
  tcc_f_mod = tcc_f * wcf0 / wcf1;

  op.SetMatrixRow(0, al_mod, ad_mod, ar); 

  // use modified boundary condition
  op.ApplyBCs(tcc_f_mod, Operators::OPERATOR_BC_DIRICHLET,
              0.0, Operators::OPERATOR_BC_NEUMANN);

  op.ApplyInverse(op.rhs(), tcc_m);  

  tcc_f = (wcf0 * tcc_f - beta * tcc_m(0)) / (wcf1 - beta);

  double h = op.mesh_cell_volume(0);
  double tmp = (flux_liquid > 0.0) ? tcc_f : tcc_m(0); 
  return flux_liquid * tmp - al / h * (tcc_f - tcc_m(0));
}

}  // namespace Transport
}  // namespace Amanzi
  
  
