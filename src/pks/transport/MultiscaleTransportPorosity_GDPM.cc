/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>

#include "DenseVector.hh"

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
  nnodes_ = sublist.get<int>("number of matrix layers", 2);

  depth_ = sublist.get<double>("matrix depth");
  geometry_ = sublist.get<std::string>("pore space geometry", "planar");
  std::vector<double> mol_diff = plist.get<Teuchos::Array<double> >("molecular diffusion").toVector();

  // make uniform mesh inside matrix
  auto mesh = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(nnodes_ + 1));
  double h = depth_ / nnodes_;
  for (int i = 0; i < nnodes_ + 1; ++i) (*mesh)(i) = h * i;

  // initialize diffusion operators for each species
  int ncomp = mol_diff.size();
  op_diff_.resize(ncomp);
  for (int i = 0; i < ncomp; ++i) {
    op_diff_[i].Init(mesh, geometry_, 1.0, 1.0);
    op_diff_[i].Setup(mol_diff[i]);
  }
}


/* ******************************************************************
* It should be called only once; otherwise, create an evaluator.
****************************************************************** */
double MultiscaleTransportPorosity_GDPM::ComputeSoluteFlux(
    double flux_liquid, double tcc_f, double tcc_m,
    int icomp, double phi, std::vector<double>* tcc_m_aux)
{
  double tmp = (flux_liquid > 0.0) ? tcc_f : tcc_m; 
  return flux_liquid;  // * tmp + omega_ * (tcc_f - tcc_m);
}


/* ******************************************************************
* It should be called only once; otherwise, create an evaluator.
****************************************************************** */
void MultiscaleTransportPorosity_GDPM::UpdateStabilityOutflux(
    double flux_liquid, double* outflux)
{ 
  *outflux += std::max(0.0, flux_liquid);  // + omega_;
}

}  // namespace Transport
}  // namespace Amanzi
  
  
