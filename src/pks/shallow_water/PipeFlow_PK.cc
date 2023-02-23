/*
  Pipe Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Pipe flow model inherited from shallow water model.

  Author: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

#include "PipeFlow_PK.hh"
#include "NumericalFluxFactory.hh"

namespace Amanzi {
namespace ShallowWater {

using CV_t = CompositeVector;
double Pi = 3.14159265359;
double TwoPi = 6.28318530718;

PipeFlow_PK::PipeFlow_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln)
                  : ShallowWater_PK(pk_tree, glist, S, soln), 
                  PK(pk_tree, glist, S, soln){

  pipe_diameter_ = sw_list_->get<double>("pipe diameter", 1.0);

  pipe_cross_section_ = Pi * 0.25 * pipe_diameter_ * pipe_diameter_;

  Manning_coeff_ = sw_list_->get<double>("Manning coefficient", 0.005);

  celerity_ = sw_list_->get<double>("celerity", 100); // m/s

}

//--------------------------------------------------------------------
// Discretization of the friction source term
//--------------------------------------------------------------------
double PipeFlow_PK::NumericalSourceFriction(double h, double qx, double WettedAngle)
{

  double S1 = 0.0;
  if (std::fabs(h) > 1.e-10) { //we have to raise this to the power of 7/3 below so the tolerance needs to be stricter
     double WettedPerimeter = 0.5 * pipe_diameter_ * WettedAngle;
     double num = - g_ * Manning_coeff_ * Manning_coeff_ * pow(WettedPerimeter, 4.0/3.0) * std::fabs(qx) * qx;
     double denom = pow( h, 7.0/3.0);
     S1 = num / denom;
  }

  return S1;
}

//--------------------------------------------------------------------
// Newton solve to compute wetted angle given wetted area
//--------------------------------------------------------------------
void PipeFlow_PK::UpdateWettedAngle(){ 

   auto& h_c = *S_->GetW<CV_t>(ponded_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& WettedAngle_c = *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);//TODO should this be all instead of owned?

   for (int c = 0; c < ncells_owned; ++c) {

       WettedAngle_c[0][c] = ComputeWettedAngleNewton(h_c[0][c]);;

   }

}

//--------------------------------------------------------------------
// Compute wetted area and wetted angle at edge location
//--------------------------------------------------------------------
std::vector<double> PipeFlow_PK::ComputeWettedQuantitiesEdge(int c, int e, double WettedAreaCell, double WettedAngleCell,
                                                          double Bc, double Bmax, const Epetra_MultiVector& B_n)
{

  std::vector <double> W(2,0.0); //vector to return 

  const auto& xc = mesh_->cell_centroid(c);
  const auto& xf = mesh_->face_centroid(e);
  Amanzi::AmanziMesh::Entity_ID_List cfaces;

  bool cell_is_dry, cell_is_fully_flooded, cell_is_partially_wet;
  cell_is_partially_wet = false;
  cell_is_dry = false;
  cell_is_fully_flooded = false;

  if (WettedAreaCell >= pipe_cross_section_) cell_is_fully_flooded = true;
  else if (std::fabs(WettedAreaCell) < 1.e-12)  cell_is_dry = true;
  else  cell_is_partially_wet = true;
  
  if (cell_is_fully_flooded) {
    W[0] = WettedAreaCell; //TODO: think about this decision
    W[1] = TwoPi;
  } else if (cell_is_dry) {
    W[0] = 0.0;
    W[1] = 0.0;
  } else if (cell_is_partially_wet) {
    Amanzi::AmanziMesh::Entity_ID_List cfaces;
    mesh_->cell_get_faces(c, &cfaces);

    double mu_eps_sum = 0.0;

    double htc = ComputePondedDepth(WettedAngleCell) + Bc;

    for (int f = 0; f < cfaces.size(); ++f) {
      Amanzi::AmanziGeometry::Point x0, x1;
      int edge = cfaces[f];

      Amanzi::AmanziMesh::Entity_ID_List face_nodes;
      mesh_->face_get_nodes(edge, &face_nodes);
      int n0 = face_nodes[0], n1 = face_nodes[1];

      mesh_->node_get_coordinates(n0, &x0);
      mesh_->node_get_coordinates(n1, &x1);

      double area = norm((xc - x0) ^ (xc - x1)) / 2.0;

      double epsilon;

      if ((htc < B_n[0][n0]) && (htc < B_n[0][n1])) {
        epsilon = 0.0;
      } else if ((htc >= B_n[0][n0]) && (htc >= B_n[0][n1])) {
        epsilon = 1.0;
      } else {
        epsilon = 0.5;
      }

      mu_eps_sum += (area / mesh_->cell_volume(c)) * (epsilon);
    }

    Amanzi::AmanziMesh::Entity_ID_List face_nodes;
    mesh_->face_get_nodes(e, &face_nodes);

    double ht_edge = 0.0;
    for (int i = 0; i < face_nodes.size(); ++i) {
      if (htc < B_n[0][face_nodes[i]]) {
        ht_edge += B_n[0][face_nodes[i]];
      } else {
        ht_edge += B_n[0][face_nodes[i]] + ((htc - Bc) / mu_eps_sum);
      }
    }
    ht_edge /= 2.0;
    double B_edge = BathymetryEdgeValue(e, B_n);
    W[1] = ComputeWettedAngle(ht_edge - B_edge);
    W[0] = ComputeWettedArea(W[1]);
  }

  return W;

}


//--------------------------------------------------------------------
// Compute wetted angle given ponded depth
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWettedAngle(double PondedDepth){

   if (PondedDepth >= pipe_diameter_) return TwoPi; //if pipe is filled wetted anlge is 2pi

   else return 2.0 * acos(1.0 - 2.0 * PondedDepth / pipe_diameter_);

}

//--------------------------------------------------------------------
// Compute wetted area given wetted angle
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWettedArea(double WettedAngle){

   if (WettedAngle >= TwoPi) return Pi * 0.25 * pipe_diameter_ * pipe_diameter_; //if pipe is filled wetter area is the full cross section

   else return pipe_diameter_ * pipe_diameter_ * 0.125 * (WettedAngle - sin(WettedAngle));

}

//--------------------------------------------------------------------
// Compute ponded depth given wetted angle
//--------------------------------------------------------------------
double PipeFlow_PK::ComputePondedDepth(double WettedAngle){

   if(WettedAngle >= TwoPi) return pipe_diameter_;

   else return pipe_diameter_ * 0.5 * (1.0 - cos(WettedAngle * 0.5)); 

}

//--------------------------------------------------------------------
// Compute wetted angle given wetted area with Newton's method
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWettedAngleNewton(double WettedArea){

   double tol = 1.e-12;
   unsigned max_iter = 1000;
   double WettedAngle;
   double PipeCrossSection = Pi * 0.25 * pipe_diameter_ * pipe_diameter_;

   if (std::fabs(WettedArea) < 1.e-12) { //cell is dry
      WettedAngle = 0.0;
   }
   else if (std::fabs(WettedArea - PipeCrossSection) < 1.e-12){ //cell is fully flooded
      WettedAngle = TwoPi;   
   }
   else { //cell is partially flooded
      unsigned iter = 0;
      if (std::fabs(WettedAngle) < 1.e-12) WettedAngle = Pi; // change initial guess to pi if was zero
       double err = WettedAngle - sin(WettedAngle) - 8.0 * WettedArea / (pipe_diameter_ * pipe_diameter_);
       while(iter < max_iter && std::fabs(err) > tol){
          WettedAngle =  WettedAngle - err / (1.0 - cos(WettedAngle));
          err = WettedAngle - sin(WettedAngle) - 8.0 * WettedArea / (pipe_diameter_ * pipe_diameter_);
          iter++;
       }
   }

   return WettedAngle;
}


//--------------------------------------------------------------------
// Initialize PK 
//--------------------------------------------------------------------
void PipeFlow_PK::Initialize(){

  ShallowWater_PK::Initialize();

  UpdateWettedAngle(); 

}


}  // namespace PipeFlow 
}  // namespace Amanzi

