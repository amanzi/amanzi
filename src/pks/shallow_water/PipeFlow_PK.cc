/*
  Pipe Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Pipe flow model inherited from shallow water model. 
*/

#include "PipeFlow_PK.hh"
#include "NumericalFluxFactory.hh"

namespace Amanzi {
namespace ShallowWater {

using CV_t = CompositeVector;

PipeFlow_PK::PipeFlow_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln)
                  : ShallowWater_PK(pk_tree, glist, S, soln), 
                  PK(pk_tree, glist, S, soln){

  pipe_diameter_ = sw_list_->get<double>("pipe diameter", 1.0);

  Manning_coeff_ = sw_list_->get<double>("Manning coefficient", 0.005);

  celerity_ = sw_list_->get<double>("celerity", 100); // m/s

}

//--------------------------------------------------------------------
// Discretization of the friction source term
//--------------------------------------------------------------------
double PipeFlow_PK::NumericalSourceFriction(double h, double qx, double WettedAngle)
{

  double S1 = 0.0;
  if (fabs(h) > 1.e-10) { //we have to raise this to the power of 7/3 below so the tolerance needs to be stricter
     double WettedPerimeter = 0.5 * pipe_diameter_ * WettedAngle;
     double num = - g_ * Manning_coeff_ * Manning_coeff_ * pow(WettedPerimeter, 4.0/3.0) * fabs(qx) * qx;
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

   double tol = 1.e-12;
   unsigned max_iter = 1000;

   //TODO: replace the body of the for loop below with calls of ComputeWettedAngleIterative

   for (int c = 0; c < ncells_owned; ++c) {

      if (fabs(h_c[0][c]) < 1.e-12) {
         WettedAngle_c[0][c] = 0.0;
      }
      else {
         unsigned iter = 0;
         if (fabs(WettedAngle_c[0][c]) < 1.e-12) WettedAngle_c[0][c] = 3.14159265359; // change initial guess to pi if was zero
         double err = WettedAngle_c[0][c] - sin(WettedAngle_c[0][c]) - 8.0 * h_c[0][c] / (pipe_diameter_ * pipe_diameter_);
         while(iter < max_iter && fabs(err) > tol){
            WettedAngle_c[0][c] =  WettedAngle_c[0][c] - err / (1.0 - cos(WettedAngle_c[0][c]));
            err = WettedAngle_c[0][c] - sin(WettedAngle_c[0][c]) - 8.0 * h_c[0][c] / (pipe_diameter_ * pipe_diameter_);
            iter++;
         }
      }

   }

}

//--------------------------------------------------------------------
// Compute wetted angle given ponded depth
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWettedAngle(double PondedDepth){

   double WettedAngle = 2.0 * acos(1.0 - 2.0 * PondedDepth / pipe_diameter_);

   return WettedAngle;

}

//--------------------------------------------------------------------
// Compute wetted area given wetted angle
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWettedArea(double WettedAngle){

   double WettedArea = pipe_diameter_ * pipe_diameter_ * 0.25 * (WettedAngle - sin(WettedAngle));

   return WettedArea;

}

//--------------------------------------------------------------------
// Compute ponded depth given wetted angle
//--------------------------------------------------------------------
double PipeFlow_PK::ComputePondedDepth(double WettedAngle){

   double PondedDepth = pipe_diameter_ * 0.5 * (1.0 - cos(WettedAngle * 0.5)); 

   return PondedDepth;

}

//--------------------------------------------------------------------
// Compute wetted angle given wetted area with Newton's method
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWettedAngleNewton(double WettedArea){

   double tol = 1.e-12;
   unsigned max_iter = 1000;
   double WettedAngle;

   if (fabs(WettedArea) < 1.e-12) {
         WettedAngle = 0.0;
   }
   else {
      unsigned iter = 0;
      if (fabs(WettedAngle) < 1.e-12) WettedAngle = 3.14159265359; // change initial guess to pi if was zero
       double err = WettedAngle - sin(WettedAngle) - 8.0 * WettedArea / (pipe_diameter_ * pipe_diameter_);
       while(iter < max_iter && fabs(err) > tol){
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


}  // namespace ShallowWater
}  // namespace Amanzi

