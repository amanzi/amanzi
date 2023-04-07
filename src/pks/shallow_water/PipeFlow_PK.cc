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

PipeFlow_PK::PipeFlow_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln)
                  : ShallowWater_PK(pk_tree, glist, S, soln), 
                  PK(pk_tree, glist, S, soln){

  pipe_cross_section_ = Pi * 0.25 * pipe_diameter_ * pipe_diameter_;

  Manning_coeff_ = sw_list_->get<double>("Manning coefficient", 0.005);

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
// Discretization of the bed slope source term
//--------------------------------------------------------------------
std::vector<double>
PipeFlow_PK::NumericalSourceBedSlope(int c, double htc, double Bc, double Bmax, const Epetra_MultiVector& B_n, std::vector<int> bc_model, std::vector<double> bc_value_h)
{

  std::vector<double> S(3, 0.0);
  std::vector<double> W(2, 0.0);
  std::vector<double> UL(2), UR(2);

  if (std::fabs(htc - Bc) >= 1.e-15) { //cell is not dry
  
     AmanziMesh::Entity_ID_List cfaces;
     mesh_->cell_get_faces(c, &cfaces);
     double vol = mesh_->cell_volume(c);
     int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
     int orientation;

     double BGrad = 0.0;
     double OtherTermLeft = 0.0;
     double OtherTermRight = 0.0;
     double FaceAreaL = 0.0;
     double FaceAreaR = 0.0;
     double denomL = 0.0;
     double denomR = 0.0;

     for (int n = 0; n < cfaces.size(); ++n) {
        int f = cfaces[n];
        double farea = mesh_->face_area(f);
        const auto& normal = mesh_->face_normal(f, false, c, &orientation);
        if (normal[0] > 0.0){ //this identifies the j+1/2 face

           AmanziMesh::Entity_ID_List cells;
           mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
           int c1 = cells[0];
           int c2 = (cells.size() == 2) ? cells[1] : -1;
           if (c1 > ncells_owned && c2 == -1) continue;
           if (c2 > ncells_owned) std::swap(c1, c2);

           BGrad = BathymetryEdgeValue(f, B_n); //B_(j+1/2)

           W = ComputeWettedQuantitiesEdge(c1, f, htc, Bc, Bmax, B_n);

           UL[0] = W[0]; //wetted area
           UL[1] = W[1]; //wetted angle

           FaceAreaL = mesh_->face_area(f);
           denomL = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

        }

        else if (normal[0] < 0.0) { //this identifies the j-1/2 face

           BGrad -= BathymetryEdgeValue(f, B_n); //B_(j+1/2) - //B_(j-1/2)

           AmanziMesh::Entity_ID_List cells;
           mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
           int c1 = cells[0];
           int c2 = (cells.size() == 2) ? cells[1] : -1;
           if (c1 > ncells_owned && c2 == -1) continue;
           if (c2 > ncells_owned) std::swap(c1, c2);

           if (c2 == -1) {
               if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
               UR[0] = bc_value_h[f];
               UR[1] = ComputeWettedAngleNewton(bc_value_h[f]);
           } else {

             // default outflow BC
             W = ComputeWettedQuantitiesEdge(c1, f, htc, Bc, Bmax, B_n);

             UR[0] = W[0]; //wetted area
             UR[1] = W[1]; //wetted angle

           }
              denomR = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

           } else {

              W = ComputeWettedQuantitiesEdge(c2, f, htc, Bc, Bmax, B_n);

              UR[0] = W[0];
              UR[1] = W[1];

              denomR = TotalDepthEdgeValue(c2, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

           }

           FaceAreaR = mesh_->face_area(f);

        }

      }

      if ((0.0 < UL[0] && UL[0] < pipe_cross_section_)){ //flow is ventilated (free-surface)

         OtherTermLeft = 3.0 * sin(UL[1] * 0.5) - pow(sin(UL[1] * 0.5),3) - 3.0 * (UL[1] * 0.5) * cos(UL[1] * 0.5);
         OtherTermLeft = OtherTermLeft * g_ * pow(pipe_diameter_,3) / 24.0;

      }

      else if (UL[0] >= pipe_cross_section_) { //flow is pressurized

         double PressurizedHead = (celerity_ * celerity_ * (UL[0] - pipe_cross_section_)) / (g_ * pipe_cross_section_);
         OtherTermLeft = g_ * UL[0] * (PressurizedHead + sqrt(UL[0]/Pi));

      }

      if ((0.0 < UR[0] && UR[0] < pipe_cross_section_)){ //flow is ventilated (free-surface)

         OtherTermRight = 3.0 * sin(UR[1] * 0.5) - pow(sin(UR[1] * 0.5),3) - 3.0 * (UR[1] * 0.5) * cos(UR[1] * 0.5);
         OtherTermRight = OtherTermRight * g_ * pow(pipe_diameter_,3) / 24.0;

      }

      else if (UR[0] >= pipe_cross_section_) { //flow is pressurized

         double PressurizedHead = (celerity_ * celerity_ * (UR[0] - pipe_cross_section_)) / (g_ * pipe_cross_section_);
         OtherTermRight = g_ * UR[0] * (PressurizedHead + sqrt(UR[0]/Pi));

      }

      BGrad /= vol; //(B_(j+1/2) - //B_(j-1/2)) / dx

      double denom = denomL - denomR;
      S[1] = std::fabs(BGrad)<1.e-14 ? 0.0 : - (FaceAreaL * OtherTermLeft - FaceAreaR * OtherTermRight) * BGrad / denom;

 } // closes cell is not dry 

  return S;
}

//--------------------------------------------------------------------
// Newton solve to compute wetted angle given wetted area
//--------------------------------------------------------------------
void PipeFlow_PK::UpdateWettedQuantities(){ 

   auto& WettedArea_c = *S_->GetW<CV_t>(ponded_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& WettedAngle_c = *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& TotalDepth_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

   for (int c = 0; c < ncells_owned; ++c) {

       WettedAngle_c[0][c] = ComputeWettedAngleNewton(WettedArea_c[0][c]);
       TotalDepth_c[0][c] = ComputeTotalDepth(WettedArea_c[0][c], WettedAngle_c[0][c], B_c[0][c]);

   }

}


//--------------------------------------------------------------------
// Compute total depth
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeTotalDepth(double WettedArea, double WettedAngle, double Bathymetry){

   double TotalDepth = 0.0;

   if( WettedArea >= 0.0 && WettedArea < pipe_cross_section_){

      TotalDepth = ComputePondedDepth(WettedAngle) + Bathymetry;

   }

   else if ( WettedArea >= pipe_cross_section_) {

      double PressurizedHead = (celerity_ * celerity_ * (WettedArea - pipe_cross_section_)) / (g_ * pipe_cross_section_);
      TotalDepth = pipe_diameter_ + Bathymetry + PressurizedHead;

   }

   else {

      std::cout << " wetted area is negative in ComputeTotalDepth " << std::endl;
      abort();

   }

   return TotalDepth;

}


//--------------------------------------------------------------------
// Compute wetted area and wetted angle at edge location
//--------------------------------------------------------------------
std::vector<double> PipeFlow_PK::ComputeWettedQuantitiesEdge(int c, int e, double htc, double Bc, double Bmax, const Epetra_MultiVector& B_n)
{

  std::vector <double> W(2,0.0); //vector to return 

  double ht_edge = TotalDepthEdgeValue(c, e, htc, Bc, Bmax, B_n);
  double B_edge = BathymetryEdgeValue(e, B_n);
  double h_edge = std::max( (ht_edge - B_edge), 0.0);
  W[1] = ComputeWettedAngle(h_edge);
  W[0] = ComputeWettedArea(W[1]);

  return W;

}


//--------------------------------------------------------------------
// Compute wetted angle given ponded depth
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWettedAngle(double PondedDepth){

   if (PondedDepth >= pipe_diameter_) return TwoPi; //if pipe is filled wetted angle is TwoPi

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

   double tol = 1.e-15;
   unsigned max_iter = 1000;
   double WettedAngle;
   double PipeCrossSection = Pi * 0.25 * pipe_diameter_ * pipe_diameter_;

   if (std::fabs(WettedArea) < 1.e-15) { //cell is dry
      WettedAngle = 0.0;
   }
   else if (WettedArea >=  PipeCrossSection){ //cell is fully flooded
      WettedAngle = TwoPi;   
   }
   else { //cell is partially flooded
      unsigned iter = 0;
      if (std::fabs(WettedAngle) < 1.e-15) WettedAngle = Pi; // change initial guess to pi if was zero
       double err = WettedAngle - sin(WettedAngle) - 8.0 * WettedArea / (pipe_diameter_ * pipe_diameter_);
       while(iter < max_iter && std::fabs(err) > tol){
          WettedAngle =  WettedAngle - err / (1.0 - cos(WettedAngle));
          err = WettedAngle - sin(WettedAngle) - 8.0 * WettedArea / (pipe_diameter_ * pipe_diameter_);
          iter++;
       }
   }

   return WettedAngle;
}


}  // namespace PipeFlow 
}  // namespace Amanzi

