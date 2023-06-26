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
#include "WaterDepthEvaluator.hh"
#include "PressureHeadEvaluator.hh"
#include "NumericalFluxFactory.hh"
#include "PK_DomainFunctionFactory.hh"

namespace Amanzi {
namespace ShallowWater {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

PipeFlow_PK::PipeFlow_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln)
                  : ShallowWater_PK(pk_tree, glist, S, soln), 
                  PK(pk_tree, glist, S, soln){

  pipe_cross_section_ = Pi * 0.25 * pipe_diameter_ * pipe_diameter_;

  Manning_coeff_ = sw_list_->get<double>("Manning coefficient", 0.005);

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = sw_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("PipeFlow", vlist));
  
}

//--------------------------------------------------------------
// Register fields and field evaluators with the state
// Conservative variables: (A, Au)
//--------------------------------------------------------------
void PipeFlow_PK::Setup()
{
 
   ShallowWater_PK::Setup();

   // -- water depth
   if (!S_->HasRecord(water_depth_key_)) {
      S_->Require<CV_t, CVS_t>(water_depth_key_, Tags::DEFAULT, water_depth_key_)
         .SetMesh(mesh_) ->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist(water_depth_key_);
      elist.set<std::string>("my key", water_depth_key_).set<std::string>("tag", Tags::DEFAULT.get())
         .set<double>("pipe diameter", pipe_diameter_);
      auto eval = Teuchos::rcp(new WaterDepthEvaluator(elist));
      S_->SetEvaluator(water_depth_key_, Tags::DEFAULT, eval);
   }

   // -- pressure head
   if (!S_->HasRecord(pressure_head_key_)) {
      S_->Require<CV_t, CVS_t>(pressure_head_key_, Tags::DEFAULT, pressure_head_key_)
         .SetMesh(mesh_) ->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist(pressure_head_key_);
      elist.set<std::string>("my key", pressure_head_key_).set<std::string>("tag", Tags::DEFAULT.get())
         .set<double>("pipe diameter", pipe_diameter_)
         .set<double>("celerity", celerity_);
      auto eval = Teuchos::rcp(new PressureHeadEvaluator(elist));
      S_->SetEvaluator(pressure_head_key_, Tags::DEFAULT, eval);
   }
   
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
PipeFlow_PK::NumericalSourceBedSlope(int c, double htc, double Bc, double Bmax, const Epetra_MultiVector& B_n, 
                                     std::vector<int> bc_model, std::vector<double> bc_value_h)
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
        const auto& normal = mesh_->face_normal(f, false, c, &orientation);
        if (normal[0] > 0.0){ //this identifies the j+1/2 face

           AmanziMesh::Entity_ID_List cells;
           mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
           int c1 = cells[0];
           int c2 = (cells.size() == 2) ? cells[1] : -1;
           if (c1 > ncells_owned && c2 == -1) continue;
           if (c2 > ncells_owned) std::swap(c1, c2);

           BGrad = BathymetryEdgeValue(f, B_n); //B_(j+1/2)

           W = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n);

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
             W = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n);

             UR[0] = W[0]; //wetted area
             UR[1] = W[1]; //wetted angle

           }
              denomR = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

           } else {

              W = ComputeFieldsOnEdge(c2, f, htc, Bc, Bmax, B_n);

              UR[0] = W[0];
              UR[1] = W[1];

              denomR = TotalDepthEdgeValue(c2, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

           }

           FaceAreaR = mesh_->face_area(f);

        }

      }

      OtherTermLeft = ComputeHydrostaticPressureForce(UL);
      OtherTermRight = ComputeHydrostaticPressureForce(UR);

      BGrad /= vol; //(B_(j+1/2) - B_(j-1/2)) / cellVol

      double denom = denomL - denomR;
      S[1] = std::fabs(BGrad)<1.e-14 ? 0.0 : - (FaceAreaL * OtherTermLeft - FaceAreaR * OtherTermRight) * BGrad / denom;

 } // closes cell is not dry 

  return S;
}

//--------------------------------------------------------------------
// Compute pressure head
//--------------------------------------------------------------------
double PipeFlow_PK::ComputePressureHead (double WettedArea){

   return (celerity_ * celerity_ * (WettedArea - pipe_cross_section_)) / (g_ * pipe_cross_section_);

}

//--------------------------------------------------------------------
// Compute hydrostatic pressure force
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeHydrostaticPressureForce (std::vector<double> SolArray){

      double I = 0.0;

      if ((0.0 < SolArray[0] && SolArray[0] < pipe_cross_section_)){ //flow is ventilated (free-surface)

         I = 3.0 * sin(SolArray[1] * 0.5) - pow(sin(SolArray[1] * 0.5),3) - 3.0 * (SolArray[1] * 0.5) * cos(SolArray[1] * 0.5);
         I = I * g_ * pow(pipe_diameter_,3) / 24.0;

      }

      else if (SolArray[0] >= pipe_cross_section_) { //flow is pressurized

         I = g_ * SolArray[0] * (ComputePressureHead(SolArray[0]) + sqrt(SolArray[0]/Pi));

      }

      return I;

}

//--------------------------------------------------------------------
// Update wetted angle and total depth for pipe model
//--------------------------------------------------------------------
void PipeFlow_PK::UpdateSecondaryFields(){ 

   auto& WettedArea_c = *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& WettedAngle_c = *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& TotalDepth_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

   for (int c = 0; c < ncells_owned; ++c) {

       WettedAngle_c[0][c] = ComputeWettedAngleNewton(WettedArea_c[0][c]);
       TotalDepth_c[0][c] = ComputeTotalDepth(WettedArea_c[0][c], B_c[0][c], WettedAngle_c[0][c]);

   }

}


//--------------------------------------------------------------------
// Compute total depth
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeTotalDepth(double PrimaryVar, double Bathymetry, double WettedAngle){

   double TotalDepth = 0.0;

   if( PrimaryVar >= 0.0 && PrimaryVar < pipe_cross_section_){

      TotalDepth = ComputeWaterDepth(WettedAngle) + Bathymetry;

   }

   else if ( PrimaryVar >= pipe_cross_section_) {

      TotalDepth = pipe_diameter_ + Bathymetry + ComputePressureHead(PrimaryVar);

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
std::vector<double> PipeFlow_PK::ComputeFieldsOnEdge(int c, int e, double htc, double Bc, double Bmax, const Epetra_MultiVector& B_n)
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
// Compute wetted angle given water depth
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWettedAngle(double WaterDepth){

   if (WaterDepth >= pipe_diameter_) return TwoPi; //if pipe is filled wetted angle is TwoPi

   else return 2.0 * acos(1.0 - 2.0 * WaterDepth / pipe_diameter_);

}

//--------------------------------------------------------------------
// Compute wetted area given wetted angle
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWettedArea(double WettedAngle){

   if (WettedAngle >= TwoPi) return Pi * 0.25 * pipe_diameter_ * pipe_diameter_; //if pipe is filled, wetted area is the full cross section

   else return pipe_diameter_ * pipe_diameter_ * 0.125 * (WettedAngle - sin(WettedAngle));

}

//--------------------------------------------------------------------
// Compute water depth given wetted angle
//--------------------------------------------------------------------
double PipeFlow_PK::ComputeWaterDepth(double WettedAngle){

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

//--------------------------------------------------------------
// Setup Primary Variable Keys
//--------------------------------------------------------------
void PipeFlow_PK::SetupPrimaryVariableKeys(){

  primary_variable_key_ = Keys::getKey(domain_, "wetted_area");
  prev_primary_variable_key_ = Keys::getKey(domain_, "prev_wetted_area");

}

//--------------------------------------------------------------
// Setup Extra Evaluators Keys
//--------------------------------------------------------------
void PipeFlow_PK::SetupExtraEvaluatorsKeys(){

  water_depth_key_ = Keys::getKey(domain_, "water_depth");
  pressure_head_key_ = Keys::getKey(domain_, "pressure_head");

}

//--------------------------------------------------------------
// Scatter Master To Ghosted Extra Evaluators 
//--------------------------------------------------------------
void PipeFlow_PK::ScatterMasterToGhostedExtraEvaluators(){

     S_->Get<CV_t>(water_depth_key_).ScatterMasterToGhosted("cell");
     S_->Get<CV_t>(pressure_head_key_).ScatterMasterToGhosted("cell");

}

//--------------------------------------------------------------
// Update Extra Evaluators
//--------------------------------------------------------------
void PipeFlow_PK::UpdateExtraEvaluators(){

     Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(
       S_->GetEvaluatorPtr(wetted_angle_key_, Tags::DEFAULT))->SetChanged();

     S_->GetEvaluator(water_depth_key_).Update(*S_, passwd_);
     S_->GetEvaluator(pressure_head_key_).Update(*S_, passwd_);
}

//--------------------------------------------------------------
// Set Primary Variable BC
//--------------------------------------------------------------
void PipeFlow_PK::SetPrimaryVariableBC(Teuchos::RCP<Teuchos::ParameterList> &bc_list){

  // -- wetted area BC
  if (bc_list->isSublist("wetted area")) {
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("wetted area");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        Teuchos::RCP<ShallowWaterBoundaryFunction> bc = bc_factory.Create(spec, "wetted area", AmanziMesh::NODE, Teuchos::null);
        bc->set_bc_name("wetted area");
        bc->set_type(WhetStone::DOF_Type::SCALAR);
        PushBackBC(bc);
      }
    }
  }

}

//--------------------------------------------------------------
// Initialize Variables
//--------------------------------------------------------------
void PipeFlow_PK::InitializeFields(){

     int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

     auto& WettedArea_c = *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
     auto& WettedAngle_c = *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
     auto& WaterDepth_c = *S_->GetW<CV_t>(water_depth_key_, Tags::DEFAULT, water_depth_key_).ViewComponent("cell");
     auto& PressureHead_c = *S_->GetW<CV_t>(pressure_head_key_, Tags::DEFAULT, pressure_head_key_).ViewComponent("cell");
     auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
     auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

     for (int c = 0; c < ncells_owned; c++) {

        double maxDepth = B_c[0][c] + pipe_diameter_;
        if (ht_c[0][c] >= maxDepth){ // cell is pressurized

            double PipeCrossSection = Pi * 0.25 * pipe_diameter_ * pipe_diameter_;
            PressureHead_c[0][c] = ht_c[0][c] - pipe_diameter_ - B_c[0][c];
            WettedArea_c[0][c] = (g_ * PipeCrossSection * PressureHead_c[0][c]) / (celerity_ * celerity_) + PipeCrossSection;
            WettedAngle_c[0][c] = TwoPi;
            WaterDepth_c[0][c] = pipe_diameter_;
        }
        else if ((std::fabs(ht_c[0][c] - B_c[0][c]) < 1.e-15) || (ht_c[0][c] < B_c[0][c])){ //cell is dry

            WettedArea_c[0][c] = 0.0;
            WettedAngle_c[0][c] = 0.0;
            WaterDepth_c[0][c] = 0.0;

        }

        else if (ht_c[0][c] < maxDepth && B_c[0][c] < ht_c[0][c]) { //cell is ventilated

           WaterDepth_c[0][c] = ht_c[0][c] - B_c[0][c];
           WettedAngle_c[0][c] = ComputeWettedAngle(WaterDepth_c[0][c]);
           WettedArea_c[0][c] = ComputeWettedArea(WettedAngle_c[0][c]);

        }

     }

     if (S_->GetRecord(discharge_key_).initialized()) {

        auto& q_c = *S_->GetW<CV_t>(discharge_key_, Tags::DEFAULT, discharge_key_).ViewComponent("cell", true);
        auto& u_c = *S_->GetW<CV_t>(velocity_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

        for (int c = 0; c < ncells_owned; c++) {
            for (int i = 0; i < 2; ++i) {
             if(std::fabs(WettedArea_c[0][c]) < 1.0e-12){
                u_c[i][c] = 0.0;
             }
             else{
                u_c[i][c] = q_c[i][c] / WettedArea_c[0][c];
             }
            }
        }

        S_->GetRecordW(velocity_key_, Tags::DEFAULT, passwd_).set_initialized();

     }

     S_->GetRecordW(primary_variable_key_, Tags::DEFAULT, passwd_).set_initialized();
     S_->GetRecordW(wetted_angle_key_, Tags::DEFAULT, passwd_).set_initialized();
     S_->GetRecordW(water_depth_key_, Tags::DEFAULT, water_depth_key_).set_initialized();
     S_->GetRecordW(pressure_head_key_, Tags::DEFAULT, pressure_head_key_).set_initialized();

}

//--------------------------------------------------------------
// Compute external forcing on cells
//--------------------------------------------------------------
void PipeFlow_PK::ComputeExternalForcingOnCells(std::vector<double> &forcing){

     for (int i = 0; i < srcs_.size(); ++i) {
         for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
             int c = it->first;
             double dx;
             GetDx_(c, dx);
             forcing[c] = it->second[0] / dx; // [m^2 / s] for pipe
         }
     }

}

}  // namespace PipeFlow 
}  // namespace Amanzi

