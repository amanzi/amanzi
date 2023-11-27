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

   // -- direction
   if (!S_->HasRecord(direction_key_)) {

      S_->Require<CV_t, CVS_t>(direction_key_, Tags::DEFAULT, direction_key_)
         .SetMesh(mesh_)
         ->SetGhosted(true)
         ->SetComponent("cell", AmanziMesh::CELL, 2);

      if(direction_key_.empty()){
         AddDefaultPrimaryEvaluator_(direction_key_);
      }
      else{
        S_->RequireEvaluator(direction_key_, Tags::DEFAULT);  
      }

   }  

}

//--------------------------------------------------------------------
// Discretization of the friction source term
//--------------------------------------------------------------------
double PipeFlow_PK::NumericalSourceFriction(double h, double qx, double qy, double WettedAngle, int component)
{

  double S1 = 0.0;
  double num = 0.0;

  if (std::fabs(h) > 1.e-10) { //we have to raise this to the power of 7/3 below so the tolerance needs to be stricter
     if (WettedAngle >= 0.0){  
        double WettedPerimeter = 0.5 * pipe_diameter_ * WettedAngle;
        num = - g_ * Manning_coeff_ * Manning_coeff_ * pow(WettedPerimeter, 4.0/3.0) * std::fabs(qx) * qx;
        double denom = pow( h, 7.0/3.0);
        S1 = num / denom;
     }
     else{ //junction
        if(component == 0){ 
           num = - g_ * Manning_coeff_ * Manning_coeff_ * qx * sqrt(qx*qx + qy*qy);
        }
        else{
           num = - g_ * Manning_coeff_ * Manning_coeff_ * qy * sqrt(qx*qx + qy*qy);
        }
        double denom = pow( h, 7.0/3.0);
        S1 = num / denom;
     }
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
  std::vector<double> V_rec(2, 0.0);
  std::vector<double> UL(2), UR(2);

  if (std::fabs(htc - Bc) >= 1.e-15) { //cell is not dry

     int orientation;
     AmanziMesh::Entity_ID_List cfaces;
     mesh_->cell_get_faces(c, &cfaces);
     double vol = mesh_->cell_volume(c);
     int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
     auto& dir_c = *S_->GetW<CV_t>(direction_key_, Tags::DEFAULT, direction_key_).ViewComponent("cell", true);

     double BGrad = 0.0;
     double OtherTermLeft = 0.0;
     double OtherTermRight = 0.0;
     double FaceAreaL = 0.0;
     double FaceAreaR = 0.0;
     double denomL = 0.0;
     double denomR = 0.0;

     for (int n = 0; n < cfaces.size(); ++n) {
        int f = cfaces[n];
        AmanziGeometry::Point normal = mesh_->face_normal(f, false, c, &orientation);
        ProjectNormalOntoMeshDirection(c, normal);

        if (normal[0] > 1.e-08){ //this identifies the j+1/2 face

           AmanziMesh::Entity_ID_List cells;
           mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
           int c1 = cells[0];
           int c2 = (cells.size() == 2) ? cells[1] : -1;
           if (c1 > ncells_owned && c2 == -1) continue;
           if (c2 > ncells_owned) std::swap(c1, c2);

           BGrad += BathymetryEdgeValue(f, B_n); //B_(j+1/2)

           V_rec = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n);

           UL[0] = V_rec[0]; //wetted area
           UL[1] = V_rec[1]; //wetted angle

           FaceAreaL = mesh_->face_area(f);
           denomL = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

        }

        else if (normal[0] < -1.e-08) { //this identifies the j-1/2 face

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
             V_rec = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n);

             UR[0] = V_rec[0]; //wetted area
             UR[1] = V_rec[1]; //wetted angle

           }
              denomR = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

           } else {

              V_rec = ComputeFieldsOnEdge(c2, f, htc, Bc, Bmax, B_n);

              UR[0] = V_rec[0];
              UR[1] = V_rec[1];

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
// Discretization of the bed slope source term for junction
//--------------------------------------------------------------------
std::vector<double>
PipeFlow_PK::NumericalSourceBedSlope(int c, double htc, double Bc, double Bmax, const Epetra_MultiVector& B_n) 
{

  std::vector<double> S(3, 0.0);
  std::vector<double> V_rec(2, 0.0);
  std::vector<double> UL(2), UR(2);

  if (std::fabs(htc - Bc) >= 1.e-15) { //cell is not dry

     AmanziMesh::Entity_ID_List cfaces;
     mesh_->cell_get_faces(c, &cfaces);
     double vol = mesh_->cell_volume(c);
     int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

     std::vector<int> xMaxFace(cfaces.size(),0);
     std::vector<int> xMinFace(cfaces.size(),0);
     std::vector<int> yMaxFace(cfaces.size(),0);
     std::vector<int> yMinFace(cfaces.size(),0);

     double tmp = -1.0e+10;
     int indexToSave = 0;
     for (int n = 0; n < cfaces.size(); ++n) {
        const auto& xf = mesh_->face_centroid(cfaces[n]);
        if(xf[0] >= tmp) {
            tmp=xf[0];
            indexToSave=n;
        }
     }
     xMaxFace[indexToSave]=1;

     tmp = 1.e+10;
     indexToSave = 0;
     for (int n = 0; n < cfaces.size(); ++n) {
        const auto& xf = mesh_->face_centroid(cfaces[n]);
        if(xf[0] <= tmp) {
            tmp=xf[0];
            indexToSave=n;
        }
     }
     xMinFace[indexToSave]=1;

     tmp = -1.0e+10;
     indexToSave = 0;
     for (int n = 0; n < cfaces.size(); ++n) {
        const auto& xf = mesh_->face_centroid(cfaces[n]);
        if(xf[1] >= tmp) {
            tmp=xf[1];
            indexToSave=n;
        }
     }
     yMaxFace[indexToSave]=1;

     tmp = 1.e+10;
     indexToSave = 0;
     for (int n = 0; n < cfaces.size(); ++n) {
        const auto& xf = mesh_->face_centroid(cfaces[n]);
        if(xf[1] <= tmp) {
            tmp=xf[1];
            indexToSave=n;
        }
     }
     yMinFace[indexToSave]=1;

     double BGrad = 0.0;
     double OtherTermLeft = 0.0;
     double OtherTermRight = 0.0;
     double FaceAreaL = 0.0;
     double FaceAreaR = 0.0;
     double denomL = 0.0;
     double denomR = 0.0;
     int c1, c2;
     for (int n = 0; n < cfaces.size(); ++n) {
        int f = cfaces[n];
        const auto& xf = mesh_->face_centroid(f);
        AmanziMesh::Entity_ID_List cells;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        c1 = cells[0];
        c2 = (cells.size() == 2) ? cells[1] : -1;
        if (c1 > ncells_owned && c2 == -1) continue;
        if (c2 > ncells_owned) std::swap(c1, c2);

        if (xMaxFace[n] == 1 && c2 != -1){ //this identifies the j+1/2 face

           c1 = cells[0];
           c2 = (cells.size() == 2) ? cells[1] : -1;
           if (c1 > ncells_owned && c2 == -1) continue;
           if (c2 > ncells_owned) std::swap(c1, c2);

           BGrad += BathymetryEdgeValue(f, B_n); //B_(j+1/2)

           V_rec = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n);

           UL[0] = V_rec[0]; //wetted area
           UL[1] = V_rec[1]; //wetted angle

           denomL = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

        }

        else if (xMinFace[n] == 1 && c2 != -1) { //this identifies the j-1/2 face

           BGrad -= BathymetryEdgeValue(f, B_n); //B_(j+1/2) - //B_(j-1/2)

           V_rec = ComputeFieldsOnEdge(c2, f, htc, Bc, Bmax, B_n);

           UR[0] = V_rec[0];
           UR[1] = V_rec[1];

           denomR = TotalDepthEdgeValue(c2, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

        }

    }

     OtherTermLeft = ComputeHydrostaticPressureForce(UL);
     OtherTermRight = ComputeHydrostaticPressureForce(UR);

     BGrad /= vol; //(B_(j+1/2) - B_(j-1/2)) / cellVol

     double denom = denomL - denomR;
     S[1] = std::fabs(BGrad)<1.e-14 ? 0.0 : - (OtherTermLeft - OtherTermRight) * BGrad / denom;

     BGrad = 0.0;
     OtherTermLeft = 0.0;
     OtherTermRight = 0.0;
     denomL = 0.0;
     denomR = 0.0;

     for (int n = 0; n < cfaces.size(); ++n) {
        int f = cfaces[n];
        const auto& xf = mesh_->face_centroid(f);
        AmanziMesh::Entity_ID_List cells;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        c1 = cells[0];
        c2 = (cells.size() == 2) ? cells[1] : -1;
        if (c1 > ncells_owned && c2 == -1) continue;
        if (c2 > ncells_owned) std::swap(c1, c2);

        if (yMaxFace[n] == 1 && c2 != -1){ //this identifies the j+1/2 face

           BGrad += BathymetryEdgeValue(f, B_n); //B_(j+1/2)

           V_rec = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n);

           UL[0] = V_rec[0]; //wetted area
           UL[1] = V_rec[1]; //wetted angle

           denomL = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

        }

        else if (yMinFace[n] == 1 && c2 !=-1) { //this identifies the j-1/2 face

           BGrad -= BathymetryEdgeValue(f, B_n); //B_(j+1/2) - //B_(j-1/2)

           V_rec = ComputeFieldsOnEdge(c2, f, htc, Bc, Bmax, B_n);

           UR[0] = V_rec[0];
           UR[1] = V_rec[1];

           denomR = TotalDepthEdgeValue(c2, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

        }

      }

      OtherTermLeft = ComputeHydrostaticPressureForce(UL);
      OtherTermRight = ComputeHydrostaticPressureForce(UR);

      BGrad /= vol; //(B_(j+1/2) - B_(j-1/2)) / cellVol

      denom = denomL - denomR;
      S[2] = std::fabs(BGrad)<1.e-14 ? 0.0 : - (OtherTermLeft - OtherTermRight) * BGrad / denom;

 } // closes cell is not dry

  return S;
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

   auto& PrimaryVar_c = *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& WettedAngle_c = *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& TotalDepth_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
   auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);

   for (int c = 0; c < model_cells_owned_.size(); ++c) {
       int cell = model_cells_owned_[c];
       WettedAngle_c[0][cell] = ComputeWettedAngleNewton(PrimaryVar_c[0][cell]);
       TotalDepth_c[0][cell] = ComputeTotalDepth(PrimaryVar_c[0][cell], B_c[0][cell], WettedAngle_c[0][cell]);
   }
   for (int c = 0; c < junction_cells_owned_.size(); ++c) {
       int cell = junction_cells_owned_[c];
       WettedAngle_c[0][cell] = -1.0;
       // NOTE: the primary variable for junctions store the water depth
       TotalDepth_c[0][cell] = ShallowWater_PK::ComputeTotalDepth(PrimaryVar_c[0][cell], 
                                                   B_c[0][cell], WettedAngle_c[0][cell]);
   }

   Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(
       S_->GetEvaluatorPtr(wetted_angle_key_, Tags::DEFAULT))->SetChanged();

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

  std::vector <double> V_rec(2,0.0); //vector to return

  double ht_edge = TotalDepthEdgeValue(c, e, htc, Bc, Bmax, B_n);
  double B_edge = BathymetryEdgeValue(e, B_n);
  double h_edge = std::max( (ht_edge - B_edge), 0.0);

  V_rec[1] = ComputeWettedAngle(h_edge);
  V_rec[0] = ComputeWettedArea(V_rec[1]);

  return V_rec;

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
  direction_key_ = sw_list_->get<std::string>("direction key", "");

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
     auto& PrimaryVar_c = *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
     auto& WettedAngle_c = *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
     auto& WaterDepth_c = *S_->GetW<CV_t>(water_depth_key_, Tags::DEFAULT, water_depth_key_).ViewComponent("cell");
     auto& PressureHead_c = *S_->GetW<CV_t>(pressure_head_key_, Tags::DEFAULT, pressure_head_key_).ViewComponent("cell");
     auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
     auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

     for (int cell = 0; cell < ncells_owned; cell++) {

        double maxDepth = B_c[0][cell] + pipe_diameter_;
        if (ht_c[0][cell] >= maxDepth){ // cell is pressurized

            double PipeCrossSection = Pi * 0.25 * pipe_diameter_ * pipe_diameter_;
            PressureHead_c[0][cell] = ht_c[0][cell] - pipe_diameter_ - B_c[0][cell];
            PrimaryVar_c[0][cell] = (g_ * PipeCrossSection * PressureHead_c[0][cell]) / (celerity_ * celerity_) + PipeCrossSection;
            WettedAngle_c[0][cell] = TwoPi;
            WaterDepth_c[0][cell] = pipe_diameter_;
        }
        else if ((std::fabs(ht_c[0][cell] - B_c[0][cell]) < 1.e-15) || (ht_c[0][cell] < B_c[0][cell])){ //cell is dry

            PrimaryVar_c[0][cell] = 0.0;
            WettedAngle_c[0][cell] = 0.0;
            WaterDepth_c[0][cell] = 0.0;

        }

        else if (ht_c[0][cell] < maxDepth && B_c[0][cell] < ht_c[0][cell]) { //cell is ventilated

           WaterDepth_c[0][cell] = ht_c[0][cell] - B_c[0][cell];
           WettedAngle_c[0][cell] = ComputeWettedAngle(WaterDepth_c[0][cell]);
           PrimaryVar_c[0][cell] = ComputeWettedArea(WettedAngle_c[0][cell]);

        }

     }

     S_->GetRecordW(primary_variable_key_, Tags::DEFAULT, passwd_).set_initialized();
     S_->GetRecordW(wetted_angle_key_, Tags::DEFAULT, passwd_).set_initialized();
     S_->GetRecordW(water_depth_key_, Tags::DEFAULT, water_depth_key_).set_initialized();
     S_->GetRecordW(pressure_head_key_, Tags::DEFAULT, pressure_head_key_).set_initialized();

}

//---------------------------------------------------------------
// Initialize cell array of pipe cells (model cells) and junction
//---------------------------------------------------------------
void PipeFlow_PK::ComputeCellArrays(){

    if(!cellArraysInitDone_){
       int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
       int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
       auto& dir_c = *S_->GetW<CV_t>(direction_key_, Tags::DEFAULT, direction_key_).ViewComponent("cell", true);
       auto& PrimaryVar_c = *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
       auto& WettedAngle_c = *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
       auto& WaterDepth_c = *S_->GetW<CV_t>(water_depth_key_, Tags::DEFAULT, water_depth_key_).ViewComponent("cell");
       auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
       auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
       // write pipe direction if not
       // specified from file
       // default is positive x-axis
       if(direction_key_.empty()){
          for (int c = 0; c < ncells_wghost; c++) {
             dir_c[0][c] = 1.0;
             dir_c[1][c] = 0.0;
          }
       }

       junction_cells_owned_.resize(0);
       model_cells_owned_.resize(0);
       junction_cells_wghost_.resize(0);
       model_cells_wghost_.resize(0);

       for (int c = 0; c < ncells_owned; c++) {
          // both components of pipe direction zero
          // is code for junction cell
          if(IsJunction(c)){ 
             junction_cells_owned_.push_back(c); 
          }
          else {
             model_cells_owned_.push_back(c);
          }
       }
    
       for (int c = 0; c < ncells_wghost; c++) {
          // both components of pipe direction zero
          // is code for junction cell
          if(IsJunction(c)){
             junction_cells_wghost_.push_back(c);
          }
          else {
             model_cells_wghost_.push_back(c);
          }
       }

       // here we are overwriting the initialization for 
       // junction cells. This should have been done within
       // InitializeFields() but at init we cannot read
       // the dir_c field from file so we need to do it here

       for (int c = 0; c < junction_cells_owned_.size(); c++) {
           int cell = junction_cells_owned_[c];
           PrimaryVar_c[0][cell] = ht_c[0][cell] - B_c[0][cell];
           WettedAngle_c[0][cell] = - 1.0;
           WaterDepth_c[0][cell] = ht_c[0][cell] - B_c[0][cell];
       }

       cellArraysInitDone_ = true;

    }   

}

//--------------------------------------------------------------
// Checks if a pipe cell is a junction
//--------------------------------------------------------------
bool PipeFlow_PK::IsJunction(const int &cell)
{
    bool isJunction = 0;

    auto& dir_c = *S_->GetW<CV_t>(direction_key_, Tags::DEFAULT, direction_key_).ViewComponent("cell", true);
    // both components of pipe direction equal to zero is the definition of junction cell
    if(std::fabs(dir_c[0][cell]) < 1.e-14 && std::fabs(dir_c[1][cell]) < 1.e-14){ 
       isJunction = 1;
    }

    return isJunction;

}

//--------------------------------------------------------------
// Compute dx
//--------------------------------------------------------------
void PipeFlow_PK::GetDx(const int &cell, double &dx)
{
  dx = 0.0;
  AmanziMesh::Entity_ID_List cell_faces;
  mesh_->cell_get_faces(cell, &cell_faces);
  AmanziGeometry::Point x1;
  AmanziGeometry::Point x2;
  auto& dir_c = *S_->GetW<CV_t>(direction_key_, Tags::DEFAULT, direction_key_).ViewComponent("cell", true);
  for (int n = 0; n < cell_faces.size(); ++n) {
     int f = cell_faces[n];
     int orient;
     AmanziGeometry::Point f_normal = mesh_->face_normal(f, false, cell, &orient);
     ProjectNormalOntoMeshDirection(cell, f_normal);
     if (f_normal[0] > 1.e-08){
        x1 = mesh_->face_centroid(f);
     }
     else if (f_normal[0] < -1.e-08){
        x2 = mesh_->face_centroid(f);
     }
  }
  for (int iSize =0; iSize < x1.dim(); iSize++){
     dx += (x1[iSize] - x2[iSize]) * (x1[iSize] - x2[iSize]);
  }
  dx = sqrt(dx);

}

//--------------------------------------------------------------
// Compute external forcing on cells
//--------------------------------------------------------------
void PipeFlow_PK::ComputeExternalForcingOnCells(std::vector<double> &forcing){

     for (int i = 0; i < srcs_.size(); ++i) {
         for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
             int c = it->first;
             double dx;
             GetDx(c, dx);
             forcing[c] = it->second[0] / dx; // [m^2 / s] for pipe
         }
     }

}

//--------------------------------------------------------------
// Project normal onto mesh direction
//--------------------------------------------------------------
void PipeFlow_PK::ProjectNormalOntoMeshDirection(int c, AmanziGeometry::Point &normal){

   //the mesh direction is the pipe direction
   auto& dir_c = *S_->GetW<CV_t>(direction_key_, Tags::DEFAULT, direction_key_).ViewComponent("cell", true);
   // first find the angle between pipe direction and x-axis
   double angle = acos(dir_c[0][c]/sqrt(dir_c[0][c]*dir_c[0][c]+dir_c[1][c]*dir_c[1][c]));
   // then rotate the unit vectors of the reference frame according to this angle
   std::vector<double> e1(2);
   std::vector<double> e2(2);
   double e1Tmp1 = 1;
   double e1Tmp2 = 0;
   double e2Tmp1 = 0;
   double e2Tmp2 = 1;
   e1[0] = e1Tmp1 * cos(angle) - e1Tmp2 * sin(angle);
   e1[1] = e1Tmp1 * sin(angle) + e1Tmp2 * cos(angle);
   e2[0] = e2Tmp1 * cos(angle) - e2Tmp2 * sin(angle);
   e2[1] = e2Tmp1 * sin(angle) + e2Tmp2 * cos(angle);
   // finally, project the normal on rotated frame
   double nTmp1 = normal[0];
   double nTmp2 = normal[1];
   normal[0] = nTmp1 * e1[0] + nTmp2 * e1[1];
   normal[1] = nTmp1 * e2[0] + nTmp2 * e2[1];

}

}  // namespace PipeFlow 
}  // namespace Amanzi

