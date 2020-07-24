/*
 Shallow water PK
 
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.
 
 Author: Svetlana Tokareva (tokareva@lanl.gov)
*/

#include <vector>

#include "PK_DomainFunctionFactory.hh"

// Amanzi::ShallowWater
#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {
    
//--------------------------------------------------------------
// Standard constructor
//--------------------------------------------------------------
ShallowWater_PK::ShallowWater_PK(Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                 const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln) :
  S_(S),
  soln_(soln),
  glist_(glist),
  passwd_("state")
{
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // Create miscellaneous lists.
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  sw_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // domain name
  domain_ = sw_list_->template get<std::string>("domain name", "surface");

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = sw_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("ShallowWater", vlist)); 
}


//--------------------------------------------------------------
// Register fields and field evaluators with the state
//--------------------------------------------------------------
void ShallowWater_PK::Setup(const Teuchos::Ptr<State>& S)
{
  // SW conservative variables: (h, hu, hv)

  mesh_ = S->GetMesh(domain_);
  dim_ = mesh_->space_dimension();

  // domain name
  domain_ = "surface";
  velocity_x_key_       = Keys::getKey(domain_, "velocity-x");
  velocity_y_key_       = Keys::getKey(domain_, "velocity-y");
  discharge_x_key_      = Keys::getKey(domain_, "discharge-x");
  discharge_y_key_      = Keys::getKey(domain_, "discharge-y");
  ponded_depth_key_     = Keys::getKey(domain_, "ponded_depth");
  total_depth_key_      = Keys::getKey(domain_, "total_depth");
  bathymetry_key_       = Keys::getKey(domain_, "bathymetry");
  velocity_x_grad_key_  = Keys::getKey(domain_, "velocity-x_grad");
  velocity_y_grad_key_  = Keys::getKey(domain_, "velocity-y_grad");
  discharge_x_grad_key_ = Keys::getKey(domain_, "discharge-x_grad");
  discharge_y_grad_key_ = Keys::getKey(domain_, "discharge-y_grad");
  total_depth_grad_key_ = Keys::getKey(domain_, "total_depth_grad");
  bathymetry_grad_key_  = Keys::getKey(domain_, "bathymetry_grad");
  myPID_                = Keys::getKey(domain_, "PID");

  //-------------------------------
  // primary fields
  //-------------------------------

  // ponded_depth_key_
  if (!S->HasField(ponded_depth_key_)) {
    S->RequireField(ponded_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // total_depth_key_
  if (!S->HasField(total_depth_key_)) {
    S->RequireField(total_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // x velocity
  if (!S->HasField(velocity_x_key_)) {
    S->RequireField(velocity_x_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // y velocity
  if (!S->HasField(velocity_y_key_)) {
    S->RequireField(velocity_y_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // x discharge
  if (!S->HasField(discharge_x_key_)) {
    S->RequireField(discharge_x_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // y discharge
  if (!S->HasField(discharge_y_key_)) {
    S->RequireField(discharge_y_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // bathymetry
  if (!S->HasField(bathymetry_key_)) {
    S->RequireField(bathymetry_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // PID
  if (!S->HasField(myPID_)) {
    S->RequireField(myPID_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  //-------------------------------
  // gradients of the primary fields
  //-------------------------------

  // total_depth_key_
  if (!S->HasField(total_depth_grad_key_)) {
    S->RequireField(total_depth_grad_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
  }

  // x velocity
  if (!S->HasField(velocity_x_grad_key_)) {
    S->RequireField(velocity_x_grad_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
  }

  // y velocity
  if (!S->HasField(velocity_y_grad_key_)) {
    S->RequireField(velocity_y_grad_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
  }

  // x discharge
  if (!S->HasField(discharge_x_grad_key_)) {
    S->RequireField(discharge_x_grad_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
  }

  // y discharge
  if (!S->HasField(discharge_y_grad_key_)) {
    S->RequireField(discharge_y_grad_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
  }

  // bathymetry
  if (!S->HasField(bathymetry_grad_key_)) {
    S->RequireField(bathymetry_grad_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
  }
}


//--------------------------------------------------------------
// Initilize internal data
//--------------------------------------------------------------
void ShallowWater_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Create BC objects
  Teuchos::RCP<ShallowWaterBoundaryFunction> bc;
  Teuchos::RCP<Teuchos::ParameterList>
      bc_list = Teuchos::rcp(new Teuchos::ParameterList(sw_list_->sublist("boundary conditions", true)));

  bcs_.clear();

  // -- velocity
  if (bc_list->isSublist("velocity")) {
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction > bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("velocity");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        bc = bc_factory.Create(spec, "no slip", AmanziMesh::CELL, Teuchos::null);
        bc->set_bc_name("velocity");
        bcs_.push_back(bc);
      }
    }
  }

  // default
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  InitializeField(S_.ptr(), passwd_, bathymetry_key_, 0.0);

  if (!S_->GetField(ponded_depth_key_, passwd_)->initialized()) {

    Epetra_MultiVector& h_vec_c = *S_->GetFieldData(ponded_depth_key_,passwd_)->ViewComponent("cell",true);
    Epetra_MultiVector& ht_vec_c = *S_->GetFieldData(total_depth_key_,passwd_)->ViewComponent("cell",true);
    Epetra_MultiVector& B_vec_c = *S_->GetFieldData(bathymetry_key_,passwd_)->ViewComponent("cell",true);

    S_->GetFieldData(ponded_depth_key_, passwd_)->PutScalar(1.0);

    for (int c = 0; c < ncells_owned; c++) {
      ht_vec_c[0][c] = h_vec_c[0][c] + B_vec_c[0][c];
    }

    S_->GetField(ponded_depth_key_, passwd_)->set_initialized();
    S_->GetField(total_depth_key_, passwd_)->set_initialized();
  }

  InitializeField(S_.ptr(), passwd_, velocity_x_key_, 0.0);
  InitializeField(S_.ptr(), passwd_, velocity_y_key_, 0.0);
  InitializeField(S_.ptr(), passwd_, discharge_x_key_, 0.0);
  InitializeField(S_.ptr(), passwd_, discharge_y_key_, 0.0);

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (!S_->GetField("surface-PID", passwd_)->initialized()) {

    Epetra_MultiVector& PID_c = *S_->GetFieldData("surface-PID",passwd_)->ViewComponent("cell",true);

    for (int c = 0; c < ncells_owned; c++) {
      PID_c[0][c] = double(MyPID);
    }

    S_->GetField("surface-PID", passwd_)->set_initialized();
  }

  // summary of initialization
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Shallow water PK was initialized." << std::endl;
  }
}


//--------------------------------------------------------------
// Upwind scheme
//--------------------------------------------------------------
bool ShallowWater_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;

  bool failed = false;

  double eps = 1.e-6;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // save a copy of primary and conservative fields
  Epetra_MultiVector& B_vec_c = *S_->GetFieldData(bathymetry_key_,passwd_)->ViewComponent("cell",true);

  Epetra_MultiVector& h_vec_c = *S_->GetFieldData(ponded_depth_key_,passwd_)->ViewComponent("cell",true);
  Epetra_MultiVector h_vec_c_tmp(h_vec_c);

  Epetra_MultiVector& ht_vec_c = *S_->GetFieldData(total_depth_key_,passwd_)->ViewComponent("cell",true);
  Epetra_MultiVector ht_vec_c_tmp(ht_vec_c);

  Epetra_MultiVector& vx_vec_c = *S_->GetFieldData(velocity_x_key_,passwd_)->ViewComponent("cell",true);
  Epetra_MultiVector vx_vec_c_tmp(vx_vec_c);

  Epetra_MultiVector& vy_vec_c = *S_->GetFieldData(velocity_y_key_,passwd_)->ViewComponent("cell",true);
  Epetra_MultiVector vy_vec_c_tmp(vy_vec_c);

  Epetra_MultiVector& qx_vec_c = *S_->GetFieldData(discharge_x_key_,passwd_)->ViewComponent("cell",true);
  Epetra_MultiVector qx_vec_c_tmp(qx_vec_c);

  Epetra_MultiVector& qy_vec_c = *S_->GetFieldData(discharge_y_key_,passwd_)->ViewComponent("cell",true);
  Epetra_MultiVector qy_vec_c_tmp(qy_vec_c);

  // distribute data to ghost cells
  S_->GetFieldData(bathymetry_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(total_depth_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(ponded_depth_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(velocity_x_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(velocity_y_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(discharge_x_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(discharge_y_key_)->ScatterMasterToGhosted("cell");

  ComputeGradients(bathymetry_key_,bathymetry_grad_key_);
  ComputeGradients(total_depth_key_,total_depth_grad_key_);
  ComputeGradients(velocity_x_key_,velocity_x_grad_key_);
  ComputeGradients(velocity_y_key_,velocity_y_grad_key_);
  ComputeGradients(discharge_x_key_,discharge_x_grad_key_);
  ComputeGradients(discharge_y_key_,discharge_y_grad_key_);

  // distribute data to ghost cells
  S_->GetFieldData(bathymetry_grad_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(total_depth_grad_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(velocity_x_grad_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(velocity_y_grad_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(discharge_x_grad_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(discharge_y_grad_key_)->ScatterMasterToGhosted("cell");

  // Shallow water equations have the form
  // U_t + F_x(U) + G_y(U) = S(U)

  std::vector<double> U_new;

  // Simplest first-order form
  // U_i^{n+1} = U_i^n - dt/vol * (F_{i+1/2}^n - F_{i-1/2}^n) + dt * S_i

  U_new.resize(3);

  for (int c = 0; c < ncells_owned; c++) {

    Amanzi::AmanziMesh::Entity_ID_List cfaces, fedges, edcells, fcells;

    mesh_->cell_get_faces(c,&cfaces);

    Amanzi::AmanziMesh::Entity_ID_List adjcells;
    mesh_->cell_get_face_adj_cells(c, Amanzi::AmanziMesh::Parallel_type::OWNED,&adjcells);

    Amanzi::AmanziGeometry::Point normal(2);
    double farea;

    // cell volume
    double vol = mesh_->cell_volume(c);

    std::vector<double> FL, FR, FNum, FNum_rot, FS;  // fluxes
    std::vector<double> S;                           // source term
    std::vector<double> UL, UR, U;                   // data for the fluxes

    UL.resize(3);
    UR.resize(3);

    FS.resize(3);

    FNum.resize(3);

    for (int i = 0; i < 3; i++) FS[i] = 0.;

    for (int f = 0; f < cfaces.size(); f++) {

      int orientation;
      normal = mesh_->face_normal(cfaces[f],false,c,&orientation);
      mesh_->face_get_cells(cfaces[f],Amanzi::AmanziMesh::Parallel_type::OWNED,&fcells);
      farea = mesh_->face_area(cfaces[f]);
      normal /= farea;

      double vn, vt;

      const AmanziGeometry::Point& xcf = mesh_->face_centroid(cfaces[f]);

      double ht_rec = Reconstruction(xcf,c,total_depth_key_,total_depth_grad_key_);
      double B_rec  = Reconstruction(xcf,c,bathymetry_key_,bathymetry_grad_key_);

      if (ht_rec < B_rec) {
        ht_rec = ht_vec_c[0][c];
        B_rec  = B_vec_c[0][c];
      }
      double h_rec = ht_rec - B_rec;

      if (h_rec < 0.) {
        std::cout << "c = " << c << std::endl;
        std::cout << "ht_rec = " << ht_rec << std::endl;
        std::cout << "B_rec  = " << B_rec << std::endl;
        std::cout << "h_rec  = " << h_rec << std::endl;
        Errors::Message msg;
        msg << "Shallow water PK: negative h.\n";
        Exceptions::amanzi_throw(msg);
      }

      double vx_rec = Reconstruction(xcf,c,velocity_x_key_,velocity_x_grad_key_);
      double vy_rec = Reconstruction(xcf,c,velocity_y_key_,velocity_y_grad_key_);
      double qx_rec = Reconstruction(xcf,c,discharge_x_key_,discharge_x_grad_key_);
      double qy_rec = Reconstruction(xcf,c,discharge_y_key_,discharge_y_grad_key_);

      vx_rec = 2.*h_rec*qx_rec/(h_rec*h_rec + fmax(h_rec*h_rec,eps*eps));
      vy_rec = 2.*h_rec*qy_rec/(h_rec*h_rec + fmax(h_rec*h_rec,eps*eps));

      vn =  vx_rec*normal[0] + vy_rec*normal[1];
      vt = -vx_rec*normal[1] + vy_rec*normal[0];

      UL[0] = h_rec;
      UL[1] = h_rec*vn;
      UL[2] = h_rec*vt;

      int cn = WhetStone::cell_get_face_adj_cell(*mesh_, c, cfaces[f]);

      if (cn == -1) {
        UR[0] = UL[0];
        UR[1] = UL[1];
        UR[2] = UL[2];
      } else {

        ht_rec = Reconstruction(xcf,cn,total_depth_key_,total_depth_grad_key_);
        B_rec  = Reconstruction(xcf,cn,bathymetry_key_,bathymetry_grad_key_);

        if (ht_rec < B_rec) {
          ht_rec = ht_vec_c[0][cn];
          B_rec  = B_vec_c[0][cn];
        }
        h_rec = ht_rec - B_rec;

        if (h_rec < 0.) {
          std::cout << "cn = " << cn << std::endl;
          std::cout << "ht_rec = " << ht_rec << std::endl;
          std::cout << "B_rec  = " << B_rec << std::endl;
          std::cout << "h_rec  = " << h_rec << std::endl;
          Errors::Message msg;
          msg << "Shallow water PK: negative h.\n";
          Exceptions::amanzi_throw(msg);
        }

        vx_rec = Reconstruction(xcf,cn,velocity_x_key_,velocity_x_grad_key_);
        vy_rec = Reconstruction(xcf,cn,velocity_y_key_,velocity_y_grad_key_);
        qx_rec = Reconstruction(xcf,cn,discharge_x_key_,discharge_x_grad_key_);
        qy_rec = Reconstruction(xcf,cn,discharge_y_key_,discharge_y_grad_key_);

        vx_rec = 2.*h_rec*qx_rec/(h_rec*h_rec + fmax(h_rec*h_rec,eps*eps));
        vy_rec = 2.*h_rec*qy_rec/(h_rec*h_rec + fmax(h_rec*h_rec,eps*eps));

        vn =  vx_rec*normal[0] + vy_rec*normal[1];
        vt = -vx_rec*normal[1] + vy_rec*normal[0];

        UR[0] = h_rec;
        UR[1] = h_rec*vn;
        UR[2] = h_rec*vt;
      }

      FNum_rot = NumFlux_x(UL,UR);

      FNum[0] = FNum_rot[0];
      FNum[1] = FNum_rot[1]*normal[0] - FNum_rot[2]*normal[1];
      FNum[2] = FNum_rot[1]*normal[1] + FNum_rot[2]*normal[0];

      for (int i = 0; i < 3; i++) {
        FS[i] += FNum[i]*farea;
      }

    } // faces

    double h, u, v, qx, qy;

    U.resize(3);

    h  = h_vec_c[0][c];
    qx = qx_vec_c[0][c];
    qy = qy_vec_c[0][c];
    u = 2.*h*qx/(h*h + fmax(h*h,eps*eps));
    v = 2.*h*qy/(h*h + fmax(h*h,eps*eps));

    U[0] = h;
    U[1] = h*u;
    U[2] = h*v;

    S = NumSrc(U,c);

    for (int i = 0; i < 3; i++) {
      U_new[i] = U[i] - dt/vol*FS[i] + dt*S[i];
    }

    h  = U_new[0];
    qx = U_new[1];
    qy = U_new[2];

    h_vec_c_tmp[0][c] = h;
    u = 2.*h*qx/(h*h + fmax(h*h,eps*eps));
    v = 2.*h*qy/(h*h + fmax(h*h,eps*eps));
    vx_vec_c_tmp[0][c] = u;
    vy_vec_c_tmp[0][c] = v;
    qx_vec_c_tmp[0][c] = h*u;
    qy_vec_c_tmp[0][c] = h*v;

    ht_vec_c_tmp[0][c] = h_vec_c_tmp[0][c] + B_vec_c[0][c];

  } // c

  h_vec_c  = h_vec_c_tmp;
  ht_vec_c = ht_vec_c_tmp;
  vx_vec_c = vx_vec_c_tmp;
  vy_vec_c = vy_vec_c_tmp;
  qx_vec_c = qx_vec_c_tmp;
  qy_vec_c = qy_vec_c_tmp;

  return failed;
}


//==========================================================================
//
// Discretization: numerical fluxes, source terms, etc
//
//==========================================================================

//--------------------------------------------------------------
// minmod function
//--------------------------------------------------------------
double minmod(double a, double b)
{
  double m;

  if (a*b > 0) {
    if (std::fabs(a) < std::fabs(b)) {
      m = a;
    } else {
      m = b;
    }
  } else {
    m = 0.;
  }

  return m;
}


//--------------------------------------------------------------
// Barth-Jespersen limiter of the gradient
//--------------------------------------------------------------
void ShallowWater_PK::BJ_lim(
    const WhetStone::DenseMatrix& grad, WhetStone::DenseMatrix& grad_lim,
    int c, const Epetra_MultiVector& field)
{
  std::vector<double> Phi_k;
  std::vector<double> u_av;
  double u_av0;
  double Phi;
  double umin, umax, uk;

  AmanziGeometry::Point xv(2);

  AmanziMesh::Entity_ID_List cfaces, cnodes, fedges, edcells, fcells;

  mesh_->cell_get_faces(c,&cfaces);
  mesh_->cell_get_nodes(c,&cnodes);

  u_av0 = field[0][c];

  u_av.resize(cfaces.size());

  for (int f = 0; f < cfaces.size(); f++) {
    int cn = WhetStone::cell_get_face_adj_cell(*mesh_, c, cfaces[f]);
    if (cn == -1) cn = c;
    u_av[f] = field[0][cn];
  }

  umin = *min_element(u_av.begin(),u_av.end());
  umax = *max_element(u_av.begin(),u_av.end());

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  Phi_k.resize(cnodes.size());

  for (int k = 0; k < cnodes.size(); k++) {

    mesh_->node_get_coordinates(cnodes[k],&xv);

    uk = u_av0 + grad(0,0)*(xv[0]-xc[0]) + grad(1,0)*(xv[1]-xc[1]);

    if (uk - u_av0 > 0.) {
      Phi_k[k] = std::min(1.,(umax-u_av0)/(uk-u_av0));
    } else {
      if (uk - u_av0 < 0.) {
        Phi_k[k] = std::min(1.,(umin-u_av0)/(uk-u_av0));
      } else {
        Phi_k[k] = 1.;
      }
    }
  } // k

  Phi = *min_element(Phi_k.begin(),Phi_k.end());

  // temporary fix to reduce the gradient
  // for better robustness on wet/dry front
  Phi = 0.1*Phi;

  grad_lim = Phi*grad;
}


//--------------------------------------------------------------
// compute gradients
//--------------------------------------------------------------
void ShallowWater_PK::ComputeGradients(const Key& field_key_, const Key& field_grad_key_)
{
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  WhetStone::DenseMatrix AtA(2,2);
  WhetStone::DenseMatrix dudx(2,1), dudx_lim(2,1);
  WhetStone::DenseMatrix Atb(2,1);

  Epetra_MultiVector& u_c = *S_->GetFieldData(field_key_,passwd_)->ViewComponent("cell",true);
  Epetra_MultiVector& gradu_c = *S_->GetFieldData(field_grad_key_,passwd_)->ViewComponent("cell",true);

  for (int c = 0; c < ncells_owned; c++) {

    Amanzi::AmanziMesh::Entity_ID_List cfaces, fedges, edcells, fcells;

    mesh_->cell_get_faces(c,&cfaces);

    Amanzi::AmanziMesh::Entity_ID_List adjcells;
    mesh_->cell_get_face_adj_cells(c, Amanzi::AmanziMesh::Parallel_type::OWNED,&adjcells);

    WhetStone::DenseMatrix A(cfaces.size(),2), At(2,cfaces.size());
    WhetStone::DenseMatrix b(cfaces.size(),1);

    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    for (int f = 0; f < cfaces.size(); f++) {

      int cn = WhetStone::cell_get_face_adj_cell(*mesh_, c, cfaces[f]);

      AmanziGeometry::Point xcn;

      if (cn == -1) {

        // formally, assign cn = c to pick correct solution
        cn = c;

        // face centroid
        const AmanziGeometry::Point& xcf = mesh_->face_centroid(cfaces[f]);
        AmanziGeometry::Point dx = xcf - xc;

        // face area
        double farea = mesh_->face_area(cfaces[f]);

        // normal
        int orientation;
        AmanziGeometry::Point normal = mesh_->face_normal(cfaces[f],false,c,&orientation);
        normal /= farea;

        double l = dx * normal;

        // reflect the centroid from the internal cell to ghost cell
        xcn = xc + 2.*l*normal;

      } else {
        xcn = mesh_->cell_centroid(cn);
      }

      A(f,0) = xcn[0] - xc[0];
      A(f,1) = xcn[1] - xc[1];

      b(f,0) = u_c[0][cn] - u_c[0][c];
    }

    At.Transpose(A);

    AtA.Multiply(At,A,false);

    Atb.Multiply(At,b,false);

    AtA.Inverse();

    dudx.Multiply(AtA,Atb,false);

    BJ_lim(dudx,dudx_lim,c,u_c);

    gradu_c[0][c] = dudx_lim(0,0);
    gradu_c[1][c] = dudx_lim(1,0);
  } // c
}


//--------------------------------------------------------------
// piecewise-linear reconstruction of the function
// from pre-computed gradients
//--------------------------------------------------------------
double ShallowWater_PK::Reconstruction(
    const AmanziGeometry::Point& xf, int c,
    const Key& field_key_, const Key& field_grad_key_)
{
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  Epetra_MultiVector& u_c = *S_->GetFieldData(field_key_,passwd_)->ViewComponent("cell",true);
  Epetra_MultiVector& grad_c = *S_->GetFieldData(field_grad_key_,passwd_)->ViewComponent("cell",true);

  return u_c[0][c] + grad_c[0][c]*(xf[0]-xc[0]) + grad_c[0][c]*(xf[1]-xc[1]);
}


//--------------------------------------------------------------
// physical flux in x-direction
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::PhysFlux_x(std::vector<double> U)
{
  std::vector<double> F;

  F.resize(3);

  double h, u, v, qx, qy, g = 9.81;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  h  = U[0];
  qx = U[1];
  qy = U[2];
  u  = 2.*h*qx/(h*h + fmax(h*h,eps*eps));
  v  = 2.*h*qy/(h*h + fmax(h*h,eps*eps));

  // Form vector of x-fluxes F(U) = (hu, hu^2 + 1/2 gh^2, huv)

  F[0] = h*u;
  F[1] = h*u*u+0.5*g*h*h;
  F[2] = h*u*v;

  return F;
}


//--------------------------------------------------------------
// physical flux in y-direction
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::PhysFlux_y(std::vector<double> U)
{
  std::vector<double> G;

  G.resize(3);

  double h, u, v, qx, qy, g = 9.81;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  h  = U[0];
  qx = U[1];
  qy = U[2];
  u  = 2.*h*qx/(h*h + fmax(h*h,eps*eps));
  v  = 2.*h*qy/(h*h + fmax(h*h,eps*eps));

  // Form vector of y-fluxes G(U) = (hv, huv, hv^2 + 1/2 gh^2)

  G[0] = h*v;
  G[1] = h*u*v;
  G[2] = h*v*v+0.5*g*h*h;

  return G;
}


//--------------------------------------------------------------
// physical source term
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::PhysSrc(std::vector<double> U)
{
  std::vector<double> S;

  S.resize(3);

  double h, u, v, qx, qy, g = 9.81;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  h  = U[0];
  qx = U[1];
  qy = U[2];
  u  = 2.*h*qx/(h*h + fmax(h*h,eps*eps));
  v  = 2.*h*qy/(h*h + fmax(h*h,eps*eps));

  // Form vector of sources Sr(U) = (0, -ghB_x, -ghB_y)

  double dBathx = 0.0, dBathy = 0.;

  S[0] = 0.;
  S[1] = -g*h*dBathx;
  S[2] = -g*h*dBathy;

  return S;
}


//--------------------------------------------------------------
// numerical flux in x-direction
// note that the SW system has a rotational invariance property
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::NumFlux_x(std::vector<double> UL,std::vector<double> UR)
{
  return NumFlux_x_Rus(UL,UR);
  // return NumFlux_x_central_upwind(UL,UR);
}


//--------------------------------------------------------------
// Rusanov numerical flux -- very simple but very diffusive
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::NumFlux_x_Rus(std::vector<double> UL,std::vector<double> UR)
{
  std::vector<double> FL, FR, F;

  double hL, uL, vL, hR, uR, vR, qxL, qyL, qxR, qyR, g = 9.81;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  hL  = UL[0];
  qxL = UL[1];
  qyL = UL[2];
  uL  = 2.*hL*qxL/(hL*hL + fmax(hL*hL,eps*eps));
  vL  = 2.*hL*qyL/(hL*hL + fmax(hL*hL,eps*eps));

  hR  = UR[0];
  qxR = UR[1];
  qyR = UR[2];
  uR  = 2.*hR*qxR/(hR*hR + fmax(hR*hR,eps*eps));
  vR  = 2.*hR*qyR/(hR*hR + fmax(hR*hR,eps*eps));

  F.resize(3);

  FL = PhysFlux_x(UL);
  FR = PhysFlux_x(UR);

  double SL, SR, Smax;

  SL = std::max(std::fabs(uL) + std::sqrt(g*hL),std::fabs(vL) + std::sqrt(g*hL));
  SR = std::max(std::fabs(uR) + std::sqrt(g*hR),std::fabs(vR) + std::sqrt(g*hR));

  Smax = std::max(SL,SR);

  for (int i = 0; i < 3; i++) {
    F[i] = 0.5*(FL[i]+FR[i]) - 0.5*Smax*(UR[i]-UL[i]);
  }

  return F;
}


//--------------------------------------------------------------
// Central-upwind numerical flux (Kurganov, Acta Numerica 2018)
//--------------------------------------------------------------
std::vector<double> ShallowWater_PK::NumFlux_x_central_upwind(std::vector<double> UL,std::vector<double> UR)
{
  std::vector<double> FL, FR, F, U_star, dU;

  double hL, uL, vL, hR, uR, vR, qxL, qyL, qxR, qyR, g = 9.81;
  double apx, amx, apy, amy;
  double ap, am;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  hL  = UL[0];
  qxL = UL[1];
  qyL = UL[2];
  uL  = 2.*hL*qxL/(hL*hL + fmax(hL*hL,eps*eps));
  vL  = 2.*hL*qyL/(hL*hL + fmax(hL*hL,eps*eps));

  hR  = UR[0];
  qxR = UR[1];
  qyR = UR[2];
  uR  = 2.*hR*qxR/(hR*hR + fmax(hR*hR,eps*eps));
  vR  = 2.*hR*qyR/(hR*hR + fmax(hR*hR,eps*eps));

  apx = std::max(std::max(std::fabs(uL)+std::sqrt(g*hL),std::fabs(uR)+std::sqrt(g*hR)),0.);
  apy = std::max(std::max(std::fabs(vL)+std::sqrt(g*hL),std::fabs(vR)+std::sqrt(g*hR)),0.);
  ap  = std::max(apx,apy);

  amx = std::min(std::min(std::fabs(uL)-std::sqrt(g*hL),std::fabs(uR)-std::sqrt(g*hR)),0.);
  amy = std::min(std::min(std::fabs(vL)-std::sqrt(g*hL),std::fabs(vR)-std::sqrt(g*hR)),0.);
  am  = std::min(amx,amy);

  F.resize(3);

  FL = PhysFlux_x(UL);
  FR = PhysFlux_x(UR);

  U_star.resize(3);

  for (int i = 0; i < 3; i++) {
    U_star[i] = ( ap*UR[i] - am*UL[i] - (FR[i] - FL[i]) ) / (ap - am + eps);
  }

  dU.resize(3);

  for (int i = 0; i < 3; i++) {
    dU[i] = minmod(UR[i]-U_star[i],U_star[i]-UL[i]);
  }

  for (int i = 0; i < 3; i++) {
    F[i] = ( ap*FL[i] - am*FR[i] + ap*am*(UR[i] - UL[i] - dU[i]) ) / (ap - am + eps);
  }

  return F;
}


//--------------------------------------------------------------------
// discretization of the source term (well-balanced for lake at rest)
//--------------------------------------------------------------------
std::vector<double> ShallowWater_PK::NumSrc(std::vector<double> U, int c)
{
  std::vector<double> S(3);

  Epetra_MultiVector& B_vec_c  = *S_->GetFieldData(bathymetry_key_,passwd_)->ViewComponent("cell",true);
  Epetra_MultiVector& ht_vec_c = *S_->GetFieldData(total_depth_key_,passwd_)->ViewComponent("cell",true);

  AmanziMesh::Entity_ID_List cfaces;

  mesh_->cell_get_faces(c,&cfaces);

  double S1, S2;

  S1 = 0.;
  S2 = 0.;

  for (int f = 0; f < cfaces.size(); f++) {

    // normal
    int orientation;
    const AmanziGeometry::Point& normal = mesh_->face_normal(cfaces[f],false,c,&orientation);

    // face centroid
    const AmanziGeometry::Point& xcf = mesh_->face_centroid(cfaces[f]);

    double ht_rec = Reconstruction(xcf,c,total_depth_key_,total_depth_grad_key_);
    double B_rec  = Reconstruction(xcf,c,bathymetry_key_,bathymetry_grad_key_);

    if (ht_rec < B_rec) {
      ht_rec = ht_vec_c[0][c];
      B_rec  = B_vec_c[0][c];
    }

    S1 += (-2.*B_rec*ht_rec + B_rec*B_rec)*normal[0];
    S2 += (-2.*B_rec*ht_rec + B_rec*B_rec)*normal[1];
  }

  S[0] = 0.;
  S[1] = 0.5*g/vol*S1;
  S[2] = 0.5*g/vol*S2;

  return S;
}


//--------------------------------------------------------------
// calculation of time step from CFL condition
//--------------------------------------------------------------
double ShallowWater_PK::get_dt()
{
  double h, u, v, g = 9.81;
  double S;
  double vol, dx;
  double dt;
  double CFL = 0.1;
  double eps = 1.e-6;

  Epetra_MultiVector& h_vec_c = *S_->GetFieldData(ponded_depth_key_,passwd_)->ViewComponent("cell");
  Epetra_MultiVector& vx_vec_c = *S_->GetFieldData(velocity_x_key_,passwd_)->ViewComponent("cell");
  Epetra_MultiVector& vy_vec_c = *S_->GetFieldData(velocity_y_key_,passwd_)->ViewComponent("cell");

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  dt = 1.e10;
  for (int c = 0; c < ncells_owned; c++) {
    h = h_vec_c[0][c];
    u = vx_vec_c[0][c];
    v = vy_vec_c[0][c];
    S = std::max(std::fabs(u) + std::sqrt(g*h),std::fabs(v) + std::sqrt(g*h)) + eps;
    vol = mesh_->cell_volume(c);
    dx = std::sqrt(vol);
    dt = std::min(dt,dx/S);
  }

  double dt_min;

  mesh_->get_comm()->MinAll(&dt, &dt_min, 1);

  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "dt = " << dt << ", dt_min = " << dt_min << std::endl;
  }

  return CFL*dt_min;
}

}  // namespace ShallowWater
}  // namespace Amanzi

