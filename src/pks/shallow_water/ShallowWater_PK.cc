/*
 Shallow water PK
 
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.
 
 Author: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include <vector>

// Amanzi::ShallowWater
#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {
    
    ShallowWater_PK::ShallowWater_PK(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& glist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln) :
    S_(S),
    soln_(soln)
    {
        std::string pk_name = pk_tree.name();
        auto found = pk_name.rfind("->");
        if (found != std::string::npos) pk_name.erase(0, found + 2);
        
        // Create miscellaneous lists.
        Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
        sw_list_ = Teuchos::sublist(pk_list, pk_name, true);
        
        // domain name
        domain_ = sw_list_->template get<std::string>("domain name", "surface");
        
        vo_ = Teuchos::null;
    }
    
    bool ShallowWater_PK::AdvanceStep(double t_old, double t_new, bool reinit)
    {

    	Comm_ptr_type comm = Amanzi::getDefaultComm();
    	int MyPID = comm->MyPID();

    	double dt = t_new - t_old;

    	bool failed = false;

    	// save a copy of primary and conservative fields
    	CompositeVector ponded_depth_tmp(*S_->GetFieldData("surface-ponded_depth", passwd_));
    	CompositeVector velocity_x_tmp(*S_->GetFieldData("surface-velocity-x", passwd_));
    	CompositeVector velocity_y_tmp(*S_->GetFieldData("surface-velocity-y", passwd_));

    	Epetra_MultiVector& h_vec_c = *S_->GetFieldData("surface-ponded_depth",passwd_)->ViewComponent("cell");
    	Epetra_MultiVector h_vec_c_tmp(h_vec_c);
        
        Epetra_MultiVector& vx_vec_c = *S_->GetFieldData("surface-velocity-x",passwd_)->ViewComponent("cell");
        Epetra_MultiVector vx_vec_c_tmp(vx_vec_c);
        
        Epetra_MultiVector& vy_vec_c = *S_->GetFieldData("surface-velocity-y",passwd_)->ViewComponent("cell");
        Epetra_MultiVector vy_vec_c_tmp(vy_vec_c);

        // distribute data to ghost cells
        S_->GetFieldData("surface-ponded_depth")->ScatterMasterToGhosted("cell");
        S_->GetFieldData("surface-velocity-x")->ScatterMasterToGhosted("cell");
        S_->GetFieldData("surface-velocity-y")->ScatterMasterToGhosted("cell");

        h_vec_c_tmp  = h_vec_c;
        vx_vec_c_tmp = vx_vec_c;
        vy_vec_c_tmp = vy_vec_c;

        // Shallow water equations have the form
        // U_t + F_x(U) + G_y(U) = Sr(U)

         std::vector<double> U, U_new;

         // Simplest first-order form
         // U_i^{n+1} = U_i^n - dt/vol * (F_{i+1/2}^n - F_{i-1/2}^n) + dt * S_i

         int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
         std::cout << "ncells_owned = " << ncells_owned << std::endl;

         std::cout << "dt = " << dt << std::endl;
         std::cout << "t_new = " << t_new << std::endl;

         U_new.resize(3);

         int c_debug = -10;

      	 for (int c = 0; c < ncells_owned; c++) {

      		 // mesh sizes
      		 double dx, dy;

      		 Amanzi::AmanziMesh::Entity_ID_List cedges, cfaces, fedges, edcells, fcells;

      		 mesh_->cell_get_edges(c,&cedges);
      		 mesh_->cell_get_faces(c,&cfaces,true);

      		 Amanzi::AmanziMesh::Entity_ID_List adjcells;
			 mesh_->cell_get_face_adj_cells(c, Amanzi::AmanziMesh::Parallel_type::OWNED,&adjcells);
			 unsigned int nadj = adjcells.size();

      		 Amanzi::AmanziGeometry::Point evec(2), normal(2);
      		 double farea;

//      		 std::cout << "------------------" << std::endl;
//      		 std::cout << "c = " << c << std::endl;

      		 dx = mesh_->edge_length(0);
      		 dy = mesh_->edge_length(1);

      		 // cell volume
      		 double vol = mesh_->cell_volume(c);

             std::vector<double> FL, FR, FNum, FNum_rot, FS;  // fluxes
             std::vector<double> S;       // source term
             std::vector<double> UL, UR, U;  // data for the fluxes
             
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

				 double vn, vt;

				 vn =  vx_vec_c[0][c]*normal[0] + vy_vec_c[0][c]*normal[1];
				 vt = -vx_vec_c[0][c]*normal[1] + vy_vec_c[0][c]*normal[0];

				 UL[0] = h_vec_c[0][c];
				 UL[1] = h_vec_c[0][c]*vn; //vx_vec_c[0][c];
				 UL[2] = h_vec_c[0][c]*vt; //vy_vec_c[0][c];

				 int cn = WhetStone::cell_get_face_adj_cell(*mesh_, c, cfaces[f]);

//				 std::cout << "face = " << f << std::endl;

				 if (cn == -1) {
 					 UR[0] = UL[0];
 					 UR[1] = UL[1];
 					 UR[2] = UL[2];
// 					 int c_GID = mesh_->GID(c, AmanziMesh::CELL);
//// 					 int cn_GID = mesh_->GID(cn, AmanziMesh::CELL);
// 					 std::cout << "MyPID = " << MyPID << ", c = " << c << ", cn = " << cn << std::endl;
// 					 std::cout << "MyPID = " << MyPID << ", c_GID = " << c << ", cn_GID = " << cn << std::endl;

				 }
				 else {
					 int c_GID = mesh_->GID(c, AmanziMesh::CELL);
 			         int cn_GID = mesh_->GID(cn, AmanziMesh::CELL);
// 			         cn = cn_GID;
// 			         if (MyPID == 0) {
//						 std::cout << "MyPID = " << MyPID << ", c = " << c << ", cn = " << cn << std::endl;
//						 std::cout << "MyPID = " << MyPID << ", c_GID = " << c << ", cn_GID = " << cn << std::endl;
// 			         }
					 vn =  vx_vec_c[0][cn]*normal[0] + vy_vec_c[0][cn]*normal[1];
					 vt = -vx_vec_c[0][cn]*normal[1] + vy_vec_c[0][cn]*normal[0];
 					 UR[0] = h_vec_c[0][cn];
 					 UR[1] = h_vec_c[0][cn]*vn; //vx_vec_c[0][cn];
 					 UR[2] = h_vec_c[0][cn]*vt; //vy_vec_c[0][cn];
				 }

				 normal[0] /= farea; normal[1] /= farea;

				 if (MyPID == 0) {
				 if (c == c_debug) {

			     std::cout << "MyPID = " << MyPID << ", c = " << c << ", cn = " << cn << std::endl;

				 std::cout << "normal = " << normal << std::endl;

				 for (int i = 0; i < 3; i++) {
					 std::cout << "UL[i] = " << UL[i] << std::endl;
				 }

				 for (int i = 0; i < 3; i++) {
					 std::cout << "UR[i] = " << UR[i] << std::endl;
				 }

				 }
				 }

				 FNum_rot = NumFlux_x(UL,UR);

//				 for (int i = 0; i < 3; i++) {
//					 std::cout << "FNum_rot[i] = " << FNum_rot[i] << std::endl;
//				 }

				 FNum[0] = FNum_rot[0];
				 FNum[1] = FNum_rot[1]*normal[0] - FNum_rot[2]*normal[1];
				 FNum[2] = FNum_rot[1]*normal[1] + FNum_rot[2]*normal[0];

//				 for (int i = 0; i < 3; i++) {
//					 std::cout << "FNum[i] = " << FNum[i] << std::endl;
//				 }

//				 std::cout << "farea = " << farea << std::endl;

				 for (int i = 0; i < 3; i++) {
					 FS[i] += FNum[i]*farea;
				 }

             } // faces

//             UR[0] = h_vec_c[0][c];
//             UR[1] = h_vec_c[0][c]*vx_vec_c[0][c];
//             UR[2] = h_vec_c[0][c]*vy_vec_c[0][c];
//             if (c == 0) {
////            	 UL[0] = 10.;
////            	 UL[1] = 10.;
////            	 UL[2] = 0.;
//            	 UL[0] = UR[0];
//            	 UL[1] = UR[1];
//            	 UL[2] = UR[2];
//             }
//             else {
//				 UL[0] = h_vec_c[0][c-1];
//				 UL[1] = h_vec_c[0][c-1]*vx_vec_c[0][c-1];
//				 UL[2] = h_vec_c[0][c-1]*vy_vec_c[0][c-1];
//             }
//             FL = NumFlux_x(UL,UR);
//
////             for (int i = 0; i < 3; i++) {
////				 std::cout << "UL[i] = " << UL[i] << std::endl;
////			 }
////
////			 for (int i = 0; i < 3; i++) {
////				 std::cout << "UR[i] = " << UR[i] << std::endl;
////			 }
////
////             for (int i = 0; i < 3; i++) {
////				 std::cout << "FL[i] = " << FL[i] << std::endl;
////			 }
//
//             UL[0] = h_vec_c[0][c];
//			 UL[1] = h_vec_c[0][c]*vx_vec_c[0][c];
//			 UL[2] = h_vec_c[0][c]*vy_vec_c[0][c];
//			 if (c == ncells_owned-1) {
//				 UR[0] = UL[0];
//				 UR[1] = UL[1];
//				 UR[2] = UL[2];
//			 }
//			 else {
//				 UR[0] = h_vec_c[0][c+1];
//				 UR[1] = h_vec_c[0][c+1]*vx_vec_c[0][c+1];
//				 UR[2] = h_vec_c[0][c+1]*vy_vec_c[0][c+1];
//			 }
//             FR = NumFlux_x(UL,UR);
//
////             for (int i = 0; i < 3; i++) {
////				 std::cout << "UL[i] = " << UL[i] << std::endl;
////			 }
////
////			 for (int i = 0; i < 3; i++) {
////				 std::cout << "UR[i] = " << UR[i] << std::endl;
////			 }
////
////             for (int i = 0; i < 3; i++) {
////				 std::cout << "FR[i] = " << FR[i] << std::endl;
////			 }

             U.resize(3);

             U[0] = h_vec_c[0][c];
			 U[1] = h_vec_c[0][c]*vx_vec_c[0][c];
			 U[2] = h_vec_c[0][c]*vy_vec_c[0][c];
             S = NumSrc(U);

             for (int i = 0; i < 3; i++) {
//            	 U_new[i] = U[i] - dt/dx*(FR[i]-FL[i]) + dt*S[i];
            	 U_new[i] = U[i] - dt/vol*FS[i] + dt*S[i];
//                 std::cout << "FS[i] = " << dt/vol*FS[i] << ", FR[i]-FL[i] = " << dt/dx*(FR[i]-FL[i]) << std::endl;
             }

             if (MyPID == 0) {
            	 if (c == c_debug) {
            		 for (int i = 0; i < 3; i++) {
            			 std::cout << "U_new[i] = " << U_new[i] << std::endl;
            		 }
            	 }
             }

             h_vec_c_tmp[0][c]  = U_new[0];
             vx_vec_c_tmp[0][c] = U_new[1]/U_new[0];
             vy_vec_c_tmp[0][c] = U_new[2]/U_new[0];
      	 }

      	h_vec_c  = h_vec_c_tmp;
      	vx_vec_c = vx_vec_c_tmp;
      	vy_vec_c = vy_vec_c_tmp;

//       	 *S_->GetFieldData("surface-ponded_depth",passwd_)->ViewComponent("cell") = h_vec_c;

    	 return failed;

    }

    void ShallowWater_PK::Setup(const Teuchos::Ptr<State>& S)
    {

        // SW conservative variables: (h, hu, hv)
        
        passwd_ = "state";  // owner's password
        
        mesh_ = S->GetMesh(domain_);
        dim_ = mesh_->space_dimension();
        
        // domain name
        domain_ = "surface";
        pressure_key_       = Keys::getKey(domain_, "pressure");
        velocity_x_key_     = Keys::getKey(domain_, "velocity-x");
        velocity_y_key_     = Keys::getKey(domain_, "velocity-y");
        ponded_depth_key_   = Keys::getKey(domain_, "ponded_depth");
        myPID_  		    = Keys::getKey(domain_, "PID");

        // primary fields
        
        // -- pressure
        if (!S->HasField(pressure_key_)) {
            S->RequireField(pressure_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
            ->SetComponent("cell", AmanziMesh::CELL, 1);

//            Teuchos::ParameterList elist;
//            elist.set<std::string>("evaluator name", "pressure");
//            pressure_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
//            S->SetFieldEvaluator("pressure", pressure_eval_);
        }
        
        
        // ponded_depth_key_
        if (!S->HasField(ponded_depth_key_)) {
            S->RequireField(ponded_depth_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
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

        // PID
	    if (!S->HasField(myPID_)) {
		    S->RequireField(myPID_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
		    ->SetComponent("cell", AmanziMesh::CELL, 1);
	    }

    }
    
    void ShallowWater_PK::Initialize(const Teuchos::Ptr<State>& S) {

    	int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
        
        if (!S_->GetField("surface-ponded_depth", passwd_)->initialized()) {

//            S_->GetFieldData("surface-ponded_depth", passwd_)->PutScalar(2.0);

        	Epetra_MultiVector& h_vec_c = *S_->GetFieldData("surface-ponded_depth",passwd_)->ViewComponent("cell");

        	for (int c = 0; c < ncells_owned; c++) {
        		AmanziGeometry::Point xc = mesh_->cell_centroid(c);
        		if (xc[0] < 0.3) {
        		   h_vec_c[0][c] = 4.;
        		}
			    else {
			       h_vec_c[0][c] = 1.;
			    }
        	}

        	S_->GetField("surface-ponded_depth", passwd_)->set_initialized();
        }
        
        if (!S_->GetField("surface-velocity-x", passwd_)->initialized()) {
            S_->GetFieldData("surface-velocity-x", passwd_)->PutScalar(0.0);
            S_->GetField("surface-velocity-x", passwd_)->set_initialized();
        }
        
        if (!S_->GetField("surface-velocity-y", passwd_)->initialized()) {
            S_->GetFieldData("surface-velocity-y", passwd_)->PutScalar(0.0);
            S_->GetField("surface-velocity-y", passwd_)->set_initialized();
        }
        
        Comm_ptr_type comm = Amanzi::getDefaultComm();
        int MyPID = comm->MyPID();

        if (!S_->GetField("surface-PID", passwd_)->initialized()) {

//            S_->GetFieldData("surface-ponded_depth", passwd_)->PutScalar(2.0);

        	Epetra_MultiVector& PID_c = *S_->GetFieldData("surface-PID",passwd_)->ViewComponent("cell");

        	for (int c = 0; c < ncells_owned; c++) {
			    PID_c[0][c] = MyPID;
        	}

        	S_->GetField("surface-PID", passwd_)->set_initialized();
        }

    }
    
    std::vector<double> ShallowWater_PK::PhysFlux_x(std::vector<double> U) {
        std::vector<double> F;
        
        F.resize(3);
        
        double h, u, v, g = 9.81;
        
        // SW conservative variables: (h, hu, hv)
        h = U[0];
        u = U[1]/U[0];
        v = U[2]/U[0];
        
        // Form vector of x-fluxes F(U) = (hu, hu^2 + 1/2 gh^2, huv)
        
        F[0] = h*u;
        F[1] = h*u*u+0.5*g*h*h;
        F[2] = h*u*v;

        return F;
    }
    
    std::vector<double> ShallowWater_PK::PhysFlux_y(std::vector<double> U) {
        std::vector<double> G;

        G.resize(3);

        double h, u, v, g = 9.81;

        // SW conservative variables: (h, hu, hv)
        h = U[0];
        u = U[1]/U[0];
        v = U[2]/U[0];

        // Form vector of y-fluxes G(U) = (hv, huv, hv^2 + 1/2 gh^2)

        G[0] = h*v;
        G[1] = h*u*v;
        G[2] = h*v*v+0.5*g*h*h;

        return G;
    }

    std::vector<double> ShallowWater_PK::PhysSrc(std::vector<double> U) {
            std::vector<double> S;

            S.resize(3);

            double h, u, v, g = 9.81;

            // SW conservative variables: (h, hu, hv)
            h = U[0];
            u = U[1]/U[0];
            v = U[2]/U[0];

            // Form vector of sources Sr(U) = (0, -ghB_x, -ghB_y)

            double dBathx = 0., dBathy = 0.;

            S[0] = 0.;
            S[1] = g*h*dBathx;
            S[2] = g*h*dBathy;

            return S;
        }

    std::vector<double> ShallowWater_PK::NumFlux_x(std::vector<double> UL,std::vector<double> UR) {
        std::vector<double> FL, FR, F;
        
        double hL, uL, vL, hR, uR, vR, g = 9.81;

		// SW conservative variables: (h, hu, hv)
		hL = UL[0];
		uL = UL[1]/UL[0];
		vL = UL[2]/UL[0];

		hR = UR[0];
		uR = UR[1]/UR[0];
		vR = UR[2]/UR[0];

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
    
    std::vector<double> ShallowWater_PK::NumFlux_xn(std::vector<double> UL,std::vector<double> UR, Amanzi::AmanziGeometry::Point normal) {
        std::vector<double> FL, FR, F;

        double hL, uL, vL, hR, uR, vR, g = 9.81;

		// SW conservative variables: (h, hu, hv)
		hL = UL[0];
		uL = UL[1]/UL[0];
		vL = UL[2]/UL[0];

		hR = UR[0];
		uR = UR[1]/UR[0];
		vR = UR[2]/UR[0];

        F.resize(3);

        FL = PhysFlux_x(UL);
        FR = PhysFlux_x(UR);

        double SL, SR, Smax;

        SL = std::max(uL + std::sqrt(g*hL),vL + std::sqrt(g*hL));
        SR = std::max(uR + std::sqrt(g*hR),vR + std::sqrt(g*hR));

        Smax = std::max(SL,SR);

        for (int i = 0; i < 3; i++) {
            F[i] = 0.5*(FL[i]+FR[i]) - 0.5*Smax*(UR[i]-UL[i]);
        }

        return F;
    }

    std::vector<double> ShallowWater_PK::NumFlux_y(std::vector<double> UL,std::vector<double> UR) {
        std::vector<double> GL, GR, G;

        double hL, uL, vL, hR, uR, vR, g = 9.81;

		// SW conservative variables: (h, hu, hv)
		hL = UL[0];
		uL = UL[1]/UL[0];
		vL = UL[2]/UL[0];

		hR = UR[0];
		uR = UR[1]/UR[0];
		vR = UR[2]/UR[0];

        G.resize(3);

        GL = PhysFlux_y(UL);
        GR = PhysFlux_y(UR);

        double SL, SR, Smax;

        SL = std::max(uL + std::sqrt(g*hL),vL + std::sqrt(g*hL));
        SR = std::max(uR + std::sqrt(g*hR),vR + std::sqrt(g*hR));

        Smax = std::max(SL,SR);

        for (int i = 0; i < 3; i++) {
            G[i] = 0.5*(GL[i]+GR[i]) - 0.5*Smax*(UR[i]-UL[i]);
        }

        return G;
    }

    std::vector<double> ShallowWater_PK::NumSrc(std::vector<double> U) {
        std::vector<double> S;

        S.resize(3);

        S = PhysSrc(U);

        return S;

        }

    double ShallowWater_PK::get_dt() {
    	return 0.0005;
    }

}  // namespace ShallowWater
}  // namespace Amanzi

