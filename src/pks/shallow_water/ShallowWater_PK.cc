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
        // U_t + F_x(U) + G_y(U) = S(U)

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

      		 dx = mesh_->edge_length(0);
      		 dy = mesh_->edge_length(1);

      		 // cell volume
      		 double vol = mesh_->cell_volume(c);

             std::vector<double> FL, FR, FNum, FNum_rot, FS;  // fluxes
             std::vector<double> S;                           // source term
             std::vector<double> UL, UR, U;  			      // data for the fluxes
             
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

				 Amanzi::AmanziGeometry::Point xcf = mesh_->face_centroid(cfaces[f]);

				 double h_rec = Reconstruction(xcf[0],xcf[1],c);
//				 double h_rec = h_vec_c[0][c];

				 std::cout << "c = " << c << "/" << ncells_owned << ", h_rec = " << h_rec << std::endl;
				 std::cout << "xcf = " << xcf << std::endl;

				 vn =  vx_vec_c[0][c]*normal[0] + vy_vec_c[0][c]*normal[1];
				 vt = -vx_vec_c[0][c]*normal[1] + vy_vec_c[0][c]*normal[0];

//				 UL[0] = h_vec_c[0][c];
//				 UL[1] = h_vec_c[0][c]*vn;
//				 UL[2] = h_vec_c[0][c]*vt;

				 UL[0] = h_rec;
				 UL[1] = h_rec*vn;
				 UL[2] = h_rec*vt;

				 int cn = WhetStone::cell_get_face_adj_cell(*mesh_, c, cfaces[f]);

				 if (cn == -1) {
 					 UR[0] = UL[0];
 					 UR[1] = UL[1];
 					 UR[2] = UL[2];
				 }
				 else {
					 int c_GID = mesh_->GID(c, AmanziMesh::CELL);
 			         int cn_GID = mesh_->GID(cn, AmanziMesh::CELL);

 					 double h_rec = Reconstruction(xcf[0],xcf[1],cn);
 	//				 double h_rec = h_vec_c[0][cn];

					 vn =  vx_vec_c[0][cn]*normal[0] + vy_vec_c[0][cn]*normal[1];
					 vt = -vx_vec_c[0][cn]*normal[1] + vy_vec_c[0][cn]*normal[0];

// 					 UR[0] = h_vec_c[0][cn];
// 					 UR[1] = h_vec_c[0][cn]*vn;
// 					 UR[2] = h_vec_c[0][cn]*vt;

 					 UR[0] = h_rec;
				     UR[1] = h_rec*vn;
				     UR[2] = h_rec*vt;

				 }

				 normal[0] /= farea; normal[1] /= farea;

				 // debug
				 if (MyPID == 0) {
					 if (c == c_debug) {

					 for (int i = 0; i < 3; i++) {
						 std::cout << "UL[i] = " << UL[i] << std::endl;
					 }

					 for (int i = 0; i < 3; i++) {
						 std::cout << "UR[i] = " << UR[i] << std::endl;
					 }

					 }
				 }

				 FNum_rot = NumFlux_x(UL,UR);

				 FNum[0] = FNum_rot[0];
				 FNum[1] = FNum_rot[1]*normal[0] - FNum_rot[2]*normal[1];
				 FNum[2] = FNum_rot[1]*normal[1] + FNum_rot[2]*normal[0];

				 for (int i = 0; i < 3; i++) {
					 FS[i] += FNum[i]*farea;
				 }

             } // faces

             U.resize(3);

             U[0] = h_vec_c[0][c];
			 U[1] = h_vec_c[0][c]*vx_vec_c[0][c];
			 U[2] = h_vec_c[0][c]*vy_vec_c[0][c];
             S = NumSrc(U);

             for (int i = 0; i < 3; i++) {
            	 U_new[i] = U[i] - dt/vol*FS[i] + dt*S[i];
             }

             // debug
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
        		if (xc[0] < 3.) {
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
    
//==========================================================================
//
// Discretization: numerical fluxes, source terms, etc
//
//==========================================================================

    //--------------------------------------------------------------
    // minmod function
    //--------------------------------------------------------------
    double minmod(double a, double b) {
        double m;

        if (a*b > 0) {
            if (std::fabs(a) < std::fabs(b)) {
                m = a;
            }
            else {
                m = b;
            }
        }
        else {
            m = 0.;
        }

		return m;

    }

    //--------------------------------------------------------------
    // bottom topography
    //--------------------------------------------------------------
    double Bathymetry(double x, double y) {
        return 0.;
    }

    //--------------------------------------------------------------
    // Barth-Jespersen limiter of the gradient
    //--------------------------------------------------------------
    void ShallowWater_PK::BJ_lim(WhetStone::DenseMatrix grad, WhetStone::DenseMatrix& grad_lim, int c) {

    	std::vector<double> Phi_k;
    	std::vector<double> u_av;
    	double u_av0;
    	double Phi;
    	double umin, umax, uk;
    	double tol = 1.e-6;

    	Amanzi::AmanziGeometry::Point xc(2), xv(2);

    	Amanzi::AmanziMesh::Entity_ID_List cedges, cfaces, cnodes, fedges, edcells, fcells;

		mesh_->cell_get_edges(c,&cedges);
		mesh_->cell_get_faces(c,&cfaces,true);
		mesh_->cell_get_nodes(c,&cnodes);

		Epetra_MultiVector& h_vec_c = *S_->GetFieldData("surface-ponded_depth",passwd_)->ViewComponent("cell");

		u_av0 = h_vec_c[0][c];

		u_av.resize(cfaces.size());

		for (int f = 0; f < cfaces.size(); f++) {

			int cn = WhetStone::cell_get_face_adj_cell(*mesh_, c, cfaces[f]);

			if (cn == -1) cn = c;

			u_av[f] = h_vec_c[0][cn];

		} // f

		umin = *min_element(u_av.begin(),u_av.end());
		umax = *max_element(u_av.begin(),u_av.end());

		xc = mesh_->cell_centroid(c);

		Phi_k.resize(cnodes.size());

		for (int k = 0; k < cnodes.size(); k++) {

			mesh_->node_get_coordinates(cnodes[k],&xv);

		    uk = u_av0 + grad(0,0)*(xv[0]-xc[0]) + grad(1,0)*(xv[1]-xc[1]);

            if (uk - u_av0 > 0.) {
                Phi_k[k] = std::min(1.,(umax-u_av0)/(uk-u_av0+tol));
            }
            else {
                if (uk - u_av0 < 0.) {
                    Phi_k[k] = std::min(1.,(umin-u_av0)/(uk-u_av0+tol));
                }
                else {
                    Phi_k[k] = 1.;
                }
            }

        } // k

        Phi = *min_element(Phi_k.begin(),Phi_k.end());

        grad_lim = Phi*grad;

    }

    //--------------------------------------------------------------
	// piecewise-linear reconstruction of the function
	//--------------------------------------------------------------
    double ShallowWater_PK::Reconstruction(double x, double y, int c) {

    	// for Cartesian grids

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

		 dx = mesh_->edge_length(0);
		 dy = mesh_->edge_length(1);

		 // cell volume
		 double vol = mesh_->cell_volume(c);

		 std::vector<double> U, Un;

		 U.resize(3);
		 Un.resize(3);

		 WhetStone::DenseMatrix A(cfaces.size(),2), At(2,cfaces.size()), AtA(2,2), InvAtA(2,2);
		 WhetStone::DenseMatrix dudx(2,1), dudx_lim(2,1);
		 WhetStone::DenseMatrix b(cfaces.size(),1);
		 WhetStone::DenseMatrix Atb(2,1);

		 AmanziGeometry::Point xc = mesh_->cell_centroid(c);

		 Epetra_MultiVector& h_vec_c = *S_->GetFieldData("surface-ponded_depth",passwd_)->ViewComponent("cell");

		 for (int f = 0; f < cfaces.size(); f++) {

			 std::cout << "f = " << f << std::endl;

			 int cn = WhetStone::cell_get_face_adj_cell(*mesh_, c, cfaces[f]);

			 if (cn == -1) cn = c;

			 AmanziGeometry::Point xcn = mesh_->cell_centroid(cn);

			 std::cout << "xc = " << xc << ", xcn = " << xcn << std::endl;

			 A(f,0) = xcn[0] - xc[0];
			 A(f,1) = xcn[1] - xc[1];

			 std::cout << "A(f,0) = " << A(f,0) << ",  A(f,1) = " <<  A(f,1) << std::endl;

			 b(f,0) = h_vec_c[0][cn] -  h_vec_c[0][c];

			 std::cout << "b(f,0) = " << b(f,0) << std::endl;

		 }

		 At.Transpose(A);

		 for (int f = 0; f < cfaces.size(); f++) {
			 std::cout << "At(0,f) = " << At(0,f) << ",  At(1,f) = " <<  At(1,f) << std::endl;
		 }

		 AtA.Multiply(At,A,false);

		 for (int f = 0; f < 2; f++) {
			 std::cout << "AtA(0,f) = " << AtA(0,f) << ",  AtA(1,f) = " <<  AtA(1,f) << std::endl;
		 }

		 Atb.Multiply(At,b,false);

		 for (int f = 0; f < 2; f++) {
			 std::cout << "Atb(f,0) = " << Atb(f,0) << std::endl;
		 }

		 AtA.Inverse();

		 InvAtA = AtA;

		 for (int f = 0; f < 2; f++) {
			 std::cout << "InvAtA(0,f) = " << InvAtA(0,f) << ",  InvAtA(1,f) = " <<  InvAtA(1,f) << std::endl;
		 }

		 dudx.Multiply(InvAtA,Atb,false);

		 std::cout << "dudx     = " << dudx(0,0) << " " << dudx(1,0) << std::endl;

		 BJ_lim(dudx,dudx_lim,c);

		 std::cout << "dudx_lim = " << dudx_lim(0,0) << " " << dudx_lim(1,0) << std::endl;

		 double h_rec = h_vec_c[0][c] + dudx_lim(0,0)*(x-xc[0]) + dudx_lim(1,0)*(y-xc[1]);

         return h_rec;

    }

    //--------------------------------------------------------------
    // physical flux in x-direction
    //--------------------------------------------------------------
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
    
    //--------------------------------------------------------------
    // physical flux in y-direction
    //--------------------------------------------------------------
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

    //--------------------------------------------------------------
    // physical source term
    //--------------------------------------------------------------
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

    //--------------------------------------------------------------
    // numerical flux in x-direction
    // note that the SW system has a rotational invariance property
    //--------------------------------------------------------------
    std::vector<double> ShallowWater_PK::NumFlux_x(std::vector<double> UL,std::vector<double> UR) {

//        return NumFlux_x_Rus(UL,UR);
    	return NumFlux_x_central_upwind(UL,UR);

    }

//    //--------------------------------------------------------------
//    // numerical flux in y-direction
//    //--------------------------------------------------------------
//    std::vector<double> ShallowWater_PK::NumFlux_y(std::vector<double> UL,std::vector<double> UR) {
//        std::vector<double> GL, GR, G;
//
//        double hL, uL, vL, hR, uR, vR, g = 9.81;
//
//		// SW conservative variables: (h, hu, hv)
//		hL = UL[0];
//		uL = UL[1]/UL[0];
//		vL = UL[2]/UL[0];
//
//		hR = UR[0];
//		uR = UR[1]/UR[0];
//		vR = UR[2]/UR[0];
//
//        G.resize(3);
//
//        GL = PhysFlux_y(UL);
//        GR = PhysFlux_y(UR);
//
//        double SL, SR, Smax;
//
//        SL = std::max(uL + std::sqrt(g*hL),vL + std::sqrt(g*hL));
//        SR = std::max(uR + std::sqrt(g*hR),vR + std::sqrt(g*hR));
//
//        Smax = std::max(SL,SR);
//
//        for (int i = 0; i < 3; i++) {
//            G[i] = 0.5*(GL[i]+GR[i]) - 0.5*Smax*(UR[i]-UL[i]);
//        }
//
//        return G;
//    }

    //--------------------------------------------------------------
    // Rusanov numerical flux -- very simple but very diffusive
    //--------------------------------------------------------------
    std::vector<double> ShallowWater_PK::NumFlux_x_Rus(std::vector<double> UL,std::vector<double> UR) {
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

    //--------------------------------------------------------------
    // Central-upwind numerical flux (Kurganov, Acta Numerica 2018)
    //--------------------------------------------------------------
    std::vector<double> ShallowWater_PK::NumFlux_x_central_upwind(std::vector<double> UL,std::vector<double> UR) {
        std::vector<double> FL, FR, F, U_star, dU;

        double hL, uL, vL, hR, uR, vR, g = 9.81;
        double apx, amx, apy, amy;
        double ap, am;

		// SW conservative variables: (h, hu, hv)
		hL = UL[0];
		uL = UL[1]/UL[0];
		vL = UL[2]/UL[0];

		hR = UR[0];
		uR = UR[1]/UR[0];
		vR = UR[2]/UR[0];

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
        	U_star[i] = ( ap*UR[i] - am*UL[i] - (FR[i] - FL[i]) ) / (ap - am);
        }

        dU.resize(3);

        for (int i = 0; i < 3; i++) {
        	dU[i] = minmod(UR[i]-U_star[i],U_star[i]-UL[i]);
        }

        for (int i = 0; i < 3; i++) {
            F[i] = ( ap*FL[i] - am*FR[i] + ap*am*(UR[i] - UL[i] - dU[i]) ) / (ap - am);
        }

        return F;
    }

    //--------------------------------------------------------------
    // discretization of the source term (not well-balanced)
    //--------------------------------------------------------------
    std::vector<double> ShallowWater_PK::NumSrc(std::vector<double> U) {
        std::vector<double> S;

        double g = 9.81;

        Epetra_MultiVector& h_vec_c = *S_->GetFieldData("surface-ponded_depth",passwd_)->ViewComponent("cell");

        double BRx, BLx, BRy, BLy;

        S.resize(3);

        S = PhysSrc(U);

        BRx = 0.;
        BLx = 0.;
        BRy = 0.;
        BLy = 0.;

        int c = 1; // for now

        double h = h_vec_c[0][c];

        Amanzi::AmanziMesh::Entity_ID_List cedges;

        mesh_->cell_get_edges(c,&cedges);
        double dx = mesh_->edge_length(0);
        double dy = mesh_->edge_length(1);

        S[0] = 0.;
        S[1] = -g*h*(BRx - BLx)/dx;
        S[2] = -g*h*(BRy - BLy)/dy;

        return S;

    }

    //--------------------------------------------------------------
    // calculation of time step from CFL condition
    //--------------------------------------------------------------
    double ShallowWater_PK::get_dt() {

    	double h, u, v, g = 9.81;
		double S, Smax;
		double vol, dx, dx_min;
		double dt;
		double CFL = 0.45;

    	Epetra_MultiVector& h_vec_c = *S_->GetFieldData("surface-ponded_depth",passwd_)->ViewComponent("cell");
		Epetra_MultiVector& vx_vec_c = *S_->GetFieldData("surface-velocity-x",passwd_)->ViewComponent("cell");
		Epetra_MultiVector& vy_vec_c = *S_->GetFieldData("surface-velocity-y",passwd_)->ViewComponent("cell");

		int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

		dt = 1.e10;
		for (int c = 0; c < ncells_owned; c++) {
			h = h_vec_c[0][c];
			u = vx_vec_c[0][c];
			v = vy_vec_c[0][c];
			S = std::max(std::fabs(u) + std::sqrt(g*h),std::fabs(v) + std::sqrt(g*h));
			vol = mesh_->cell_volume(c);
			dx = std::sqrt(vol);
			dt = std::min(dt,dx/S);
		}

		Comm_ptr_type comm = Amanzi::getDefaultComm();
		int MyPID = comm->MyPID();

		double dt_min;

		mesh_->get_comm()->MinAll(&dt, &dt_min, 1);

		dt_min = 0.0005;

//		std::cout << "dt = " << dt << ", dt_min = " << dt_min << std::endl;

    	return CFL*dt_min;
    }

}  // namespace ShallowWater
}  // namespace Amanzi

