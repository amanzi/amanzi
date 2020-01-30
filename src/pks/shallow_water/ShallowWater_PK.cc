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
        
        // Create miscaleneous lists.
        Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
        sw_list_ = Teuchos::sublist(pk_list, pk_name, true);
        
        // domain name
        domain_ = sw_list_->template get<std::string>("domain name", "surface");
        
        vo_ = Teuchos::null;
    }
    
    bool ShallowWater_PK::AdvanceStep(double t_old, double t_new, bool reinit)
    {


    	double dt = t_new - t_old;

    	bool failed = false;

    	// save a copy of primary and conservative fields
    	CompositeVector ponded_depth_copy(*S_->GetFieldData("surface-ponded_depth", passwd_));
    	CompositeVector velocity_x_copy(*S_->GetFieldData("surface-velocity-x", passwd_));
    	CompositeVector velocity_y_copy(*S_->GetFieldData("surface-velocity-y", passwd_));

    	Epetra_MultiVector& h_vec_c = *S_->GetFieldData("surface-ponded_depth",passwd_)->ViewComponent("cell");
    	Epetra_MultiVector& h_vec_c_old(h_vec_c);


        // Shallow water equations have the form
         // U_t + F_x(U) + G_y(U) = Sr(U)

         std::vector<double> U, F, G, Sr, U_new;

         double h, u, v, g;

         // Form vector of conservative variables U = (h,hu,hv)

         U.push_back(h);
         U.push_back(h*u);
         U.push_back(h*v);

         // Form vector of x-fluxes F(U) = (hu, hu^2 + 1/2 gh^2, huv)

         F.push_back(h*u);
         F.push_back(h*u*u+0.5*g*h*h);
         F.push_back(h*u*v);

         // Form vector of y-fluxes G(U) = (hv, huv, hv^2 + 1/2 gh^2)

         G.push_back(h*v);
         G.push_back(h*u*v);
         G.push_back(h*v*v+0.5*g*h*h);

         // Form vector of sources Sr(U) = (0, -ghB_x, -ghB_y)

         // Simplest first-order form
         // U_i^{n+1} = U_i^n + dt/vol * (F_{i+1/2}^n - F_{i-1/2}^n)

         int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
         std::cout << "ncells_owned = " << ncells_owned << std::endl;

         std::cout << "dt = " << dt << std::endl;
         std::cout << "t_new = " << t_new << std::endl;

         U_new.resize(3);

      	 for (int c = 0; c < ncells_owned; c++) {
      		// cell volume
      		 double vol = mesh_->cell_volume(c);
     		 h_vec_c[0][c] = h_vec_c_old[0][c] + dt/vol*c*t_new;
      	 }

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

    }
    
    void ShallowWater_PK::Initialize(const Teuchos::Ptr<State>& S) {
        
        if (!S_->GetField("surface-ponded_depth", passwd_)->initialized()) {
            S_->GetFieldData("surface-ponded_depth", passwd_)->PutScalar(2.0);
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
        
    }
    
}  // namespace ShallowWater
}  // namespace Amanzi

