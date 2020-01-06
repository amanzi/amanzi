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
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln) :
        S_(S),
        soln_(soln)
        {
            // Expect a single global list with sublist Flow
            glist_ = Teuchos::rcp(new Teuchos::ParameterList(*global_list));
            
            ti_list_ = glist_->sublist("cycle driver").sublist("time intervals").sublist("TI 0");
            
        }
        
        bool ShallowWater_PK::AdvanceStep(double t_old, double t_new, bool reinit)
        {
            bool failed = false;
            
            if ((step_count + 2)%3 == 0) {
                failed = true;
                dummy_dt = 0.8*dummy_dt;
                std::cout<<"Step failed\n";
            }
            else {
                failed = false;
                dummy_dt = 1.2*dummy_dt;
                std::cout<<"Step succeed. New time "<<t_new<<"\n";
            }
            
            double dt = t_new - t_old;
            step_count++;
            
            return failed;
        }
        
    }  // namespace ShallowWater
}  // namespace Amanzi

