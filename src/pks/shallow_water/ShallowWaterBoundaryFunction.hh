/*
 Shallow water PK
 
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.
 
 Author: Svetlana Tokareva (tokareva@lanl.gov)
 */

#ifndef AMANZI_SHALLOW_WATER_BOUNDARY_FUNCTION_HH_
#define AMANZI_SHALLOW_WATER_BOUNDARY_FUNCTION_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "WhetStoneDefs.hh"

#include "PK_DomainFunction.hh"

namespace Amanzi {
    namespace ShallowWater {
        
        class ShallowWaterBoundaryFunction : public PK_DomainFunction {
            public:
            ShallowWaterBoundaryFunction() : bc_name_("undefined") {};
            ShallowWaterBoundaryFunction(const Teuchos::ParameterList& plist);
            
            void ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
            
            // modifiers and access
            void set_bc_name(const std::string& name) { bc_name_ = name; }
            std::string bc_name() { return bc_name_; }
            
            void set_type(WhetStone::DOF_Type type) { type_ = type; }
            WhetStone::DOF_Type type() { return type_; }
            
            private:
            std::string bc_name_;
            WhetStone::DOF_Type type_;  // type of dofs related to this bc
            
            std::vector<std::string> regions_;
        };
        
    }  // namespace ShallowWater
}  // namespace Amanzi

#endif

