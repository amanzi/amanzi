/* -*-  mode: c++; indent-tabs-mode: nil -*- */
// Boundary conditions base classes.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_BC_FACTORY_HH_
#define AMANZI_BC_FACTORY_HH_

/*!

In general, boundary conditions are provided in a heirarchical list by
boundary condition type, then functional form.  Boundary condition specs are
split between two types -- those which require a user-provided function
(i.e. Dirichlet data, etc) and those which do not (i.e. zero gradient
conditions).

A list of conditions might pull in both Dirichlet and Neumann data on
different regions, or use different functions on different regions.  The
following example illustrates how boundary conditions are prescribed across
the domain for a typical PK:

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="DIRICHLET_TYPE">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="DIRICHLET_FUNCTION_NAME">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
     <ParameterList name="BC east">
       <Parameter name="regions" type="Array(string)" value="{east}"/>
       <ParameterList name="DIRICHLET_FUNCTION_NAME">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="102325."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
   <ParameterList name="mass flux">
     <ParameterList name="BC north">
       <Parameter name="regions" type="Array(string)" value="{north}"/>
       <ParameterList name="outward mass flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
   <ParameterList name="zero gradient">
     <ParameterList name="BC south">
       <Parameter name="regions" type="Array(string)" value="{south}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Different PKs populate this general format with different names, replacing
DIRICHLET_TYPE and DIRICHLET_FUNCTION_NAME.
  
 */


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Point.hh"
#include "Mesh.hh"
#include "BoundaryFunction.hh"

namespace Amanzi {
namespace Functions {

class BoundaryFunctionFactory {

public:
  BoundaryFunctionFactory() {}
  BoundaryFunctionFactory(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh,
                const Teuchos::ParameterList& plist)
     : mesh_(mesh), plist_(plist) {}

  Teuchos::RCP<Functions::BoundaryFunction>
  CreateWithFunction(std::string list_name, std::string function_name) const;

  Teuchos::RCP<Functions::BoundaryFunction>
  CreateWithoutFunction(std::string list_name) const;

  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) { mesh_ = mesh; }
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() const { return mesh_; }
  void set_names(const std::string& list_name, const std::string& function_name="") {
    list_name_ = list_name;
    function_name_ = function_name;
  }
  void set_parameterlist(const Teuchos::ParameterList& plist) { plist_ = plist; }

  Teuchos::RCP<Functions::BoundaryFunction>
  Create() const {
    // ERROR CHECKING! --etc
    if (function_name_.empty()) return CreateWithoutFunction(list_name_);
    else return CreateWithFunction(list_name_, function_name_);
  } 
  
 private:

  void ProcessListWithFunction_(const Teuchos::ParameterList&,
          std::string function_name,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessListWithoutFunction_(const Teuchos::ParameterList&,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessSpecWithFunction_(const Teuchos::ParameterList&,
          std::string function_name,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessSpecWithoutFunction_(const Teuchos::ParameterList&,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string list_name_;
  std::string function_name_;
  Teuchos::ParameterList plist_;
};

}  // namespace
}  // namespace

#endif // AMANZI_BC_FACTORY_HH_
