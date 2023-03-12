/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

This list contains data collected by the input parser of a higher-level spec. 

* `"used boundary condition regions`" [Array(string)] provides list of boundary regions 
  for analysis. The simulator will print number of faces and total area of these regions
  if verbosity level is equal to or above *high*.

* `"used source and sink regions`" [Array(string)] provides list of source and sink regions
  for analysis. The simulator will print number of cells and the total volume of these regions
  if verbosity level is equal to or above *high*.

* `"used observation regions`" [Array(string)] provides list of observation regions
  for analysis. The simulator will print number of faces(or cells) and the total area 
  (or volume) of these regions if verbosity level is equal to or above *high*.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="analysis">
    <Parameter name="used boundary condition regions" type="Array(string)" value="{_REG1,_REG2}"/>
    <Parameter name="used source and sink regions" type="Array(string)" value="{_REG3,_REG4}"/>
    <Parameter name="used observation regions" type="Array(string)" value="{_REG5}"/>
    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/
  
#ifndef AMANZI_INPUT_ANALYSIS_HH_
#define AMANZI_INPUT_ANALYSIS_HH_

#include <string>

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "InputAnalysis.hh"
#include "Mesh.hh"
#include "VerboseObject.hh"

namespace Amanzi {

class InputAnalysis {
 public:
  InputAnalysis(Teuchos::RCP<const AmanziMesh::Mesh> mesh, const std::string& domain)
    : mesh_(mesh), domain_(domain), vo_(NULL){};
  ~InputAnalysis()
  {
    if (vo_ != NULL) delete vo_;
  };

  // main members
  void Init(Teuchos::ParameterList& plist);
  void RegionAnalysis();
  void OutputBCs();

  // supporting members
  template <class Iterator>
  Iterator SelectUniqueEntries(Iterator first, Iterator last);

 private:
  Teuchos::ParameterList* plist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string domain_;
  VerboseObject* vo_;
};

} // namespace Amanzi

#endif
