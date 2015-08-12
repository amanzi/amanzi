/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson (jnjohnson@lbl.gov)
*/

#ifndef AMANZI_INPUT_CONVERTER_STRUCTURED_HH_
#define AMANZI_INPUT_CONVERTER_STRUCTURED_HH_

#include "ParmParse.H"
#include "InputConverter.hh"

namespace Amanzi {
namespace AmanziInput {

class InputConverterS : public InputConverter {
 public:

  InputConverterS();
  ~InputConverterS();

  // main members
  void Translate();

 private:

  void ParseUnits();
  void ParseDefinitions();
  void ParseExecutionControls();
  void ParseNumericalControls();
  void ParseMesh();
  void ParseRegions();
  void ParseGeochemistry();
  void ParseMaterials();
  void ParseProcessKernels();
  void ParsePhases();
  void ParseInitialConditions();
  void ParseBoundaryConditions();
  void ParseOutput();
  void ParseMisc();

  void ParseMechProperty_(xercesc::DOMElement* mech_prop_node, 
                          const std::string& material_name, 
                          const std::string& property_name,
                          std::list<ParmParse::PP_entry>& table,
                          bool required);
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
