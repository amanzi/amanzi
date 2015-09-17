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

#ifdef ENABLE_Structured

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

  void ParseUnits_();
  void ParseDefinitions_();
  void ParseExecutionControls_();
  void ParseNumericalControls_();
  void ParseMesh_();
  void ParseRegions_();
  void ParseGeochemistry_();
  void ParseMaterials_();
  void ParseProcessKernels_();
  void ParsePhases_();
  void ParseInitialConditions_();
  void ParseBoundaryConditions_();
  void ParseOutput_();
  void ParseMisc_();

  void ParseMechProperty_(xercesc::DOMElement* mech_prop_node, 
                          const std::string& material_name, 
                          const std::string& property_name,
                          std::list<ParmParse::PP_entry>& table,
                          bool required);

  // Private data.
  int dim_;
  int nx_, ny_, nz_;
  std::vector<double> lo_coords_, hi_coords_;
  std::string chemistry_engine_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif

#endif
