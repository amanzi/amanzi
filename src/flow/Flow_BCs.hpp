#ifndef __FlowBCs_hpp__
#define __FlowBCs_hpp__

#include "Mesh_maps.hh"
#include "Teuchos_ParameterList.hpp"

// we assume that the BCs have already read into
// the internal Teuchos parameter list database
// 
// we also assume that the flow BCs are named
// FlowBC01 ... FlowBCXX, where XX is the number
// of flow BCs (really, this is the number of 
// side sets


struct flow_bc {
  double value;
  unsigned int side_set;
  std::string bc_type;
  int num_faces;
  std::vector<int> faces;
  std::vector<double> aux;
};

enum Flow_BC_Kind
{
  PRESS_CONST,
  NO_FLOW
};

class Flow_BCs {

public:
  Flow_BCs(Teuchos::RCP<const STK_mesh::Mesh_maps> mesh_maps_, Teuchos::ParameterList &parameterList):
    mesh(mesh_maps_), paramList(parameterList) {};
  ~Flow_BCs() {};


  void read_Flow_BCs();
  const int get_num_BCs() { return num_BCs; };
  const std::vector<flow_bc> get_BCs() { return BCs; };
  
private:
  const Teuchos::RCP<const STK_mesh::Mesh_maps> mesh;
  Teuchos::ParameterList paramList; 
  int num_BCs;
  std::vector<flow_bc> BCs;
};

#endif
