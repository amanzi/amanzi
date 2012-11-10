#ifndef DEFORM_HH_
#define DEFORM_HH_

#include<iostream>
using namespace std;

#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_MSTK.hh"
#include "Mesh.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

#define LINE(A)                                     \
  do {                                              \
    for ( int i=0; i<20; ++i ) { std::cout<<(#A); } \
    std::cout<<std::endl ;                          \
  } while(0)

class DeformMesh {
private:
  Teuchos::ParameterList& plist_;
  const Teuchos::RCP<Mesh>& mesh_;

  // aux methods
  void loop_monitor(int k, int kmax);

public:
  DeformMesh( Teuchos::ParameterList& plist, 
	      const Teuchos::RCP<AmanziMesh::Mesh>& mesh );
 ~DeformMesh() {}

  void check_mesh_nodes();

  void move_a_single_node();
  void move_a_node_column();
  void parabolic_profile();
  void parabolic_profile( double ss );

  double stretch( double old_val );

  void print_VTK_structured_mesh( string fname );
  void print_VTK_unstructured_mesh( string fname );

  void print_goodbye();
};

#endif
