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

#include "composite_vector.hh"
#include "tree_vector.hh"
#include "state.hh"

#include "PK.hh"
#include "pk_factory.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

#define LINE(A)                                     \
  do {                                              \
    for ( int i=0; i<20; ++i ) { std::cout<<(#A); } \
    std::cout<<std::endl ;                          \
  } while(0)

namespace Amanzi {
namespace Deform {


//temporary class, useful to store info about the vertical columns of
//nodes when building the starting mesh from the last one for the
//final movie (it will be probably removed in next versions)
class PCol {
private:
  vector<int> vlist;
  vector<double> zlist;
  double zmin, zmax;
public:
  PCol(): zmin(+1e+20), zmax(-1e+20) {}
  ~PCol(){}
  void insert( int _iV, double _zV ) {
    vlist.push_back(_iV);
    zlist.push_back(_zV);
    zmin = min( zmin, _zV ) ;
    zmax = max( zmax, _zV ) ;
  }
  int    get_size() { return vlist.size(); }
  double get_zmin() { return zmin; }
  double get_zmax() { return zmax; }
  int    get_iV( int i ) { return vlist[i]; }
  double get_zV( int i ) { return zlist[i]; }
};


class DeformMesh {
private:
  Teuchos::ParameterList& plist_;
  const Teuchos::RCP<Mesh>& mesh0_; // first mesh
  const Teuchos::RCP<Mesh>& mesh1_; // final mesh
  const Teuchos::RCP<Mesh>& mesh_;  // current mesh

  // aux methods
  void loop_monitor(int k, int kmax);

public:
  // constructor with three meshes (useful for the movie test)
  DeformMesh( Teuchos::ParameterList& plist, 
	      const Teuchos::RCP<AmanziMesh::Mesh>& mesh0,
	      const Teuchos::RCP<AmanziMesh::Mesh>& mesh1,
	      const Teuchos::RCP<AmanziMesh::Mesh>& mesh );
  // usual pk constructor
  DeformMesh( Teuchos::ParameterList& plist, 
	      const Teuchos::RCP<AmanziMesh::Mesh>& mesh );
  // default destructor
 ~DeformMesh() {}

  // main PK methods
  // -- Initialize owned (dependent) variables.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);

  // impose a final bell_shaped profile
  void bell_shaped_profile();
  void bell_shaped_profile( double ss );
  void bell_shaped_profile( double ss, Amanzi::AmanziGeometry::Point &P0 );

  // impose a final layer profile 
  void layer_profile( double ss );
  void layer_profile();

  void set_mesh_name( int k, string &fname );

  // VTK output
  void print_VTK_domain_boundary( string fname );
  void print_VTK_unstructured_mesh( string fname );
  void print_VTK_submesh(string fname);

  // mesh deformation from the final mesh
  void build_the_starting_mesh( Entity_ID_List & newnod );
  void mesh_deformation();
  void mesh_deformation_top_nodes();

  // remap the column of nodes (it uses PCol to store the node information)
  void analyze_final_mesh( vector<PCol> & pcol );
  void analyze_final_mesh();

  // final print
  void print_goodbye();
};

}
}

#endif // end of DEFORM_HH_
