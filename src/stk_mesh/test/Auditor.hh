// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: Auditor.hh
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December  9, 2010 by William A. Perkins
// Last Change: Tue May 17 11:46:33 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _Auditor_hh_
#define _Auditor_hh_

#include <iostream>
#include <sstream>
#include <Teuchos_RCPDecl.hpp>
#include "MeshAudit.hh"
#include "../Mesh_STK.hh"

// -------------------------------------------------------------
//  class Auditor
// -------------------------------------------------------------
class Auditor {
 protected:

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_map;
  Teuchos::RCP<MeshAudit> audit;
  std::ofstream ofs;

  /// Protected copy constructor to avoid unwanted copies.
  Auditor(const Auditor& old);

 public:

  /// Default constructor.
  Auditor(std::string oname, Amanzi::AmanziMesh::STK::Mesh_STK_Impl_p mesh)
      : mesh_map(new Amanzi::AmanziMesh::Mesh_STK(mesh)) {
    std::ostringstream ofile;
    ofile << oname
          << std::setfill('0') << std::setw(4) 
          << mesh->communicator().MyPID() << ".tst";
    ofs.open(ofile.str().c_str());
    if (mesh->communicator().MyPID() == 0)
      std::cout << "Writing results to " << ofile.str() << ", etc." << std::endl;
    audit.reset(new MeshAudit(mesh_map, ofs));
  }

  Auditor(std::string oname, Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mm)
      : mesh_map(mm) {
    std::ostringstream ofile;
    ofile << oname
          << std::setfill('0') << std::setw(4) 
          << mesh_map->get_comm()->MyPID() << ".tst";
    ofs.open(ofile.str().c_str());
    if (mesh_map->get_comm()->MyPID() == 0)
      std::cout << "Writing results to " << ofile.str() << ", etc." << std::endl;
    audit.reset(new MeshAudit(mesh_map, ofs));
  }

  /// Destructor
  ~Auditor(void) {
    ofs.close();
  }

  void 
  operator() (void) {

    mesh_map->node_epetra_map(true).Print(std::cerr);
    mesh_map->face_epetra_map(true).Print(std::cerr);
    mesh_map->cell_epetra_map(true).Print(std::cerr);
    CHECK(audit->Verify() == 0);
    // CHECK(audit->check_entity_counts() == 0);
    // CHECK(audit->check_cell_to_nodes_refs() == 0);
    // CHECK(audit->check_cell_to_faces_refs() == 0);
    // CHECK(audit->check_face_to_nodes_refs() == 0);
    // CHECK(audit->check_cell_to_nodes_consistency() == 0);
    // CHECK(audit->check_cell_to_faces_consistency() == 0);
    // CHECK(audit->check_face_to_nodes_consistency() == 0);
    // CHECK(audit->check_cell_to_face_dirs_basic() == 0);
    // CHECK(audit->check_node_to_coordinates() == 0);
    // CHECK(audit->check_cell_degeneracy() == 0);
    // CHECK(audit->check_cell_to_faces() == 0);
    // CHECK(audit->check_cell_to_coordinates() == 0);
    // CHECK(audit->check_face_to_coordinates() == 0);
    // CHECK(audit->check_cell_topology() == 0);
    // CHECK(audit->check_node_maps() == 0);
    // CHECK(audit->check_face_maps() == 0);
    // CHECK(audit->check_cell_maps() == 0);
    // CHECK(audit->check_node_to_coordinates_ghost_data() == 0);
    // CHECK(audit->check_face_to_nodes_ghost_data() == 0);
    // CHECK(audit->check_cell_to_nodes_ghost_data() == 0);
    // CHECK(audit->check_cell_to_faces_ghost_data() == 0);
        
    // CHECK(audit->check_node_partition() == 0);
    // CHECK(audit->check_face_partition() == 0);
    // CHECK(audit->check_node_sets() == 0);
    // CHECK(audit->check_face_sets() == 0);
    // CHECK(audit->check_cell_sets() == 0);
  }

};

#endif
