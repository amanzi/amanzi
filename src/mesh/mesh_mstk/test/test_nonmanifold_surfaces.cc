#include <iostream>

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

#include "../Mesh_MSTK.hh"


TEST(MSTK_READ_NONMANIFOLD_SURFACES)
{
  std::string expcsetnames[2] = {"FRACTURE 1", "FRACTURE 2"};
  
  Teuchos::RCP<Epetra_MpiComm> comm_(new Epetra_MpiComm(MPI_COMM_WORLD));

  Teuchos::ParameterList parameterlist;
 
  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions"); 
  
  Teuchos::ParameterList& frac1 = reg_spec.sublist("FRACTURE 1");
  Teuchos::ParameterList& frac1_def = frac1.sublist("region: plane");
  Teuchos::Array<double> loc1 = Teuchos::tuple(0.0,0.0,0.5);
  Teuchos::Array<double> dir1 = Teuchos::tuple(0.0,0.0,1.0);
  frac1_def.set< Teuchos::Array<double> >("point",loc1);
  frac1_def.set< Teuchos::Array<double> >("normal",dir1);
  frac1_def.set< std::string >("entity","cell");

  Teuchos::ParameterList& frac2 = reg_spec.sublist("FRACTURE 2");
  Teuchos::ParameterList& frac2_def = frac2.sublist("region: plane");
  Teuchos::Array<double> loc2 = Teuchos::tuple(0.0,0.5,0.0);
  Teuchos::Array<double> dir2 = Teuchos::tuple(0.0,-1.0,0.0);
  frac2_def.set< Teuchos::Array<double> >("point",loc2);
  frac2_def.set< Teuchos::Array<double> >("normal",dir2);
  frac2_def.set< std::string >("entity","cell");


  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, comm_.get()));


  // Load a mesh consisting of quads representing two perpendicular
  // intersecting surfaces in 3 space. Each surface is meshed with 100
  // regular quads

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_MSTK("test/fractures.exo",comm_.get(),3,gm));


  CHECK_EQUAL(2, mesh->manifold_dimension());
  CHECK_EQUAL(3, mesh->space_dimension());

  int nc = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);

  Teuchos::ParameterList::ConstIterator i;
  for (i = reg_spec.begin(); i != reg_spec.end(); i++) {
    const std::string reg_name = reg_spec.name(i);     

    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // See if the geometric model has a region by this name
  
    Teuchos::RCP<const Amanzi::AmanziGeometry::Region> reg = gm->FindRegion(reg_name);

    CHECK(reg.get());

    // Do their names match ?

    CHECK_EQUAL(reg->name(),reg_name);

    // Get the region info directly from the XML and compare
  
    Teuchos::ParameterList::ConstIterator j = reg_params.begin(); 

    std::string shape = reg_params.name(j);

    if (shape == "region: plane") {

      if (reg_name == "FRACTURE 1" || reg_name == "FRACTURE 2") {

        // Do we have a valid cellset by this name
        
        CHECK(mesh->valid_set_name(reg_name,Amanzi::AmanziMesh::CELL));

        if (expcsetnames[0] != reg_name && expcsetnames[1] != reg_name) break;
        
        // Verify that we can get the number of entities in the set
        
        int set_size = mesh->get_set_size(reg_name, Amanzi::AmanziMesh::CELL,
                                          Amanzi::AmanziMesh::OWNED);

        CHECK_EQUAL(nc/2, set_size);
     
        // Verify that we can retrieve the set entities
        
        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh->get_set_entities(reg_name, Amanzi::AmanziMesh::CELL,
                               Amanzi::AmanziMesh::OWNED, &setents);

        if (reg_name == "FRACTURE 1") {  
          // VERIFY THAT EACH 2D CELL HAS A NORMAL THAT POINTS IN [0,0,1] DIR
          std::vector<Amanzi::AmanziGeometry::Point> ccoords;
          for (auto const& cellid : setents) {
            mesh->cell_get_coordinates(cellid, &ccoords);

            // Assumes that for 2D cells the coordinates are returned in 
            // ccw direction around the cell
            Amanzi::AmanziGeometry::Point vec0 = ccoords[2]-ccoords[1];
            Amanzi::AmanziGeometry::Point vec1 = ccoords[0]-ccoords[1];
            Amanzi::AmanziGeometry::Point cpvec = vec0^vec1;
            double len = norm(cpvec);
            cpvec /= len;

            CHECK_CLOSE(0.0, cpvec[0], 1.0e-08);
            CHECK_CLOSE(0.0, cpvec[1], 1.0e-08);
            CHECK_CLOSE(1.0, cpvec[2], 1.0e-08);
          }
        } else {
          // VERIFY THAT EACH 2D CELL HAS A NORMAL THAT POINTS IN [0,1,0] DIR
          std::vector<Amanzi::AmanziGeometry::Point> ccoords;
          for (auto const& cellid : setents) {
            mesh->cell_get_coordinates(cellid, &ccoords);

            // Assumes that for 2D cells the coordinates are returned in 
            // ccw direction around the cell
            Amanzi::AmanziGeometry::Point vec0 = ccoords[2]-ccoords[1];
            Amanzi::AmanziGeometry::Point vec1 = ccoords[0]-ccoords[1];
            Amanzi::AmanziGeometry::Point cpvec = vec0^vec1;
            double len = norm(cpvec);
            cpvec /= len;

            CHECK_CLOSE(0.0, cpvec[0], 1.0e-08);
            CHECK_CLOSE(-1.0, cpvec[1], 1.0e-08);
            CHECK_CLOSE(0.0, cpvec[2], 1.0e-08);
          }
        }
      }
    }
  }
}

