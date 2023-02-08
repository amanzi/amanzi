/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include <UnitTest++.h>

#include <iostream>

#include "AmanziComm.hh"

#include "Geometry.hh"
#include "MeshException.hh"
#include "MeshFrameworkTraits.hh"
#include "MeshFrameworkFactory.hh"

TEST(MESH_MAPS)
{
  auto comm = Amanzi::getDefaultComm();

  // Create the mesh
  Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);
  auto mesh = meshfactory.create(0.0,0.0,1.0,1.0,2,2);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm->SumAll(&ierr, &aerr, 1);

    CHECK_EQUAL(aerr,0);


    double exp_getCellVolume[4] = {0.25,0.25,0.25,0.25};
    double exp_getCellCentroid[4][2] = {{0.25,0.25},
                                      {0.25,0.75},
                                      {0.75,0.25},
                                      {0.75,0.75}};
    double exp_getFaceArea[12] = {0.5,0.5,0.5,0.5,
                                0.5,0.5,0.5,0.5,
                                0.5,0.5,0.5,0.5};
    double exp_getFaceCentroid[12][2] = {{0.25,0.0},{0.5,0.25},
                                       {0.25,0.5},{0.0,0.25},
                                       {0.5,0.75},{0.25,1.0},
                                       {0.0,0.75},{0.75,0.0},
                                       {1.0,0.25},{0.75,0.5},
                                       {1.0,0.75},{0.75,1.0}};

    int ncells = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,Amanzi::AmanziMesh::Parallel_kind::OWNED);
    int nfaces = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE,Amanzi::AmanziMesh::Parallel_kind::ALL);

    int space_dim_ = 2;

    for (int i = 0; i < ncells; i++) {

      Amanzi::AmanziGeometry::Point centroid = mesh->getCellCentroid(i);

      // Search for a cell with the same centroid in the 
      // expected list of centroid

      bool found = false;

      for (int j = 0; j < ncells; j++) {
        if (fabs(exp_getCellCentroid[j][0]-centroid[0]) < 1.0e-10 &&
            fabs(exp_getCellCentroid[j][1]-centroid[1]) < 1.0e-10) {

          found = true;
          CHECK_EQUAL(exp_getCellVolume[j],mesh->getCellVolume(i));
          break;

        }
      }

      CHECK_EQUAL(found,true);

      Amanzi::AmanziMesh::Entity_ID_View cfaces;
      Amanzi::AmanziGeometry::Point normal_sum(2), normal(2);      

      mesh->getCellFaces(i,cfaces);
      normal_sum.set(0.0);

      for (int j = 0; j < cfaces.size(); j++) {
        normal = mesh->getFaceNormal(cfaces[j],i);
        normal_sum += normal;
      }

      double val = L22(normal_sum);
      CHECK_CLOSE(val,0.0,1.0e-20);                
    }

    for (int i = 0; i < nfaces; i++) {
      Amanzi::AmanziGeometry::Point centroid = mesh->getFaceCentroid(i);

      bool found = false;

      for (int j = 0; j < nfaces; j++) {
        if (fabs(exp_getFaceCentroid[j][0]-centroid[0]) < 1.0e-10 &&
            fabs(exp_getFaceCentroid[j][1]-centroid[1]) < 1.0e-10) {

          found = true;

          CHECK_EQUAL(exp_getFaceArea[j],mesh->getFaceArea(i));

          // Check the natural normal

          Amanzi::AmanziGeometry::Point normal = mesh->getFaceNormal(i);
            
      
          // Check the normal with respect to each connected cell
          
          Amanzi::AmanziMesh::Entity_ID_View cellids;
          mesh->getFaceCells(i,Amanzi::AmanziMesh::Parallel_kind::ALL,cellids);
          
          for (int k = 0; k < cellids.size(); k++) {
            int dir;
            Amanzi::AmanziGeometry::Point normal_wrt_cell =
              mesh->getFaceNormal(i,cellids[k],&dir);

            //            Amanzi::AmanziMesh::Entity_ID_View cellfaces;
            //            std::vector<int> cellfacedirs;
            //            mesh->getCellFacesAndDirs(cellids[k],cellfaces,&cellfacedirs);
            //
            //            bool found2 = false;
            //            int dir = 1;
            //            for (int m = 0; m < cellfaces.size(); m++) {
            //              if (cellfaces[m] == i) {
            //                found2 = true;
            //                dir = cellfacedirs[m];
            //                break;
            //              }
            //            }
            //
            //            CHECK_EQUAL(found2,true);

            Amanzi::AmanziGeometry::Point normal1(normal);
            normal1 *= dir;

            CHECK_ARRAY_EQUAL(&(normal1[0]),&(normal_wrt_cell[0]),space_dim_);

            
            Amanzi::AmanziGeometry::Point cellcentroid = mesh->getCellCentroid(cellids[k]);
            Amanzi::AmanziGeometry::Point facecentroid = mesh->getFaceCentroid(i);

            Amanzi::AmanziGeometry::Point outvec = facecentroid-cellcentroid;


            double dp = outvec*normal_wrt_cell;
            dp /= (norm(outvec)*norm(normal_wrt_cell));         


            CHECK_CLOSE(dp,1.0,1e-10);
          }

          break;
        }
      }

      CHECK_EQUAL(found,true);
    }

  } // for each framework i

}


TEST(MESH_GEOMETRY_SURFACE)
{
  // DISABLED FOR NOW
  return;

  auto comm = Amanzi::getDefaultComm();

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh
  std::vector<Amanzi::AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

#ifdef HAVE_MSTK_MESH
  frameworks.push_back(Amanzi::AmanziMesh::Framework::MSTK);
  framework_names.push_back("MSTK");
#endif

  for (int frm = 0; frm < frameworks.size(); ++frm) {
    // Set the framework
    std::cerr << "Testing geometry operators with " << framework_names[frm] << std::endl;

    // Create the mesh
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);
    Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::Preference prefs(meshfactory.get_preference());
      prefs.clear(); 
      prefs.push_back(frameworks[frm]);

      meshfactory.set_preference(prefs);

      mesh = meshfactory.create("test/surfquad.exo");

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm->SumAll(&ierr, &aerr, 1);

    CHECK_EQUAL(aerr,0);


    double exp_getCellVolume[4] = {0.25,0.25,0.25,0.25};
    double exp_getCellCentroid[4][3] = {{0.25,0.25,0.0},
                                      {0.25,0.75,0.0},
                                      {0.5,0.25,0.25},
                                      {0.5,0.75,0.25}};
    double exp_getFaceArea[12] = {0.5,0.5,0.5,0.5,
                                0.5,0.5,0.5,0.5,
                                0.5,0.5,0.5,0.5};
    double exp_getFaceCentroid[12][3] = {{0.25,0.0,0.0},{0.5,0.25,0.0},
                                       {0.25,0.5,0.0},{0.0,0.25,0.0},
                                       {0.5,0.75,0.0},{0.25,1.0,0.0},
                                       {0.0,0.75,0.0},{0.5,0.0,0.25},
                                       {0.5,0.25,0.5},{0.5,0.5,0.25},
                                       {0.5,0.75,0.5},{0.5,1.0,0.25}};

    int ncells = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,Amanzi::AmanziMesh::Parallel_kind::OWNED);
    int nfaces = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE,Amanzi::AmanziMesh::Parallel_kind::ALL);

    for (int i = 0; i < ncells; i++) {

      Amanzi::AmanziGeometry::Point centroid = mesh->getCellCentroid(i);

      // Search for a cell with the same centroid in the 
      // expected list of centroid

      bool found = false;

      for (int j = 0; j < ncells; j++) {
        if (fabs(exp_getCellCentroid[j][0]-centroid[0]) < 1.0e-10 &&
            fabs(exp_getCellCentroid[j][1]-centroid[1]) < 1.0e-10 &&
            fabs(exp_getCellCentroid[j][2]-centroid[2]) < 1.0e-10) {

          found = true;
          CHECK_EQUAL(exp_getCellVolume[j],mesh->getCellVolume(i));
          break;

        }
      }

      CHECK_EQUAL(found,true);
    }

    for (int i = 0; i < nfaces; i++) {
      Amanzi::AmanziGeometry::Point centroid = mesh->getFaceCentroid(i);

      bool found = false;

      for (int j = 0; j < nfaces; j++) {
        if (fabs(exp_getFaceCentroid[j][0]-centroid[0]) < 1.0e-10 &&
            fabs(exp_getFaceCentroid[j][1]-centroid[1]) < 1.0e-10 &&
            fabs(exp_getFaceCentroid[j][2]-centroid[2]) < 1.0e-10) {

          found = true;

          CHECK_EQUAL(exp_getFaceArea[j],mesh->getFaceArea(i));

      
          // Check the normal with respect to each connected cell
          
          Amanzi::AmanziMesh::Entity_ID_View cellids;
          mesh->getFaceCells(i,Amanzi::AmanziMesh::Parallel_kind::ALL,cellids);


          Amanzi::AmanziGeometry::Point facecentroid = mesh->getFaceCentroid(i);
          
          for (int k = 0; k < cellids.size(); k++) {
            int dir;
            Amanzi::AmanziGeometry::Point normal_wrt_cell = 
              mesh->getFaceNormal(i,cellids[k],&dir);

            //            Amanzi::AmanziMesh::Entity_ID_View cellfaces;
            //            std::vector<int> cellfacedirs;
            //            mesh->getCellFacesAndDirs(cellids[k],&cellfaces,&cellfacedirs);

            //            bool found2 = false;
            //            int dir = 1;
            //            for (int m = 0; m < cellfaces.size(); m++) {
            //              if (cellfaces[m] == i) {
            //                found2 = true;
            //                dir = cellfacedirs[m];
            //                break;
            //              }
            //            }

            //            CHECK_EQUAL(found2,true);


            Amanzi::AmanziGeometry::Point cellcentroid = mesh->getCellCentroid(cellids[k]);

            Amanzi::AmanziGeometry::Point outvec = facecentroid-cellcentroid;

            double dp = outvec*normal_wrt_cell;
            dp /= (norm(outvec)*norm(normal_wrt_cell));         

            CHECK_CLOSE(dp,1.0,1e-10);

          }


          if (cellids.size() == 2 && 
              ((fabs(facecentroid[0]-0.5) < 1.0e-16) &&
               (fabs(facecentroid[2]) < 1.0e-16))) {

            // An edge on the crease. The two normals should be different 
   
            Amanzi::AmanziGeometry::Point n0 = mesh->getFaceNormal(i,cellids[0]);
            Amanzi::AmanziGeometry::Point n1 = mesh->getFaceNormal(i,cellids[1]);

            double dp = n0*n1/(norm(n0)*norm(n1));

            CHECK_CLOSE(dp,0.0,1e-10);

          }

          break;
        }
      }

      CHECK_EQUAL(found,true);
    }

  } // for each framework i
}


TEST(MESH_GEOMETRY_SOLID)
{
  auto comm = Amanzi::getDefaultComm();

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh
  std::vector<Amanzi::AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

  frameworks.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
  framework_names.push_back("Simple");
  
#ifdef HAVE_MSTK_MESH
  frameworks.push_back(Amanzi::AmanziMesh::Framework::MSTK);
  framework_names.push_back("MSTK");
#endif
  
  for (int frm = 0; frm < frameworks.size(); ++frm) {
    // Set the framework
    std::cerr << "Testing geometry operators with " << framework_names[frm] << std::endl;

    // Create the mesh
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);
    Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::Preference prefs(meshfactory.get_preference());
      prefs.clear(); 
      prefs.push_back(frameworks[frm]);
      meshfactory.set_preference(prefs);

      mesh = meshfactory.create(0.0,0.0,0.0,1.0,1.0,1.0,2,2,2);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm->SumAll(&ierr, &aerr, 1);
    CHECK_EQUAL(0, aerr);


    double exp_getCellVolume[8] = {0.125,0.125,0.125,0.125,
                                 0.125,0.125,0.125,0.125};
    double exp_getCellCentroid[8][3] = {{0.25,0.25,0.25},
                                      {0.75,0.25,0.25},
                                      {0.25,0.75,0.25},
                                      {0.75,0.75,0.25},
                                      {0.25,0.25,0.75},
                                      {0.75,0.25,0.75},
                                      {0.25,0.75,0.75},
                                      {0.75,0.75,0.75}};
    double exp_getFaceArea[36] = {0.25,0.25,0.25,0.25,
                                0.25,0.25,0.25,0.25,
                                0.25,0.25,0.25,0.25,
                                0.25,0.25,0.25,0.25,
                                0.25,0.25,0.25,0.25,
                                0.25,0.25,0.25,0.25,
                                0.25,0.25,0.25,0.25,
                                0.25,0.25,0.25,0.25,
                                0.25,0.25,0.25,0.25};
    double exp_getFaceCentroid[36][3] = {{0.0,0.25,0.25},
                                       {0.0,0.75,0.25},
                                       {0.0,0.25,0.75},
                                       {0.0,0.75,0.75},

                                       {0.5,0.25,0.25},
                                       {0.5,0.75,0.25},
                                       {0.5,0.25,0.75},
                                       {0.5,0.75,0.75},

                                       {1.0,0.25,0.25},
                                       {1.0,0.75,0.25},
                                       {1.0,0.25,0.75},
                                       {1.0,0.75,0.75},

                                       {0.25,0.0,0.25},
                                       {0.75,0.0,0.25},
                                       {0.25,0.0,0.75},
                                       {0.75,0.0,0.75},
                                                 
                                       {0.25,0.5,0.25},
                                       {0.75,0.5,0.25},
                                       {0.25,0.5,0.75},
                                       {0.75,0.5,0.75},
                                                 
                                       {0.25,1.0,0.25},
                                       {0.75,1.0,0.25},
                                       {0.25,1.0,0.75},
                                       {0.75,1.0,0.75},

                                       {0.25,0.25,0.0},
                                       {0.75,0.25,0.0},
                                       {0.25,0.75,0.0},
                                       {0.75,0.75,0.0},
                                                     
                                       {0.25,0.25,0.5},
                                       {0.75,0.25,0.5},
                                       {0.25,0.75,0.5},
                                       {0.75,0.75,0.5},
                                                     
                                       {0.25,0.25,1.0},
                                       {0.75,0.25,1.0},
                                       {0.25,0.75,1.0},
                                       {0.75,0.75,1.0},
                                       
    };


    int ncells = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,Amanzi::AmanziMesh::Parallel_kind::OWNED);
    int nfaces = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE,Amanzi::AmanziMesh::Parallel_kind::ALL);

    int space_dim_ = 3;

    for (int i = 0; i < ncells; i++) {

      Amanzi::AmanziGeometry::Point centroid = mesh->getCellCentroid(i);

      // Search for a cell with the same centroid in the 
      // expected list of centroid

      bool found = false;

      for (int j = 0; j < ncells; j++) {
        if (fabs(exp_getCellCentroid[j][0]-centroid[0]) < 1.0e-10 &&
            fabs(exp_getCellCentroid[j][1]-centroid[1]) < 1.0e-10 &&
            fabs(exp_getCellCentroid[j][2]-centroid[2]) < 1.0e-10) {

          found = true;
          CHECK_EQUAL(exp_getCellVolume[j],mesh->getCellVolume(i));
          break;

        }
      }

      CHECK_EQUAL(found,true);

      Amanzi::AmanziMesh::Entity_ID_View cfaces;
      Amanzi::AmanziGeometry::Point normal_sum(3), normal(3);      

      mesh->getCellFaces(i,cfaces);
      normal_sum.set(0.0);

      for (int j = 0; j < cfaces.size(); j++) {
        normal = mesh->getFaceNormal(cfaces[j],i);
        normal_sum += normal;
      }

      double val = L22(normal_sum);
      CHECK_CLOSE(val,0.0,1.0e-20);                
    }

    for (int i = 0; i < nfaces; i++) {
      Amanzi::AmanziGeometry::Point centroid = mesh->getFaceCentroid(i);

      bool found = false;

      for (int j = 0; j < nfaces; j++) {
        if (fabs(exp_getFaceCentroid[j][0]-centroid[0]) < 1.0e-10 &&
            fabs(exp_getFaceCentroid[j][1]-centroid[1]) < 1.0e-10 &&
            fabs(exp_getFaceCentroid[j][2]-centroid[2]) < 1.0e-10) {

          found = true;

          CHECK_EQUAL(exp_getFaceArea[j],mesh->getFaceArea(i));

          // Check the natural normal

          Amanzi::AmanziGeometry::Point normal = mesh->getFaceNormal(i);
            
      
          // Check the normal with respect to each connected cell
          
          Amanzi::AmanziMesh::Entity_ID_View cellids;
          mesh->getFaceCells(i,Amanzi::AmanziMesh::Parallel_kind::ALL,cellids);
          
          for (int k = 0; k < cellids.size(); k++) {
            int dir;
            Amanzi::AmanziGeometry::Point normal_wrt_cell = 
              mesh->getFaceNormal(i,cellids[k],&dir);

            // Amanzi::AmanziMesh::Entity_ID_View cellfaces;
            // std::vector<int> cellfacedirs;
            // mesh->getCellFacesAndDirs(cellids[k],&cellfaces,cellfacedirs);

            // bool found2 = false;
            // int dir = 1;
            // for (int m = 0; m < cellfaces.size(); m++) {
            //   if (cellfaces[m] == i) {
            //     found2 = true;
            //     dir = cellfacedirs[m];
            //     break;
            //   }
            // }

            // CHECK_EQUAL(found2,true);

            Amanzi::AmanziGeometry::Point normal1(normal);
            normal1 *= dir;

            CHECK_ARRAY_EQUAL(&(normal1[0]),&(normal_wrt_cell[0]),space_dim_);

            
            Amanzi::AmanziGeometry::Point cellcentroid = mesh->getCellCentroid(cellids[k]);
            Amanzi::AmanziGeometry::Point facecentroid = mesh->getFaceCentroid(i);

            Amanzi::AmanziGeometry::Point outvec = facecentroid-cellcentroid;

            double dp = outvec*normal_wrt_cell;
            dp /= (norm(outvec)*norm(normal_wrt_cell));         


            CHECK_CLOSE(dp,1.0,1e-10);
          }

          break;
        }
      }

      CHECK_EQUAL(found,true);
    }

    
    // Now deform the mesh a little and verify that the sum of the
    // outward normals of all faces of cell is still zero

    Amanzi::AmanziGeometry::Point ccoords = mesh->getNodeCoordinate(13); // central node

    // Lets be sure this is the central node
    CHECK_EQUAL(ccoords[0],0.5);
    CHECK_EQUAL(ccoords[1],0.5);
    CHECK_EQUAL(ccoords[2],0.5);
    
    // Perturb it
    ccoords.set(0.7,0.7,0.7);
    mesh->setNodeCoordinate(13,ccoords);

    // Now check the normals

    for (int i = 0; i < ncells; i++) {

      Amanzi::AmanziMesh::Entity_ID_View cfaces;
      Amanzi::AmanziGeometry::Point normal_sum(3), normal(3);      

      mesh->getCellFaces(i,cfaces);
      normal_sum.set(0.0);

      for (int j = 0; j < cfaces.size(); j++) {
        normal = mesh->getFaceNormal(cfaces[j],i);
        normal_sum += normal;
      }

      double val = L22(normal_sum);
      CHECK_CLOSE(val,0.0,1.0e-20);                

    }    

  } // for each framework i

}

