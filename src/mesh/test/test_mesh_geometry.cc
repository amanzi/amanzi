/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

// -------------------------------------------------------------
/**
 * @file   test_mesh_geometry.cc
 * @author Rao V. Garimella
 * @date   Tue May 15, 2012
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <iostream>

#include "AmanziComm.hh"
#include "AmanziMap.hh"

#include "Geometry.hh"
#include "MeshException.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

struct test {};

// TEST(MESH_GEOMETRY_PLANAR)
TEST_FIXTURE(test, MESH_GEOMETRY_PLANAR)
{
  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());


  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh
  std::vector<Amanzi::AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

  if (Amanzi::AmanziMesh::framework_enabled(
        Amanzi::AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(Amanzi::AmanziMesh::Framework::MSTK);
    framework_names.push_back("MSTK");
  }

  for (int i = 0; i < frameworks.size(); i++) {
    // Set the framework
    std::cerr << "Testing geometry operators with " << framework_names[i]
              << std::endl;

    // Create the mesh
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::Preference prefs(meshfactory.preference());
      prefs.clear();
      prefs.push_back(frameworks[i]);
      meshfactory.set_preference(prefs);

      mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 2, 2);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &aerr);

    CHECK_EQUAL(aerr, 0);

    const int ncells = mesh->num_entities(
      Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    const int nfaces = mesh->num_entities(
      Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::ALL);
    const int nnodes = mesh->num_entities(
      Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::ALL);

    Kokkos::View<Amanzi::AmanziGeometry::Point*> exp_cell_centroid_view("ccv",
                                                                        4);
    exp_cell_centroid_view(0) = { 0.25, 0.25 };
    exp_cell_centroid_view(1) = { 0.25, 0.75 };
    exp_cell_centroid_view(2) = { 0.75, 0.25 };
    exp_cell_centroid_view(3) = { 0.75, 0.75 };
    Kokkos::View<double*> exp_cell_volume_view("cvv", 4);
    exp_cell_volume_view(0) = 0.25;
    exp_cell_volume_view(1) = 0.25;
    exp_cell_volume_view(2) = 0.25;
    exp_cell_volume_view(3) = 0.25;
    Kokkos::View<double*> exp_face_area_view("fav", 12);
    exp_face_area_view(0) = 0.5;
    exp_face_area_view(1) = 0.5;
    exp_face_area_view(2) = 0.5;
    exp_face_area_view(3) = 0.5;
    exp_face_area_view(4) = 0.5;
    exp_face_area_view(5) = 0.5;
    exp_face_area_view(6) = 0.5;
    exp_face_area_view(7) = 0.5;
    exp_face_area_view(8) = 0.5;
    exp_face_area_view(9) = 0.5;
    exp_face_area_view(10) = 0.5;
    exp_face_area_view(11) = 0.5;
    Kokkos::View<Amanzi::AmanziGeometry::Point*> exp_face_centroid_view("fcv",
                                                                        12);
    exp_face_centroid_view(0) = { 0.25, 0.0 };
    exp_face_centroid_view(1) = { 0.5, 0.25 };
    exp_face_centroid_view(2) = { 0.25, 0.5 };
    exp_face_centroid_view(3) = { 0.0, 0.25 };
    exp_face_centroid_view(4) = { 0.5, 0.75 };
    exp_face_centroid_view(5) = { 0.25, 1.0 };
    exp_face_centroid_view(6) = { 0.0, 0.75 };
    exp_face_centroid_view(7) = { 0.75, 0.0 };
    exp_face_centroid_view(8) = { 1.0, 0.25 };
    exp_face_centroid_view(9) = { 0.75, 0.5 };
    exp_face_centroid_view(10) = { 1.0, 0.75 };
    exp_face_centroid_view(11) = { 0.75, 1.0 };

    const int space_dim_ = 2;

    Amanzi::AmanziMesh::Mesh* m = mesh.get();

    // Perform tests on GPU
    Kokkos::View<bool*> result_cv("result cv", ncells);
    Kokkos::View<bool*> result_found("result found", 1);
    Kokkos::parallel_for(ncells, KOKKOS_LAMBDA(const Amanzi::LO& i) {
      Kokkos::View<int*> test;
      Amanzi::AmanziGeometry::Point centroid = m->cell_centroid(i);
      for (int j = 0; j < ncells; ++j) {
        if (fabs(exp_cell_centroid_view(j)[0] - centroid[0]) < 1.0e-10 &&
            fabs(exp_cell_centroid_view(j)[1] - centroid[1]) < 1.0e-10) {
          result_found(0) = true;
          assert(exp_cell_volume_view(j) == m->cell_volume(i));
          break;
        }
      }
      assert(result_found(0) == true);

      Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cfaces;
      Amanzi::AmanziGeometry::Point normal_sum(2), normal(2);
      m->cell_get_faces(i, cfaces);
      normal_sum.set(0.0);
      for (int j = 0; j < cfaces.extent(0); ++j) {
        normal = m->face_normal(cfaces(j), false, i);
        normal_sum += normal;
      }

      double val = L22(normal_sum);
      assert(val < 1.0e-20);
    });

    Kokkos::parallel_for(nfaces, KOKKOS_LAMBDA(const Amanzi::LO& i) {
      // for(int i = 0; i < nfaces; ++i){
      Amanzi::AmanziGeometry::Point centroid = m->face_centroid(i);
      bool found = false;

      for (int j = 0; j < nfaces; j++) {
        if (fabs(exp_face_centroid_view(j)[0] - centroid[0]) < 1.0e-10 &&
            fabs(exp_face_centroid_view(j)[1] - centroid[1]) < 1.0e-10) {
          found = true;

          assert(exp_face_area_view[j] == m->face_area(i));

          // Check the natural normal
          Amanzi::AmanziGeometry::Point normal = m->face_normal(i);
          // Check the normal with respect to each connected cell
          Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cellids;
          m->face_get_cells(i, Amanzi::AmanziMesh::Parallel_type::ALL, cellids);

          for (int k = 0; k < cellids.extent(0); k++) {
            int dir;
            Amanzi::AmanziGeometry::Point normal_wrt_cell =
              m->face_normal(i, false, cellids(k), &dir);
            Amanzi::AmanziGeometry::Point normal1(normal);
            normal1 *= dir;
            for (int dim = 0; dim < space_dim_; ++dim) {
              assert(normal1[dim] == normal_wrt_cell[dim]);
            }
            Amanzi::AmanziGeometry::Point cellcentroid =
              m->cell_centroid(cellids(k));
            Amanzi::AmanziGeometry::Point facecentroid = m->face_centroid(i);
            Amanzi::AmanziGeometry::Point outvec = facecentroid - cellcentroid;
            double dp = outvec * normal_wrt_cell;
            dp /= (norm(outvec) * norm(normal_wrt_cell));
            assert(dp - 1.0 < 1e-10);
          }
          break;
        }
      }
      assert(found == true);
    });
  } // for each framework i
}


#if 0 
TEST(MESH_GEOMETRY_SURFACE)
{

// DISABLED FOR NOW

 return;


  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh
  std::vector<Amanzi::AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

#  ifdef HAVE_MSTK_MESH
  frameworks.push_back(Amanzi::AmanziMesh::Framework::MSTK);
  framework_names.push_back("MSTK");
#  endif

  for (int i = 0; i < frameworks.size(); i++) {
    // Set the framework
    std::cerr << "Testing geometry operators with " << framework_names[i] << std::endl;

    // Create the mesh
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::Preference prefs(meshfactory.preference());
      prefs.clear();
      prefs.push_back(frameworks[i]);

      meshfactory.set_preference(prefs);

      mesh = meshfactory.create("test/surfquad.exo");

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &aerr);

    CHECK_EQUAL(aerr,0);


    double exp_cell_volume[4] = {0.25,0.25,0.25,0.25};
    double exp_cell_centroid[4][3] = {{0.25,0.25,0.0},
                                      {0.25,0.75,0.0},
                                      {0.5,0.25,0.25},
                                      {0.5,0.75,0.25}};
    double exp_face_area[12] = {0.5,0.5,0.5,0.5,
                                0.5,0.5,0.5,0.5,
                                0.5,0.5,0.5,0.5};
    double exp_face_centroid[12][3] = {{0.25,0.0,0.0},{0.5,0.25,0.0},
                                       {0.25,0.5,0.0},{0.0,0.25,0.0},
                                       {0.5,0.75,0.0},{0.25,1.0,0.0},
                                       {0.0,0.75,0.0},{0.5,0.0,0.25},
                                       {0.5,0.25,0.5},{0.5,0.5,0.25},
                                       {0.5,0.75,0.5},{0.5,1.0,0.25}};

    int ncells = mesh->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
    int nfaces = mesh->num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::ALL);
    int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::ALL);

    int space_dim_ = 3;


    for (int i = 0; i < ncells; i++) {

      Amanzi::AmanziGeometry::Point centroid = mesh->cell_centroid(i);

      // Search for a cell with the same centroid in the
      // expected list of centroid

      bool found = false;

      for (int j = 0; j < ncells; j++) {
        if (fabs(exp_cell_centroid[j][0]-centroid[0]) < 1.0e-10 &&
            fabs(exp_cell_centroid[j][1]-centroid[1]) < 1.0e-10 &&
            fabs(exp_cell_centroid[j][2]-centroid[2]) < 1.0e-10) {

          found = true;
          CHECK_EQUAL(exp_cell_volume[j],mesh->cell_volume(i,false));
          break;

        }
      }

      CHECK_EQUAL(found,true);
    }

    for (int i = 0; i < nfaces; i++) {
      Amanzi::AmanziGeometry::Point centroid = mesh->face_centroid(i);

      bool found = false;

      for (int j = 0; j < nfaces; j++) {
        if (fabs(exp_face_centroid[j][0]-centroid[0]) < 1.0e-10 &&
            fabs(exp_face_centroid[j][1]-centroid[1]) < 1.0e-10 &&
            fabs(exp_face_centroid[j][2]-centroid[2]) < 1.0e-10) {

          found = true;

          CHECK_EQUAL(exp_face_area[j],mesh->face_area(i));


          // Check the normal with respect to each connected cell

          Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cellids;
          mesh->face_get_cells(i,Amanzi::AmanziMesh::Parallel_type::ALL,cellids);


          Amanzi::AmanziGeometry::Point facecentroid = mesh->face_centroid(i);

          for (int k = 0; k < cellids.extent(0); k++) {
            int dir;
            Amanzi::AmanziGeometry::Point normal_wrt_cell =
              mesh->face_normal(i,false,cellids(k),&dir);

            //            Amanzi::AmanziMesh::Entity_ID_List cellfaces;
            //            std::vector<int> cellfacedirs;
            //            mesh->cell_get_faces_and_dirs(cellids[k],&cellfaces,&cellfacedirs);

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


            Amanzi::AmanziGeometry::Point cellcentroid = mesh->cell_centroid(cellids(k));

            Amanzi::AmanziGeometry::Point outvec = facecentroid-cellcentroid;

            double dp = outvec*normal_wrt_cell;
            dp /= (norm(outvec)*norm(normal_wrt_cell));

            CHECK_CLOSE(dp,1.0,1e-10);

          }


          if (cellids.extent(0) == 2 &&
              ((fabs(facecentroid[0]-0.5) < 1.0e-16) &&
               (fabs(facecentroid[2]) < 1.0e-16))) {

            // An edge on the crease. The two normals should be different

            Amanzi::AmanziGeometry::Point n0 = mesh->face_normal(i,false,cellids(0));
            Amanzi::AmanziGeometry::Point n1 = mesh->face_normal(i,false,cellids(1));

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
#endif

// TEST(MESH_GEOMETRY_SOLID)
TEST_FIXTURE(test, MESH_GEOMETRY_SOLID)
{
  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

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

  for (int i = 0; i < frameworks.size(); i++) {
    // Set the framework
    std::cerr << "Testing geometry operators with " << framework_names[i]
              << std::endl;

    // Create the mesh
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::Preference prefs(meshfactory.preference());
      prefs.clear();
      prefs.push_back(frameworks[i]);

      meshfactory.set_preference(prefs);

      mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &aerr);

    CHECK_EQUAL(aerr, 0);

    Kokkos::View<double*> exp_cell_volume_view("ecvv", 8);
    exp_cell_volume_view(0) = 0.125;
    exp_cell_volume_view(1) = 0.125;
    exp_cell_volume_view(2) = 0.125;
    exp_cell_volume_view(3) = 0.125;
    exp_cell_volume_view(4) = 0.125;
    exp_cell_volume_view(5) = 0.125;
    exp_cell_volume_view(6) = 0.125;
    exp_cell_volume_view(7) = 0.125;
    Kokkos::View<Amanzi::AmanziGeometry::Point*> exp_cell_centroid_view("eccv",
                                                                        8);
    exp_cell_centroid_view(0) = { 0.25, 0.25, 0.25 };
    exp_cell_centroid_view(1) = { 0.75, 0.25, 0.25 };
    exp_cell_centroid_view(2) = { 0.25, 0.75, 0.25 };
    exp_cell_centroid_view(3) = { 0.75, 0.75, 0.25 };
    exp_cell_centroid_view(4) = { 0.25, 0.25, 0.75 };
    exp_cell_centroid_view(5) = { 0.75, 0.25, 0.75 };
    exp_cell_centroid_view(6) = { 0.25, 0.75, 0.75 };
    exp_cell_centroid_view(7) = { 0.75, 0.75, 0.75 };
    Kokkos::View<double*> exp_face_area_view("efav", 36);
    for (int i = 0; i < 36; ++i) exp_face_area_view(i) = 0.25;
    Kokkos::View<Amanzi::AmanziGeometry::Point*> exp_face_centroid_view("efcv",
                                                                        36);
    exp_face_centroid_view(0) = { 0.0, 0.25, 0.25 };
    exp_face_centroid_view(1) = { 0.0, 0.75, 0.25 };
    exp_face_centroid_view(2) = { 0.0, 0.25, 0.75 };
    exp_face_centroid_view(3) = { 0.0, 0.75, 0.75 };

    exp_face_centroid_view(4) = { 0.5, 0.25, 0.25 };
    exp_face_centroid_view(5) = { 0.5, 0.75, 0.25 };
    exp_face_centroid_view(6) = { 0.5, 0.25, 0.75 };
    exp_face_centroid_view(7) = { 0.5, 0.75, 0.75 };

    exp_face_centroid_view(8) = { 1.0, 0.25, 0.25 };
    exp_face_centroid_view(9) = { 1.0, 0.75, 0.25 };
    exp_face_centroid_view(10) = { 1.0, 0.25, 0.75 };
    exp_face_centroid_view(11) = { 1.0, 0.75, 0.75 };

    exp_face_centroid_view(12) = { 0.25, 0.0, 0.25 };
    exp_face_centroid_view(13) = { 0.75, 0.0, 0.25 };
    exp_face_centroid_view(14) = { 0.25, 0.0, 0.75 };
    exp_face_centroid_view(15) = { 0.75, 0.0, 0.75 };

    exp_face_centroid_view(16) = { 0.25, 0.5, 0.25 };
    exp_face_centroid_view(17) = { 0.75, 0.5, 0.25 };
    exp_face_centroid_view(18) = { 0.25, 0.5, 0.75 };
    exp_face_centroid_view(19) = { 0.75, 0.5, 0.75 };

    exp_face_centroid_view(20) = { 0.25, 1.0, 0.25 };
    exp_face_centroid_view(21) = { 0.75, 1.0, 0.25 };
    exp_face_centroid_view(22) = { 0.25, 1.0, 0.75 };
    exp_face_centroid_view(23) = { 0.75, 1.0, 0.75 };

    exp_face_centroid_view(24) = { 0.25, 0.25, 0.0 };
    exp_face_centroid_view(25) = { 0.75, 0.25, 0.0 };
    exp_face_centroid_view(26) = { 0.25, 0.75, 0.0 };
    exp_face_centroid_view(27) = { 0.75, 0.75, 0.0 };

    exp_face_centroid_view(28) = { 0.25, 0.25, 0.5 };
    exp_face_centroid_view(29) = { 0.75, 0.25, 0.5 };
    exp_face_centroid_view(30) = { 0.25, 0.75, 0.5 };
    exp_face_centroid_view(31) = { 0.75, 0.75, 0.5 };

    exp_face_centroid_view(32) = { 0.25, 0.25, 1.0 };
    exp_face_centroid_view(33) = { 0.75, 0.25, 1.0 };
    exp_face_centroid_view(34) = { 0.25, 0.75, 1.0 };
    exp_face_centroid_view(35) = { 0.75, 0.75, 1.0 };

    int ncells = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED);
    int nfaces = mesh->num_entities(Amanzi::AmanziMesh::FACE,
                                    Amanzi::AmanziMesh::Parallel_type::ALL);
    int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                                    Amanzi::AmanziMesh::Parallel_type::ALL);

    int space_dim_ = 3;

    Amanzi::AmanziMesh::Mesh* m = mesh.get();

    Kokkos::parallel_for(ncells, KOKKOS_LAMBDA(const int& i) {
      Amanzi::AmanziGeometry::Point centroid = m->cell_centroid(i);
      // Search for a cell with the same centroid in the
      // expected list of centroid
      bool found = false;
      for (int j = 0; j < ncells; j++) {
        if (fabs(exp_cell_centroid_view(j)[0] - centroid[0]) < 1.0e-10 &&
            fabs(exp_cell_centroid_view(j)[1] - centroid[1]) < 1.0e-10 &&
            fabs(exp_cell_centroid_view(j)[2] - centroid[2]) < 1.0e-10) {
          found = true;
          assert(exp_cell_volume_view(j) == m->cell_volume(i));
          break;
        }
      }
      assert(found == true);

      Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cfaces;
      Amanzi::AmanziGeometry::Point normal_sum(3), normal(3);

      m->cell_get_faces(i, cfaces);
      normal_sum.set(0.0);

      for (int j = 0; j < cfaces.extent(0); j++) {
        normal = m->face_normal(cfaces(j), false, i);
        normal_sum += normal;
      }

      double val = L22(normal_sum);
      assert(val < 1.0e-20);
    });


    // for(int i = 0 ; i < nfaces; ++i){
    Kokkos::parallel_for(nfaces, KOKKOS_LAMBDA(const int& i) {
      Amanzi::AmanziGeometry::Point centroid = m->face_centroid(i);

      bool found = false;

      for (int j = 0; j < nfaces; j++) {
        if (fabs(exp_face_centroid_view(j)[0] - centroid[0]) < 1.0e-10 &&
            fabs(exp_face_centroid_view(j)[1] - centroid[1]) < 1.0e-10 &&
            fabs(exp_face_centroid_view(j)[2] - centroid[2]) < 1.0e-10) {
          found = true;
          assert(exp_face_area_view(j) == m->face_area(i));
          // Check the natural normal
          Amanzi::AmanziGeometry::Point normal = m->face_normal(i);
          // Check the normal with respect to each connected cell
          Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cellids;
          m->face_get_cells(i, Amanzi::AmanziMesh::Parallel_type::ALL, cellids);

          for (int k = 0; k < cellids.extent(0); k++) {
            int dir;
            Amanzi::AmanziGeometry::Point normal_wrt_cell =
              m->face_normal(i, false, cellids(k), &dir);
            Amanzi::AmanziGeometry::Point normal1(normal);
            normal1 *= dir;

            for (int sd = 0; sd < space_dim_; ++sd)
              assert(normal1[sd] == normal_wrt_cell[sd]);

            Amanzi::AmanziGeometry::Point cellcentroid =
              m->cell_centroid(cellids(k));
            Amanzi::AmanziGeometry::Point facecentroid = m->face_centroid(i);
            Amanzi::AmanziGeometry::Point outvec = facecentroid - cellcentroid;
            double dp = outvec * normal_wrt_cell;
            dp /= (norm(outvec) * norm(normal_wrt_cell));
            assert(fabs(dp - 1.0) < 1e-10);
          }
          break;
        }
      }
      assert(found == true);
    });

    // Now deform the mesh a little and verify that the sum of the
    // outward normals of all faces of cell is still zero

    Amanzi::AmanziGeometry::Point ccoords(3);
    m->node_get_coordinates(13, &ccoords); // central node

    // Lets be sure this is the central node
    assert(ccoords[0] == 0.5);
    assert(ccoords[1] == 0.5);
    assert(ccoords[2] == 0.5);

    // Perturb it
    ccoords.set(0.7, 0.7, 0.7);
    m->node_set_coordinates(13, ccoords);

    // Now check the normals

    Kokkos::parallel_for(ncells, KOKKOS_LAMBDA(const int& i) {
      Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cfaces;
      Amanzi::AmanziGeometry::Point normal_sum(3), normal(3);

      m->cell_get_faces(i, cfaces);
      normal_sum.set(0.0);

      for (int j = 0; j < cfaces.extent(0); j++) {
        normal = m->face_normal(cfaces(j), false, i);
        normal_sum += normal;
      }

      double val = L22(normal_sum);
      assert(val < 1.0e-20);
    });

  } // for each framework i
}
