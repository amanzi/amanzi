/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! test utilities used in a few tests
#include <string>
#include <fstream>
#include <sstream>

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "MeshFactory.hh"
#include "AmanziVector.hh"
#include "CompositeVector.hh"
#include "Output.hh"
#include "Input.hh"

//
// Reads a file and returns a string
//
std::string
read_text_file(const std::string& filename)
{
  std::ifstream t(filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  return std::move(buffer.str());
}

//
// Compares (text) two files
//
bool
compareTextFiles(const std::string& p1, const std::string& p2)
{
  std::cout << " ... comparing text files:" << std::endl
            << "       - " << p1 << std::endl
            << "       - " << p2 << std::endl;
  auto s1 = read_text_file(p1);
  auto s2 = read_text_file(p2);
  auto res = (s1 == s2);
  std::cout << "    " << (res ? "matches" : "differs") << std::endl;
  return res;
}

//
// Compares (bitwise) two arbitrary binary files.
//
bool
compareBinaryFiles(const std::string& p1, const std::string& p2)
{
  std::cout << " ... comparing binary files:" << std::endl
            << "       - " << p1 << std::endl
            << "       - " << p2 << std::endl;

  std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    std::cout << "   file problem" << std::endl;
    return false; // file problem
  }

  if (f1.tellg() != f2.tellg()) {
    std::cout << "   binary size mismatch" << std::endl;
    return false; // size mismatch
  }

  // seek back to beginning and use std::equal to compare contents
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  auto res = std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                        std::istreambuf_iterator<char>(),
                        std::istreambuf_iterator<char>(f2.rdbuf()));
  std::cout << "   " << (res ? "matches" : "differs") << std::endl;
  return res;
}


bool
compareH5Files(const std::string& p1, const std::string& p2)
{
  std::cout << " ... comparing h5 files:" << std::endl
            << "       - " << p1 << std::endl
            << "       - " << p2 << std::endl;
  std::string command = std::string("h5diff ") + p1 + " " + p2;
  int ierr = system(command.c_str());
  bool res = (ierr == 0);
  std::cout << "   " << (res ? "matches" : "differs") << std::endl;
  return res;
}


using namespace Amanzi;

struct output_test_harness {
  Comm_ptr_type comm;
  Teuchos::RCP<AmanziMesh::MeshHost> mesh;
  Teuchos::RCP<Vector_type> cell_quantity;
  Teuchos::RCP<Vector_type> node_quantity;
  Teuchos::RCP<MultiVector_type> multi;
  Teuchos::RCP<IntVector_type> indices;

  output_test_harness() { comm = getDefaultComm(); }

  void create_data()
  {
    unsigned int num_nodes = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE,
                                                  Amanzi::AmanziMesh::Parallel_kind::OWNED);
    unsigned int num_cells = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                                  Amanzi::AmanziMesh::Parallel_kind::OWNED);

    // Setup node quantity
    std::vector<int> node_index_list{ 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
    std::vector<double> node_values{ 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120 };
    auto node_map = mesh->getMap(Amanzi::AmanziMesh::Entity_kind::NODE, false);
    node_quantity = Teuchos::rcp(new Vector_type(node_map));
    for (int i = 0; i != node_index_list.size(); ++i) {
      auto lid = node_map->getLocalElement(node_index_list[i]);
      if (lid >= 0) { node_quantity->replaceLocalValue(lid, node_values[i]); }
    }

    // Setup cell quantity
    std::vector<int> cell_index_list{ 0, 1, 2, 3 };
    std::vector<double> cell_values{ 10, 20, 30, 40 };
    auto cell_map = mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false);
    cell_quantity = Teuchos::rcp(new Vector_type(cell_map));
    for (int i = 0; i != cell_index_list.size(); ++i) {
      auto lid = cell_map->getLocalElement(cell_index_list[i]);
      if (lid >= 0) { cell_quantity->replaceLocalValue(lid, cell_values[i]); }
    }

    // multivector
    multi = Teuchos::rcp(new MultiVector_type(cell_map, 2));
    multi->getVectorNonConst(0)->putScalar(1.0);
    multi->getVectorNonConst(1)->putScalar(2.0);

    // intvector
    indices = Teuchos::rcp(new IntVector_type(cell_map));
    indices->putScalar(42);
  }


  void create_mesh_structured()
  {
    AmanziMesh::MeshFactory meshfactory(comm);
    AmanziMesh::Preference pref{ AmanziMesh::Framework::SIMPLE,
                                 AmanziMesh::Framework::MSTK,
                                 AmanziMesh::Framework::MOAB };
    meshfactory.set_preference(pref);

    mesh = AmanziMesh::onMemSpace<MemSpace_kind::HOST>(
      meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 1, 1));
  }

  void create_mesh_polyhedral()
  {
    AmanziMesh::MeshFactory meshfactory(comm);
    AmanziMesh::Preference pref{ AmanziMesh::Framework::MSTK, AmanziMesh::Framework::MOAB };
    meshfactory.set_preference(pref);
    mesh =
      AmanziMesh::onMemSpace<MemSpace_kind::HOST>(meshfactory.create("./test/four_polygon.exo"));
  }

  void test_write(const Output& out)
  {
    Teuchos::ParameterList node_plist("node_quantity");
    node_plist.set("location", AmanziMesh::NODE);
    out.write(node_plist, *node_quantity);

    Teuchos::ParameterList cell_plist("cell_quantity");
    cell_plist.set("location", AmanziMesh::CELL);
    out.write(cell_plist, *cell_quantity);

    Teuchos::ParameterList multi_plist("multivector");
    multi_plist.set("location", AmanziMesh::CELL);
    std::vector<std::string> subfieldnames{ "one", "two" };
    multi_plist.set<Teuchos::Array<std::string>>("subfieldnames", subfieldnames);
    out.write(multi_plist, *multi);

    Teuchos::ParameterList int_plist("int_quantity");
    int_plist.set("location", AmanziMesh::CELL);
    out.write(int_plist, *indices);

    out.write(Teuchos::ParameterList("six_one"), 6.1);
    out.write(Teuchos::ParameterList("seven"), 7);
    out.write(Teuchos::ParameterList("name"), "my name");
  }

  void test_read(const Input& in)
  {
    std::string nm;
    in.read(Teuchos::ParameterList("name"), nm);
    CHECK_EQUAL(std::string("name"), nm);

    int sv;
    in.read(Teuchos::ParameterList("name"), sv);
    CHECK_EQUAL(7, sv);

    double db;
    in.read(Teuchos::ParameterList("name"), db);
    CHECK_CLOSE(6.1, db, 1.e-12);

    Teuchos::ParameterList cell_plist("cell_quantity");
    cell_plist.set("location", AmanziMesh::CELL);
    Vector_type cells(*cell_quantity);
    cells.putScalar(0.);
    in.read(cell_plist, cells);
    cells.update(-1, *cell_quantity, 1);
    CHECK_CLOSE(0, cells.normInf(), 1e-10);

    Teuchos::ParameterList node_plist("node_quantity");
    node_plist.set("location", AmanziMesh::NODE);
    Vector_type nodes(*node_quantity);
    nodes.putScalar(0.);
    in.read(node_plist, nodes);
    nodes.update(-1, *node_quantity, 1);
    CHECK_CLOSE(0, nodes.normInf(), 1e-10);

    cell_plist.setName("indices");
    IntVector_type inds(*indices);
    inds.putScalar(0);
    in.read(cell_plist, inds);
    inds.update(-1, *indices, 1);
    CHECK_EQUAL(0, inds.normInf());

    cell_plist.setName("multivector");
    std::vector<std::string> subfieldnames{ "one", "two" };
    cell_plist.set<Teuchos::Array<std::string>>("subfieldnames", subfieldnames);
    MultiVector_type multiv(*multi);
    multiv.putScalar(0);
    in.read(cell_plist, multiv);
    multiv.update(-1, *multi, 1);
    Teuchos::Array<double> norms(multiv.getNumVectors());
    multiv.normInf(norms);
    CHECK_CLOSE(0, norms[0], 1e-10);
    CHECK_CLOSE(0, norms[1], 1e-10);
  }

  void test_write_multiple(Output& out)
  {
    double time = 0.0;
    int NITS = 5;
    for (int i = 0; i < NITS; i++) {
      std::cout << "writing iteration... " << i << std::endl;

      // write time step data
      out.createTimestep(time, i);
      test_write(out);
      // close file
      out.finalizeTimestep();

      // advance time and values
      time += 2.0;
      double scalar = ((double)i + 2) / (i + 1);
      int scalar_int = i;
      node_quantity->scale(scalar);
      cell_quantity->scale(scalar);
      indices->scale(scalar_int);
      multi->getVectorNonConst(0)->scale(scalar);
      multi->getVectorNonConst(1)->scale(scalar + 1);
    }
  }
};
