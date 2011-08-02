// Emacs Mode Line: -*- Mode:c++; c-default-style: "google"; indent-tabs-mode: nil -*-
// -------------------------------------------------------------
/**
 * @file   verify_hex.cc
 * @author William A. Perkins
 * @date Tue Aug  2 13:27:55 2011
 * 
 * @brief  A simple test of hex-mesh generation -- serial or parallel
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May 24, 2011 by William A. Perkins
// Last Change: Tue Aug  2 13:27:55 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <boost/format.hpp>
// TODO: We are using depreciated parts of boost::filesystem
#define BOOST_FILESYSTEM_VERSION 2
#include <boost/filesystem/path.hpp>
namespace bf = boost::filesystem;
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <Teuchos_GlobalMPISession.hpp>
#include <Epetra_Vector.h>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Epetra_MpiComm.h>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "MeshException.hh"

#include "cgns_mesh_par.hh"

// -------------------------------------------------------------
// dump_cgns
// -------------------------------------------------------------
/** 
 * Dump a CGNS file using the mesh specified by @c maps. Include a
 * solution field for cells that identifies which process owns them.
 * 
 * @param me this process' id
 * @param maps mesh to output
 * @param cgnsout CGNS format file to produce
 */
void
dump_cgns(const int& me, Amanzi::AmanziMesh::Mesh &mesh, const std::string& cgnsout)
{
  Amanzi::CGNS_PAR::create_mesh_file(mesh, cgnsout);

  const Epetra_Map& cmap(mesh.cell_epetra_map(false));
  Epetra_Vector part(cmap);
  int nmycell(cmap.NumMyElements());
  std::vector<int> myidx(nmycell, 0);

  Amanzi::CGNS_PAR::open_data_file(cgnsout);
  Amanzi::CGNS_PAR::create_timestep(0.0, 0, Amanzi::AmanziMesh::CELL);

  std::vector<double> mypart(nmycell, static_cast<double>(me+1.0));
  for (unsigned int i = 0; i < nmycell; i++) myidx[i] = i;
  part.ReplaceMyValues(nmycell, &mypart[0], &myidx[0]);
  
  Amanzi::CGNS_PAR::write_field_data(part, "Partition");

  std::fill(mypart.begin(), mypart.end(), -1);

  Amanzi::AmanziMesh::Set_ID_List setids;
  mesh.get_set_ids(Amanzi::AmanziMesh::CELL, &setids);
  for (Amanzi::AmanziMesh::Set_ID_List::const_iterator i = setids.begin(); 
       i != setids.end(); ++i) {
    Amanzi::AmanziMesh::Entity_ID_List gids;
    mesh.get_set_entities(*i, Amanzi::AmanziMesh::CELL, 
                           Amanzi::AmanziMesh::OWNED, &gids);
    for (Amanzi::AmanziMesh::Entity_ID_List::const_iterator g = gids.begin();
         g != gids.end(); ++g) {
      int lidx(*g);
      mypart[lidx] = *i;
      // std::cerr << me << ": set " << *i << ", cell " << *g << " (" << lidx << ")" << std::endl;
    }
  }

  part.ReplaceMyValues(nmycell, &mypart[0], &myidx[0]);
  Amanzi::CGNS_PAR::write_field_data(part, "Block");


  Amanzi::CGNS_PAR::close_data_file();
}


// -------------------------------------------------------------
// do_the_audit
// -------------------------------------------------------------
int
do_the_audit(const int& me, 
             Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh, 
             const std::string& name)
{
  int lresult(0);

  std::string ofile =
    boost::str(boost::format("%s%04d.txt") % name % me);

  std::ofstream ofs(ofile.c_str());
  if (me == 0)
    std::cout << "Writing results to " << ofile.c_str() << ", etc." << std::endl;
  Amanzi::MeshAudit audit(mesh, ofs);
  lresult = audit.Verify();

  int gresult;

  mesh->get_comm()->MaxAll(&lresult, &gresult, 1);

  return gresult;

}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  bf::path progpath(argv[0]);
  std::string progname = progpath.leaf();

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());

  unsigned int xcells(4), ycells(4), zcells(4);
  double xdelta(1.0), ydelta(1.0), zdelta(1.0);
  double xorigin(0.0), yorigin(0.0), zorigin(0.0);
  std::string inname;
  std::string outname("verify_hex");
  std::string outcgns;

  bool dosimple(false);
  bool doaudit(true);

  po::options_description desc("Available options");

  try {

    desc.add_options()
      ("help", "produce this help message")

      ("noaudit", "do not audit the generated mesh")

      ("simple", "use the Mesh_Simple framework instead of stk::mesh")

      ("xcells", po::value<unsigned int>()->default_value(xcells), "number of cells in the x-direction")
      ("ycells", po::value<unsigned int>()->default_value(ycells), "number of cells in the y-direction")
      ("zcells", po::value<unsigned int>()->default_value(zcells), "number of cells in the z-direction")

      ("xdelta", po::value<double>()->default_value(xdelta), "cell size in the x-direction")
      ("ydelta", po::value<double>()->default_value(ydelta), "cell size in the y-direction")
      ("zdelta", po::value<double>()->default_value(zdelta), "cell size in the z-direction")

      ("xorigin", po::value<double>()->default_value(xorigin), "x origin")
      ("yorigin", po::value<double>()->default_value(yorigin), "y origin")
      ("zorigin", po::value<double>()->default_value(zorigin), "z origin")

      ("xml-file", po::value<std::string>(), "XML file from which to parameters")

      ("output", po::value<std::string>()->default_value(outname), "output file base name")
      ;

    po::positional_options_description p;
    p.add("output", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);    
  
    if (vm.count("help")) {
      std::cerr << boost::str(boost::format("Usage: %s [options]") % progname) << std::endl;
      std::cerr << desc << std::endl;
      return 3;
    }

    if (vm.count("noaudit")) doaudit = false;

    if (vm.count("xml-file") > 0) {
      inname = vm["xml-file"].as<std::string>();
    } else {
      inname.clear();

      xcells = vm["xcells"].as<unsigned int>();
      ycells = vm["ycells"].as<unsigned int>();
      zcells = vm["zcells"].as<unsigned int>();
      
      xdelta = vm["xdelta"].as<double>();
      ydelta = vm["ydelta"].as<double>();
      zdelta = vm["zdelta"].as<double>();
      
      xorigin = vm["xorigin"].as<double>();
      yorigin = vm["yorigin"].as<double>();
      zorigin = vm["zorigin"].as<double>();
    }

    dosimple = (vm.count("simple") > 0);

    outname = vm["output"].as<std::string>();
    outcgns = outname;
    outcgns += ".cgns";


  } catch (po::error& e) {
    if (me == 0) {
      std::cerr << boost::str(boost::format("%s: command line error: %s") %
                              progname % e.what()) 
                << std::endl;
      std::cerr << boost::str(boost::format("Usage: %s [options]") % progname) << std::endl;
      std::cerr << desc << std::endl;
    }
    return 3;
  } catch (boost::bad_any_cast& e) {
    if (me == 0) {
      std::cerr << boost::str(boost::format("%s: command line error: %s") %
                              progname % e.what()) 
                << std::endl;
      std::cerr << boost::str(boost::format("Usage: %s [options]") % progname) << std::endl;
      std::cerr << desc << std::endl;
    }
    return 3;
  } catch (...) {
    if (me == 0) {
      std::cerr << boost::str(boost::format("Usage: %s [options]") % progname) << std::endl;
      std::cerr << desc << std::endl;
    }
    return 3;
  }

  // generate a mesh

  Amanzi::AmanziMesh::MeshFactory factory(comm);
  Amanzi::AmanziMesh::FrameworkPreference pref;
  if (dosimple) {
    pref.push_back(Amanzi::AmanziMesh::Simple);
  } else {
    pref.push_back(Amanzi::AmanziMesh::STKMESH);
  }
  factory.preference(pref);

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  if (inname.empty()) {
    mesh = factory(xorigin, yorigin, zorigin,
                   xorigin+xdelta*xcells,
                   yorigin+ydelta*ycells,
                   zorigin+zdelta*zcells,
                   xcells, ycells, zcells);
  } else {
    Teuchos::ParameterList parameter_list;
    
    int ierr(0), aerr(0);

    try {
      Teuchos::ParameterList all_parameter_list;
      Teuchos::updateParametersFromXmlFile(inname, &all_parameter_list);
      Teuchos::ParameterList mesh_parameter_list = all_parameter_list.sublist("Mesh");
      parameter_list = mesh_parameter_list.sublist("Generate");
    } catch (const std::runtime_error& e) {
      std::cerr << me << ": error parsing xml-file: " << e.what() << std::endl;
      ierr++;
    }

    comm.SumAll(&ierr, &aerr, 1);
    if (aerr > 0) {
      return 1;
    }
   
    mesh = factory(parameter_list);
  }

  std::cout << "Generated mesh has " 
            << mesh->num_sets(Amanzi::AmanziMesh::CELL)
            << " cell sets" 
            << std::endl;

  // make sure it's OK

  if (doaudit) {
    int notok(do_the_audit(me, mesh, outname));
    if (me == 0) {
      std::cout << "Mesh \"" << outname << "\" " 
                << (notok ? "has errors" : "OK") << std::endl;
    }
  }

  // dump it out

  dump_cgns(me, *mesh, outcgns);

  return 0;
}

