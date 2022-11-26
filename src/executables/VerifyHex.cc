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
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem/path.hpp>
namespace bf = boost::filesystem;
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <Teuchos_GlobalMPISession.hpp>
#include <Epetra_Vector.h>
#include <AmanziComm.hh>

#include "Teuchos_ParameterXMLFileReader.hpp"
// DEPRECATED #include "Teuchos_XMLParameterListHelpers.hpp"

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "MeshException.hh"
#include "VerboseObject_objs.hh"

#include "HDF5_MPI.hh"

// -------------------------------------------------------------
// grab_filename
// -------------------------------------------------------------
/**
  * Simple routine to parse the filename from a path 
  * boost::filesystem object that handles the differences between
  * version 2 (path.leaf()) and 3 (path.filename())
  *
  * @param some_path a boost::filesystem path 
  * 
  * Return string that defines the filename
  */
std::string
grab_filename(const bf::path& some_path)
{
#if BOOST_FILESYSTEM_VERSION == 2
  return some_path.leaf();
#elif BOOST_FILESYSTEM_VERSION == 3
  return some_path.filename().generic_string();
#else
#  error Invalid Boost Filesystem library version
#endif
}

// -------------------------------------------------------------
// dump_output
// -------------------------------------------------------------
/** 
 * Dump a viz file using the mesh specified by @c maps. Include a
 * solution field for cells that identifies which process owns them.
 * 
 * @param me this process' id
 * @param maps mesh to output
 * @param filenameout HDF5/XDMF format file to produce
 */
void
dump_output(const int& me, Amanzi::AmanziMesh::Mesh& mesh, const std::string& filenameout)
{
  Amanzi::HDF5_MPI* viz_output = new Amanzi::HDF5_MPI(mesh.get_comm());
  viz_output->setTrackXdmf(true);
  viz_output->createMeshFile(Teuchos::rcp(&mesh), filenameout);
  viz_output->createDataFile(filenameout);

  const Epetra_Map& cmap(mesh.cell_map(false));
  Epetra_Vector part(cmap);
  int nmycell(cmap.NumMyElements());
  std::vector<int> myidx(nmycell, 0);

  viz_output->createTimestep(0.0, 0, "");

  std::vector<double> mypart(nmycell, static_cast<double>(me + 1.0));
  for (unsigned int i = 0; i < nmycell; i++) myidx[i] = i;
  part.ReplaceMyValues(nmycell, &mypart[0], &myidx[0]);

  viz_output->writeCellDataReal(part, "Partition");

  std::fill(mypart.begin(), mypart.end(), -1);

  // Amanzi::AmanziMesh::Set_ID_List setids;
  // mesh.get_set_ids(Amanzi::AmanziMesh::CELL, &setids);
  // for (Amanzi::AmanziMesh::Set_ID_List::const_iterator i = setids.begin();
  //      i != setids.end(); ++i) {
  //   Amanzi::AmanziMesh::Entity_ID_List gids;
  //   mesh.get_set_entities(*i, Amanzi::AmanziMesh::CELL,
  //                          Amanzi::AmanziMesh::Parallel_type::OWNED, &gids);
  //   for (Amanzi::AmanziMesh::Entity_ID_List::const_iterator g = gids.begin();
  //        g != gids.end(); ++g) {
  //     int lidx(*g);
  //     mypart[lidx] = *i;
  //     // std::cerr << me << ": set " << *i << ", cell " << *g << " (" << lidx << ")" << std::endl;
  //   }
  // }

  part.ReplaceMyValues(nmycell, &mypart[0], &myidx[0]);
  viz_output->writeCellDataReal(part, "Block");


  viz_output->endTimestep();

  delete viz_output;
}


// -------------------------------------------------------------
// do_the_audit
// -------------------------------------------------------------
int
do_the_audit(const int& me, Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh, const std::string& name)
{
  int lresult(0);

  std::string ofile = boost::str(boost::format("%s%04d.txt") % name % me);

  std::ofstream ofs(ofile.c_str());
  if (me == 0) std::cout << "Writing results to " << ofile.c_str() << ", etc." << std::endl;
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
main(int argc, char** argv)
{
  bf::path progpath(argv[0]);
  std::string progname = grab_filename(progpath);


  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  auto comm = Amanzi::getDefaultComm();
  const int me(comm->MyPID());

  unsigned int xcells(4), ycells(4), zcells(4);
  double xdelta(1.0), ydelta(1.0), zdelta(1.0);
  double xorigin(0.0), yorigin(0.0), zorigin(0.0);
  std::string inname;
  std::string outname("verify_hex");
  std::string outfilename;

  bool dosimple(false);
  bool doaudit(true);

  po::options_description desc("Available options");

  try {
    desc.add_options()("help", "produce this help message")

      ("noaudit", "do not audit the generated mesh")

        ("simple", "use the Mesh_Simple framework instead of stk::mesh")

          ("xcells",
           po::value<unsigned int>()->default_value(xcells),
           "number of cells in the x-direction")("ycells",
                                                 po::value<unsigned int>()->default_value(ycells),
                                                 "number of cells in the y-direction")(
            "zcells",
            po::value<unsigned int>()->default_value(zcells),
            "number of "
            "cells in "
            "the "
            "z-"
            "direction")

            ("xdelta", po::value<double>()->default_value(xdelta), "cell size in the x-direction")(
              "ydelta", po::value<double>()->default_value(ydelta), "cell size in the y-direction")(
              "zdelta",
              po::value<double>()->default_value(zdelta),
              "cell size in the "
              "z-direction")

              ("xorigin", po::value<double>()->default_value(xorigin), "x origin")(
                "yorigin", po::value<double>()->default_value(yorigin), "y origin")(
                "zorigin", po::value<double>()->default_value(zorigin), "z origin")

                ("xml-file", po::value<std::string>(), "XML file from which to parameters")

                  ("output",
                   po::value<std::string>()->default_value(outname),
                   "output file base name");

    po::positional_options_description p;
    p.add("output", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cerr << "Usage: " << progname << " [options]" << std::endl;
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
    outfilename = outname;

  } catch (po::error& e) {
    if (me == 0) {
      std::cerr << boost::str(boost::format("%s: command line error: %s") % progname % e.what())
                << std::endl;
      std::cerr << boost::str(boost::format("Usage: %s [options]") % progname) << std::endl;
      std::cerr << desc << std::endl;
    }
    return 3;
  } catch (boost::bad_any_cast& e) {
    if (me == 0) {
      std::cerr << boost::str(boost::format("%s: command line error: %s") % progname % e.what())
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

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  Amanzi::AmanziMesh::Preference pref;
  if (dosimple) {
    pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
  } else {
    pref.push_back(Amanzi::AmanziMesh::Framework::STK);
  }
  meshfactory.set_preference(pref);

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  if (inname.empty()) {
    mesh = meshfactory.create(xorigin,
                              yorigin,
                              zorigin,
                              xorigin + xdelta * xcells,
                              yorigin + ydelta * ycells,
                              zorigin + zdelta * zcells,
                              xcells,
                              ycells,
                              zcells);
    // } else {
    //   Teuchos::ParameterList parameter_list;

    //   int ierr(0), aerr(0);

    //   try {
    //     Teuchos::ParameterXMLFileReader xmlreader(inname);
    //     Teuchos::ParameterList all_parameter_list(xmlreader.getParameters());
    //     Teuchos::ParameterList mesh_parameter_list = all_parameter_list.sublist("mesh");
    //     parameter_list = mesh_parameter_list.sublist("Generate");
    //   } catch (const std::runtime_error& e) {
    //     std::cerr << me << ": error parsing xml-file: " << e.what() << std::endl;
    //     ierr++;
    //   }

    //   comm->SumAll(&ierr, &aerr, 1);
    //   if (aerr > 0) {
    //     return 1;
    //   }

    //   mesh = meshfactory(parameter_list);
  }

  //  std::cout << "Generated mesh has "
  //            << mesh->num_sets(Amanzi::AmanziMesh::CELL)
  //            << " cell sets"
  //            << std::endl;

  // make sure it's OK

  if (doaudit) {
    int notok(do_the_audit(me, mesh, outname));
    if (me == 0) {
      std::cout << "Mesh \"" << outname << "\" " << (notok ? "has errors" : "OK") << std::endl;
    }
  }

  // dump it out

  dump_output(me, *mesh, outfilename);

  return 0;
}
