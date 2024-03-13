/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William A. Perkins
*/

#include <filesystem>
#include <iostream>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Epetra_Vector.h"
#include "AmanziComm.hh"

#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "MeshException.hh"
#include "VerboseObject_objs.hh"

#include "HDF5_MPI.hh"

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
dump_output(const int& me,
            Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh,
            const std::string& filenameout)
{
  Amanzi::HDF5_MPI* viz_output = new Amanzi::HDF5_MPI(mesh->getComm());
  viz_output->setTrackXdmf(true);
  viz_output->createMeshFile(mesh, filenameout);
  viz_output->createDataFile(filenameout);

  const Epetra_Map& cmap(mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false));
  Epetra_Vector part(cmap);
  int nmycell(cmap.NumMyElements());
  std::vector<int> myidx(nmycell, 0);

  viz_output->createTimestep(0.0, 0, "");

  std::vector<double> mypart(nmycell, static_cast<double>(me + 1.0));
  for (unsigned int i = 0; i < nmycell; i++) myidx[i] = i;
  part.ReplaceMyValues(nmycell, &mypart[0], &myidx[0]);

  viz_output->open_h5file();
  viz_output->writeCellDataReal(part, "Partition");

  std::fill(mypart.begin(), mypart.end(), -1);

  part.ReplaceMyValues(nmycell, &mypart[0], &myidx[0]);
  viz_output->writeCellDataReal(part, "Block");

  viz_output->endTimestep();
  viz_output->close_h5file();

  delete viz_output;
}


// -------------------------------------------------------------
// do_the_audit
// -------------------------------------------------------------
int
do_the_audit(const int& me, Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh, const std::string& name)
{
  int lresult(0);

  std::stringstream ss;
  ss << name << std::setw(4) << std::setfill('0') << me << ".txt";
  std::string ofile = ss.str();


  std::ofstream ofs(ofile.c_str());
  if (me == 0) std::cout << "Writing results to " << ofile.c_str() << ", etc." << std::endl;
  Amanzi::AmanziMesh::MeshAudit audit(mesh, ofs);
  lresult = audit.Verify();

  int gresult;

  mesh->getComm()->MaxAll(&lresult, &gresult, 1);

  return gresult;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char** argv)
{
  std::filesystem::path progpath(argv[0]);
  std::string progname = progpath.filename();

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();

  {
    auto comm = Amanzi::getDefaultComm();
    const int me(comm->MyPID());

    int xcells(4), ycells(4), zcells(4);
    double xdelta(1.0), ydelta(1.0), zdelta(1.0);
    double xorigin(0.0), yorigin(0.0), zorigin(0.0);
    std::string inname;
    std::string outname("verify_hex");
    std::string outfilename;

    bool dosimple(false);
    bool doaudit(true);

    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString("The framework test for hex meshes.\n");

    try {
      CLP.setOption("audit", "noaudit", &doaudit, "do not audit the generated mesh");

      std::string framework("SIMPLE");
      CLP.setOption("framework", &framework, "mesh framework");
      dosimple = (framework == "SIMPLE");

      CLP.setOption("xcells", &xcells, "number of cells in the x-direction");
      CLP.setOption("ycells", &ycells, "number of cells in the y-direction");
      CLP.setOption("zcells", &zcells, "number of cells in the z-direction");

      CLP.setOption("xdelta", &xdelta, "cell size in the x-direction");
      CLP.setOption("ydelta", &ydelta, "cell size in the y-direction");
      CLP.setOption("zdelta", &zdelta, "cell size in the z-direction");

      CLP.setOption("xorigin", &xorigin, "x origin");
      CLP.setOption("yorigin", &yorigin, "y origin");
      CLP.setOption("zorigin", &zorigin, "z origin");

      // CLP.setOption("xml-file", &zorigin, "XML file from which to parameters");
      CLP.setOption("output", &outname, "output file base name");
      outfilename = outname;

      CLP.throwExceptions(false);
      CLP.recogniseAllOptions(true);

      Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = CLP.parse(argc, argv);

      if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
        throw std::string("Program not run");

      if (parseReturn == Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION)
        throw std::string("Program not run");

      if (parseReturn == Teuchos::CommandLineProcessor::PARSE_ERROR)
        throw std::string("Program not run");

    } catch (...) {
      if (me == 0) {
        std::cerr << progname << ": command line error." << std::endl;
        std::cerr << "Usage: " << progname << " --help" << std::endl;
      }
      return 1;
    }

    // generate a mesh
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
    Amanzi::AmanziMesh::Preference pref;
    if (dosimple) { pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE); }
    meshfactory.set_preference(pref);

    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
    if (inname.empty()) {
      if (me == 0) {
        std::cout << "Mesh: " << xcells << " x " << ycells << " x " << zcells << std::endl;
      }
      mesh = meshfactory.create(xorigin,
                                yorigin,
                                zorigin,
                                xorigin + xdelta * xcells,
                                yorigin + ydelta * ycells,
                                zorigin + zdelta * zcells,
                                xcells,
                                ycells,
                                zcells);
    }

    // make sure it's OK
    if (doaudit) {
      int notok(do_the_audit(me, mesh, outname));
      if (me == 0) {
        std::cout << "Mesh \"" << outname << "\" " << (notok ? "has errors" : "OK") << std::endl;
      }
    }

    // dump it out
    dump_output(me, mesh, outfilename);
  }

  Kokkos::finalize();
  return 0;
}
