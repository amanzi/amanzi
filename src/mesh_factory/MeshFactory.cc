// -------------------------------------------------------------
// file: MeshFactory.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 10, 2011 by William A. Perkins
// Last Change: Tue Apr  5 12:44:58 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <boost/format.hpp>

#include "MeshFactory.hh"
#include "MeshFileType.hh"
#include "FrameworkTraits.hh"

namespace Amanzi{
namespace AmanziMesh {
  

  // -------------------------------------------------------------
  //  class MeshFactory
  // -------------------------------------------------------------

  // -------------------------------------------------------------
  // MeshFactory:: constructors / destructor
  // -------------------------------------------------------------
  MeshFactory::MeshFactory(const Epetra_MpiComm &communicator)
    : my_comm(communicator), my_preference(default_preference())
  {
  
  }

  MeshFactory::~MeshFactory(void)
  {
  }

  // -------------------------------------------------------------
  // MeshFactory::preference
  // -------------------------------------------------------------
  /** 
   * local -- but better be the same on all processes
   *
   * This routine populates the framework preference list, but only
   * with available frameworks.  If none of the preferred frameworks
   * are available, the preference list is left empty and an exception
   * is thrown.
   * 
   * @param pref list of mesh framework preferences
   */
  void
  MeshFactory::preference(const FrameworkPreference& pref)
  {
    my_preference.clear();
    my_preference = available_preference(pref);
    if (my_preference.empty()) {
      Message e("specified framework(s) not available: ");
      for (FrameworkPreference::const_iterator i = pref.begin(); 
           i != pref.end(); i++) {
        e.add_data(framework_name(*i).c_str());
        e.add_data(" ");
        amanzi_throw(e);
      }
    }
  }

  // -------------------------------------------------------------
  // MeshFactory::create
  // -------------------------------------------------------------
  /** 
   * Collective
   *
   * This creates a mesh by reading the specified file (or file set).  
   * 
   * @param filename mesh file to read
   * 
   * @return mesh instance
   */
  Teuchos::RCP<Mesh> 
  MeshFactory::create(const std::string& filename)
  {
    // check the file format
    Format fmt = file_format(my_comm, filename);

    if (fmt == UnknownFormat) {
      FileMessage 
        e(boost::str(boost::format("%s: unknown file format") %
                     filename).c_str());
      amanzi_throw(e);
    }
      
    Message e("MeshFactory::create: error: ");
    int ierr[1];
    ierr[0] = 0;

    Teuchos::RCP<Mesh> result;
    for (FrameworkPreference::const_iterator i = my_preference.begin(); 
         i != my_preference.end(); i++) {
      if (framework_reads(*i, fmt, my_comm.NumProc() > 1)) {
        try {
          result = framework_read(my_comm, *i, filename);
          return result;
        } catch (const Message& msg) {
          ierr[0] += 1;
          e.add_data(msg.what());
        } catch (const std::exception& stde) {
          ierr[0] += 1;
          e.add_data("internal error: ");
          e.add_data(stde.what());
        }
        my_comm.SumAll(ierr, ierr, 1);
        if (ierr[0] > 0) amanzi_throw(e);
      }
    }
    e.add_data(boost::str(boost::format("%s: unable to read mesh file") %
                          filename).c_str());
    amanzi_throw(e);
  }

  /** 
   * Coellective
   *
   * This creates a mesh by generating a block of hexahedral cells.
   *
   * Hopefully, if any one process has an error, all processes will
   * throw an Mesh::Message exception.
   * 
   * @param x0 origin x-coordinate
   * @param y0 origin y-coordinate
   * @param z0 origin z-coordinate
   * @param x1 maximum x-coordinate
   * @param y1 maximum y-coordinate
   * @param z1 maximum z-coordinate
   * @param nx number of cells in the x-direction
   * @param ny number of cells in the y-direction
   * @param nz number of cells in the z-direction
   * 
   * @return mesh instance
   */
  Teuchos::RCP<Mesh> 
  MeshFactory::create(double x0, double y0, double z0,
                      double x1, double y1, double z1,
                      int nx, int ny, int nz)
  {
    Teuchos::RCP<Mesh> result;
    Message e("MeshFactory::create: error: ");
    int ierr[1];
    ierr[0] = 0;
    if (nx <= 0 || ny <= 0 || nz <= 0) {
      ierr[0] += 1;
      e.add_data(boost::str(boost::format("invalid mesh cells requested: %d x %d x %d") %
                            nx % ny % nz).c_str());
    }
    my_comm.SumAll(ierr, ierr, 1);
    if (ierr[0] > 0) amanzi_throw(e);

    if (x1 - x0 <= 0.0 || y1 - y0 <= 0.0 || z1 - z0 <= 0.0) {
      ierr[0] += 1;
      e.add_data(boost::str(boost::format("invalid mesh dimensions requested: %.6g x %.6g x %.6g") %
                            (x1 - x0) % (y1 - y0) % (z1 - z0)).c_str());
    }
    my_comm.SumAll(ierr, ierr, 1);
    if (ierr[0] > 0) amanzi_throw(e);
      
    for (FrameworkPreference::const_iterator i = my_preference.begin(); 
         i != my_preference.end(); i++) {
      if (framework_generates(*i, my_comm.NumProc() > 1)) {
        try {
          result = framework_generate(my_comm, *i, 
                                      x0, y0, z0, x1, y1, z1, 
                                      nx, ny, nz);
          return result;
        } catch (const Message& msg) {
          ierr[0] += 1;
          e.add_data(msg.what());
        } catch (const std::exception& stde) {
          ierr[0] += 1;
          e.add_data("internal error: ");
          e.add_data(stde.what());
        }
        my_comm.SumAll(ierr, ierr, 1);
        if (ierr[0] > 0) amanzi_throw(e);
      }
    }
    e.add_data("unable to generate mesh");
    amanzi_throw(e);
  }

  /** 
   * This creates a mesh by generating a block of hexahedral cells,
   * but using a parameter list with the limits and cell counts.
   * 
   * @param parameter_list 
   * 
   * @return 
   */
  Teuchos::RCP<Mesh> 
  MeshFactory::create(Teuchos::ParameterList &parameter_list)
  {
    double x0, y0, z0, x1, y1, z1;
    int nx, ny, nz, nblock;

    // read the parameters from the parameter list
    
    nx = parameter_list.get<int>("Numer of Cells in X");
    ny = parameter_list.get<int>("Numer of Cells in Y");
    nz = parameter_list.get<int>("Numer of Cells in Z");
    
    x0 = parameter_list.get<double>("X_Min");
    x1 = parameter_list.get<double>("X_Max");
    
    y0 = parameter_list.get<double>("Y_Min");
    y1 = parameter_list.get<double>("Y_Max");
    
    z0 = parameter_list.get<double>("Z_Min");
    z1 = parameter_list.get<double>("Z_Max");
    
    // Note that this is ignored.  It doesn't seem to be used by
    // Mesh_maps_simple anyway.
    nblock = parameter_list.get<int>("Number of mesh blocks", 0);
      
    return create(x0, y0, z0, x1, y1, z1, nx, ny, nz);
  }

} // namespace AmanziMesh
} // namespace Amanzi
