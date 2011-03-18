// -------------------------------------------------------------
// file: MeshFactory.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 10, 2011 by William A. Perkins
// Last Change: Fri Mar 18 15:06:54 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include "MeshFactory.hh"
#include "MeshFileType.hh"
#include "FrameworkTraits.hh"
namespace Mesh {
  

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
  void
  MeshFactory::preference(const FrameworkPreference& pref)
  {
    my_preference.clear();
    my_preference = available_preference(pref);
    if (my_preference.empty()) {
      Mesh::Message e("specified framework(s) not available: ");
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
  Teuchos::RCP<Mesh_maps_base> 
  MeshFactory::create(const std::string& filename)
  {
    Format fmt = file_format(my_comm, filename);
    Teuchos::RCP<Mesh_maps_base> result;
    for (FrameworkPreference::const_iterator i = my_preference.begin(); 
         i != my_preference.end(); i++) {
      if (framework_reads(*i, fmt, my_comm.NumProc() > 1)) {
        result = framework_read(my_comm, *i, filename);
        return result;
      }
    }
    return result;
  }

  Teuchos::RCP<Mesh_maps_base> 
  MeshFactory::create(double x0, double y0, double z0,
                      double x1, double y1, double z1,
                      int nx, int ny, int nz)
  {
    Teuchos::RCP<Mesh_maps_base> result;
    for (FrameworkPreference::const_iterator i = my_preference.begin(); 
         i != my_preference.end(); i++) {
      if (framework_generates(*i, my_comm.NumProc() > 1)) {
        result = framework_generate(my_comm, *i, 
                                    x0, y0, z0, x1, y1, z1, 
                                    nx, ny, nz);
        return result;
      }
    }
    return result;
  }

} // namespace Mesh
