/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   MeshFactory.hh
 * @author William A. Perkins
 * @date Mon Aug  1 12:12:13 2011
 * 
 * @brief  declaration of the MeshFactory class
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 10, 2011 by William A. Perkins
// Last Change: Mon Aug  1 12:12:13 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _MeshFactory_hh_
#define _MeshFactory_hh_

#include <string>
#include <vector>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_ParameterList.hpp>

#include "MeshException.hh"
#include "MeshFramework.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

  // -------------------------------------------------------------
  //  class MeshFactory
  // -------------------------------------------------------------
  class MeshFactory {
  protected:

    /// The parallel environment
    Epetra_MpiComm my_comm;

    /// A list of preferred mesh frameworks to consider
    FrameworkPreference my_preference;

  private:

    /// private, undefined copy constructor to avoid unwanted copies
    MeshFactory(MeshFactory& old);

  public:

    /// Default constructor.
    explicit MeshFactory(const Epetra_MpiComm& communicator);

    /// Destructor
    ~MeshFactory(void);

    /// Get the framework preference
    const FrameworkPreference& preference(void) const
    { return my_preference; }

    /// Set the framework preference
    void preference(const FrameworkPreference& pref);

    /// Create a mesh by reading the specified file (or set of files)
    Teuchos::RCP<Mesh> create(const std::string& filename);


    /// Create a hexahedral mesh of the specified dimensions
    Teuchos::RCP<Mesh> create(double x0, double y0, double z0,
                              double x1, double y1, double z1,
                              int nx, int ny, int nz);
    
    /// Create a hexahedral mesh using the specified parameter list
    Teuchos::RCP<Mesh> create(Teuchos::ParameterList &parameter_list);

    /// Create a mesh by reading the specified file (or set of files) -- operator
    Teuchos::RCP<Mesh> operator() (const std::string& filename)
    {
      return create(filename);
    }
  
    /// Create a hexahedral mesh of the specified dimensions -- operator
    Teuchos::RCP<Mesh> operator() (double x0, double y0, double z0,
                                             double x1, double y1, double z1,
                                             int nx, int ny, int nz)
    { 
      return create(x0, y0, z0, x1, y1, z1, nx, ny, nz);
    }

    /// Create a hexahedral mesh using the specified parameter list
    Teuchos::RCP<Mesh> operator() (Teuchos::ParameterList &parameter_list)
    {
      return create(parameter_list);
    }


  };

} // namespace AmanziMesh
} // namespace Amanzi

#endif
