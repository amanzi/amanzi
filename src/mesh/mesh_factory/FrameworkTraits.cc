/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   FrameworkTraits.cc
 * @author William A. Perkins
 * @date Tue Oct  4 06:14:05 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 14, 2011 by William A. Perkins
// Last Change: Tue Oct  4 06:14:05 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <boost/format.hpp>

#include <iostream>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>

namespace mpl = boost::mpl;

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

// How stupid is this? MOAB mesh stuff evidently cannot be included
// before Epetra_MpiComm.h

#ifdef HAVE_MOAB_MESH
#define MOAB_FLAG true
#define USE_MPI 
#include "Mesh_MOAB.hh"
#undef USE_MPI
#endif

#include "FrameworkTraits.hh"


#include "MeshFramework.hh"
#include "MeshFileType.hh"
#include "MeshException.hh"
#include "Mesh.hh"
#include "GenerationSpec.hh"

// -------------------------------------------------------------
//  class bogus_maps
// -------------------------------------------------------------
/**
 * This class is just used to throw an exception when a mesh framework
 * is used, but not available, or is available, but misused.
 * 
 */

class bogus_maps : public Amanzi::AmanziMesh::Mesh {
 public:
  
  /// Default constructor.
  bogus_maps(const char *filename, const Epetra_MpiComm *c,
             const Amanzi::AmanziGeometry::GeometricModelPtr& gm) 
      : Mesh(), bogus_map_(NULL) 
  {
    Exceptions::amanzi_throw(Errors::Message("reading not supported"));
  }
  
  bogus_maps(const char *filename, const Epetra_MpiComm *c, int dim,
             const Amanzi::AmanziGeometry::GeometricModelPtr& gm) 
      : Mesh(), bogus_map_(NULL) 
  {
    Exceptions::amanzi_throw(Errors::Message("reading not supported"));
  }
  
  bogus_maps(double x0, double y0, double z0,
	     double x1, double y1, double z1,
	     int nx, int ny, int nz, 
	     const Epetra_MpiComm *communicator,
             const Amanzi::AmanziGeometry::GeometricModelPtr& gm)
      : Amanzi::AmanziMesh::Mesh(), bogus_map_(NULL) 
  {
    Exceptions::amanzi_throw(Errors::Message("generation not supported"));
  }

  bogus_maps(double x0, double y0,
	     double x1, double y1,
	     int nx, int ny,
	     const Epetra_MpiComm *communicator,
             const Amanzi::AmanziGeometry::GeometricModelPtr& gm)
      : Amanzi::AmanziMesh::Mesh(), bogus_map_(NULL) 
  {
    Exceptions::amanzi_throw(Errors::Message("generation not supported"));
  }

  bogus_maps(const Amanzi::AmanziMesh::GenerationSpec& gspec, 
             const Epetra_MpiComm *communicator,
             const Amanzi::AmanziGeometry::GeometricModelPtr& gm)
      : Amanzi::AmanziMesh::Mesh(), bogus_map_(NULL) 
  {
    Exceptions::amanzi_throw(Errors::Message("generation not supported"));
  }

  // all of these virtual methods need to be implemented, even though
  // they'll never be used (because instantiating this class will
  // throw an exception)

  Amanzi::AmanziMesh::Parallel_type 
  entity_get_ptype(const Amanzi::AmanziMesh::Entity_kind kind, 
                   const Amanzi::AmanziMesh::Entity_ID entid) const
  { return Amanzi::AmanziMesh::OWNED; }

  Amanzi::AmanziMesh::Cell_type 
  cell_get_type(const Amanzi::AmanziMesh::Entity_ID cellid) const
  { return Amanzi::AmanziMesh::UNKNOWN; }

  unsigned int 
  num_entities (const Amanzi::AmanziMesh::Entity_kind kind,
                const Amanzi::AmanziMesh::Parallel_type ptype) const
  { return 0; }

  unsigned int 
  GID(const Amanzi::AmanziMesh::Entity_ID lid, 
      const Amanzi::AmanziMesh::Entity_kind kind) const
  { return 0; }

  void cell_get_faces_and_dirs (const Amanzi::AmanziMesh::Entity_ID cellid,
                                Amanzi::AmanziMesh::Entity_ID_List *faceids,
                                std::vector<int> *face_dirs,
                                const bool ordered=false) const
  {}


  void cell_get_nodes (const Amanzi::AmanziMesh::Entity_ID cellid, 
                       Amanzi::AmanziMesh::Entity_ID_List *nodeids) const
  {}

  void face_get_nodes (const Amanzi::AmanziMesh::Entity_ID faceid, 
                       Amanzi::AmanziMesh::Entity_ID_List *nodeids) const
  {}

  void node_get_cells (const Amanzi::AmanziMesh::Entity_ID nodeid, 
                       const Amanzi::AmanziMesh::Parallel_type ptype,
                       Amanzi::AmanziMesh::Entity_ID_List *cellids) const
  {}

  void node_get_faces (const Amanzi::AmanziMesh::Entity_ID nodeid, 
                       const Amanzi::AmanziMesh::Parallel_type ptype,
                       Amanzi::AmanziMesh::Entity_ID_List *faceids) const
  {}
    
  void node_get_cell_faces (const Amanzi::AmanziMesh::Entity_ID nodeid, 
                            const Amanzi::AmanziMesh::Entity_ID cellid,
                            const Amanzi::AmanziMesh::Parallel_type ptype,
                            Amanzi::AmanziMesh::Entity_ID_List *faceids) const
  {}
    
  void face_get_cells (const Amanzi::AmanziMesh::Entity_ID faceid, 
                       const Amanzi::AmanziMesh::Parallel_type ptype,
                       Amanzi::AmanziMesh::Entity_ID_List *cellids) const
  {}

  void cell_get_face_adj_cells(const Amanzi::AmanziMesh::Entity_ID cellid,
                               const Amanzi::AmanziMesh::Parallel_type ptype,
                               Amanzi::AmanziMesh::Entity_ID_List *fadj_cellids) const
  {}

  void cell_get_node_adj_cells(const Amanzi::AmanziMesh::Entity_ID cellid,
                               const Amanzi::AmanziMesh::Parallel_type ptype,
                               Amanzi::AmanziMesh::Entity_ID_List *nadj_cellids) const
  {}

  Amanzi::AmanziMesh::Cell_type 
  cell_get_type_4viz(const Amanzi::AmanziMesh::Entity_ID cellid) const
  { return Amanzi::AmanziMesh::UNKNOWN; }
    
  void 
  cell_get_nodes_4viz (const Amanzi::AmanziMesh::Entity_ID cellid, 
                       Amanzi::AmanziMesh::Entity_ID_List *nodeids) const
  {}

  void 
  node_get_coordinates (const Amanzi::AmanziMesh::Entity_ID nodeid, 
                        Amanzi::AmanziGeometry::Point *ncoord) const
  {}

  void face_get_coordinates (const Amanzi::AmanziMesh::Entity_ID faceid, 
			     std::vector<Amanzi::AmanziGeometry::Point> *fcoords) const
  {}

  void cell_get_coordinates (const Amanzi::AmanziMesh::Entity_ID cellid, 
			     std::vector<Amanzi::AmanziGeometry::Point> *ccoords) const
  {}

  void node_set_coordinates(const Amanzi::AmanziMesh::Entity_ID nodeid, 
                                      const double *coords) 
  {}

  void node_set_coordinates(const Amanzi::AmanziMesh::Entity_ID nodeid,
                            const Amanzi::AmanziGeometry::Point coords)
  {}

  const Epetra_Map& cell_epetra_map (const bool include_ghost) const
  { return *bogus_map_; }
    
  const Epetra_Map& face_epetra_map (const bool include_ghost) const
  { return *bogus_map_; }
    
  const Epetra_Map& node_epetra_map (const bool include_ghost) const
  { return *bogus_map_; }


  unsigned int get_set_size (const std::string setname, 
                             const Amanzi::AmanziMesh::Entity_kind kind,
                             const Amanzi::AmanziMesh::Parallel_type ptype) const
  { return 0; }

  unsigned int get_set_size (const char *setname, 
                             const Amanzi::AmanziMesh::Entity_kind kind,
                             const Amanzi::AmanziMesh::Parallel_type ptype) const
  { return 0; }

  unsigned int get_set_size (const Amanzi::AmanziMesh::Set_ID setid, 
                             const Amanzi::AmanziMesh::Entity_kind kind,
                             const Amanzi::AmanziMesh::Parallel_type ptype) const
  { return 0; }

  void get_set_entities (const std::string setname, 
                         const Amanzi::AmanziMesh::Entity_kind kind, 
                         const Amanzi::AmanziMesh::Parallel_type ptype, 
                         Amanzi::AmanziMesh::Entity_ID_List *entids) const
  {}

  void get_set_entities (const char *setid, 
                         const Amanzi::AmanziMesh::Entity_kind kind, 
                         const Amanzi::AmanziMesh::Parallel_type ptype, 
                         Amanzi::AmanziMesh::Entity_ID_List *entids) const
  {}
  void get_set_entities (const Amanzi::AmanziMesh::Set_ID setid, 
                         const Amanzi::AmanziMesh::Entity_kind kind, 
                         const Amanzi::AmanziMesh::Parallel_type ptype, 
                         Amanzi::AmanziMesh::Entity_ID_List *entids) const
  {}

 private:

  Epetra_Map *bogus_map_;

};


// Here, and in the Mesh::Framework unit test are hopefully the only
// places where there has to be ifdef's for mesh frameworks.

#ifdef HAVE_MOAB_MESH
#define MOAB_FLAG true
#define USE_MPI
#include "Mesh_MOAB.hh"
#undef USE_MPI
#else
#define MOAB_FLAG false
typedef bogus_maps Mesh_MOAB;
#endif

#ifdef HAVE_STK_MESH
#define STK_FLAG true
#include "Mesh_STK.hh"
#else
#define STK_FLAG false
typedef bogus_maps Mesh_STK;
#endif

#ifdef HAVE_MSTK_MESH
#define MSTK_FLAG true
#include "Mesh_MSTK.hh"
#else
#define MSTK_FLAG false
typedef bogus_maps Mesh_MSTK;
#endif


#include "Mesh_simple.hh"


namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------
// Mesh::FrameworkTraits
// 
// The idea here is to make as many decisions as possible at compile
// time.  This hopefully will reduce the code necessary to make
// appropriate framework choices at runtime.
//
// This is a pretty straightforward use of the boost::mpl library.
// Refer to Abrahams and Gurtovoy (2005). C++ Template
// Metaprogramming, Addison-Wesley.
//
// There are several things that need to be figured out.  
// 
//   1. is the framework available (compiled into the code) 
// 
//   2. can the framework read a file of a certain format,
//   considering whether the environment is parallel or not
//   
//   3. can the framework generate a mesh, in parallel or not
//
//   4. determine the appropriate Mesh_maps_* constructor to use to
//   make an instance
//
// The template argument M is expected to be one of the enum
// Amanzi::AmanziMesh::Framework.
// -------------------------------------------------------------
template < int M = 0 > 
struct FrameworkTraits {
  
  // a type that's used to see if a the specified mesh framework (M)
  // is available
  typedef mpl::bool_<
    M == Simple  || 
    ( M == MOAB && MOAB_FLAG ) ||
    ( M == STKMESH && STK_FLAG ) ||
    ( M == MSTK && MSTK_FLAG ) 
    > available;
  
  // this defines a type, there constructor of which is used to
  // instantiate a mesh when it's generated
  typedef mpl::eval_if<
    mpl::bool_<M == Simple>
    , mpl::identity<Mesh_simple>
    , mpl::eval_if<
        mpl::bool_<M == STKMESH>
        , mpl::identity<Mesh_STK>
        , mpl::eval_if<
            mpl::bool_<M == MSTK>
            , mpl::identity<Mesh_MSTK>
            , mpl::identity<bogus_maps>
            >
        >
    > generate_maps;

  // this defines a type, there constructor of which is used to
  // instantiate a mesh when it's read from a file or file set
  typedef mpl::eval_if<
    mpl::bool_<M == MOAB>
    , mpl::identity<Mesh_MOAB>
    , mpl::eval_if<
        mpl::bool_<M == STKMESH>
        , mpl::identity<Mesh_STK>
        , mpl::eval_if<
            mpl::bool_<M == MSTK>
            , mpl::identity<Mesh_MSTK>
            , mpl::identity<bogus_maps>
            >
        >
    > read_maps;
    
    
  // -------------------------------------------------------------
  // FrameworkTraits<M>::canread
  // -------------------------------------------------------------
  /// A type to indicate whether this framework can mesh of a specific format

  template < int FMT = 0 > 
  struct canread {

    struct parallel :
        mpl::eval_if<
      mpl::bool_< M == MOAB >
      , mpl::bool_< FMT == MOABHDF5 >
      , mpl::eval_if<
          mpl::bool_< M == STKMESH >
          , mpl::bool_< FMT == Nemesis >
          , mpl::eval_if<
              mpl::bool_< M == MSTK >
              , mpl::bool_< FMT == ExodusII >
              , mpl::false_
              >
          >
      >::type {};

    struct serial :
        mpl::eval_if<
      mpl::bool_< M == MOAB >
      , mpl::bool_< FMT == ExodusII || FMT == MOABHDF5 >
      , mpl::eval_if<
          mpl::bool_< M == STKMESH >
          , mpl::bool_< FMT == ExodusII >
          , mpl::eval_if<
              mpl::bool_< M == MSTK >
              , mpl::bool_< FMT == ExodusII >
              , mpl::false_
              >
          >
      >::type {};
  };

  /// Construct a mesh from a Exodus II file or file set
  static Teuchos::RCP<Mesh>
  read(const Epetra_MpiComm *comm, const std::string& fname,
     const AmanziGeometry::GeometricModelPtr& gm)
  {
    Teuchos::RCP<Mesh> 
    result(new typename read_maps::type(fname.c_str(), comm, gm));
    return result;
  }

  /// A type to indicate whether this framework can generate meshes
  template < unsigned int DIM = 0 >
  struct cangenerate {

    struct parallel : 
      mpl::eval_if<
      mpl::bool_< M == STKMESH >
      , mpl::bool_< DIM == 3 >
      , mpl::eval_if< 
          mpl::bool_< M == MSTK >
          , mpl::bool_< DIM == 2 || DIM == 3 >
          , mpl::false_
          >
      >::type {};
    
    struct serial :
      mpl::eval_if<
      mpl::bool_< M == Simple >
      , mpl::bool_< DIM == 3 >
      , mpl::eval_if< 
          mpl::bool_< M == STKMESH >
          , mpl::bool_< DIM == 3 >
          , mpl::eval_if<
              mpl::bool_< M == MSTK >
              , mpl::bool_< DIM == 2 || DIM == 3 >
              , mpl::false_
              >
          >
      >::type {};
  };

/// Generate a hex mesh from explicit arguments
  static Teuchos::RCP<Mesh>
  generate(const double& x0, const double& y0, const double& z0,
           const double& x1, const double& y1, const double& z1,
           const unsigned int& nx, const unsigned int& ny, const unsigned int& nz, 
         const Epetra_MpiComm *comm, 
         const AmanziGeometry::GeometricModelPtr& gm)
  {
    Teuchos::RCP<Mesh> 
    result(new typename generate_maps::type(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm, gm));
    return result;
  }

/// Generate a quad mesh from explicit arguments
  static Teuchos::RCP<Mesh>
  generate(const double& x0, const double& y0,
           const double& x1, const double& y1,
           const unsigned int& nx, const unsigned int& ny,
         const Epetra_MpiComm *comm, 
         const AmanziGeometry::GeometricModelPtr& gm)
  {
    Teuchos::RCP<Mesh> 
    result(new typename generate_maps::type(x0, y0, x1, y1, nx, ny, comm, gm));
    return result;
  }

/// Generate a hex mesh from arguments sent in through a parameter list
  static Teuchos::RCP<Mesh>
  generate(Teuchos::ParameterList &parameter_list, const Epetra_MpiComm *comm,
         const AmanziGeometry::GeometricModelPtr& gm)
  {
    GenerationSpec gspec(parameter_list);
    Teuchos::RCP<Mesh> 
    result(new typename generate_maps::type(gspec, comm, gm));
    return result;
  }
};

// -------------------------------------------------------------
// framework_available
// -------------------------------------------------------------
bool
framework_available(const Framework& f)
{
  bool result;
  switch (f) {
    case Simple:
      result = FrameworkTraits<Simple>::available::value;
      break;
    case STKMESH:
      result = FrameworkTraits<STKMESH>::available::value;
      break;
    case MOAB:
      result = FrameworkTraits<MOAB>::available::value;
      break;
    case MSTK:
      result = FrameworkTraits<MSTK>::available::value;
      break;
    default:
      {
        std::string msg = 
            boost::str(boost::format("unknown mesh framework: %d") % static_cast<int>(f));
        Exceptions::amanzi_throw(Errors::Message(msg.c_str()));
      }
  }
  return result;
}

// -------------------------------------------------------------
// parallel_test
// -------------------------------------------------------------
template < class thetest >
static bool 
parallel_test(const bool& isp)
{
  bool result;
  if (isp) {
    result = thetest::parallel::value;
  } else {
    result = thetest::serial::value;
  }
  return result;
}

// -------------------------------------------------------------
// framework_reads
// -------------------------------------------------------------

template < int F >
static bool 
framework_reads(const Format& fmt, const bool& parallel)
{
  typedef FrameworkTraits<F> traits;
  bool result = false;
  switch (fmt) {
    case ExodusII:
      result = parallel_test< typename traits::template canread<ExodusII> >(parallel);
      break;
    case MOABHDF5:
      result = parallel_test< typename traits::template canread<MOABHDF5> >(parallel);
      break;
    case Nemesis:
      result = parallel_test< typename traits::template canread<Nemesis> >(parallel);
      break;
    default:
      result = false;
  }
  return result;
}

bool
framework_reads(const Framework& f, const Format& fmt, const bool& parallel)
{
  bool result;
  switch (f) {
    case Simple:
      result = framework_reads<Simple>(fmt, parallel);
      break;
    case STKMESH:
      result = framework_reads<STKMESH>(fmt, parallel);
      break;
    case MOAB:
      result = framework_reads<MOAB>(fmt, parallel);
      break;
    case MSTK:
      result = framework_reads<MSTK>(fmt, parallel);
      break;
    default:
      {
        std::string msg = 
            boost::str(boost::format("unknown mesh framework: %d") % static_cast<int>(f));
        Exceptions::amanzi_throw(Errors::Message(msg.c_str()));
      }
  }
  return result;
}

// -------------------------------------------------------------
// framework_read
// -------------------------------------------------------------
Teuchos::RCP<Mesh> 
framework_read(const Epetra_MpiComm *comm, const Framework& f, const std::string& fname, 
               const AmanziGeometry::GeometricModelPtr& gm)
{
  Teuchos::RCP<Mesh> result;
  switch (f) {
    case Simple:
      result = FrameworkTraits<Simple>::read(comm, fname, gm);
      break;
    case STKMESH:
      result = FrameworkTraits<STKMESH>::read(comm, fname, gm);
      break;
    case MOAB:
      result = FrameworkTraits<MOAB>::read(comm, fname, gm);
      break;
    case MSTK:
      result = FrameworkTraits<MSTK>::read(comm, fname, gm);
      break;
    default:
      {
        std::string msg = 
            boost::str(boost::format("unknown mesh framework: %d") % static_cast<int>(f));
        Exceptions::amanzi_throw(Errors::Message(msg.c_str()));
      }
  }
  return result;
}

// -------------------------------------------------------------
// framework_generates
// -------------------------------------------------------------

template < int F >
static bool 
framework_generates(const bool& parallel, const unsigned int& dimension)
{
  typedef FrameworkTraits<F> traits;
  bool result = false;
  switch (dimension) {
  case 2:
    result = parallel_test< typename traits::template cangenerate<2> >(parallel);
    break;
  case 3:
    result = parallel_test< typename traits::template cangenerate<3> >(parallel);
    break;
  default:
    result = false;
  }
  return result;
}


bool
framework_generates(const Framework& f, const bool& parallel, const unsigned int& dimension)
{
  bool result;  

  switch (f) {
    case Simple:
      result = framework_generates<Simple>(parallel, dimension);
      break;
    case STKMESH:  
      result = framework_generates<STKMESH>(parallel, dimension);
      break;
    case MOAB:
      result = framework_generates<MOAB>(parallel, dimension);
      break;
    case MSTK:
      result = framework_generates<MSTK>(parallel, dimension);
      break;
    default:
      {
        std::string msg = 
            boost::str(boost::format("Cannot generate dimension %d meshes") % static_cast<int>(dimension));
        Exceptions::amanzi_throw(Errors::Message(msg.c_str()));
      }
  }
  return result;
}

// -------------------------------------------------------------
// framework_generate
// -------------------------------------------------------------
Teuchos::RCP<Mesh> 
framework_generate(const Epetra_MpiComm *comm, const Framework& f, 
                   const double& x0, const double& y0, const double& z0,
                   const double& x1, const double& y1, const double& z1,
                   const unsigned int& nx, const unsigned int& ny, const unsigned int& nz,
                   const AmanziGeometry::GeometricModelPtr& gm)
{
  Teuchos::RCP<Mesh> result;
  switch (f) {
    case Simple:
      result = FrameworkTraits<Simple>::generate(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm, gm);
      break;
    case STKMESH:
      result = FrameworkTraits<STKMESH>::generate(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm, gm);
      break;
    case MOAB:
      result = FrameworkTraits<MOAB>::generate(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm, gm);
      break;
    case MSTK:
      result = FrameworkTraits<MSTK>::generate(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm, gm);
      break;
    default:
      {
        std::string msg = 
            boost::str(boost::format("unknown mesh framework: %d") % static_cast<int>(f));
        Exceptions::amanzi_throw(Errors::Message(msg.c_str()));
      }
  }
  return result;
}

// -------------------------------------------------------------
// framework_generate
// -------------------------------------------------------------
Teuchos::RCP<Mesh> 
framework_generate(const Epetra_MpiComm *comm, const Framework& f, 
                   const double& x0, const double& y0,
                   const double& x1, const double& y1,
                   const unsigned int& nx, const unsigned int& ny, 
                   const AmanziGeometry::GeometricModelPtr& gm)
{
  Teuchos::RCP<Mesh> result;
  switch (f) {
    case Simple:
      result = FrameworkTraits<Simple>::generate(x0, y0, x1, y1, nx, ny, comm, gm);
      break;
    case STKMESH:
      result = FrameworkTraits<STKMESH>::generate(x0, y0, x1, y1, nx, ny, comm, gm);
      break;
    case MOAB:
      result = FrameworkTraits<MOAB>::generate(x0, y0, x1, y1, nx, ny, comm, gm);
      break;
    case MSTK:
      result = FrameworkTraits<MSTK>::generate(x0, y0, x1, y1, nx, ny, comm, gm);
      break;
    default:
      {
        std::string msg = 
            boost::str(boost::format("unknown mesh framework: %d") % static_cast<int>(f));
        Exceptions::amanzi_throw(Errors::Message(msg.c_str()));
      }
  }
  return result;
}

Teuchos::RCP<Mesh> 
framework_generate(const Epetra_MpiComm *comm, const Framework& f, 
                   Teuchos::ParameterList &parameter_list,
                   const AmanziGeometry::GeometricModelPtr& gm)
{
  Teuchos::RCP<Mesh> result;
  switch (f) {
    case Simple:
      result = FrameworkTraits<Simple>::generate(parameter_list, comm, gm);
      break;
    case STKMESH:
      result = FrameworkTraits<STKMESH>::generate(parameter_list, comm, gm);
      break;
    case MOAB:
      result = FrameworkTraits<MOAB>::generate(parameter_list, comm, gm);
      break;
    case MSTK:
      result = FrameworkTraits<MSTK>::generate(parameter_list, comm, gm);
      break;
    default:
      {
        std::string msg = 
            boost::str(boost::format("unknown mesh framework: %d") % static_cast<int>(f));
        Exceptions::amanzi_throw(Errors::Message(msg.c_str()));
      }
  }
  return result;
}


} // namespace AmanziMesh
} // namespace Amanzi
