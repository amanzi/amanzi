// -------------------------------------------------------------
/**
 * @file   FrameworkTraits.cc
 * @author William A. Perkins
 * @date Fri Mar 18 15:22:43 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 14, 2011 by William A. Perkins
// Last Change: Fri Mar 18 15:22:43 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <boost/format.hpp>

#include <iostream>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>

namespace mpl = boost::mpl;

#include <Teuchos_RCP.hpp>

#include "FrameworkTraits.hh"
#include "MeshFramework.hh"
#include "MeshFileType.hh"
#include "MeshException.hh"
#include "Mesh_maps_base.hh"

// -------------------------------------------------------------
//  class bogus_maps
// -------------------------------------------------------------
/**
 * This class is just used to throw an exception when a mesh framework
 * is used, but not available, or is available, but misused.
 * 
 */

class bogus_maps : public Mesh_maps_base {
public:

  /// Default constructor.
  bogus_maps(const char *filename, MPI_Comm comm) 
    : Mesh_maps_base() 
  {
    amanzi_throw(Mesh::Message("reading not supported"));
  }

  bogus_maps(double x0, double y0, double z0,
             double x1, double y1, double z1,
             int nx, int ny, int nz, 
             Epetra_MpiComm *communicator)
    : Mesh_maps_base() 
  {
    amanzi_throw(Mesh::Message("generation not supported"));
  }

};

// Here, and in the Mesh::Framework unit test are hopefully the only
// places where there has to be ifdef's for mesh frameworks.

#ifdef HAVE_MOAB_MESH_NOT
#define MOAB_FLAG true
#include "Mesh_maps_moab.hh"
#else
#define MOAB_FLAG false
typedef bogus_maps Mesh_maps_moab;
#endif

#ifdef HAVE_STK_MESH
#define STK_FLAG true
#include "Mesh_maps_stk.hh"
typedef STK_mesh::Mesh_maps_stk Mesh_maps_stk;
#else
#define STK_FLAG false
typedef bogus_maps Mesh_maps_stk;
#endif

#ifdef HAVE_MSTK_MESH
#define MSTK_FLAG true
#include "Mesh_maps_mstk.hh"
#else
#define MSTK_FLAG false
typedef bogus_maps Mesh_maps_mstk;
#endif

#include "Mesh_maps_simple.hh"

namespace Mesh {



// template <int N>
// struct foo
// {
//   foo() { std::cout << "foo<" << N << ">" << std::endl; }
// };

// template <int Index>
// struct type_at
// {
//   // this is the type container..
//   typedef typename mpl::vector<foo<0>, foo<1>, foo<2>, foo<3> > seq;
//   // this allows us to get the type at the given index
//   typedef typename mpl::at<seq, mpl::int_<Index> >::type type;
// };
  

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
  //   -------------------------------------------------------------
  template < int M = 0 > 
  struct FrameworkTraits {
    
    // a type that's used to see if a the specified mesh framework (M)
    // is available
    typedef mpl::bool_<
      M == Simple  || 
      ( M == MOAB && MOAB_FLAG ) ||
      ( M == STK && STK_FLAG ) ||
      ( M == MSTK && MSTK_FLAG ) 
      > available;
  
    // this defines a type, there constructor of which is used to
    // instantiate a mesh when it's generated
    typedef mpl::eval_if<
      mpl::bool_<M == Simple>
      , mpl::identity<Mesh_maps_simple>
      , mpl::eval_if<
          mpl::bool_<M == STK>
          , mpl::identity<Mesh_maps_stk>
          , mpl::identity<bogus_maps>
          >
      > generate_maps;
    
    // // this defines a type, there constructor of which is used to
    // // instantiate a mesh when it's read from a file or file set
    // typedef mpl::eval_if<
    //   mpl::bool_<M == MOAB>
    //   , mpl::identity<Mesh_maps_moab>
    //   , mpl::if_<
    //       mpl::bool_<M == STK>
    //       , mpl::identity<Mesh_maps_stk>
    //       , mpl::if_<
    //           mpl::bool_<M == MSTK>
    //           , mpl::identity<Mesh_maps_mstk>
    //           , mpl::identity<bogus_maps>
    //           >
    //       >
    //   > read_maps;
    
    
    // -------------------------------------------------------------
    // FrameworkTraits<M>::canread
    // -------------------------------------------------------------
    /// A type to indicate whether this framework can mesh of a specific format
    template < int FMT = 0 > 
    struct canread {

      struct parallel :
        mpl::eval_if<
        mpl::bool_< M == MOAB >
        , mpl::bool_< FMT == ExodusII || FMT == MOABHDF5 >
        , mpl::eval_if<
            mpl::bool_< M == STK >
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
            mpl::bool_< M == STK >
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
    static Teuchos::RCP<Mesh_maps_base>
    read(const Epetra_MpiComm& comm, const std::string& fname)
    {
      Teuchos::RCP<Mesh_maps_base> result;
      return result;
    }

    /// A type to indicate whether this framework can generate meshes
    struct cangenerate {
      struct parallel : 
        mpl::bool_< M == STK >::type
      {};
      struct serial :
        mpl::bool_<M == Simple || M == STK >::type
      {};
    };

    /// Generate a hex mesh 
    static Teuchos::RCP<Mesh_maps_base>
    generate(const double& x0, const double& y0, const double& z0,
         const double& x1, const double& y1, const double& z1,
         const unsigned int& nx, const unsigned int& ny, const unsigned int& nz, 
         Epetra_MpiComm& comm)
    {
      Teuchos::RCP<Mesh_maps_base> 
        result(new typename generate_maps::type(x0, y0, z0, x1, y1, z1, nx, ny, nz, &comm));
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
    case STK:
      result = FrameworkTraits<STK>::available::value;
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
        amanzi_throw(Message(msg.c_str()));
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
    case STK:
      result = framework_reads<STK>(fmt, parallel);
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
        amanzi_throw(Message(msg.c_str()));
      }
    }
    return result;
  }

  // -------------------------------------------------------------
  // framework_read
  // -------------------------------------------------------------
  Teuchos::RCP<Mesh_maps_base> 
  framework_read(Epetra_MpiComm& comm, const Framework& f, const std::string& fname)
  {
    Teuchos::RCP<Mesh_maps_base> result;
    switch (f) {
    case Simple:
      result = FrameworkTraits<Simple>::read(comm, fname);
      break;
    case STK:
      result = FrameworkTraits<STK>::read(comm, fname);
      break;
    case MOAB:
      result = FrameworkTraits<MOAB>::read(comm, fname);
      break;
    case MSTK:
      result = FrameworkTraits<MSTK>::read(comm, fname);
      break;
    default:
      {
        std::string msg = 
          boost::str(boost::format("unknown mesh framework: %d") % static_cast<int>(f));
        amanzi_throw(Message(msg.c_str()));
      }
    }
    return result;
  }

  // -------------------------------------------------------------
  // framework_generates
  // -------------------------------------------------------------
  bool
  framework_generates(const Framework& f, const bool& parallel)
  {
    bool result;
    switch (f) {
    case Simple:
      result = parallel_test< FrameworkTraits<Simple>::cangenerate >(parallel);
      break;
    case STK:
      result = parallel_test< FrameworkTraits<STK>::cangenerate >(parallel);
      break;
    case MOAB:
      result = parallel_test< FrameworkTraits<MOAB>::cangenerate >(parallel);
      break;
    case MSTK:
      result = parallel_test< FrameworkTraits<MSTK>::cangenerate >(parallel);
      break;
    default:
      {
        std::string msg = 
          boost::str(boost::format("unknown mesh framework: %d") % static_cast<int>(f));
        amanzi_throw(Message(msg.c_str()));
      }
    }
    return result;
  }

  // -------------------------------------------------------------
  // framework_generates
  // -------------------------------------------------------------
  Teuchos::RCP<Mesh_maps_base> 
  framework_generate(Epetra_MpiComm& comm, const Framework& f, 
                     const double& x0, const double& y0, const double& z0,
                     const double& x1, const double& y1, const double& z1,
                     const unsigned int& nx, const unsigned int& ny, const unsigned int& nz)
  {
    Teuchos::RCP<Mesh_maps_base> result;
    switch (f) {
    case Simple:
      result = FrameworkTraits<Simple>::generate(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm);
      break;
    case STK:
      result = FrameworkTraits<STK>::generate(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm);
      break;
    case MOAB:
      result = FrameworkTraits<MOAB>::generate(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm);
      break;
    case MSTK:
      result = FrameworkTraits<MSTK>::generate(x0, y0, z0, x1, y1, z1, nx, ny, nz, comm);
      break;
    default:
      {
        std::string msg = 
          boost::str(boost::format("unknown mesh framework: %d") % static_cast<int>(f));
        amanzi_throw(Message(msg.c_str()));
      }
    }
    return result;
  }


} // namespace Mesh

