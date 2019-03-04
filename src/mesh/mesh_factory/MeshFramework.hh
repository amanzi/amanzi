/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Mesh Framework Manager

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: William A. Perkins, Ethan Coon
*/

/*!

Implementation note: this is a stripped down traits mechanism.

Hopefully it is simple enough for us to use.

*/
  
#ifndef AMANZI_MESH_FRAMEWORK_HH_
#define AMANZI_MESH_FRAMEWORK_HH_

#include <string>
#include <vector>
#include "boost/type_traits.hpp"

namespace Amanzi {
namespace AmanziMesh {

// A type to identify available mesh frameworks
struct enum Framework { SIMPLE=0, MSTK=1, MOAB=2, STK=3 };


// A list of frameworks, maybe available, maybe not.
typedef std::vector<Framework> Preference;


// This namespace is impelemntation detail, and really only should be used by
// the mesh factory.
namespace AmanziFrameworks {

//
// Available traits -- all default to false_type, so that you must explicitly
// ENABLE traits only.
// -----------------------------------------------------------------------------

// A generic struct for mapping from Framework to the class
template <Framework f>
struct Mesh_type {};

// Statically informs whether a framework is available.
template <Framework f>
struct is_available() : public false_type{};

// Statically informs whether a framework can read from file
template <Framework f, int num_proc>
struct reads() : public false_type{};

// Statically informs whether a framework can generate
template <Framework f, int num_proc, int dim>
struct generates() : public false_type{};

// Statically informs whether a framework can extract
template <Framework f, int num_proc>
struct extracts() : public false_type{};


//
// SIMPLE traits
//
// simple only generates 3D meshes.  It is largely only for testing.
// -----------------------------------------------------------------------------
#include "Mesh_simple.hh"

template <>
struct Mesh_type<Framework::SIMPLE> : public Mesh_Simple {};

template <> struct is_available<Framework::SIMPLE>() : public true_type{};

// simple generates in serial for 3D meshes
template <> struct generates<Framework::SIMPLE, 1, 3>() : public true_type{};


//
// MSTK traits
//
// MSTK does everything
// -----------------------------------------------------------------------------
#ifdef HAVE_MSTK_MESH
#include "Mesh_MSTK.hh"

template <>
struct Mesh_type<Framework::MSTK> : public Mesh_MSTK {};

template <> struct is_available<Framework::MSTK>() : public true_type{};

// mstk reads in any parallel type
template <int num_proc> struct reads<Framework::MSTK, num_proc>() : public true_type{};

// mstk generates in serial and parallel and in 2D and 3D
template <Parallel p> struct generates<Framework::MSTK, p, 2>() : public true_type{};
template <Parallel p> struct generates<Framework::MSTK, p, 3>() : public true_type{};

// mstk extracts in any parallel type
template <Parallel p> struct extracts<Framework::MSTK, p>() : public true_type{};
#endif 


//
// STK traits
//
// STK IS OFFICIALLY DISABLED!
// -----------------------------------------------------------------------------
#ifdef HAVE_STK_MESH
#endif 


//
// MOAB traits
//
// MOAB can be used to read, in serial or parallel.
// -----------------------------------------------------------------------------
#ifdef HAVE_MOAB_MESH
#include "Mesh_MOAB.hh"

template <>
struct Mesh_type<Framework::MOAB> : public Mesh_MOAB {};

template <> struct is_available<Framework::MOAB>() : public true_type{};
// moab reads in any parallel type
template <Parallel p> struct reads<Framework::MOAB, p>() : public true_type{};
#endif 


//
// returns null if it didn't work
//
template<Mesh_type, bool b, Args... args>
Teuchos::RCP<Mesh_type> construct(const boost::integral_constant<bool, b>&,
        Args... args) {
  return Teuchos::null;
}

//
// returns a mesh if it did
template<Mesh_type, Args... args>
Teuchos::RCP<Mesh_type> construct(const boost::true_type&, Args... args) {
  return Teuchos::rcp(new Mesh_type(args...));
}

  









} // close namespace AmanziFrameworks
} // close namespace AmanziMesh
} // close namespace Amanzi

#endif
