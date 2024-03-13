/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Julien Loiseau (jloiseau@lanl.gov)
*/

/*

Developer note:

The organization of files is a bit convoluted here, due to code sharing between
MeshFrameworks and MeshCache.  The basic concepts/files are:

MeshFramework is a CPU only implementation of an unstructured mesh, that can be
slow and allows virtual functions.

MeshCache is a CPU or DEVICE implementation that must be fast, has no virtual
lookups, and supports multiple ways of getting data (see MeshCache_decl.hh).

When we refer to the Mesh API, we mean the shared API of both MeshCache and
MeshFramework.

MeshInternals are internally used functions, namespaced in an Impl namespace to
indicate that they are not used outside of this library.  They are templated on
mesh type, and are used by both MeshFramework and MeshCache objects.

MeshAlgorithms is a virtual class that implements things that we may wish to
compute and not cache.  These algorithms are not templated as they must be
virtual to allow "strange" meshes like MeshLogical and others to override them,
and are separate from the MeshFramework to allow us to destroy the Framework
and hold onto these virtual functions.  These are only used internally to the
mesh library.

Lastly MeshHelpers are functions that use the public Mesh API and may be used
externally to this library.

*/

#pragma once

#include "MeshDefs.hh"
#include "MeshInternals_decl.hh"
#include "MeshCache_decl.hh"
#include "MeshAlgorithms.hh"
#include "MeshHelpers.hh"
#include "MeshInternals_impl.hh"
#include "MeshCache_impl.hh"
