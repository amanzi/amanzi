/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Julien Loiseau (jloiseau@lanl.gov)
           Rao Garimella (rao@lanl.gov)
*/

//! Caches mesh information for fast repeated access.
// See MeshCache_decl.hh for documentation!

#pragma once

#include "MeshCacheHost_decl.hh"
#include "MeshCacheDevice_decl.hh"

namespace Amanzi::AmanziMesh 
{
  
// -----------------------------------------------------------------------------
// Memory space transfer
// -----------------------------------------------------------------------------
inline const MeshCacheDevice
onMemDevice(const MeshCacheHost& mc_in)
{
  MeshCacheHost* mc_in_nc = const_cast<MeshCacheHost*>(&mc_in);
  MeshCacheDevice mc_out(*mc_in_nc);
  return mc_out;
}

inline Teuchos::RCP<const MeshCacheDevice>
onMemDevice(const Teuchos::RCP<const MeshCacheHost>& mc_in)
{
  MeshCacheHost* mc_in_nc = const_cast<MeshCacheHost*>(mc_in.ptr().get());
  Teuchos::RCP<MeshCacheDevice> mc_on_out = Teuchos::rcp(new MeshCacheDevice(*mc_in_nc));
  return (Teuchos::RCP<const MeshCacheDevice>)mc_on_out;
}

inline Teuchos::RCP<MeshCacheDevice>
onMemDevice(const Teuchos::RCP<MeshCacheHost>& mc_in)
{
  return Teuchos::rcp(new MeshCacheDevice(*mc_in));
}

inline const MeshCacheHost
onMemHost(const MeshCacheDevice& mc_in)
{
  MeshCacheDevice* mc_in_nc = const_cast<MeshCacheDevice*>(&mc_in);
  MeshCacheHost mc_out(*mc_in_nc);
  return mc_out;
} 
  
inline Teuchos::RCP<const MeshCacheHost>
onMemHost(const Teuchos::RCP<const MeshCacheDevice>& mc_in)
{
  MeshCacheDevice* mc_in_nc = const_cast<MeshCacheDevice*>(mc_in.ptr().get());
  Teuchos::RCP<MeshCacheHost> mc_on_out = Teuchos::rcp(new MeshCacheHost(*mc_in_nc));
  return (Teuchos::RCP<const MeshCacheHost>)mc_on_out;
}

inline Teuchos::RCP<MeshCacheHost>
onMemHost(const Teuchos::RCP<MeshCacheDevice>& mc_in)
{
  return Teuchos::rcp(new MeshCacheHost(*mc_in));
}
  
}
