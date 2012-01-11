/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   RegionFactory.hh
 * @author Rao V Garimella
 * @date 
 * 
 * @brief  Factory to create specific types of regions
 * 
 * 
 */

#ifndef _RegionFactory_hh_
#define _RegionFactory_hh_

#include "Epetra_MpiComm.h"
#include <Teuchos_ParameterList.hpp>


#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {


// Factory for creating specific types of regions based on the shape
// specification. We cannot use a constructor because the we have to 
// create derived region classes based on the shape parameter

RegionPtr RegionFactory(const std::string reg_name, const unsigned int reg_id,
                        const Teuchos::ParameterList& reg_spec,
                        const Epetra_MpiComm *comm);

} // namespace AmanziGeometry
} // namespace Amanzi

#endif

