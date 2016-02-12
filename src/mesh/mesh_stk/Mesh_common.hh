// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   Mesh_common.hh
 * @author William A. Perkins
 * @date Wed Dec 29 11:11:11 2010
 * 
 * @brief  Declarations of some stk::mesh utility routines
 * 
 * 
 */

// -------------------------------------------------------------
// Created December 29, 2010 by William A. Perkins
// Last Change: Wed Dec 29 11:11:11 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _Mesh_common_hh_
#define _Mesh_common_hh_

#include <stk_mesh/base/BulkData.hpp>

namespace stk {
  namespace mesh {

    /// Make sure faces are owned by the same process as a connecting cell
    void check_face_ownership(BulkData& bulk_data);

    /// Make sure nodes are owned by the same process as a connecting cell
    void check_node_ownership(BulkData& bulk_data);
    
  } // close namespace mesh
} // close namespace stk
#endif
