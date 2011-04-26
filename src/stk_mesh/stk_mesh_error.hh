// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: stk_mesh_error.hh
// -------------------------------------------------------------
/**
 * @file   stk_mesh_error.hh
 * @author William A. Perkins
 * @date Thu Apr 21 12:56:54 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// Created April 21, 2011 by William A. Perkins
// Last Change: Thu Apr 21 12:56:54 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _stk_mesh_error_hh_
#define _stk_mesh_error_hh_

#include "errors.hh"

namespace STK_mesh
{
    class STKMeshError : public Errors::Message 
    {
    public:
        explicit STKMeshError(void) : Errors::Message() {};
        explicit STKMeshError(const char* message) : Errors::Message(message) {};
        ~STKMeshError(void) throw() {};
    };
} // namespace STK_mesh
#endif
