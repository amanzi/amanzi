/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * ATS C Language interface prototype.
 *----------------------------------------------------------------------------*/

#if !defined(_ats_source) && !defined(_include_ats_h)
#error "Error: do not include this file directly, use #include <ats.h>"
#endif

#ifndef ats_interface_h
#define ats_interface_h

#include <stdint.h>
#include <mpi.h>

#include <ats_defines.h>

#if defined(__cplusplus)
extern "C" {
#endif

/*
 */
int32_t ats_init(MPI_Comm comm);

int32_t ats_finalize();

#if defined(__cplusplus)
} // extern
#endif

#endif // ats_interface_h
