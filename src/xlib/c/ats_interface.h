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
int32_t ats_init(MPI_Comm mpi_comm, int32_t * type_ids,
	int32_t num_type_ids, int32_t num_types);

int32_t ats_finalize();

int32_t ats_set_init_clm_data(double * T, double * Sl, double * Si);
int32_t ats_set_clm_data(double * e_flux, double * w_flux);
int32_t ats_get_clm_data(double * T, double * Sl, double * Si);

int32_t ats_advance(double dt, int32_t force_viz);

#if defined(__cplusplus)
} // extern
#endif

#endif // ats_interface_h
