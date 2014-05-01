/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * ATS C Language interface implementation.
 *----------------------------------------------------------------------------*/

#define _ats_source

#include <ats_interface.h>
#include <mpi.h>

int32_t ats_init_f90(MPI_Fint * comm, int32_t * type_ids,
	int32_t * num_type_ids, int32_t * num_types) {
	MPI_Comm _comm = MPI_Comm_f2c(*comm);
	return ats_init(_comm, type_ids, *num_type_ids, *num_types);
} // ats_init_f90

int32_t ats_finalize_f90() {
	return ats_finalize();
} // ats_finalize_f90

int32_t ats_set_init_clm_data_f90(double * T, double * Sl, double * Si) {
	return ats_set_init_clm_data(T, Sl, Si);
} // ats_set_init_clm_data_f90

int32_t ats_set_clm_data_f90(double * e_flux, double * w_flux) {
	return ats_set_clm_data(e_flux, w_flux);
} // ats_set_clm_data_f90

int32_t ats_get_clm_data_f90(double * T, double * Sl, double * Si) {
	return ats_get_clm_data(T, Sl, Si);
} // ats_get_clm_data_f90

int32_t ats_advance_f90(double * dt, int32_t * force_viz) {
	return ats_advance(*dt, *force_viz);
} // ats_advance_f90
