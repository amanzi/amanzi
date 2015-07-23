/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * ATS C Language interface implementation.
 *----------------------------------------------------------------------------*/

#define _ats_source

#include <iostream>

#include <ats_interface.h>
#include <ats_state.hh>

ats_state_t & _state = ats_state_t::instance();

int32_t ats_init(MPI_Comm mpi_comm, int32_t * type_ids,
	int32_t num_type_ids, int32_t num_types) {
	return _state.clm_driver().Initialize(mpi_comm, type_ids,
		num_type_ids, num_types);
} // ats_init

int32_t ats_finalize() {
	return _state.clm_driver().Finalize();
} // ats_finalize

int32_t ats_set_init_clm_data(double * T, double * Sl, double * Si) {
	return _state.clm_driver().SetInitCLMData(T, Sl, Si);
} // ats_set_init_clm_data

int32_t ats_set_clm_data(double * e_flux, double * w_flux) {
	return _state.clm_driver().SetCLMData(e_flux, w_flux);
} // ats_set_clm_data

int32_t ats_get_clm_data(double * T, double * Sl, double * Si) {
	return _state.clm_driver().GetCLMData(T, Sl, Si);
} // ats_get_clm_data

int32_t ats_advance(double dt, int32_t force_viz) {
	return _state.clm_driver().Advance(dt, force_viz == 1 ? true : false);
} // ats_advance
