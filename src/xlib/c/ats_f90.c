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

int32_t ats_init_f90(MPI_Fint * comm) {
	MPI_Comm _comm = MPI_Comm_f2c(*comm);
	return ats_init(_comm);
} // ats_init

int32_t ats_finalize_f90() {
	return ats_finalize();
} // ats_finalize
