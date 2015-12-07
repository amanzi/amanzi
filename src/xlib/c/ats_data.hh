/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * Data structure prototypes for ATS C Language interface.
 *----------------------------------------------------------------------------*/

#ifndef ats_data_hh
#define ats_data_hh

#if !defined(_ats_source) && !defined(_include_ats_h)
#error "Error: do not include this file directly, use #include <ats.h>"
#endif

#include <mpi.h>
#include <exceptions.hh>

#include <ats_defines.h>

struct mpi_data_t {
	MPI_Comm comm;
	int32_t cached_rank;
	int32_t cached_size;

	int32_t rehash(MPI_Comm & _comm) {
		comm = _comm;
		std::cerr << "Here" << std::endl;
		Exceptions::amanzi_throw("Help me");

		if(MPI_Comm_rank(comm, &cached_rank) != MPI_SUCCESS) {
			return ATS_MPI_ERROR;
		} // if

		return 0;
	} // rehash

}; // mpi_data_t

#endif // ats_data_hh
