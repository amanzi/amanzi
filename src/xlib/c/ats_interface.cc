/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * ATS C Language interface implementation.
 *----------------------------------------------------------------------------*/

#define _ats_source

#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include <ats_interface.h>
#include <ats_state.hh>

ats_state_t & _state = ats_state_t::instance();

int32_t ats_init(MPI_Comm comm) {

#ifdef HAVE_MPI
	Epetra_MpiComm * ecomm = new Epetra_MpiComm(comm);
#else
	Epetra_SerialComm * ecomm = new Epetra_SerialComm();
#endif

	return _state.mpi_rehash(comm);
} // ats_init

int32_t ats_finalize() {
} // ats_finalize
