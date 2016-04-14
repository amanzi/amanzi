/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * ATS state singleton to hold C++ state for C Language interface.
 *----------------------------------------------------------------------------*/

#if !defined(_ats_source) && !defined(_include_ats_h)
#error "Error: do not include this file directly, use #include <ats.h>"
#endif

#ifndef ats_state_hh
#define ats_state_hh

#include <ats_data.hh>
#include <ats_clm_driver.hh>

/*----------------------------------------------------------------------------*
 * ats_state_t
 *----------------------------------------------------------------------------*/

class ats_state_t
{
public:

	/*-------------------------------------------------------------------------*
	 * Meyer's singleton instance
	 *-------------------------------------------------------------------------*/

	static ats_state_t & instance() {
		static ats_state_t _s;
		return _s;
	} // instance

	/*-------------------------------------------------------------------------*
	 * MPI data
	 *-------------------------------------------------------------------------*/

	int32_t mpi_rehash(MPI_Comm & comm) { return mpi_data_.rehash(comm); }

	const MPI_Comm & mpi_comm() const { return mpi_data_.comm; }
	MPI_Comm & mpi_comm() { return mpi_data_.comm; }

	const int32_t & mpi_rank() const { return mpi_data_.cached_rank; }
	int32_t & mpi_rank() { return mpi_data_.cached_rank; }

	const int32_t & mpi_size() const { return mpi_data_.cached_size; }
	int32_t & mpi_size() { return mpi_data_.cached_size; }

	Amanzi::ATSCLMDriver & clm_driver() { return clm_driver_; }

private:

	/*-------------------------------------------------------------------------*
	 * These are private to keep the singleton single...
	 *-------------------------------------------------------------------------*/

	ats_state_t() {}
	ats_state_t(const ats_state_t &) {}
	~ats_state_t() {}
	ats_state_t & operator = (const ats_state_t &);

	/*-------------------------------------------------------------------------*
	 * Data members
	 *-------------------------------------------------------------------------*/

	mpi_data_t mpi_data_;	
	Amanzi::ATSCLMDriver clm_driver_;

}; // class ats_state_t

#endif // ats_state_hh
