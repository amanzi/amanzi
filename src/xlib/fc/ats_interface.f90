!------------------------------------------------------------------------------!
! ATS
!
! License: see $ATS_DIR/COPYRIGHT
!
! ATS Fortran interface.
!------------------------------------------------------------------------------!

module ats_interface

   use, intrinsic :: ISO_C_BINDING
   use :: ats_bindings

   contains

   !---------------------------------------------------------------------------!
   ! ats_init
   !---------------------------------------------------------------------------!

   subroutine ats_init(mpi_comm, type_ids, num_type_ids, num_types, ierr)
      use :: ats_data
      implicit none
      integer(int32_t) :: mpi_comm
      type(c_ptr) :: type_ids
      integer(int32_t) :: num_type_ids
      integer(int32_t) :: num_types
      integer(int32_t) :: ierr

      ierr = ats_init_f90(mpi_comm, type_ids, num_type_ids, num_types)
   end subroutine ats_init

   !---------------------------------------------------------------------------!
   ! ats_finalize
   !---------------------------------------------------------------------------!

   subroutine ats_finalize(ierr)
      use :: ats_data
      implicit none
      integer(int32_t) :: ierr

      ierr = ats_finalize_f90()
   end subroutine ats_finalize

   !---------------------------------------------------------------------------!
   ! ats_init_clm_data
   !---------------------------------------------------------------------------!

   subroutine ats_set_init_clm_data(T, Sl, Si, ierr)
      use :: ats_data
      implicit none
      type(c_ptr), value :: T
      type(c_ptr), value :: Sl
      type(c_ptr), value :: Si
      integer(int32_t) :: ierr

      ierr = ats_set_init_clm_data_f90(T, Sl, Si)
   end subroutine ats_set_init_clm_data

   !---------------------------------------------------------------------------!
   ! ats_set_clm_data
   !---------------------------------------------------------------------------!

   subroutine ats_set_clm_data(e_flux, w_flux, ierr)
      use :: ats_data
      implicit none
      type(c_ptr), value :: e_flux
      type(c_ptr), value :: w_flux
      integer(int32_t) :: ierr

      ierr = ats_set_clm_data_f90(e_flux, w_flux)
   end subroutine ats_set_clm_data

   !---------------------------------------------------------------------------!
   ! ats_get_clm_data
   !---------------------------------------------------------------------------!

   subroutine ats_get_clm_data(T, Sl, Si, ierr)
      use :: ats_data
      implicit none
      type(c_ptr), value :: T
      type(c_ptr), value :: Sl
      type(c_ptr), value :: Si
      integer(int32_t) :: ierr

      ierr = ats_get_clm_data_f90(T, Sl, Si)
   end subroutine ats_get_clm_data

   !---------------------------------------------------------------------------!
   ! ats_advance
   !---------------------------------------------------------------------------!

   subroutine ats_advance(dt, force_viz)
      use :: ats_data
      implicit none
      real(c_double) :: dt
      integer(int32_t) :: force_viz
      integer(int32_t) :: ierr

      ierr = ats_advance_f90(dt, force_viz)
   end subroutine ats_advance

end module ats_interface
