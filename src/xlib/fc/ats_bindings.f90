!------------------------------------------------------------------------------!
! ATS
!
! License: see $ATS_DIR/COPYRIGHT
!
! C Language binding layer using Fortran 2003 interoperability features.
!------------------------------------------------------------------------------!

module ats_bindings

interface

   !---------------------------------------------------------------------------!
   ! ats_init_f90
   !---------------------------------------------------------------------------!

   function ats_init_f90(mpi_comm, type_ids, num_type_ids, num_types) &
      result(ierr) bind(C, name="ats_init_f90")
      use, intrinsic :: ISO_C_BINDING
      use :: ats_data
      implicit none
      integer(int32_t) :: mpi_comm
      type(c_ptr) :: type_ids
      integer(int32_t) :: num_type_ids
      integer(int32_t) :: num_types
      integer(int32_t) :: ierr
   end function ats_init_f90

   !---------------------------------------------------------------------------!
   ! ats_finalize_f90
   !---------------------------------------------------------------------------!

   function ats_finalize_f90() &
      result(ierr) bind(C, name="ats_finalize_f90")
      use, intrinsic :: ISO_C_BINDING
      use :: ats_data
      implicit none
      integer(int32_t) :: ierr
   end function ats_finalize_f90

   !---------------------------------------------------------------------------!
   ! ats_set_init_clm_data_f90
   !---------------------------------------------------------------------------!

   function ats_set_init_clm_data_f90(T, Sl, Si) &
      result(ierr) bind(C, name="ats_set_init_clm_data_f90")
      use, intrinsic :: ISO_C_BINDING
      use :: ats_data
      implicit none
      type(c_ptr), value :: T
      type(c_ptr), value :: Sl
      type(c_ptr), value :: Si
      integer(int32_t) :: ierr
   end function ats_set_init_clm_data_f90

   !---------------------------------------------------------------------------!
   ! ats_set_clm_data_f90
   !---------------------------------------------------------------------------!

   function ats_set_clm_data_f90(e_flux, w_flux) &
      result(ierr) bind(C, name="ats_set_clm_data_f90")
      use, intrinsic :: ISO_C_BINDING
      use :: ats_data
      implicit none
      type(c_ptr), value :: e_flux
      type(c_ptr), value :: w_flux
      integer(int32_t) :: ierr
   end function ats_set_clm_data_f90

   !---------------------------------------------------------------------------!
   ! ats_get_clm_data_f90
   !---------------------------------------------------------------------------!

   function ats_get_clm_data_f90(T, Sl, Si) &
      result(ierr) bind(C, name="ats_get_clm_data_f90")
      use, intrinsic :: ISO_C_BINDING
      use :: ats_data
      implicit none
      type(c_ptr), value :: T
      type(c_ptr), value :: Sl
      type(c_ptr), value :: Si
      integer(int32_t) :: ierr
   end function ats_get_clm_data_f90

   !---------------------------------------------------------------------------!
   ! ats_advance_f90
   !---------------------------------------------------------------------------!

   function ats_advance_f90(dt, force_viz) &
      result(ierr) bind(C, name="ats_advance_f90")
      use, intrinsic :: ISO_C_BINDING
      use :: ats_data
      implicit none
      real(c_double) :: dt
      integer(int32_t) :: force_viz
      integer(int32_t) :: ierr
   end function ats_advance_f90

end interface

end module
