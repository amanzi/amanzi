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

   function ats_init_f90(comm) &
      result(ierr) bind(C, name="ats_init_f90")
      use, intrinsic :: ISO_C_BINDING
      use :: ats_data
      implicit none
      integer(int32_t) :: comm
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

end interface

end module
