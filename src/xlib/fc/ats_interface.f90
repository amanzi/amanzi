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

   subroutine ats_init(comm, ierr)
      use :: ats_data
      implicit none
      integer(int32_t) :: comm
      integer(int32_t) :: ierr

      ierr = ats_init_f90(comm)
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

end module ats_interface
