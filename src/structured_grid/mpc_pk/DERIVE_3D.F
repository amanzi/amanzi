#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H" 
#include "CONSTANTS.H" 
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H" 
#include "DERIVE_F.H" 
#include "PROB_PM_F.H"
#include "ArrayLim.H"

#define SDIM 3

      subroutine FORT_DERGRDPX (grdpx,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     $                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
c
c     This routine computes pressure gradient in X direction
c
      integer lo(3), hi(3)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(3), domhi(3)
      integer nv, ncomp
      integer bc(3,2,ncomp)
      REAL_T  delta(3), xlo(3), time, dt
      REAL_T  grdpx(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)

      integer level, grid_no
      integer i,j,k
      REAL_T  dtmp

      dtmp = 1.0d0 / delta(1)

!$omp parallel do private(i,j,k)
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               grdpx(i,j,k,1) = 0.5d0*(dat(i+1,j,k,1) - dat(i-1,j,k,1))*dtmp
            end do
         end do 
      end do 
!$omp end parallel do

      if (lo(1) .eq. domlo(1)) then
         i = lo(1)
         if (bc(1,1,1) .eq. FOEXTRAP .or. bc(1,1,1) .eq. EXT_DIR) then
            do k = lo(3),hi(3)
               do j = lo(2), hi(2)
                 grdpx(i,j,k,1) = (-1.5d0 * dat(i,j,k,1) + 2.d0 * dat(i+1,j,k,1) - 0.5d0 * dat(i+2,j,k,1))*dtmp
              end do
           end do
        end if
      end if

      if (hi(1) .eq. domhi(1)) then
         i = hi(1)
         if (bc(1,2,1) .eq. FOEXTRAP .or. bc(1,2,1) .eq. EXT_DIR) then
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                 grdpx(i,j,k,1) = (1.5d0 * dat(i,j,k,1) - 2.d0 * dat(i-1,j,k,1) + 0.5d0 * dat(i-2,j,k,1))*dtmp
              end do
           end do
         end if
      end if

      end

      subroutine FORT_DERGRDPY (grdpy,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     $                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
c
c     This routine computes pressure gradient in X direction
c
      integer lo(3), hi(3)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(3), domhi(3)
      integer nv, ncomp
      integer bc(3,2,ncomp)
      REAL_T  delta(3), xlo(3), time, dt
      REAL_T  grdpy(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)

      integer level, grid_no
      integer i,j,k
      REAL_T  dtmp

      dtmp = 1.0d0 / delta(2)

!$omp parallel do private(i,j,k)
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               grdpy(i,j,k,1) = 0.5d0*(dat(i,j+1,k,1) - dat(i,j-1,k,1))*dtmp
            end do
         end do 
      end do 
!$omp end parallel do

      if (lo(2) .eq. domlo(2) .and. bc(2,1,1) .eq. FOEXTRAP) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1), hi(1)
               grdpy(i,j,k,1) = (-1.5d0 * dat(i,j,k,1) + 2.d0 * dat(i,j+1,k,1) - 0.5d0 * dat(i,j+2,k,1))*dtmp
            end do
         end do
      end if

      if (hi(2) .eq. domhi(2) .and. bc(2,2,1) .eq. FOEXTRAP) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1), hi(1)
               grdpy(i,j,k,1) = (1.5d0 * dat(i,j,k,1) - 2.d0 * dat(i,j-1,k,1) + 0.5d0 * dat(i,j-2,k,1))*dtmp
            end do
         end do
      end if

      end

      subroutine FORT_DERGRDPZ (grdpz,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     $                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
c
c     This routine computes pressure gradient in X direction
c
      integer lo(3), hi(3)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(3), domhi(3)
      integer nv, ncomp
      integer bc(3,2,ncomp)
      REAL_T  delta(3), xlo(3), time, dt
      REAL_T  grdpz(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)
      integer level, grid_no

      integer i,j,k
      REAL_T  dtmp

      dtmp = 1.0d0 / delta(3)

!$omp parallel do private(i,j,k)
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               grdpz(i,j,k,1) = 0.5d0*(dat(i,j,k+1,1) - dat(i,j,k-1,1))*dtmp
            end do
         end do 
      end do 
!$omp end parallel do

      if (lo(3) .eq. domlo(3) .and. bc(3,1,1) .eq. FOEXTRAP) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1), hi(1)
               grdpz(i,j,k,1) = (-1.5d0 * dat(i,j,k,1) + 2.d0 * dat(i,j,k+1,1) - 0.5d0 * dat(i,j,k+2,1))*dtmp
            end do
         end do
      end if

      if (hi(3) .eq. domhi(3) .and. bc(3,2,1) .eq. FOEXTRAP) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1), hi(1)
               grdpz(i,j,k,1) = (1.5d0 * dat(i,j,k,1) - 2.d0 * dat(i,j,k-1,1) + 0.5d0 * dat(i,j,k-2,1))*dtmp
            end do
         end do
      end if

      end


      subroutine FORT_DERGRDN  (grdn,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     $                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
c
c     This routine computes pressure gradient in dn/dx + dn/dy direction
c
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(SDIM), domhi(SDIM)
      integer nv, ncomp
      integer bc(SDIM,2,ncomp)
      REAL_T  delta(SDIM), xlo(SDIM), time, dt
      REAL_T  grdn(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)
      integer level, grid_no

      integer i,j,k
      REAL_T  xd1, xd2, xd3

      xd1 = 1.d0 / delta(1)
      xd2 = 1.d0 / delta(2)
      xd3 = 1.d0 / delta(3)

!$omp parallel do private(i,j,k)
      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               grdn(i,j,k,1) = half * abs(dat(i,j,k+1,1) - dat(i,j,k-1,1)) * xd3 +
     &              half * abs(dat(i,j+1,k,1) - dat(i,j-1,k,1)) * xd2 +
     &              half * abs(dat(i+1,j,k,1) - dat(i-1,j,k,1)) * xd1
            end do
         end do
      end do
!$omp end parallel do
 
      if (lo(3) .eq. domlo(3) .and. bc(3,1,1) .eq. FOEXTRAP) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1), hi(1)
               grdn(i,j,k,1) = abs(-1.5d0 * dat(i,j,k,1) 
     &              + 2.d0 * dat(i,j,k+1,1) - 0.5d0 * dat(i,j,k+2,1)) * xd3 +
     &              half * abs(dat(i+1,j,k,1) - dat(i-1,j,k,1)) * xd1 +
     &              half * abs(dat(i,j+1,k,1) - dat(i,j-1,k,1)) * xd2
            end do
         end do
      end if

      if (hi(3) .eq. domhi(3) .and. bc(3,2,1) .eq. FOEXTRAP) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1), hi(1)
               grdn(i,j,k,1) = abs(1.5d0 * dat(i,j,k,1) 
     &              - 2.d0 * dat(i,j,k-1,1) + 0.5d0 * dat(i,j,k-2,1)) * xd3 +
     &              half *abs(dat(i+1,j,k,1) - dat(i-1,j,k,1)) * xd1 +
     &              half *abs(dat(i,j+1,k,1) - dat(i,j-1,k,1)) * xd2 
            end do
         end do
      end if

      if (lo(2) .eq. domlo(2) .and. bc(2,1,1) .eq. FOEXTRAP) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1), hi(1)
               grdn(i,j,k,1) = abs(-1.5d0 * dat(i,j,k,1) 
     &              + 2.d0 * dat(i,j+1,k,1) - 0.5d0 * dat(i,j+2,k,1))* xd2 +
     &              half * abs(dat(i+1,j,k,1) - dat(i-1,j,k,1)) * xd1 +
     &              half * abs(dat(i,j,k+1,1) - dat(i,j,k-1,1)) * xd3
            end do
         end do
      end if

      if (hi(2) .eq. domhi(2) .and. bc(2,2,1) .eq. FOEXTRAP) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1), hi(1)
               grdn(i,j,k,1) = abs(1.5d0 * dat(i,j,k,1) 
     &              - 2.d0 * dat(i,j-1,k,1) + 0.5d0 * dat(i,j-2,k,1)) * xd2 +
     &              half *abs(dat(i+1,j,k,1) - dat(i-1,j,k,1)) * xd1 +
     &              half *abs(dat(i,j,k+1,1) - dat(i,j,k-1,1)) * xd3 
            end do
         end do
      end if

      if (lo(1) .eq. domlo(1)) then
         i = lo(1)
         if (bc(1,1,1) .eq. FOEXTRAP .or. bc(1,1,1) .eq. EXT_DIR) then
            do k = lo(3),hi(3)
               do j = lo(2), hi(2)
                 grdn(i,j,k,1) = abs(-1.5d0 * dat(i,j,k,1) 
     &                + 2.d0 * dat(i+1,j,k,1) - 0.5d0 * dat(i+2,j,k,1)) * xd1 +
     &                half * abs(dat(i,j+1,k,1) - dat(i,j-1,k,1)) * xd2 +
     &                half * abs(dat(i,j,k+1,1) - dat(i,j,k-1,1)) * xd3
              end do
           end do
         end if
      end if

      if (hi(1) .eq. domhi(1)) then
         i = hi(1)
         if (bc(1,2,1) .eq. FOEXTRAP .or. bc(1,2,1) .eq. EXT_DIR) then
            do k = lo(3),hi(3)
               do j = lo(2), hi(2)
                 grdn(i,j,k,1) = abs(1.5d0 * dat(i,j,k,1) 
     &                - 2.d0 * dat(i-1,j,k,1) + 0.5d0 * dat(i-2,j,k,1)) * xd1 +
     &                half * abs(dat(i,j+1,k,1) - dat(i,j-1,k,1)) * xd2 +
     &                half * abs(dat(i,j,k+1,1) - dat(i,j,k-1,1)) * xd3
              end do
           end do
         end if
      end if

      end
