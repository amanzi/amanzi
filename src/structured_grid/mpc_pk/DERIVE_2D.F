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

#define SDIM 2

      subroutine FORT_DERGRDPX (grdpx,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     $                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
c
c     This routine computes pressure gradient in X direction
c
      integer lo(2), hi(2)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(2), domhi(2)
      integer nv, ncomp
      integer bc(2,2,ncomp)
      REAL_T  delta(2), xlo(2), time, dt
      REAL_T  grdpx(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)
      integer level, grid_no
c
      integer i,j
c
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            grdpx(i,j,1) = half * (dat(i+1,j,1) - dat(i-1,j,1)) / delta(1)
         end do
      end do 

      if (lo(1) .eq. domlo(1)) then
         i = lo(1)
         if (bc(1,1,1) .eq. FOEXTRAP .or. bc(1,1,1) .eq. EXT_DIR) then
           do j = lo(2), hi(2)
              grdpx(i,j,1) = (-1.5d0 * dat(i,j,1) + 2.d0 * dat(i+1,j,1) 
     &             - 0.5d0 * dat(i+2,j,1))/delta(1)
           end do
         end if
      end if

      if (hi(1) .eq. domhi(1)) then
         i = hi(1)
         if (bc(1,2,1) .eq. FOEXTRAP .or. bc(1,2,1) .eq. EXT_DIR) then
           do j = lo(2), hi(2)
              grdpx(i,j,1) = (1.5d0 * dat(i,j,1) - 2.d0 * dat(i-1,j,1) 
     &             + 0.5d0 * dat(i-2,j,1))/delta(1)
           end do
         end if
      end if

      end

      subroutine FORT_DERGRDPY (grdpy,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     $                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
c
c     This routine computes pressure gradient in Y direction
c
      integer lo(2), hi(2)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(2), domhi(2)
      integer nv, ncomp
      integer bc(2,2,ncomp)
      REAL_T  delta(2), xlo(2), time, dt
      REAL_T  grdpy(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)
      integer level, grid_no
c
      integer i,j
c
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            grdpy(i,j,1) = half * (dat(i,j+1,1) - dat(i,j-1,1)) / delta(2)
         end do
      end do 

      if (lo(2) .eq. domlo(2) .and. bc(2,1,1) .eq. FOEXTRAP) then
         j = lo(2)
         do i = lo(1), hi(1)
            grdpy(i,j,1) = (-1.5d0 * dat(i,j,1) + 2.d0 * dat(i,j+1,1) 
     &           - 0.5d0 * dat(i,j+2,1))/delta(2)
         end do
      end if

      if (hi(2) .eq. domhi(2) .and. bc(2,2,1) .eq. FOEXTRAP) then
         j = hi(2)
         do i = lo(1), hi(1)
            grdpy(i,j,1) = (1.5d0 * dat(i,j,1) - 2.d0 * dat(i,j-1,1) 
     &           + 0.5d0 * dat(i,j-2,1))/delta(2)
         end do
      end if

      end

      subroutine FORT_DERGRDN  (grdn,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     $                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
c
c     This routine computes pressure gradient in dn/dx + dn/dy direction
c
      integer lo(2), hi(2)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(2), domhi(2)
      integer nv, ncomp
      integer bc(2,2,ncomp)
      REAL_T  delta(2), xlo(2), time, dt,x
      REAL_T  grdn(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)
      integer level, grid_no
      integer i,j,n,idx
      grdn = zero
      do n = 1,ncomp
         idx = n
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               grdn(i,j,1) = grdn(i,j,1) 
     &              + half*abs(dat(i,j+1,idx)-dat(i,j-1,idx))/delta(2) 
     &              + half*abs(dat(i+1,j,idx)-dat(i-1,j,idx))/delta(1)
            end do
         end do 

         if (lo(2) .eq. domlo(2) .and. bc(2,1,1) .eq. FOEXTRAP) then
            j = lo(2)
            do i = lo(1), hi(1)
               grdn(i,j,1) = grdn(i,j,1)+abs(-1.5d0 * dat(i,j,idx) 
     &              + 2.d0 * dat(i,j+1,idx) 
     &              - 0.5d0 * dat(i,j+2,idx))/delta(2)
     &              + half *abs(dat(i+1,j,idx)-dat(i-1,j,idx))/delta(1)
            end do
         end if


         if (hi(2) .eq. domhi(2) .and. bc(2,2,1) .eq. FOEXTRAP) then
            j = hi(2)
            do i = lo(1), hi(1)
               grdn(i,j,1) = grdn(i,j,1)+abs(1.5d0 * dat(i,j,idx) 
     &              - 2.d0 * dat(i,j-1,idx) 
     &              + 0.5d0 * dat(i,j-2,idx))/delta(2)
     &              + half *abs(dat(i+1,j,idx)-dat(i-1,j,idx))/delta(1)
            end do
         end if

         if (lo(1) .eq. domlo(1)) then
            i = lo(1)
            if (bc(1,1,1) .eq. FOEXTRAP .or. bc(1,1,1) .eq. EXT_DIR) then
               do j = lo(2), hi(2)
                  grdn(i,j,1) = grdn(i,j,1)+abs(-1.5d0 * dat(i,j,idx) 
     &                 + 2.d0 * dat(i+1,j,idx)
     &                 - 0.5d0 * dat(i+2,j,idx))/delta(1)
     &                 + half*abs(dat(i,j+1,idx)-dat(i,j-1,idx))/delta(2)
               end do
            end if
         end if

         if (hi(1) .eq. domhi(1)) then
            i = hi(1)
            if (bc(1,2,1) .eq. FOEXTRAP .or. bc(1,2,1) .eq. EXT_DIR) then
               do j = lo(2), hi(2)
                  grdn(i,j,1) = grdn(i,j,1)+abs(1.5d0 * dat(i,j,idx) 
     &                 - 2.d0 * dat(i-1,j,idx) 
     &                 + 0.5d0* dat(i-2,j,idx))/delta(1)
     &                 + half*abs(dat(i,j+1,idx)-dat(i,j-1,idx))/delta(2)
               end do
            end if
         end if  

      end do

      end
