#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "MFTOWER_F.H"
#include "ArrayLim.H"

#define SDIM 3

#if 0

 NEED to edit this for 3D!

      subroutine FORT_CC_TO_EC_UPW(coef,DIMS(coef),
     &     lam,DIMS(lam),vel,DIMS(vel),lo,hi,mult,dir,nc)
      implicit none      
      integer DIMDEC(coef)
      integer DIMDEC(lam)
      integer DIMDEC(vel)
      integer lo(SDIM),hi(SDIM)
      REAL_T  coef(DIMV(coef),0:nc-1)
      REAL_T  lam(DIMV(lam),0:nc-1)
      REAL_T  vel(DIMV(vel))
      REAL_T  mult
      integer dir, nc
      
      integer i, j, n, do_upwind

      do_upwind = 1

      if (dir.eq.0) then
         if (do_upwind.eq.1) then
            do n=0,nc-1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     coef(i,j,n)
     &                    = merge(lam(i-1,j,n),lam(i,j,n),vel(i,j)>=zero)
                  enddo
               enddo
            enddo
         else
            do n=0,nc-1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     coef(i,j,n) = 0.5d0 * (lam(i-1,j,n)+lam(i,j,n))
                  enddo
               enddo
            enddo
         endif
      else if (dir.eq.1) then
         if (do_upwind.eq.1) then
            do n=0,nc-1
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     coef(i,j,n) 
     &                    = merge(lam(i,j-1,n),lam(i,j,n),vel(i,j)>=zero)
                  enddo
               enddo
            enddo
         else
            do n=0,nc-1
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     coef(i,j,n) = 0.5d0 * (lam(i,j-1,n)+lam(i,j,n))
                  enddo
               enddo
            enddo
         endif
      else
         print *, 'FORT_CALC_RICHARD_H::bad dir: ',dir
         call bl_pd_abort()
      endif
      end     
#endif

