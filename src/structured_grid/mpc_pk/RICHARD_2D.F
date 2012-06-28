#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "RICHARD_F.H"
#include "ArrayLim.H"

#define SDIM 2

      subroutine FORT_RICHARD_GXPA(e,DIMS(e),
     &     c,DIMS(c),lo,hi,dx,a,dir,nc)
      implicit none      
      integer dir, nc
      integer DIMDEC(e)
      integer DIMDEC(c)
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  e(DIMV(e),nc)
      REAL_T  c(DIMV(c),nc)
      REAL_T  a(nc)
      
      integer i, j, n
      REAL_T tmult
      if (dir.eq.0) then
         do n=1,nc
            tmult = a(n)/dx(dir)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  e(i,j,n) = tmult*(c(i,j,n) - c(i-1,j,n))
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  e(i,j,n) = tmult*(c(i,j,n) - c(i,j-1,n))
               enddo
            enddo
         enddo
      else
         print *, 'FORT_RICHARD_GXPA::bad dir: ',dir
         call bl_pd_abort()
      endif
      end     


      subroutine FORT_RICHARD_YMULTX(Y,DIMS(Y),X,DIMS(X),lo,hi,nc)
      implicit none      
      integer DIMDEC(Y)
      integer DIMDEC(X)
      integer lo(SDIM),hi(SDIM)
      REAL_T  Y(DIMV(Y),nc)
      REAL_T  X(DIMV(X),nc)
      integer nc
      
      integer i, j, n

      do n=1,nc
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Y(i,j,n) = Y(i,j,n) * X(i,j,n)
            enddo
         enddo
      enddo
      end     

      subroutine FORT_RICHARD_CTE_UPW(coef,DIMS(coef),
     &     lam,DIMS(lam),vel,DIMS(vel),lo,hi,dir,nc)
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

