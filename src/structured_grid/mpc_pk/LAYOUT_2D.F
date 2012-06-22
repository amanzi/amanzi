#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "LAYOUT_F.H"
#include "ArrayLim.H"

#define SDIM 2

      subroutine FORT_CC_TO_EC_GRAD(grad,DIMS(grad),
     &     c,DIMS(c),lo,hi,dx,mult,dir,nc)
      implicit none      
      integer DIMDEC(grad)
      integer DIMDEC(c)
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(0:SDIM-1)
      REAL_T  grad(DIMV(grad),0:nc-1)
      REAL_T  c(DIMV(c),0:nc-1)
      REAL_T  mult
      integer dir, nc
      
      integer i, j, n
      REAL_T tmult
      tmult = mult/dx(dir)
      if (dir.eq.0) then
         do n=0,nc-1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  grad(i,j,n) = tmult*(c(i,j,n) - c(i-1,j,n))
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=0,nc-1
            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  grad(i,j,n) = tmult*(c(i,j,n) - c(i,j-1,n))
               enddo
            enddo
         enddo
      else
         print *, 'FORT_CC_TO_EC_GRAD::bad dir: ',dir
         call bl_pd_abort()
      endif
      end     

      subroutine FORT_EC_TO_CC_DIV(div,DIMS(div),
     &     e,DIMS(e),lo,hi,dx,mult,dir,nc)
      implicit none      
      integer DIMDEC(div)
      integer DIMDEC(e)
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(0:SDIM-1)
      REAL_T  div(DIMV(div),0:nc-1)
      REAL_T  e(DIMV(e),0:nc-1)
      REAL_T  mult
      integer dir, nc
      
      integer i, j, n
      REAL_T tmult
      tmult = mult/dx(dir)
      if (dir.eq.0) then
         do n=0,nc-1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  div(i,j,n) = div(i,j,n) + tmult*(e(i+1,j,n) - e(i,j,n))
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=0,nc-1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  div(i,j,n) = div(i,j,n) + tmult*(e(i,j+1,n) - e(i,j,n))
               enddo
            enddo
         enddo
      else
         print *, 'FORT_EC_TO_CC_DIR::bad dir: ',dir
         call bl_pd_abort()
      endif
      end     
