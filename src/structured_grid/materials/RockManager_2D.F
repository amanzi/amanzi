#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <ArrayLim.H>
#include <RockManager_F.H> 

#define SDIM 2

      subroutine FORT_FILLPMAT (p,DIMS(p),id,DIMS(id),vals,DIMS(vals),
     &     myid,lo,hi,nc)
      implicit none
      integer nc
      integer DIMDEC(p)
      integer DIMDEC(id)
      integer DIMDEC(vals)
      integer lo(SDIM), hi(SDIM)
      REAL_T   p(DIMV(p),nc)
      integer id(DIMV(id))
      REAL_T  vals(DIMV(vals),nc)
      integer myid
      integer i,j,n
      do n=1,nc
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (id(i,j).eq.myid) then
                  p(i,j,n) = vals(i,j,n)
               endif
            enddo
         enddo
      enddo
      end

