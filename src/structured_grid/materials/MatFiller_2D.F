#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <ArrayLim.H>
#include <MatFiller_F.H> 

#define SDIM 2

      subroutine FORT_FILLP (p, DIMS(p),id,DIMS(id),lo,hi,vals)
      implicit none
      integer DIMDEC(p)
      integer DIMDEC(id)
      integer lo(SDIM), hi(SDIM)
      REAL_T   p(DIMV(p))
      integer id(DIMV(id))
      REAL_T vals(*)
      integer i,j
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            p(i,j) = vals(id(i,j)+1)
         end do
      end do
      end

      subroutine FORT_CRSNHARM (fdat, DIMS(fdat),
     &     cdat,DIMS(cdat),lo,hi,ratio,nc)
      implicit none
      integer DIMDEC(fdat)
      integer DIMDEC(cdat)
      integer lo(SDIM), hi(SDIM),nc
      integer ratio(SDIM)
      REAL_T  fdat(DIMV(fdat),nc)
      REAL_T  cdat(DIMV(cdat),nc)

      integer i,j,ii,jj,fi,fj,comp(SDIM)
      REAL_T ktmp, fac1, fac2

      if (nc.eq.SDIM) then
         do i=1,SDIM
            comp(i) = i
         enddo
      else if (nc.eq.1) then
         do i=1,SDIM
            comp(i) = 1
         enddo
      else
         call bl_abort("Invalid nc in MatFiller::FORT_CRSNHARM")
      endif

      cdat(lo(1):hi(1),lo(2):hi(2),1:SDIM) = zero
      fac1 = DBLE(ratio(1)) / DBLE(ratio(2))
      fac2 = DBLE(ratio(2)) / DBLE(ratio(1))

c     ::: x-dir  Harmonic in x, arithmetic in y
      do jj = 0,ratio(2)-1
         do j = lo(2),hi(2)
            fj = ratio(2)*j
            do i = lo(1),hi(1)
               fi = ratio(1)*i
               ktmp = zero
               do ii = 0,ratio(1)-1
                  ktmp = ktmp + one/fdat(fi+ii,fj+jj,comp(1))
               enddo
               cdat(i,j,1) = cdat(i,j,1) + fac1 / ktmp
            end do
         end do
      end do

c     ::: y-dir  Harmonic in y, arithmetic in x
      do ii = 0,ratio(1)-1
         do j = lo(2),hi(2)
            fj = ratio(2)*j
            do i = lo(1),hi(1)
               fi = ratio(1)*i
               ktmp = zero
               do jj = 0,ratio(2)-1
                  ktmp = ktmp + one/fdat(fi+ii,fj+jj,comp(2))
               enddo
               cdat(i,j,2) = cdat(i,j,2) + fac2 / ktmp
            end do
         end do
      end do
      end


      subroutine FORT_CRSNARITH (fdat,DIMS(fdat),
     &     cdat,DIMS(cdat),lo,hi,ratio,nc)

      implicit none
      integer DIMDEC(cdat)
      integer DIMDEC(fdat)
      integer lo(SDIM), hi(SDIM),nc
      integer ratio(SDIM)
      REAL_T  cdat(DIMV(cdat),nc)
      REAL_T  fdat(DIMV(fdat),nc)

      integer i,j,ii,jj,fi,fj
      REAL_T  fac

      cdat(lo(1):hi(1),lo(2):hi(2),1:nc) = zero
      fac = 1.d0 / dble(ratio(1)*ratio(2))
      do ii = 0,ratio(1)-1
         do jj = 0,ratio(2)-1
            do j = lo(2), hi(2)
               fj = ratio(2)*j
               do i = lo(1), hi(1)
                  fi = ratio(1)*i
                  cdat(i,j,1:nc) = cdat(i,j,1:nc)  +  fac * fdat(fi+ii,fj+jj,1:nc)
               end do
            end do
         end do
      end do
      end


      subroutine FORT_REFINEPC (cdat,DIMS(cdat),
     &     fdat,DIMS(fdat),lo,hi,ratio,nc)

      implicit none
      integer DIMDEC(fdat)
      integer DIMDEC(cdat)
      integer lo(SDIM), hi(SDIM), nc
      integer ratio(SDIM)
      REAL_T  fdat(DIMV(fdat),nc)
      REAL_T  cdat(DIMV(cdat),nc)

      integer i,j,ii,jj,fi,fj
      do ii = 0,ratio(1)-1
         do jj = 0,ratio(2)-1
            do j = lo(2), hi(2)
               fj = ratio(2)*j
               do i = lo(1), hi(1)
                  fi = ratio(1)*i
                  fdat(fi+ii,fj+jj,1:nc) = cdat(i,j,1:nc)
               end do
            end do
         end do
      end do
      end

