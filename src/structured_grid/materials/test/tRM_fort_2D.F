#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <ArrayLim.H>
#include <RockManager_F.H> 
#include <tRM_fort_F.H> 

#define SDIM 2

      subroutine FORT_TEST_PCAP (rmID, time, lo, hi,
     &     s, DIMS(s),
     &     m, DIMS(m),
     &     p, DIMS(p))
      implicit none
      integer rmID
      REAL_T  time
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(s)
      integer DIMDEC(m)
      integer DIMDEC(p)
      REAL_T  s(DIMV(s))
      integer m(DIMV(m))
      REAL_T  p(DIMV(p))
      integer j,Nx

      Nx = hi(1) - lo(1) + 1
      do j = lo(2),hi(2)
         call ROCK_MANAGER_PCAP(rmID,s(ARG_L1(s),j),m(ARG_L1(m),j),
     &        time,p(ARG_L1(p),j),Nx)
      end do
      end

      subroutine FORT_TEST_INVPCAP (rmID, time, lo, hi,
     &     p, DIMS(p),
     &     m, DIMS(m),
     &     s, DIMS(s))
      implicit none
      integer rmID
      REAL_T  time
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(p)
      integer DIMDEC(m)
      integer DIMDEC(s)
      REAL_T  p(DIMV(p))
      integer m(DIMV(m))
      REAL_T  s(DIMV(s))
      integer j,Nx

      Nx = hi(1) - lo(1) + 1
      do j = lo(2),hi(2)
         call ROCK_MANAGER_INVPCAP(rmID,p(ARG_L1(p),j),m(ARG_L1(m),j),
     &        time,s(ARG_L1(s),j),Nx)
      end do
      end

      subroutine FORT_TEST_RELPERM (rmID, time, lo, hi,
     &     s,  DIMS(s),
     &     m,  DIMS(m),
     &     kr, DIMS(kr))
      implicit none
      integer rmID
      REAL_T  time
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(s)
      integer DIMDEC(m)
      integer DIMDEC(kr)
      REAL_T  s(DIMV(s))
      integer m(DIMV(m))
      REAL_T  kr(DIMV(kr))
      integer j,Nx

      Nx = hi(1) - lo(1) + 1
      do j = lo(2),hi(2)
         call ROCK_MANAGER_RELPERM(rmID,s(ARG_L1(s),j),m(ARG_L1(m),j),
     &        time,kr(ARG_L1(kr),j),Nx)
      end do
      end

      subroutine FORT_TEST_DSDPCAP (rmID, time, lo, hi,
     &     s,  DIMS(s),
     &     m,  DIMS(m),
     &     dsdp, DIMS(dsdp))
      implicit none
      integer rmID
      REAL_T  time
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(s)
      integer DIMDEC(m)
      integer DIMDEC(dsdp)
      REAL_T  s(DIMV(s))
      integer m(DIMV(m))
      REAL_T  dsdp(DIMV(dsdp))
      integer j,Nx

      Nx = hi(1) - lo(1) + 1
      do j = lo(2),hi(2)
         call ROCK_MANAGER_DSDPCAP(rmID,s(ARG_L1(s),j),m(ARG_L1(m),j),
     &        time,dsdp(ARG_L1(dsdp),j),Nx)
      end do
      end

