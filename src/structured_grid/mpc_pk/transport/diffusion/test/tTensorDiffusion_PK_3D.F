#include <CONSTANTS.H>
#include <REAL.H>

#include <tTensorDiffusion_PK_F.H>

      subroutine FORT_FILLCRSE (
     $     crse, DIMS(crse),
     $     h, nc, dlo, dhi, plo,
     $     flag)
      implicit none
      integer nc
      integer DIMDEC(crse)
      integer dlo(BL_SPACEDIM), dhi(BL_SPACEDIM)
      REAL_T crse(DIMV(crse),nc)
      REAL_T h(BL_SPACEDIM)
      REAL_T plo(BL_SPACEDIM)
      integer flag

      print *,'tTensorDiffusion_PK_3D not ready yet'
      end
      

      
      subroutine FORT_FILLFINE (
     $     fine, DIMS(fine),
     $     h, nc, dlo, dhi, plo,
     $     lo, hi, flag)
      implicit none
      integer nc
      integer DIMDEC(fine)
      integer dlo(BL_SPACEDIM), dhi(BL_SPACEDIM)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T fine(DIMV(fine),nc)
      REAL_T h(BL_SPACEDIM)
      REAL_T plo(BL_SPACEDIM)
      integer flag

      integer  i, j, n
      REAL_T x, y

      print *,'tTensorDiffusion_PK_3D not ready yet'
      end
      
      subroutine FORT_FILLRHS (
     $     rhs, DIMS(rhs),
     $     h, nc, dlo, dhi, plo
     $     )
      implicit none
      integer nc
      integer DIMDEC(rhs)
      integer dlo(BL_SPACEDIM), dhi(BL_SPACEDIM)
      REAL_T rhs(DIMV(rhs),nc)
      REAL_T h(BL_SPACEDIM)
      REAL_T plo(BL_SPACEDIM)
      print *,'tTensorDiffusion_PK_3D not ready yet'
      end

      subroutine FORT_MAKEMU(
     $     mu,DIMS(mu),
     $     h, idir, aT, aL, u, v )
      implicit none
      integer idir
      integer DIMDEC(mu)
      REAL_T mu(DIMV(mu))
      REAL_T h(BL_SPACEDIM)
      REAL_T aT, aL, u, v

      integer i,j
      REAL_T x,y
      REAL_T vel, velInv
      print *,'tTensorDiffusion_PK_3D not ready yet'
      end

c     \begin{array}{cc}
c     \a_T v + (\a_L-\a_T)\dfrac{v_x^2}{v} & (\a_L-\a_T)\dfrac{v_x v_y}{v} \\
c     (\a_L-\a_T)\dfrac{v_y v_x}{v} & \a_T v + (\a_L-\a_T)\dfrac{v_y^2}{v}
c     \end{array}

      subroutine FORT_MAKEMU1(
     $     mu1,DIMS(mu1),
     $     h, idir, aT, aL, u, v )
      implicit none
      integer idir
      integer DIMDEC(mu1)
      REAL_T mu1(DIMV(mu1))
      REAL_T h(BL_SPACEDIM)
      REAL_T aT, aL, u, v

      integer i,j
      REAL_T x,y
      print *,'tTensorDiffusion_PK_3D not ready yet'
      end
