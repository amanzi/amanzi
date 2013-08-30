#include <CONSTANTS.H>
#include <REAL.H>

#include <tTensorDiffusion_PK_F.H>

      subroutine FORT_FILLCRSE (
     $     crse, DIMS(crse),
     $     h, nc, dlo, dhi, plo
     $     )
      implicit none
      integer nc
      integer DIMDEC(crse)
      integer dlo(BL_SPACEDIM), dhi(BL_SPACEDIM)
      REAL_T crse(DIMV(crse),nc)
      REAL_T h(BL_SPACEDIM)
      REAL_T plo(BL_SPACEDIM)
c
      integer  i, j, n
      REAL_T x, y

      do j=ARG_L2(crse), ARG_H2(crse)
         y = (j + 0.5)*h(2)
         if (j.lt.dlo(2)) y = plo(2)
         if (j.gt.dhi(2)) y = plo(2)+(dhi(2)-dlo(2)+1)*h(2)
         do i=ARG_L1(crse), ARG_H1(crse)
            x = (i + 0.5)*h(1)
            if (i.lt.dlo(1)) x = plo(1)
            if (i.gt.dhi(1)) x = plo(1)+(dhi(1)-dlo(1)+1)*h(1)
            do n=1,nc
               crse(i,j,n)=-(x-0.5d0)**2
c               crse(i,j,n)=-(y-0.5d0)**2
c               crse(i,j,n)=y
            enddo
         enddo
      enddo
      
      end
      

      
      subroutine FORT_FILLFINE (
     $     fine, DIMS(fine),
     $     h, nc, dlo, dhi, plo
     $     )
      implicit none
      integer nc
      integer DIMDEC(fine)
      integer dlo(BL_SPACEDIM), dhi(BL_SPACEDIM)
      REAL_T fine(DIMV(fine),nc)
      REAL_T h(BL_SPACEDIM)
      REAL_T plo(BL_SPACEDIM)

      integer  i, j, n
      REAL_T x, y

      do j=ARG_L2(fine), ARG_H2(fine)
         y = plo(2) + (j + 0.5)*h(2)
         if (j.lt.dlo(2)) y = plo(2)
         if (j.gt.dhi(2)) y = plo(2)+(dhi(2)-dlo(2)+1)*h(2)
         do i=ARG_L1(fine), ARG_H1(fine)
            x = (i + 0.5)*h(1)
            if (i.lt.dlo(1)) x = plo(1)
            if (i.gt.dhi(1)) x = plo(1)+(dhi(1)-dlo(1)+1)*h(1)
            do n=1,nc
               fine(i,j,n)=-(x-0.5d0)**2
c               fine(i,j,n)=-(y-0.5d0)**2
c               fine(i,j,n)=y**2
c               fine(i,j,n)=y
            enddo
         enddo
      enddo
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

      integer imid, jmid
      imid = 0.5*(ARG_L1(rhs)+ARG_H1(rhs))
      jmid = 0.5*(ARG_L2(rhs)+ARG_H2(rhs))

      rhs(:,:,:) = 0.d0
c      rhs(imid,jmid,:) = 1.d0
c      rhs(imid+1,jmid+1,:) = -1.d0
      
      end

      subroutine FORT_MAKEMU(
     $     mu,DIMS(mu),
     $     h, idir )
      implicit none
      integer idir
      integer DIMDEC(mu)
      REAL_T mu(DIMV(mu))
      REAL_T h(BL_SPACEDIM)

      integer i,j
      REAL_T x,y

      REAL_T aT, aL, u, v, vel, velInv
      namelist /fortin/ aT, aL, u, v
      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(9)
      vel = SQRT(u*u + v*v)
      if (vel .ne. 0.d0) then
         velInv = 1.d0 / vel
      else
         velInv = 0.d0
      endif

      if( idir .eq. 0 ) then
         do j=ARG_L2(mu),ARG_H2(mu)
            y = (j+0.5)*h(2)
            do i=ARG_L1(mu),ARG_H1(mu)
               x = (i)*h(1)
               mu(i,j) = aT*vel + (aL - aT)*u*u*velInv
            enddo
         enddo
      else if( idir .eq. 1 ) then
         do j=ARG_L2(mu),ARG_H2(mu)
            y = (j)*h(2)
            do i=ARG_L1(mu),ARG_H1(mu)
               x = (i+0.5)*h(1)
               mu(i,j) = aT*vel + (aL - aT)*v*v*velInv
            enddo
         enddo
      else
         write(6,*)'FORT_MAKEMU: bad idir'
         stop
      endif
      end

c     \begin{array}{cc}
c     \a_T v + (\a_L-\a_T)\dfrac{v_x^2}{v} & (\a_L-\a_T)\dfrac{v_x v_y}{v} \\
c     (\a_L-\a_T)\dfrac{v_y v_x}{v} & \a_T v + (\a_L-\a_T)\dfrac{v_y^2}{v}
c     \end{array}

      subroutine FORT_MAKEMU1(
     $     mu1,DIMS(mu1),
     $     h, idir )
      implicit none
      integer idir
      integer DIMDEC(mu1)
      REAL_T mu1(DIMV(mu1))
      REAL_T h(BL_SPACEDIM)

      integer i,j
      REAL_T x,y

      REAL_T aT, aL, u, v, vel, velInv
      namelist /fortin/ aT, aL, u, v
      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(9)
      vel = SQRT(u*u + v*v)
      if (vel .ne. 0.d0) then
         velInv = 1.d0 / vel
      else
         velInv = 0.d0
      endif

      if( idir .eq. 0 ) then
         do j=ARG_L2(mu1),ARG_H2(mu1)
            y = (j+0.5)*h(2)
            do i=ARG_L1(mu1),ARG_H1(mu1)
               x = (i)*h(1)
               mu1(i,j) = (aL - aT)*u*v*velInv
            enddo
         enddo
      else if( idir .eq. 1 ) then
         do j=ARG_L2(mu1),ARG_H2(mu1)
            y = (j)*h(2)
            do i=ARG_L1(mu1),ARG_H1(mu1)
               x = (i+0.5)*h(1)
               mu1(i,j) = (aL - aT)*u*v*velInv
            enddo
         enddo
      else
         write(6,*)'FORT_MAKEMU1: bad idir'
         stop
      endif
      
      return
      end
