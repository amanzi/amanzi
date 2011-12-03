#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "MACOPERATOR_F.H"
#include "ArrayLim.H"

#define SDIM 2

c :: ----------------------------------------------------------
c :: MACCOEF
c ::             Compute the coefficents for MAC solv
c ::             Changes were made to reflect dependence on
c ::             pressure bc.
c ::
c :: INPUTS / OUTPUTS:
c ::  cx,cy       <=  edge coef arrays
c ::  DIMS(cx)     => index limits for cx
c ::  DIMS(cy)     => index limits for cy
c ::  lo,hi        => index limits for rhs
c ::  ax,ay        => edge based area arrays
c ::  DIMS(ax)     => index limits for ax
c ::  DIMS(ay)     => index limits for ay
c ::  coef         => cell centered density array
c ::  DIMS(coef)   => index limits for coef
c ::  dx           => cell size
c :: ----------------------------------------------------------
c ::
      subroutine FORT_MACCOEF (cx,DIMS(cx),cy,DIMS(cy),
     &                         lx,DIMS(lx),ly,DIMS(ly),
     &                         ax,DIMS(ax),ay,DIMS(ay),
     &                         lo,hi,domlo,domhi,dx,
     &                         lo_bc,hi_bc)
      integer DIMDEC(cx)
      integer DIMDEC(cy)
      integer DIMDEC(lx)
      integer DIMDEC(ly)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer lo(SDIM), hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  cx(DIMV(cx))
      REAL_T  cy(DIMV(cy))
      REAL_T  lx(DIMV(lx))
      REAL_T  ly(DIMV(ly))
      REAL_T  ax(DIMV(ax))
      REAL_T  ay(DIMV(ay))
      integer lo_bc(SDIM)
      integer hi_bc(SDIM)
      
      integer i, j
      REAL_T coef_avg
c     
c     ::::: make coef in X direction 
c     
      
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            coef_avg = lx(i,j)
            cx(i,j) = dx(1) * ax(i,j) * coef_avg
         end do
      end do
       
      if (lo(1) .eq. domlo(1)) then
         i = lo(1)
         if (lo_bc(1).eq.FOEXTRAP) then
            do j = lo(2), hi(2)
               cx(i,j) = zero
            end do
         end if
      end if
       
      if (hi(1) .eq. domhi(1)) then
         i = hi(1)+1
         if (hi_bc(1).eq.FOEXTRAP) then
            do j = lo(2), hi(2)
               cx(i,j) = zero
            end do
         end if
      end if
      
c
c      ::::: make coef in Y direction 
c
      do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            coef_avg = ly(i,j)
            cy(i,j) = dx(2) * ay(i,j) * coef_avg
         end do
      end do
      
      if (lo(2) .eq. domlo(2)) then
         j = lo(2)
         if (lo_bc(2).eq.FOEXTRAP) then
            do i = lo(1), hi(1)
               cy(i,j) = zero
            end do
         end if
      end if

      if (hi(2) .eq. domhi(2)) then
         j = hi(2)+1
         if (hi_bc(2).eq.FOEXTRAP) then
            do i = lo(1), hi(1)
               cy(i,j) = zero
            end do
         end if
      end if

      end

c :: ----------------------------------------------------------
c :: MACRHS
c ::             Compute the RHS for MAC solve
c ::             Changes were made to reflect dependence on
c ::             pressure bc.
c ::
c :: INPUTS / OUTPUTS:
c ::  ux,uy       <=  edge velocity arrays
c ::  DIMS(ux)     => index limits for ux
c ::  DIMS(uy)     => index limits for uy
c ::  lo,hi        => index limits for rhs
c ::  ax,ay        => edge based area arrays
c ::  DIMS(ax)     => index limits for ax
c ::  DIMS(ay)     => index limits for ay
c ::  vol          => cell centered volume array
c ::  vlo,vhi      => index limits of vol array
c ::  DIMS(vol)    => index limits for vol
c ::  rhs         <=> cell centered rhs array
c ::  DIMS(rhs)    => index limits for rhs
c ::  scale        => scale factor
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACRHS (ux,DIMS(ux),uy,DIMS(uy),
     &                         rx,DIMS(rx),ry,DIMS(ry),
     &                         ax,DIMS(ax),ay,DIMS(ay),
     &                         vol,DIMS(vol),rhs,DIMS(rhs),
     &                         lo,hi,domlo,domhi,lo_bc,hi_bc,scale)

       implicit none

       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(rx)
       integer DIMDEC(ry)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(vol)
       integer DIMDEC(rhs)
       integer lo(SDIM), hi(SDIM)
       integer domlo(SDIM), domhi(SDIM)
       REAL_T  scale
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  rx(DIMV(rx))
       REAL_T  ry(DIMV(ry))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  vol(DIMV(vol))
       REAL_T  rhs(DIMV(rhs))
       integer lo_bc(SDIM)
       integer hi_bc(SDIM)

       integer i, j
       REAL_T  divu,coef_l,coef_r,coef_b,coef_t
       REAL_T  rcoef_l,rcoef_r,rcoef_b,rcoef_t

c
c      ::::: rhs holds the divergence condition (possibly zero)
c
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rcoef_l = one
             rcoef_r = one
             rcoef_b = one
             rcoef_t = one
             coef_t  = one
             coef_b  = one
             coef_l  = one
             coef_r  = one 
             if (i.eq.domlo(1)) then 
                if  (lo_bc(1) .eq. FOEXTRAP) then
                   coef_l  = zero
                   rcoef_l = zero
                end if
             end if

             if (i.eq.domhi(1)) then 
                if  (hi_bc(1) .eq. FOEXTRAP) then
                   coef_r  = zero
                   rcoef_r = zero
                end if
             end if

             if (j.eq.domlo(2)) then 
                if  (lo_bc(2) .eq. FOEXTRAP) then
                   coef_b  = zero
                   rcoef_b = zero
                end if
             end if

            if (j.eq.domhi(2)) then 
                if  (hi_bc(2) .eq. FOEXTRAP) then
                   coef_t  = zero
                   rcoef_t = zero
                end if
             end if

             divu = ax(i+1,j)*(ux(i+1,j)*coef_r-rx(i+1,j)*rcoef_r)
     &            - ax(i  ,j)*(ux(i  ,j)*coef_l-rx(i  ,j)*rcoef_l)
     &            + ay(i,j+1)*(uy(i,j+1)*coef_t-ry(i,j+1)*rcoef_t)
     &            - ay(i,j  )*(uy(i,j  )*coef_b-ry(i,j  )*rcoef_b)
             rhs(i,j) = scale*(divu + vol(i,j)*rhs(i,j))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: MACUPDATE
c ::             Compute the update to velocity field 
c ::             Changes were made to reflect dependence on
c ::             pressure bc.
c ::
c :: INPUTS / OUTPUTS:
c ::  ux,uy       <=  edge based velocity arrays
c ::  DIMS(ux)     => index limits for ux
c ::  DIMS(uy)     => index limits for uy
c ::  phi          => soln from MAC project
c ::  DIMS(phi)    => index limits for phi
c ::  coef         => MAC coefficient at time N
c ::  DIMS(coef)   => index limits for coef
c ::  dx           => cell size
c ::  mult         => scalar multiplier
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACUPDATE(init,ux,DIMS(ux),uy,DIMS(uy),
     &     lx,DIMS(lx),ly,DIMS(ly),
     &     phi,DIMS(phi),
     &     lo,hi,domlo,domhi,dx,lo_bc,hi_bc,mult)

       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(lx)
       integer DIMDEC(ly)
       integer DIMDEC(phi)
       integer lo(SDIM), hi(SDIM)
       integer domlo(SDIM), domhi(SDIM)
       integer lo_bc(SDIM), hi_bc(SDIM)
       integer init
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  lx(DIMV(lx))
       REAL_T  ly(DIMV(ly))
       REAL_T  phi(DIMV(phi))
       REAL_T  dx(SDIM), mult

       integer i, j
       REAL_T  coef_avg, gp
c
c     set values to 0.0 if initializing
c
       if ( init .eq. 1 ) then
          ux = zero
          uy = zero
       end if
c
c     compute the x mac gradient
c

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             coef_avg = lx(i,j)
             gp = (phi(i,j)-phi(i-1,j))/dx(1)
             if ((i .eq.domlo(1)   .and. lo_bc(1) .eq. FOEXTRAP) .or.
     &           (i .eq.domhi(1)+1 .and. hi_bc(1) .eq. FOEXTRAP)) then
                ux(i,j) = ux(i,j) + coef_avg*mult*gp
             else
c                ux(i,j) = coef_avg * (ux(i,j) + mult * gp )
                ux(i,j) = ux(i,j) + coef_avg*mult*gp;
             end if
          end do
       end do
c
c      compute the y mac gradient
c
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             coef_avg = ly(i,j)
             gp = (phi(i,j)-phi(i,j-1))/dx(2)
             if ((j .eq.domlo(2)   .and. lo_bc(2) .eq. FOEXTRAP) .or.
     &           (j .eq.domhi(2)+1 .and. hi_bc(2) .eq. FOEXTRAP)) then
                uy(i,j) = uy(i,j) + coef_avg*mult*gp
             else
c                uy(i,j) = coef_avg * (uy(i,j) + mult * gp )
                uy(i,j) = uy(i,j) + coef_avg*mult*gp
             end if
          end do
       end do

       end

       subroutine FORT_MAC_SYNC( ux,DIMS(ux),uy,DIMS(uy),
     &                           rhodx,DIMS(rx),rhody,DIMS(ry), 
     &                           ax,DIMS(ax),ay,DIMS(ay),
     &                           lx,DIMS(lx),ly,DIMS(ly),
     &                           phi,DIMS(phi),lo,hi,dx,mult)

       implicit none

       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(rx)
       integer DIMDEC(ry)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(lx)
       integer DIMDEC(ly)
       integer DIMDEC(phi)
       integer lo(SDIM), hi(SDIM)
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  rhodx(DIMV(rx))
       REAL_T  rhody(DIMV(ry))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  lx(DIMV(lx))
       REAL_T  ly(DIMV(ly))
       REAL_T  phi(DIMV(phi))
       REAL_T  dx(SDIM), mult

       integer i, j
       REAL_T  coef_avg, gp

c
c     compute the x mac gradient
c

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             coef_avg = lx(i,j)
             gp = (phi(i,j)-phi(i-1,j))/dx(1)
             rhodx(i,j) = - rhodx(i,j)/ax(i,j)
             ux(i,j) = ux(i,j) + coef_avg*mult*gp + rhodx(i,j)
          end do
       end do
c
c      compute the y mac gradient
c
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             coef_avg = ly(i,j)
             gp = (phi(i,j)-phi(i,j-1))/dx(2)
             rhody(i,j) = -rhody(i,j)/ay(i,j)
             uy(i,j) = uy(i,j) + coef_avg*mult*gp + rhody(i,j)
             
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: MACNORML1
c :: ----------------------------------------------------------
       subroutine FORT_MACNORML1(normval,ux,DIMS(ux),uy,DIMS(uy),
     &                           lo,hi,dx)

       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer lo(SDIM), hi(SDIM)
       REAL_T  normval
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  dx(SDIM)

       integer i, j

       normval = 0
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             normval = normval + abs(ux(i,j))*dx(1)*dx(1) 
     &            + abs(uy(i,j))*dx(2)*dx(2)
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: MASYNCRHS
c ::        Modify the RHS for M SYNC solve
c ::
c :: INPUTS / OUTPUTS:
c ::  rhs         <=  right hand side array
c ::  DIMS(rhs)    => index limits for rhs
c ::  vol          => cell centered volume array
c ::  DIMS(vol)    => index limits of vol array
c ::  rhsscale     => const multiplier to rhs
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACSYNCRHS(rhs,DIMS(rhs),
     &                            ux,DIMS(ux),uy,DIMS(uy),
     &                            rx,DIMS(rx),ry,DIMS(ry),
     &                            ax,DIMS(ax),ay,DIMS(ay),
     &                            vol,DIMS(vol),
     &                            lo,hi,domlo,domhi,
     &                            lo_bc,hi_bc,rhsscale)
       integer DIMDEC(rhs)
       integer DIMDEC(vol)
       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(rx)
       integer DIMDEC(ry)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer lo(SDIM), hi(SDIM)
       integer domlo(SDIM), domhi(SDIM)
       integer lo_bc(SDIM), hi_bc(SDIM)
       REAL_T  rhsscale
       REAL_T  rhs(DIMV(rhs))
       REAL_T  vol(DIMV(vol))
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  rx(DIMV(rx))
       REAL_T  ry(DIMV(ry))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))

       integer i, j
c
c      ::::: multiply by volume since reflux step (which computed rhs)
c      ::::: divided by volume.
c
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhs(i,j) = rhsscale*( vol(i,j)*rhs(i,j))
          end do
       end do
       
       end

c :: -------------------------------------------------------------
c :: RHOGBC: Define the hydrostatic values of phi along boundaries
c :: -------------------------------------------------------------

      subroutine FORT_RHOGBC(rhog,DIMS(rhog),phi,DIMS(phi),
     &     face,gravity,dx,domlo,domhi,
     &     lo_bc,hi_bc,press_lo,press_hi)
c
c    Compute the contribution of gravity to the boundary conditions
c      for phi at outflow faces only.
c
      implicit none

      integer DIMDEC(rhog)
      integer DIMDEC(phi)
      integer face
      integer domlo(SDIM)
      integer domhi(SDIM)
      integer lo_bc(SDIM)
      integer hi_bc(SDIM)
      REAL_T  rhog(DIMV(rhog))
      REAL_T  phi(DIMV(phi))
      REAL_T  dx(SDIM)
      REAL_T  gravity
      REAL_T  press_lo(SDIM)
      REAL_T  press_hi(SDIM)

c     Local variables
      integer i,j
      REAL_T rhog_step, rhoG_ext, rhosum

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3


      if (face .eq. YLO) then
         phi(:,ARG_L2(phi)) = press_lo(2)
      elseif (face .eq. YHI) then
         phi(:,ARG_H2(phi)) = press_hi(2)
      end if

c     phi is only one node wide in i-direction.
      i = ARG_L1(phi)
      if (face .eq. XLO) then
         
         j = ARG_L2(phi)
         rhoG_ext = rhog(i+1,j)
         rhog_step = rhoG_ext*half*dx(2) 
         phi(i,j)  = phi(i,j) + rhog_step  
         do j = ARG_L2(phi)+1,ARG_H2(phi)
            rhoG_ext  = rhog(i+1,j)
            rhog_step = rhog_step + rhoG_ext * dx(2)
            phi(i,j)  = phi(i,j) + rhog_step 
         end do
         rhosum = phi(i,ARG_H2(phi))+half*rhog(i+1,ARG_H2(phi)+1)
         do j = ARG_L2(phi),ARG_H2(phi)
            phi(i,j) = phi(i,j) - rhosum + press_lo(1)
         end do

      else if (face .eq. XHI) then
         
         j = ARG_L2(phi)
         rhoG_ext = rhog(i-1,j)
         rhog_step = rhoG_ext*half*dx(2) 
         phi(i,j)  = phi(i,j) + rhog_step  
         
         do j = ARG_L2(phi)+1,ARG_H2(phi)
            rhoG_ext = rhog(i-1,j)
            rhog_step = rhog_step + rhoG_ext * dx(2)
            phi(i,j) = phi(i,j) + rhog_step 
         end do
         rhosum = phi(i,ARG_H2(phi))+half*rhog(i-1,ARG_H2(phi)+1)
         do j = ARG_L2(phi),ARG_H2(phi)
            phi(i,j) = phi(i,j) - rhosum + press_hi(1)
         end do
         
      endif

#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end
