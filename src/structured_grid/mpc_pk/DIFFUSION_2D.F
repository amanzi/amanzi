c
c $Id: DIFFUSION_2D.F,v 1.16 2011-08-30 19:00:04 gpau Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "DIFFUSION_F.H"
#include "ArrayLim.H"

#define SDIM 2

c ::
c :: ----------------------------------------------------------
c :: compute the time averaged viscous flux at the given edge
c ::  for variable transport coefficients
c ::
c :: NOTE:
c ::   These are only computed for a reflux operation, and require
c ::   that the ghost cells have been filled according to the 
c ::   diffusion operator (maybe this isn't such a good idea...)
c :: ----------------------------------------------------------
c ::

      subroutine FORT_VISCFLUX_VC (s_o,s_n,DIMS(s),lo,hi,
     $                             flux, DIMS(flux), area, DIMS(area),
     $                             bn, bnp1, DIMS(b),
     $                             dx,mult,theta,dir)

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(s)
      integer DIMDEC(flux)
      integer DIMDEC(area)
      integer DIMDEC(b)
      REAL_T  s_o(DIMV(s))
      REAL_T  s_n(DIMV(s))
      REAL_T  flux(DIMV(flux))
      REAL_T  area(DIMV(area))
      REAL_T  bn(DIMV(b))
      REAL_T  bnp1(DIMV(b))
      REAL_T  dx
      REAL_T  mult,theta
      integer dir

      integer i, j
      REAL_T  sx1, sx2, sy1, sy2
      REAL_T  onemintheta

      onemintheta = one-theta
      if (dir .EQ. 0) then
c     ::::: compute X fluxes
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            sx1 = (s_o(i,j) - s_o(i-1,j))*bn(i,j)
            sx2 = (s_n(i,j) - s_n(i-1,j))*bnp1(i,j)
            flux(i,j) = mult*(onemintheta*sx1 + theta*sx2)*area(i,j)/dx
         end do
         end do
      else
c     ::::: compute Y fluxes
         do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            sy1 = (s_o(i,j) - s_o(i,j-1))*bn(i,j)
            sy2 = (s_n(i,j) - s_n(i,j-1))*bnp1(i,j)
            flux(i,j) = mult*(onemintheta*sy1 + theta*sy2)*area(i,j)/dx
         end do
         end do
      end if
      end

c ::
c :: ----------------------------------------------------------
c :: compute the time averaged viscous flux at the given edge
c ::  for constant transport coefficients
c ::
c :: NOTE:
c ::   These are only computed for a reflux operation, and require
c ::   that the ghost cells have been filled according to the 
c ::   diffusion operator (maybe this isn't such a good idea...)
c :: ----------------------------------------------------------
c ::

      subroutine FORT_VISCFLUX_CC (s_o,s_n,DIMS(s),lo,hi,
     $                             flux, DIMS(flux), area, DIMS(area),
     $                             dx,mult,theta,dir)

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(s)
      integer DIMDEC(flux)
      integer DIMDEC(area)
      REAL_T  s_o(DIMV(s))
      REAL_T  s_n(DIMV(s))
      REAL_T  flux(DIMV(flux))
      REAL_T  area(DIMV(area))
      REAL_T  dx
      REAL_T  mult,theta
      integer dir

      integer i, j
      REAL_T  sx1, sx2, sy1, sy2
      REAL_T  onemintheta

      onemintheta = one-theta
      if (dir .EQ. 0) then
c     ::::: compute X fluxes
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            sx1 = s_o(i,j) - s_o(i-1,j)
            sx2 = s_n(i,j) - s_n(i-1,j)
            flux(i,j) = mult*(onemintheta*sx1 + theta*sx2)*area(i,j)/dx
         end do
         end do
      else
c     ::::: compute Y fluxes
         do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            sy1 = s_o(i,j) - s_o(i,j-1)
            sy2 = s_n(i,j) - s_n(i,j-1)
            flux(i,j) = mult*(onemintheta*sy1 + theta*sy2)*area(i,j)/dx
         end do
         end do
      end if
      end

      subroutine FORT_VISCSYNCFLUX (ssync,DIMS(ssync),lo,hi,
     $                              xflux,DIMS(xf),yflux,DIMS(yf),
     $                              xarea,DIMS(ax),yarea,DIMS(ay),dx,mult)

      integer lo(2), hi(2)
      integer DIMDEC(ssync)
      integer DIMDEC(xf)
      integer DIMDEC(yf)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T  ssync(DIMV(ssync))
      REAL_T  xflux(DIMV(xf))
      REAL_T  yflux(DIMV(yf))
      REAL_T  xarea(DIMV(ax))
      REAL_T  yarea(DIMV(ay))
      REAL_T  dx(2)
      REAL_T  mult

      REAL_T  sx
      REAL_T  sy
      integer i, j
c
c     ::::: compute X fluxes
c
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
	       sx = ssync(i,j) - ssync(i-1,j)
               xflux(i,j) = mult*sx*xarea(i,j)/dx(1)
            end do
         end do
c
c     ::::: compute Y fluxes
c
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
	       sy = ssync(i,j) - ssync(i,j-1)
               yflux(i,j) = mult*sy*yarea(i,j)/dx(2)
	    end do
         end do

      end

c :: ----------------------------------------------------------
c :: HOOPSRC
c ::             fab(i,j) = fab(i,j) - mu*u/(r(i)^2)
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab)  => index limits of fab
c ::  mu         => viscous coefficient
c :: ----------------------------------------------------------
c ::
       subroutine FORT_HOOPSRC (DIMS(grid), fab, DIMS(fab), u, DIMS(u), r, mu)

       integer DIMDEC(grid)
       integer DIMDEC(fab)
       integer DIMDEC(u)
       REAL_T  fab(DIMV(fab))
       REAL_T  u(DIMV(u))
       REAL_T  r(DIM1(grid))
       REAL_T  mu

       integer i, j

c      if (ARG_L1(u) .lt. ARG_L1(fab) .or. ARG_H1(u) .gt. ARG_H1(fab)) then
c         write(6,*) "FORT_HOOPSRC: bad index limits"
c         stop
c      end if

       do j = ARG_L2(grid), ARG_H2(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             fab(i,j) = fab(i,j) - mu*u(i,j)/(r(i)*r(i))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: HOOPRHS
c ::             rhs(i,j) = rhs(i,j) - (one-theta)*dt*u*mu*vol/(r(i)^2)
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab)  => index limits of fab
c ::  u          => array to be modified
c ::  DIMS(u)    => index limits of u
c ::  r          => 1-D r array (in first coordinate direction)
c ::  mu         => scalar viscosity
c ::  dt         => time step
c ::  vol        => volume array
c ::  DIMS(vol)  => index limits of vol
c ::  b          => (one-theta)*dt
c :: ----------------------------------------------------------
c ::
       subroutine FORT_HOOPRHS (fab, DIMS(fab), u, DIMS(u), r, b,
     &                          vol, DIMS(vol))
       integer DIMDEC(fab)
       integer DIMDEC(u)
       integer DIMDEC(vol)
       REAL_T  fab(DIMV(fab))
       REAL_T  u(DIMV(u))
       REAL_T  vol(DIMV(vol))
       REAL_T  r(DIM1(fab))
       REAL_T  b

       integer i, j

       do j = ARG_L2(fab), ARG_H2(fab)
          do i = ARG_L1(fab), ARG_H1(fab)
             fab(i,j) = fab(i,j) - b*vol(i,j)*u(i,j)/(r(i)*r(i))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: TENSOR_HOOPRHS
c ::             rhs(i,j) = rhs(i,j) - (1-theta)*dt*u*two*mu_cen*vol/(r(i)^2)
c ::                                                  ^^^ yes, that is correct
c ::                                                      for variable mu
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab)  => index limits of fab
c ::  u          => array to be modified
c ::  DIMS(u)    => index limits of u
c ::  r          => 1-D r array (in first coordinate direction)
c ::  mu         => scalar viscosity
c ::  dt         => time step
c ::  vol        => volume array
c ::  DIMS(vol)  => index limits of vol
c ::  b          => (1-theta)*dt
c :: ----------------------------------------------------------
c ::
       subroutine FORT_TENSOR_HOOPRHS (xvelcomp, fab, DIMS(fab), u, DIMS(u), 
     &                          r, b,
     &                          vol, DIMS(vol), betax, DIMS(betax),
     &                          betay, DIMS(betay))
       integer xvelcomp
       integer DIMDEC(fab)
       integer DIMDEC(u)
       integer DIMDEC(vol)
       integer DIMDEC(betax)
       integer DIMDEC(betay)
       REAL_T  fab(DIMV(fab),2)
       REAL_T  u(DIMV(u),2)
       REAL_T  vol(DIMV(vol))
       REAL_T  r(DIM1(fab))
       REAL_T  betax(DIMV(betax))
       REAL_T  betay(DIMV(betay))
       REAL_T  b

       REAL_T  betacen
       integer i, j

       do j = ARG_L2(fab), ARG_H2(fab)
          do i = ARG_L1(fab), ARG_H1(fab)
             betacen = fourth*(betax(i,j)+betax(i+1,j)+
     &                           betay(i,j)+betay(i,j+1))
             fab(i,j,xvelcomp) = fab(i,j,xvelcomp) - 
     &           b*two*betacen*vol(i,j)*u(i,j,xvelcomp)/(r(i)*r(i))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: TENSOR_HOOPSRC
c ::             fab(i,j) = fab(i,j) - two*mu*u/(r(i)^2)
c ::                                   ^^^ yes, that is correct
c ::                                       for variable mu
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab)  => index limits of fab
c ::  mu         => viscous coefficient
c :: ----------------------------------------------------------
c ::
       subroutine FORT_TENSOR_HOOPSRC (comp, DIMS(grid), fab, DIMS(fab), u, 
     &      DIMS(u), r, betax, DIMS(betax), betay, DIMS(betay))

       integer comp
       integer DIMDEC(grid)
       integer DIMDEC(fab)
       integer DIMDEC(u)
       integer DIMDEC(betax)
       integer DIMDEC(betay)
       REAL_T  fab(DIMV(fab),2)
       REAL_T  u(DIMV(u),2)
       REAL_T  r(DIM1(grid))
       REAL_T  betax(DIMV(betax))
       REAL_T  betay(DIMV(betay))

       integer i, j
       REAL_T  betacen

       do j = ARG_L2(grid), ARG_H2(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             betacen  = fourth*(betax(i,j)+betax(i+1,j)+
     &                           betay(i,j)+betay(i,j+1))
             fab(i,j,comp) = fab(i,j,comp) - two*betacen*u(i,j,comp)/(r(i)*r(i))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: SETALPHA
c ::             alpha(i,j) = vol(i,j)*rho(i,j)
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab) => index limits of fab
c ::  lo,hi     => index limits of box
c ::  vol       =>  volume array
c ::  DIMS(vol) => index limits of fab
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SETALPHA (fab, DIMS(fab), 
     &                           rho, DIMS(rho), 
     &                           lo, hi, 
     &                           vol, DIMS(vol))

       integer DIMDEC(fab)
       integer DIMDEC(rho)
       integer lo(SDIM), hi(SDIM)
       integer DIMDEC(vol)
       REAL_T  fab(DIMV(fab))
       REAL_T  rho(DIMV(rho))
       REAL_T  vol(DIMV(vol))

       integer i, j

       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             fab(i,j) = vol(i,j)*rho(i,j)
          end do
       end do

       end

      
      subroutine FORT_CCTOEDGE(lo, hi, DIMS(beta), beta,
     &   DIMS(betax), betax, DIMS(betay), betay)

      implicit none
c
c ... inputs
c
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(beta)
      REAL_T  beta(DIMV(beta))
      integer DIMDEC(betax)
      integer DIMDEC(betay)
c
c ... outputs
c
      REAL_T  betax(DIMV(betax))
      REAL_T  betay(DIMV(betay))
c
c ... local
c
      integer i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            betax(i,j) = .5D0*(beta(i-1,j)+beta(i,j))
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            betay(i,j) = .5D0*(beta(i,j-1)+beta(i,j))
         end do
      end do

      end

      subroutine FORT_DRHOG_RICHARD(rg,DIMS(rg),
     &     cx,DIMS(cx),cy,DIMS(cy),rho,
     &     lo,hi,dx,grav,mult)
c
c    ----------------------------------------------------------
c    RHOPRESS
c      compute gravity term in the pressure equation
c      \sum_a 1/\rho_a \nabla \rho_a^2 \lambda_a g 
c    ----------------------------------------------------------
c
      implicit none
      
      integer DIMDEC(rg)
      integer DIMDEC(cx)
      integer DIMDEC(cy)
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  rg(DIMV(rg))
      REAL_T  cx(DIMV(cx))
      REAL_T  cy(DIMV(cy))
      REAL_T  rho
      REAL_T  mult,grav
      
      integer i, j
c
c     gravity is nonzero only in the y direction
c
      do i = lo(1),hi(1)
         do j = lo(2),hi(2)
            rg(i,j) = rg(i,j) - mult*rho*grav*(cy(i,j+1) - cy(i,j))/dx(2)
         end do
      end do

      end     

      subroutine FORT_FRHOG_RICHARD(fx, DIMS(fx), cx, DIMS(cx),
     &     fy, DIMS(fy),cy,DIMS(cy),rho,
     &     lo,hi,dx,grav,mult)
c
c    ----------------------------------------------------------
c    RHOPRESS
c      compute gravity term in the pressure equation
c      \sum_a 1/\rho_a \nabla \rho_a^2 \lambda_a g 
c    ----------------------------------------------------------
c
      implicit none
      
      integer DIMDEC(fx)
      integer DIMDEC(fy)
      integer DIMDEC(cx)
      integer DIMDEC(cy)
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  fx(DIMV(fx))
      REAL_T  fy(DIMV(fy))
      REAL_T  cx(DIMV(cx))
      REAL_T  cy(DIMV(cy))
      REAL_T  rho
      REAL_T  mult,grav
      
      integer i, j
c
c     gravity is nonzero only in the y direction
c
      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)
            fy(i,j) = fy(i,j) + mult*rho*grav*cy(i,j) 
         end do
      end do

      end     


c ----------------------------------------------------------------------
c      
c     DFXDOTY: Dot product of X and Y
c     
c ----------------------------------------------------------------------
      subroutine FORT_DFXDOTY(
     &     pw,
     &     pp, DIMS(pp),np,idxp,
     &     ww, DIMS(ww),nw,idxw,
     &     lo, hi )
      implicit none
      integer np,nw,idxp,idxw
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(ww)
      REAL_T ww(DIMV(ww),nw)
      integer DIMDEC(pp)
      REAL_T pp(DIMV(pp),np)
      REAL_T pw
c
      integer i, j
c
      pw = 0.0D0
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            pw = pw + pp(i,j,idxp)*ww(i,j,idxw)
         end do
      end do
c
      end

c ----------------------------------------------------------------------
c      
c     DFSXAY: S = X + A*Y
c     
c ----------------------------------------------------------------------
      subroutine FORT_DFSXAY(
     &     ss, DIMS(ss),
     &     xx, DIMS(xx),
     &     a,
     &     yy, DIMS(yy),
     &     lo, hi, nc
     &     )
      implicit none
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(ss)
      REAL_T ss(DIMV(ss),nc)
      integer DIMDEC(xx)
      REAL_T xx(DIMV(xx),nc)
      integer DIMDEC(yy)
      REAL_T yy(DIMV(yy),nc)
      REAL_T a
c
      integer i, j, n
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ss(i,j,n) = xx(i,j,n) + a*yy(i,j,n)
            end do
         end do
      end do
c
      end


      subroutine FORT_ARNOLDI(soln,H,norm_rhs,nc)

      implicit none
      integer nc
      REAL_T  norm_rhs,soln(nc),H(nc,nc)
      
      integer IPIV(nc), INFO,i

            
c     set up right hand side: norm_rhs*e_1
      do i = 1,nc
         soln(i) = 0
      end do
      soln(1) = norm_rhs

      call DGETRF_PM( nc, nc, H, nc, IPIV, INFO)
      if (INFO .eq. 0) then
         call DGETRS_PM( 'N', nc, 1, H, nc, IPIV, soln, nc, INFO )
         if (INFO .ne. 0) then
            call bl_abort('FORT_ARNOLDI::DGETRS failed')
         end if
      else
         call bl_abort('FORT_ARNOLDI::DGETRF failed')
      end if

      end
