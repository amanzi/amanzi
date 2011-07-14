c
c $Id: DIFFUSION_3D.F,v 1.10 2011-07-13 00:16:16 gpau Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES_PM.H"
#include "DIFFUSION_F.H"
#include "ArrayLim.H"

#define SDIM 3

c ::
c :: ----------------------------------------------------------
c :: compute the time averaged viscous flux at the given edge
c ::  for variable transport coefficients
c ::
c :: NOTE:
c ::   These are only computed for a reflux operation, and require
c ::   that the ghost cells have been filled according to the 
c ::   diffusion operator (maybe this isnt such a good idea...)
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

      integer i, j, k
      REAL_T  sx1, sx2, sy1, sy2, sz1, sz2
      REAL_T  onemintheta

      onemintheta = one-theta
      if (dir .EQ. 0) then
c     ::::: compute X fluxes
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
	    sx1 = (s_o(i,j,k) - s_o(i-1,j,k))*bn(i,j,k)
	    sx2 = (s_n(i,j,k) - s_n(i-1,j,k))*bnp1(i,j,k)
	    flux(i,j,k) = mult*(onemintheta*sx1 + theta*sx2)*area(i,j,k)/dx
         end do
         end do
         end do
      else if (dir .EQ. 1) then
c     ::::: compute Y fluxes
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            sy1 = (s_o(i,j,k) - s_o(i,j-1,k))*bn(i,j,k)
            sy2 = (s_n(i,j,k) - s_n(i,j-1,k))*bnp1(i,j,k)
            flux(i,j,k) = mult*(onemintheta*sy1 + theta*sy2)*area(i,j,k)/dx
	 end do
         end do
         end do
      else if (dir .EQ. 2) then
c     ::::: compute Z fluxes
         do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            sz1 = (s_o(i,j,k) - s_o(i,j,k-1))*bn(i,j,k)
            sz2 = (s_n(i,j,k) - s_n(i,j,k-1))*bnp1(i,j,k)
            flux(i,j,k) = mult*(onemintheta*sz1 + theta*sz2)*area(i,j,k)/dx
	 end do
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
c ::   diffusion operator (maybe this isnt such a good idea...)
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

      integer i, j, k
      REAL_T  sx1, sx2, sy1, sy2, sz1, sz2
      REAL_T  onemintheta

      onemintheta = one-theta
      if (dir .EQ. 0) then
c     ::::: compute X fluxes
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            sx1 = s_o(i,j,k) - s_o(i-1,j,k)
            sx2 = s_n(i,j,k) - s_n(i-1,j,k)
            flux(i,j,k) = mult*(onemintheta*sx1 + theta*sx2)*area(i,j,k)/dx
         end do
         end do
         end do
      else if (dir .EQ. 1) then
c     ::::: compute Y fluxes
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            sy1 = s_o(i,j,k) - s_o(i,j-1,k)
            sy2 = s_n(i,j,k) - s_n(i,j-1,k)
            flux(i,j,k) = mult*(onemintheta*sy1 + theta*sy2)*area(i,j,k)/dx
         end do
         end do
         end do
      else
c     ::::: compute Z fluxes
         do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            sz1 = s_o(i,j,k) - s_o(i,j,k-1)
            sz2 = s_n(i,j,k) - s_n(i,j,k-1)
            flux(i,j,k) = mult*(onemintheta*sz1 + theta*sz2)*area(i,j,k)/dx
	 end do
         end do
         end do
      end if
      end

      subroutine FORT_VISCSYNCFLUX (ssync,DIMS(ssync),lo,hi,
     $                              xflux,DIMS(xf),yflux,DIMS(yf),
     $                              zflux,DIMS(zf),
     $                              xarea,DIMS(ax),yarea,DIMS(ay),
     $                              zarea,DIMS(az),dx,mult)

      integer lo(3), hi(3)
      integer DIMDEC(ssync)
      integer DIMDEC(xf)
      integer DIMDEC(yf)
      integer DIMDEC(zf)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      REAL_T  ssync(DIMV(ssync))
      REAL_T  xflux(DIMV(xf))
      REAL_T  yflux(DIMV(yf))
      REAL_T  zflux(DIMV(zf))
      REAL_T  xarea(DIMV(ax))
      REAL_T  yarea(DIMV(ay))
      REAL_T  zarea(DIMV(az))
      REAL_T  dx(3)
      REAL_T  mult

      integer i, j, k
      REAL_T  sx
      REAL_T  sy
      REAL_T  sz
c
c     ::::: compute X fluxes
c
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               sx = ssync(i,j,k) - ssync(i-1,j,k)
               xflux(i,j,k) = mult*sx*xarea(i,j,k)/dx(1)
            end do
         end do
      end do
c
c     ::::: compute Y fluxes
c
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
               sy = ssync(i,j,k) - ssync(i,j-1,k)
               yflux(i,j,k) = mult*sy*yarea(i,j,k)/dx(2)
            end do
         end do
      end do
c
c     ::::: compute Z fluxes
c
      do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sz = ssync(i,j,k) - ssync(i,j,k-1)
               zflux(i,j,k) = mult*sz*zarea(i,j,k)/dx(3)
            end do
         end do
      end do
      
      end

c :: ----------------------------------------------------------
c :: SETALPHA
c ::             alpha(i,j,k) = vol*(1+b/(r(i)^2)) / density
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab) => index limits of fab
c ::  lo,hi     => index limits of box
c ::  r         =>  1-d array of radius
c ::  DIMS(r)   => index limits of r
c ::  b         =>  viscous coefficient
c ::  vol       =>  volume array
c ::  DIMS(vol) => index limits of fab
c ::  denfab    => array of density at time n+1/2
c ::  DIMS(den) => index limits of fab
c ::  usehoop   => do we add hoop stress?   NOT IN 3-D
c ::  useden    => do we divide by density? (only if velocity component)
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SETALPHA (fab, DIMS(fab), 
     &                           rho, DIMS(rho),  
     &                           lo, hi,
     &                           vol, DIMS(vol))

       integer DIMDEC(fab)
       integer DIMDEC(rho)
       integer DIMDEC(vol)
       integer lo(SDIM), hi(SDIM)
       REAL_T  fab(DIMV(fab))
       REAL_T  rho(DIMV(rho))
       REAL_T  vol(DIMV(vol))

       integer i, j, k

!$omp parallel do private(i,j,k)
       do k = lo(3)-1, hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                fab(i,j,k) = vol(i,j,k)*rho(i,j,k)
             end do
          end do
       end do
!$omp end parallel do

       end

      subroutine FORT_CCTOEDGE(lo, hi, DIMS(beta), beta,
     &   DIMS(betax), betax, DIMS(betay), betay, DIMS(betaz), betaz)

      implicit none
c
c ... inputs
c
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(beta)
      REAL_T  beta(DIMV(beta))
      integer DIMDEC(betax)
      integer DIMDEC(betay)
      integer DIMDEC(betaz)
c
c ... outputs
c
      REAL_T  betax(DIMV(betax))
      REAL_T  betay(DIMV(betay))
      REAL_T  betaz(DIMV(betaz))
c
c ... local
c
      integer i,j,k
 
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               betax(i,j,k) = .5D0*(beta(i-1,j,k)+beta(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               betay(i,j,k) = .5D0*(beta(i,j-1,k)+beta(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               betaz(i,j,k) = .5D0*(beta(i,j,k-1)+beta(i,j,k))
            end do
         end do
      end do

      end




      subroutine FORT_DRHOG_RICHARD(rg,DIMS(rg),
     &     cx,DIMS(cx),cy,DIMS(cy),cz,DIMS(cz),rho,
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
      integer DIMDEC(cz)
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  rg(DIMV(rg))
      REAL_T  cx(DIMV(cx))
      REAL_T  cy(DIMV(cy))
      REAL_T  cz(DIMV(cz))
      REAL_T  rho
      REAL_T  mult,grav
      
      integer i, j, k
c
c     gravity is nonzero only in the z direction
c
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rg(i,j,k) = rg(i,j,k) - mult*rho*grav*(cz(i,j,k+1) - cz(i,j,k))/dx(3)
            end do
         end do
      end do

      end     

      subroutine FORT_FRHOG_RICHARD(fx, DIMS(fx), cx, DIMS(cx),
     &     fy, DIMS(fy),cy,DIMS(cy), fz, DIMS(fz), cz, DIMS(cz), rho,
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
      integer DIMDEC(fz)
      integer DIMDEC(cx)
      integer DIMDEC(cy)
      integer DIMDEC(cz)
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  fx(DIMV(fx))
      REAL_T  fy(DIMV(fy))
      REAL_T  fz(DIMV(fz))
      REAL_T  cx(DIMV(cx))
      REAL_T  cy(DIMV(cy))
      REAL_T  cz(DIMV(cz))
      REAL_T  rho
      REAL_T  mult,grav
      
      integer i, j, k
c
c     gravity is nonzero only in the z direction
c
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               fz(i,j,k) = fz(i,j,k) - mult*rho*grav*cz(i,j,k) 
            end do
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
      integer i, j, k
c
      pw = 0.0D0
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               pw = pw + pp(i,j,k,idxp)*ww(i,j,k,idxw)
            end do
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
      integer i, j, k, n
c
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ss(i,j,k,n) = xx(i,j,k,n) + a*yy(i,j,k,n)
               end do
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

      call DGETRF( nc, nc, H, nc, IPIV, INFO)
      if (INFO .eq. 0) then
         call DGETRS( 'N', nc, 1, H, nc, IPIV, soln, nc, INFO )
         if (INFO .ne. 0) then
            call bl_abort('FORT_ARNOLDI::DGETRS failed')
         end if
      else
         call bl_abort('FORT_ARNOLDI::DGETRF failed')
      end if

      end
