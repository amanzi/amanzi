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

#define SDIM 3

c :: ----------------------------------------------------------
c :: MACCOEF
c ::             Compute the coefficents for MAC solve
c ::
c :: INPUTS / OUTPUTS:
c ::  cx,cy,cz    <=  edge coef arrays
c ::  DIMS(cx)     => index limits for cx
c ::  DIMS(cy)     => index limits for cy
c ::  DIMS(cz)     => index limits for cz
c ::  ax,ay,az     => edge based area arrays
c ::  DIMS(ax)     => index limits for ax
c ::  DIMS(ay)     => index limits for ay
c ::  DIMS(az)     => index limits for az
c ::  coef         => cell centered density array
c ::  DIMS(coef)   => index limits of coef array
c ::  lo,hi        => index limits for rhs
c ::  dx           => cell size
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACCOEF (cx,DIMS(cx),cy,DIMS(cy),cz,DIMS(cz),
     &                          lx,DIMS(lx),ly,DIMS(ly),lz,DIMS(lz),
     &                          ax,DIMS(ax),ay,DIMS(ay),az,DIMS(az),
     &                          lo,hi,domlo,domhi,dx,
     &                          lo_bc,hi_bc)

       implicit none

       integer DIMDEC(cx)
       integer DIMDEC(cy)
       integer DIMDEC(cz)
       integer DIMDEC(lx)
       integer DIMDEC(ly)
       integer DIMDEC(lz)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(az)
       integer lo(SDIM), hi(SDIM)
       integer domlo(SDIM), domhi(SDIM)
       REAL_T  dx(SDIM)
       REAL_T  cx(DIMV(cx))
       REAL_T  cy(DIMV(cy))
       REAL_T  cz(DIMV(cz))
       REAL_T  lx(DIMV(lx))
       REAL_T  ly(DIMV(ly))
       REAL_T  lz(DIMV(lz))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  az(DIMV(az))
       integer lo_bc(SDIM)
       integer hi_bc(SDIM)

       integer i, j, k
c
c      ::::: make coef in X direction
c
!$omp parallel do private(i,j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1
                cx(i,j,k) = dx(1) * ax(i,j,k) * lx(i,j,k)
             end do
          end do
       end do
!$omp end parallel do

       if (lo(1) .eq. domlo(1)) then 
         i = lo(1)
         if (lo_bc(1).eq.FOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  cx(i,j,k) = zero
               end do
            end do
         end if
       end if
 
       if (hi(1) .eq. domhi(1)) then 
         i = hi(1)+1
         if (hi_bc(1).eq.FOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  cx(i,j,k) = zero
               end do
            end do
         end if
       end if
c
c      ::::: make coef in Y direction
c
!$omp parallel do private(i,j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)
                cy(i,j,k) = dx(2) * ay(i,j,k) * ly(i,j,k)
             end do
          end do
       end do
!$omp end parallel do

       if (lo(2) .eq. domlo(2)) then 
         j = lo(2)
         if (lo_bc(2) .eq. FOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  cy(i,j,k) = zero
               end do
            end do
         end if
       end if
 
       if (hi(2) .eq. domhi(2)) then 
         j = hi(2)+1
         if (hi_bc(2) .eq. FOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  cy(i,j,k) = zero
               end do
            end do
         end if
       end if
c
c      ::::: make coef in Z direction
c
!$omp parallel do private(i,j,k)
       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                cz(i,j,k) = dx(3) * az(i,j,k) * lz(i,j,k)
             end do
          end do
       end do
!$omp end parallel do

       if (lo(3) .eq. domlo(3)) then
         k = lo(3)
         if  (lo_bc(3).eq.FOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  cz(i,j,k) = zero
               end do
            end do
         end if
       end if
 
       if (hi(3) .eq. domhi(3)) then
         k = hi(3)+1
         if (hi_bc(3).eq.FOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  cz(i,j,k) = zero
               end do
            end do
         end if
       end if

       end

c :: ----------------------------------------------------------
c :: Averages: harmonic or otherwise
c :: ----------------------------------------------------------
c ::
      REAL_T function cave(cf1,cf2,flag)
      REAL_T cf1,cf2
      logical flag

      if (flag) then
         cave = one/cf1+one/cf2
         cave = two/cave
      else
         cave = half * (cf1 + cf2)
      end if

      end 

c :: ----------------------------------------------------------
c :: MACRHS
c ::             Compute the RHS for MAC solve
c ::
c :: INPUTS / OUTPUTS:
c ::  ux,uy,uz    <=  edge velocity arrays
c ::  DIMS(ux)     => index limits for ux
c ::  DIMS(uy)     => index limits for uy
c ::  DIMS(uz)     => index limits for uz
c ::  ax,ay,az     => edge based area arrays
c ::  DIMS(ax)     => index limits for ax
c ::  DIMS(ay)     => index limits for ay
c ::  DIMS(az)     => index limits for az
c ::  vol          => cell centered volume array
c ::  DIMS(vol)    => index limits of vol array
c ::  rhs         <=> cell centered rhs array
c ::  DIMS(rhs)    => index limits of rhs array
c ::  lo,hi        => index limits for rhs
c ::  scale        => scale factor
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACRHS (ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz),
     &                         rx,DIMS(rx),ry,DIMS(ry),rz,DIMS(rz),
     &                         ax,DIMS(ax),ay,DIMS(ay),az,DIMS(az),
     &                         vol,DIMS(vol),rhs,DIMS(rhs),
     &                         lo,hi,domlo,domhi,lo_bc,hi_bc,scale)

       implicit none

       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(uz)
       integer DIMDEC(rx)
       integer DIMDEC(ry)
       integer DIMDEC(rz)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(az)
       integer DIMDEC(vol)
       integer DIMDEC(rhs)
       integer lo(SDIM), hi(SDIM)
       integer domlo(SDIM), domhi(SDIM)
       REAL_T  scale
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  uz(DIMV(uz))
       REAL_T  rx(DIMV(rx))
       REAL_T  ry(DIMV(ry))
       REAL_T  rz(DIMV(rz))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  az(DIMV(az))
       REAL_T  vol(DIMV(vol))
       REAL_T  rhs(DIMV(rhs))
       integer lo_bc(SDIM)
       integer hi_bc(SDIM)

       integer i, j, k
       REAL_T  divu
       REAL_T  coef_l ,coef_r ,coef_b ,coef_t ,coef_d ,coef_u
       REAL_T  rcoef_l,rcoef_r,rcoef_b,rcoef_t,rcoef_d,rcoef_u
c
c      ::::: rhs holds the divergence condition (possibly zero)
c
!$omp parallel do private(divu,coef_l,coef_r,coef_b,coef_t,coef_d)
!$omp& private(coef_u,rcoef_l,rcoef_r,rcoef_b,rcoef_t,rcoef_d,rcoef_u)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rcoef_l = one
                rcoef_r = one
                rcoef_t = one
                rcoef_b = one
                rcoef_u = one
                rcoef_d = one
                coef_l  = one
                coef_r  = one
                coef_t  = one
                coef_b  = one
                coef_u  = one
                coef_d  = one
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

                if (k.eq.domlo(3)) then 
                   if  (lo_bc(3) .eq. FOEXTRAP) then
                      coef_d  = zero
                      rcoef_d = zero
                   end if
                end if

                if (k.eq.domhi(3)) then 
                  if  (hi_bc(3) .eq. FOEXTRAP) then
                      coef_u  = zero
                      rcoef_u = zero
                   end if
                end if

                divu = ax(i+1,j,k)*(ux(i+1,j,k)*coef_r-rx(i+1,j,k)*rcoef_r) 
     &               - ax(i  ,j,k)*(ux(i  ,j,k)*coef_l-rx(i  ,j,k)*rcoef_l)
     &               + ay(i,j+1,k)*(uy(i,j+1,k)*coef_t-ry(i,j+1,k)*rcoef_t) 
     &               - ay(i,j  ,k)*(uy(i,j  ,k)*coef_b-ry(i,j  ,k)*rcoef_b)
     &               + az(i,j,k+1)*(uz(i,j,k+1)*coef_u-rz(i,j,k+1)*rcoef_u) 
     &               - az(i,j,k  )*(uz(i,j,k  )*coef_d-rz(i,j,k  )*rcoef_d)
                rhs(i,j,k) = scale*(divu + vol(i,j,k)*rhs(i,j,k))
             end do
          end do
       end do
!$omp end parallel do

       end

c :: ----------------------------------------------------------
c :: MACUPDATE
c ::             Compute the update to velocity field to
c ::             make it divergence free
c ::
c :: INPUTS / OUTPUTS:
c ::  ux,uy,uz    <=  edge based velocity arrays
c ::  DIMS(ux)     => index limits for ux
c ::  DIMS(uy)     => index limits for uy
c ::  DIMS(uz)     => index limits for uz
c ::  phi          => soln from MAC project
c ::  DIMS(phi)    => index limits for phi
c ::  coef         => density at time N
c ::  DIMS(coef)   => index limits for coef
c ::  dx           => cell size
c ::  mult         => scalar multiplier
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACUPDATE(
     &     init,
     &     ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz),
     &     lx,DIMS(lx),ly,DIMS(ly),lz,DIMS(lz),
     &     phi,DIMS(phi),
     &     lo,hi,domlo,domhi,dx,lo_bc,hi_bc,mult)


       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(uz)
       integer DIMDEC(lx)
       integer DIMDEC(ly)
       integer DIMDEC(lz)
       integer DIMDEC(phi)
       integer lo(SDIM), hi(SDIM)
       integer domlo(SDIM), domhi(SDIM)
       REAL_T  dx(SDIM), mult
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  uz(DIMV(uz))
       REAL_T  lx(DIMV(lx))
       REAL_T  ly(DIMV(ly))
       REAL_T  lz(DIMV(lz))
       REAL_T  phi(DIMV(phi))
       integer lo_bc(SDIM), hi_bc(SDIM)
       integer init
       
       integer i, j, k
       REAL_T  coef_avg, gp, xdx1, xdx2, xdx3
c
c     set gradient to zero if initializing
c
       if ( init .eq. 1 ) then
          ux = zero
          uy = zero
          uz = zero
       end if

       xdx1 = one/dx(1)
       xdx2 = one/dx(2)
       xdx3 = one/dx(3)

!$omp parallel private(i,j,k,coef_avg,gp)
c
c      compute x MAC gradient
c
!$omp do
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                coef_avg = lx(i,j,k)
                gp = (phi(i,j,k)-phi(i-1,j,k))*xdx1
                if ((i .eq.domlo(1) .and. lo_bc(1) .eq. FOEXTRAP) .or.
     &            (i .eq.domhi(1)+1 .and. hi_bc(1) .eq. FOEXTRAP)) then
                   ux(i,j,k) = ux(i,j,k) + coef_avg*mult*gp
                else
                   ux(i,j,k) = ux(i,j,k) + coef_avg*mult*gp
                end if
             end do
          end do
       end do
!$omp end do nowait
c     
c     compute y mac gradient
c
!$omp do
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                coef_avg = ly(i,j,k)
                gp = (phi(i,j,k)-phi(i,j-1,k))*xdx2
                if ((j .eq.domlo(2) .and. lo_bc(2) .eq. FOEXTRAP) .or.
     &              (j .eq.domhi(2)+1 .and. hi_bc(2) .eq. FOEXTRAP)) then
                   uy(i,j,k) = uy(i,j,k) + coef_avg*mult*gp
                else
                   uy(i,j,k) = uy(i,j,k) + coef_avg*mult*gp
                end if
             end do
          end do
       end do
!$omp end do nowait
c       
c     compute z mac gradient
c
!$omp do
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                coef_avg = lz(i,j,k)
                gp = (phi(i,j,k)-phi(i,j,k-1))*xdx3
                if ((k .eq.domlo(3) .and. lo_bc(3) .eq. FOEXTRAP) .or.
     &              (k .eq.domhi(3)+1 .and. hi_bc(3) .eq. FOEXTRAP)) then
                   uz(i,j,k) = uz(i,j,k) + coef_avg*mult*gp
                else
                   uz(i,j,k) = uz(i,j,k) + coef_avg*mult*gp
                end if
             end do
          end do
       end do
!$omp end do
!$omp end parallel

       end

       subroutine FORT_MAC_SYNC( ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz),
     &                           rhodx,DIMS(rx),rhody,DIMS(ry),rhodz,DIMS(rz),
     &                           ax,DIMS(ax),ay,DIMS(ay),az,DIMS(az),
     &                           lx,DIMS(lx),ly,DIMS(ly),lz,DIMS(lz),
     &                           phi,DIMS(phi),lo,hi,dx,mult)

       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(uz)
       integer DIMDEC(rx)
       integer DIMDEC(ry)
       integer DIMDEC(rz)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(az)
       integer DIMDEC(lx)
       integer DIMDEC(ly)
       integer DIMDEC(lz)
       integer DIMDEC(phi)
       integer lo(SDIM), hi(SDIM)
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  uz(DIMV(uz))
       REAL_T  rhodx(DIMV(rx))
       REAL_T  rhody(DIMV(ry))
       REAL_T  rhodz(DIMV(rz))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  az(DIMV(az))
       REAL_T  lx(DIMV(lx))
       REAL_T  ly(DIMV(ly))
       REAL_T  lz(DIMV(lz))
       REAL_T  phi(DIMV(phi))
       REAL_T  dx(SDIM), mult

       integer i, j,k
       REAL_T  coef_avg, gp, xdx1, xdx2, xdx3

       xdx1 = one/dx(1)
       xdx2 = one/dx(2)
       xdx3 = one/dx(3)

!$omp parallel private(i,j,k,coef_avg,gp)
c
c     compute the x mac gradient
c
!$omp do
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                coef_avg = lx(i,j,k)
                gp = (phi(i,j,k)-phi(i-1,j,k))*xdx1
                rhodx(i,j,k) = - rhodx(i,j,k)/ax(i,j,k)
                ux(i,j,k) = ux(i,j,k) + coef_avg*mult*gp + rhodx(i,j,k)
             end do
          end do
       end do
!$omp end do nowait
c
c      compute the y mac gradient
c
!$omp do
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                coef_avg = ly(i,j,k)
                gp = (phi(i,j,k)-phi(i,j-1,k))*xdx2
                rhody(i,j,k) = -rhody(i,j,k)/ay(i,j,k)
                uy(i,j,k) = uy(i,j,k) + coef_avg*mult*gp  + rhody(i,j,k)
             end do
          end do
       end do
!$omp end do nowait
c
c      compute the z mac gradient
c
!$omp do
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                coef_avg = lz(i,j,k)
                gp = (phi(i,j,k)-phi(i,j,k-1))*xdx3
                rhodz(i,j,k) = -rhodz(i,j,k)/az(i,j,k)
                uz(i,j,k) = uz(i,j,k) + coef_avg*mult*gp + rhodz(i,j,k)
             end do
          end do
       end do
!$omp end do
!$omp end parallel

       end

c :: ----------------------------------------------------------
c :: MACNORML1
c :: ----------------------------------------------------------
c :: 
       subroutine FORT_MACNORML1(normval,ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz),
     &                           lo,hi,dx)

       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(uz)
       integer lo(SDIM), hi(SDIM)
       REAL_T  normval
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  uz(DIMV(uz))
       REAL_T  dx(SDIM)

       integer i, j,k

       normval = 0
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                normval = normval + abs(ux(i,j,k))*dx(1)*dx(1) + 
     &               abs(uy(i,j,k))*dx(2)*dx(2) + abs(uz(i,j,k))*dx(3)*dx(3)
             end do
          end do
       end do

       end
c :: ----------------------------------------------------------
c :: MACSYNCRHS
c ::        Modify the RHS for MAC SYNC solve
c ::
c :: INPUTS / OUTPUTS:
c ::  rhs         <=  right hand side array
c ::  lo,hi        => index limits for rhs
c ::  vol          => cell centered volume array
c ::  vlo,vhi      => index limits of vol array
c ::  rhsscale     => const multiplier to rhs
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACSYNCRHS(rhs,DIMS(rhs),
     &                            ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz),
     &                            rx,DIMS(rx),ry,DIMS(ry),rz,DIMS(rz),
     &                            ax,DIMS(ax),ay,DIMS(ay),az,DIMS(az),
     &                            vol,DIMS(vol),
     &                            lo,hi,domlo,domhi,
     &                            lo_bc,hi_bc,rhsscale)
       integer DIMDEC(rhs)
       integer DIMDEC(vol)
       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(uz)
       integer DIMDEC(rx)
       integer DIMDEC(ry)
       integer DIMDEC(rz)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(az)
       integer lo(SDIM), hi(SDIM)
       integer domlo(SDIM), domhi(SDIM)
       integer lo_bc(SDIM), hi_bc(SDIM)
       REAL_T  rhsscale
       REAL_T  rhs(DIMV(rhs))
       REAL_T  vol(DIMV(vol))
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  uz(DIMV(uz))
       REAL_T  rx(DIMV(rx))
       REAL_T  ry(DIMV(ry))
       REAL_T  rz(DIMV(rz))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  az(DIMV(az))

       integer i, j, k
c
c      ::::: multiply by volume since reflux step (which computed rhs)
c      ::::: divided by volume.
c
!$omp parallel do private(i,j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rhs(i,j,k) = rhsscale*vol(i,j,k)*rhs(i,j,k)
             end do
          end do
       end do
!$omp end parallel do

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
      integer i,j,k
      REAL_T rhog_step
      REAL_T rhoG_ext, rhosum

#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5

      if (face .eq. ZLO) then
         phi(:,:,ARG_L3(phi)) = press_lo(3)
      else if (face .eq. ZHI) then
         phi(:,:,ARG_H3(phi)) = press_hi(3)     
      end if 

c     phi is only one node wide in i-direction.
      i = ARG_L1(phi)
      
      if (face .eq. XLO ) then
         do j = ARG_L2(phi),ARG_H2(phi)
            k = ARG_L3(phi)
            rhoG_ext = half * rhog(i+1,j,k)  
            rhog_step = rhoG_ext * dx(3) 
            phi(i,j,k)  = phi(i,j,k) + rhog_step  

            do k = ARG_L3(phi)+1,ARG_H3(phi)
               rhoG_ext   = rhog(i+1,j,k) 
               rhog_step  = rhog_step + rhoG_ext * dx(3)
               phi(i,j,k) = phi(i,j,k) + rhog_step 
            end do

            rhosum = phi(i,j,ARG_H3(phi))+half*rhog(i+1,j,ARG_H3(phi)+1)
            do k = ARG_L3(phi),ARG_H3(phi)
               phi(i,j,k) = phi(i,j,k) - rhosum + press_lo(1);
            end do

         end do

      elseif ( face .eq. XHI) then
         do j = ARG_L2(phi),ARG_H2(phi)
            k = ARG_L3(phi)
            rhoG_ext   = half * rhog(i-1,j,k) 
            rhog_step  = rhoG_ext * dx(3) 
            phi(i,j,k) = phi(i,j,k) + rhog_step  

            do k = ARG_L3(phi)+1,ARG_H3(phi)
               rhoG_ext   = rhog(i-1,j,k) 
               rhog_step  = rhog_step + rhoG_ext * dx(3)
               phi(i,j,k) = phi(i,j,k) + rhog_step 
            end do

            rhosum = phi(i,j,ARG_H3(phi))+half*rhog(i-1,j,ARG_H3(phi)+1)
            do k = ARG_L3(phi),ARG_H3(phi)
               phi(i,j,k) = phi(i,j,k) - rhosum + press_hi(1)
            end do
         end do
      end if

      j = ARG_L2(phi)
      if (face .eq. YLO ) then
         do i = ARG_L1(phi),ARG_H1(phi)
            k = ARG_L3(phi)
            rhoG_ext = half * rhog(i,j+1,k) 
            rhog_step = rhoG_ext * dx(3) 
            phi(i,j,k)  = phi(i,j,k) + rhog_step  

            do k = ARG_L3(phi)+1,ARG_H3(phi)
               rhoG_ext   = rhog(i,j+1,k) 
               rhog_step  = rhog_step + rhoG_ext * dx(3)
               phi(i,j,k) = phi(i,j,k) + rhog_step 
            end do

            rhosum = phi(i,j,ARG_H3(phi))+half*rhog(i,j+1,ARG_H3(phi)+1)
            do k = ARG_L3(phi),ARG_H3(phi)
               phi(i,j,k) = phi(i,j,k) - rhosum + press_lo(2)
            end do
         end do

      elseif ( face .eq. YHI) then
         do i = ARG_L1(phi),ARG_H1(phi)
            k = ARG_L3(phi)
            rhoG_ext = half * rhog(i,j-1,k)
            rhog_step = rhoG_ext * dx(3) 
            phi(i,j,k)  = phi(i,j,k) + rhog_step  

            do k = ARG_L3(phi)+1,ARG_H3(phi)
               rhoG_ext   = rhog(i,j-1,k) 
               rhog_step  = rhog_step + rhoG_ext * dx(3)
               phi(i,j,k) = phi(i,j,k) + rhog_step 
            end do

            rhosum = phi(i,j,ARG_H3(phi))+half*rhog(i,j-1,ARG_H3(phi)+1)
            do k = ARG_L3(phi),ARG_H3(phi)
               phi(i,j,k) = phi(i,j,k) - rhosum + press_hi(2)
            end do

         end do
      end if

#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI


      end
