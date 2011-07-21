#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "GODUNOV_F.H"
#include "ArrayLim.H"

#define SDIM 3
#define XVEL 1
#define YVEL 2
#define ZVEL 3

#define ALL  999

      subroutine FORT_EST_EIG (s, DIMS(s), lam, dlam, d2lam, DIMS(lam),
     &     rphi, DIMS(rphi), krcoef, DIMS(krcoef), nkrcoef,
     &     uedge, DIMS(uedge), kappax, DIMS(kpx),
     &     vedge, DIMS(vedge), kappay, DIMS(kpy),
     &     wedge, DIMS(wedge), kappaz, DIMS(kpz),
     &     lo,hi,dx,bc,eigmax)
c
c ::: This subroutine estimates the eigenvalues for timestep prediction
c
      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(lam)
      integer DIMDEC(rphi)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM)
      integer bc(SDIM,2)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  lam(DIMV(lam),2)
      REAL_T  dlam(DIMV(lam),3)
      REAL_T  d2lam(DIMV(lam),3)
      REAL_T  uedge(DIMV(uedge))
      REAL_T  vedge(DIMV(vedge))
      REAL_T  wedge(DIMV(wedge))
      REAL_T  kappax(DIMV(kpx)) 
      REAL_T  kappay(DIMV(kpy))
      REAL_T  kappaz(DIMV(kpz))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  eigmax(SDIM), dx(SDIM)

c ::: internal variables
      REAL_T  gravx, gravy, gravz
      
      gravx = zero
      gravy = zero
      gravz = -gravity

      call get_dlbd_cc(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     ncomps,rhoval,muval,model)
      call riemann_est_eig(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     uedge,DIMS(uedge),kappax,DIMS(kpx),rphi,DIMS(rphi),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,0,lo,hi,ncomps,gravx,rhoval,muval,model,dx,bc,eigmax(1))
      call riemann_est_eig(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     vedge,DIMS(vedge),kappay,DIMS(kpy),rphi,DIMS(rphi),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,0,lo,hi,ncomps,gravy,rhoval,muval,dx,bc,eigmax(2))
      call riemann_est_eig(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     wedge,DIMS(wedge),kappaz,DIMS(kpz),rphi,DIMS(rphi),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,0,1,lo,hi,ncomps,gravz,rhoval,muval,dx,bc,eigmax(3))

      end

      subroutine FORT_EST_EIG_LIN (rphi,  DIMS(rphi), 
     &     uedge, DIMS(uedge),
     &     vedge, DIMS(vedge), 
     &     wedge, DIMS(wedge), 
     &     lo,hi,eigmax)
c
c ::: This subroutine estimates the eigenvalues for timestep prediction
c
      implicit none
      
c ::: arguments
      integer DIMDEC(rphi)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer lo(SDIM),hi(SDIM)
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge)), wedge(DIMV(wedge))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  eigmax(SDIM)

c ::: internal variables
      integer i,j,k
      REAL_T  eigtmp
      
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               eigtmp = uedge(i,j,k)/(half*(rphi(i,j,k)+rphi(i-1,j,k)))
               eigmax(1) = max(eigmax(1),dabs(eigtmp))
            end do
         end do
      end do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               eigtmp = vedge(i,j,k)/(half*(rphi(i,j,k)+rphi(i,j-1,k)))
               eigmax(2) = max(eigmax(2),dabs(eigtmp))
            end do
         end do
      end do
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               eigtmp = wedge(i,j,k)/(half*(rphi(i,j,k)+rphi(i,j,k-1)))
               eigmax(3) = max(eigmax(3),dabs(eigtmp))
            end do
         end do
      end do

      end 

      subroutine FORT_EST_EIG_TRC (s, DIMS(s), nc,
     &     rphi,  DIMS(rphi), 
     &     uedge, DIMS(uedge),
     &     vedge, DIMS(vedge), 
     &     wedge, DIMS(wedge), 
     &     lo,hi,eigmax)
c
c ::: This subroutine estimates the eigenvalues for timestep prediction
c
      implicit none

#include "probdata.H"
      
c ::: arguments
      integer nc
      integer DIMDEC(s)
      integer DIMDEC(rphi)
      integer DIMDEC(uedge),DIMDEC(vedge), DIMDEC(wedge)
      integer lo(SDIM),hi(SDIM)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge)), wedge(DIMV(wedge))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  eigmax(SDIM)

c ::: internal variables
      integer i,j,k
      REAL_T  eigtmp,sat
      
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               sat = half*(s(i,j,k,nc)+s(i-1,j,k,nc))/rhoval(nc)
               eigtmp = uedge(i,j,k)/(half*(rphi(i,j,k)+rphi(i-1,j,k)))/sat
               eigmax(1) = max(eigmax(1),dabs(eigtmp))
            end do
         end do
      end do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               sat = half*(s(i,j,k,nc)+s(i,j-1,k,nc))/rhoval(nc)
               eigtmp = vedge(i,j,k)/(half*(rphi(i,j,k)+rphi(i,j-1,k)))/sat
               eigmax(2) = max(eigmax(2),dabs(eigtmp))
            end do
         end do
      end do
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               sat = half*(s(i,j,k,nc)+s(i,j,k-1,nc))/rhoval(nc)
               eigtmp = wedge(i,j,k)/(half*(rphi(i,j,k)+rphi(i,j,k-1)))/sat
               eigmax(3) = max(eigmax(3),dabs(eigtmp))
            end do
         end do
      end do

      end 

      subroutine FORT_EST_EIG_CPL (s, DIMS(s), lam, dlam, d2lam, DIMS(lam),
     &     rphi, DIMS(rphi), krcoef, DIMS(krcoef), nkrcoef,
     &     pc,DIMS(pc),uedge, DIMS(uedge), kappax, DIMS(kpx),
     &     vedge, DIMS(vedge), kappay, DIMS(kpy),
     &     wedge, DIMS(wedge), kappaz, DIMS(kpz),
     &     lo,hi,dx,bc,eigmax)
c
c ::: This subroutine estimates the eigenvalues for timestep prediction. Dummy for CPL scheme
c
      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(pc)
      integer DIMDEC(lam)
      integer DIMDEC(rphi)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM)
      integer bc(SDIM,2)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  pc(DIMV(pc))
      REAL_T  lam(DIMV(lam),2)
      REAL_T  dlam(DIMV(lam),3)
      REAL_T  d2lam(DIMV(lam),3)
      REAL_T  uedge(DIMV(uedge))
      REAL_T  vedge(DIMV(vedge))
      REAL_T  wedge(DIMV(wedge))
      REAL_T  kappax(DIMV(kpx)) 
      REAL_T  kappay(DIMV(kpy))
      REAL_T  kappaz(DIMV(kpz))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  eigmax(SDIM), dx(SDIM)

      call FORT_EST_EIG (s, DIMS(s), lam, dlam, d2lam, DIMS(lam),
     &     rphi, DIMS(rphi), krcoef, DIMS(krcoef), nkrcoef,
     &     uedge, DIMS(uedge), kappax, DIMS(kpx),
     &     vedge, DIMS(vedge), kappay, DIMS(kpy),
     &     wedge, DIMS(wedge), kappaz, DIMS(kpz),
     &     lo,hi,dx,bc,eigmax)
      
      end 

      subroutine FORT_ESTATE_FPU(s, tforces, divu, DIMS(s),
     &     xlo, xhi, sx, slxscr, stxlo, stxhi,
     &     uedge, DIMS(uedge), xstate, DIMS(xstate),

     &     ylo, yhi, sy, slyscr, stylo, styhi,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate),

     &     zlo, zhi, sz, slzscr, stzlo, stzhi,
     &     wedge, DIMS(wedge), zstate, DIMS(zstate),

     &     DIMS(work),
     &     bc,lo,hi,dt,dx,n,use_minion,iconserv)
c
c     This subroutine computes edges states, right now it uses
c     a lot of memory, but there becomes a trade off between
c     simplicity-efficiency in the new way of computing states
c     and complexity in the old way.  By eliminating loops over
c     state components though, the new way uses much less memory.
c
c     This routine differs from the default ESTATE function above in that
c     it assumes that the edge velocities are valid in a grow cell outside
c     the box, and no *ad (unprojected) velocities are used.  This routine
c     will fail if the UMAC coming in hasn't been "fillpatched"
c
      integer i,j,k,n,inc
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer imin,jmin,kmin,imax,jmax,kmax
      REAL_T hx, hy, hz, dt, dth, dthx, dthy, dthz
      REAL_T tr,tr1,tr2,ubar,vbar,wbar,stx,sty,stz,fu,fv,fw,dx(SDIM)
      REAL_T eps,eps_for_bc
      parameter (eps        = 1.d-6 )
      parameter (eps_for_bc = 1.d-10)

      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(xstate)
      integer DIMDEC(vedge)
      integer DIMDEC(ystate)
      integer DIMDEC(wedge)
      integer DIMDEC(zstate)

      REAL_T s(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)
      REAL_T stzlo(DIM3(s)),stzhi(DIM3(s)),slzscr(DIM3(s),4)

      REAL_T uedge(DIMV(uedge)), xstate(DIMV(xstate))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(ystate))
      REAL_T wedge(DIMV(wedge)), zstate(DIMV(zstate))

      REAL_T xlo(DIMV(work)), xhi(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work))
      REAL_T zlo(DIMV(work)), zhi(DIMV(work))
      REAL_T  sx(DIMV(work))
      REAL_T  sy(DIMV(work))
      REAL_T  sz(DIMV(work))
      REAL_T tforces(DIMV(work))
      REAL_T    divu(DIMV(work))

      integer use_minion, iconserv
      REAL_T spx,smx,spy,smy,spz,smz,st

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     compute the slopes
c
      call FORT_SLOPES(ALL,
     &     s,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,bc)
c
c     trace the state to the cell edges
c
      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j,k) = s(i-1,j,k) + (half  - dthx*uedge(i,j,k))*sx(i-1,j,k)
               xhi(i,j,k) = s(i,  j,k) + (-half - dthx*uedge(i,j,k))*sx(i,  j,k)
            end do
         end do
      end do

      if(use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = xlo(i,j,k) + dth*tforces(i-1,j,k)
                  xhi(i,j,k) = xhi(i,j,k) + dth*tforces(i,  j,k)
               end do
            end do
         end do
         if (iconserv .eq. 1) then
           do k = kmin-1,kmax+1
             do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = xlo(i,j,k) - dth*s(i-1,j,k)*divu(i-1,j,k)
                  xhi(i,j,k) = xhi(i,j,k) - dth*s(i  ,j,k)*divu(i,  j,k)
               end do
             end do
           end do
         end if

      end if

      call trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(work),uedge,DIMS(uedge),
     &     lo,hi,n,bc,eps_for_bc)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               fu  = cvmgt(zero,one,abs(uedge(i,j,k)).lt.eps)
               stx = cvmgp(xlo(i,j,k),xhi(i,j,k),uedge(i,j,k))
               xlo(i,j,k) = fu*stx + (one - fu)*half*(xhi(i,j,k)+xlo(i,j,k))
            end do
         end do
      end do

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               ylo(i,j,k) = s(i,j-1,k) + (half  - dthy*vedge(i,j,k))*sy(i,j-1,k)
               yhi(i,j,k) = s(i,j, k)  + (-half - dthy*vedge(i,j,k))*sy(i,j, k)
            end do
         end do
      end do

      if (use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin, jmax+1
               do i = imin-1,  imax+1
                  ylo(i,j,k) = ylo(i,j,k) + dth*tforces(i,j-1,k)
                  yhi(i,j,k) = yhi(i,j,k) + dth*tforces(i,j,  k)
               end do
            end do
         end do
         if (iconserv .eq. 1) then
           do k = kmin-1,kmax+1
             do j = jmin, jmax+1
               do i = imin-1,  imax+1
                  ylo(i,j,k) = ylo(i,j,k) - dth*s(i,j-1,k)*divu(i,j-1,k)
                  yhi(i,j,k) = yhi(i,j,k) - dth*s(i,j  ,k)*divu(i,j,  k)
               end do
             end do
           end do
         end if
      end if

      call trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(work),vedge,DIMS(vedge),
     &     lo,hi,n,bc,eps_for_bc)

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               fv  = cvmgt(zero,one,abs(vedge(i,j,k)).lt.eps)
               sty = cvmgp(ylo(i,j,k),yhi(i,j,k),vedge(i,j,k))
               ylo(i,j,k) = fv*sty + (one - fv)*half*(yhi(i,j,k)+ylo(i,j,k))
            end do
         end do
      end do

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               zlo(i,j,k) = s(i,j,k-1) + (half  - dthz*wedge(i,j,k))*sz(i,j,k-1)
               zhi(i,j,k) = s(i,j,k  ) + (-half - dthz*wedge(i,j,k))*sz(i,j,k  )
            end do
         end do
      end do

      if (use_minion.eq.1)then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  zlo(i,j,k) = zlo(i,j,k) + dth*tforces(i,j,k-1)
                  zhi(i,j,k) = zhi(i,j,k) + dth*tforces(i,j,k)
               end do
            end do
         end do
         if (iconserv .eq. 1) then
           do k = kmin,kmax+1
             do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  zlo(i,j,k) = zlo(i,j,k) - dth*s(i,j,k-1)*divu(i,j,k-1)
                  zhi(i,j,k) = zhi(i,j,k) - dth*s(i,j,k  )*divu(i,j,k  )
               end do
             end do
           end do
         end if
      end if

      call trans_zbc(
     &     s,DIMS(s),
     &     zlo,zhi,DIMS(work),wedge,DIMS(wedge),
     &     lo,hi,n,bc,eps_for_bc)

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               fw  = cvmgt(zero,one,abs(wedge(i,j,k)).lt.eps)
               stz = cvmgp(zlo(i,j,k),zhi(i,j,k),wedge(i,j,k))
               zlo(i,j,k) = fw*stz + (one-fw)*half*(zhi(i,j,k)+zlo(i,j,k))
            end do
         end do
      end do
c
c     compute the xedge states
c
      do k = kmin,kmax
            do j = jmin,jmax
               do i = imin-1,imax+1

                  spx = s(i,j,k) + half*sx(i,j,k)
                  smx = s(i,j,k) - half*sx(i,j,k)

                  if (iconserv.eq.1) then

                     tr =
     &                    (vedge(i,j+1,k)*ylo(i,j+1,k) - vedge(i,j,k)*ylo(i,j,k))/hy +   
     &                    (wedge(i,j,k+1)*zlo(i,j,k+1) - wedge(i,j,k)*zlo(i,j,k))/hz   

                     st = -dth*tr + dth*(tforces(i,j,k) - s(i,j,k)*divu(i,j,k))
     &                    + dth*s(i,j,k)*(vedge(i,j+1,k)-vedge(i,j,k))/hy
     &                    + dth*s(i,j,k)*(wedge(i,j,k+1)-wedge(i,j,k))/hz
c    &                    - dth*s(i,j,k)*(uedge(i+1,j,k)-uedge(i,j,k))/hx

                  else

                     if (vedge(i,j,k)*vedge(i,j+1,k).le.zero) then
                        vbar = half*(vedge(i,j,k)+vedge(i,j+1,k))
                        if (vbar.lt.zero) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr1 = vbar*(s(i,j+inc,k)-s(i,j+inc-1,k))/hy
                     else
                        tr1 = half*(vedge(i,j+1,k) + vedge(i,j,k)) *
     &                               (ylo(i,j+1,k) -   ylo(i,j,k)  ) / hy
                     endif
                     if (wedge(i,j,k)*wedge(i,j,k+1).lt.zero) then
                        wbar = half*(wedge(i,j,k)+wedge(i,j,k+1))
                        if (wbar.lt.zero) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr2 = wbar*(s(i,j,k+inc)-s(i,j,k+inc-1))/hz
                     else
                        tr2 = half*(wedge(i,j,k+1) + wedge(i,j,k)) *
     &                               (zlo(i,j,k+1) -   zlo(i,j,k)  ) / hz
                     endif

                     st = -dth*(tr1 + tr2) + dth*tforces(i,j,k)
                  endif

                  stxlo(i+1)= spx - dthx*uedge(i+1,j,k)*sx(i,j,k) + st
                  stxhi(i  )= smx - dthx*uedge(i  ,j,k)*sx(i,j,k) + st

               end do

               if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j,k).ge.zero) then
                  stxhi(imin) = s(imin-1,j,k)
                  stxlo(imin) = s(imin-1,j,k)
               else if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j,k).lt.zero) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP
     &                 .or.bc(1,1).eq.REFLECT_EVEN) then
                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.REFLECT_ODD) then
                  stxhi(imin) = zero
                  stxlo(imin) = zero
               end if
               if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j,k).le.zero) then
                  stxlo(imax+1) = s(imax+1,j,k)
                  stxhi(imax+1) = s(imax+1,j,k)
               else if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j,k).gt.zero) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP
     &                 .or.bc(1,2).eq.REFLECT_EVEN) then
                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.REFLECT_ODD) then
                  stxlo(imax+1) = zero
                  stxhi(imax+1) = zero
               end if

               do i = imin, imax+1
                  xstate(i,j,k) = cvmgp(stxlo(i),stxhi(i),uedge(i,j,k))
                  xstate(i,j,k) = cvmgt(half*(stxlo(i)+stxhi(i)),xstate(i,j,k)
     &                 ,abs(uedge(i,j,k)).lt.eps)
               end do
            end do
      end do
c
c     compute the yedge states
c
      do k = kmin,kmax
         do i = imin,imax
            do j = jmin-1,jmax+1

                  spy = s(i,j,k) + half*sy(i,j,k)
                  smy = s(i,j,k) - half*sy(i,j,k)

                  if (iconserv.eq.1) then

                     tr =
     &                    (uedge(i+1,j,k)*xlo(i+1,j,k) - uedge(i,j,k)*xlo(i,j,k))/hx +   
     &                    (wedge(i,j,k+1)*zlo(i,j,k+1) - wedge(i,j,k)*zlo(i,j,k))/hz   

                     st = -dth*tr + dth*(tforces(i,j,k) - s(i,j,k)*divu(i,j,k))
     &                    + dth*s(i,j,k)*(uedge(i+1,j,k)-uedge(i,j,k))/hx
     &                    + dth*s(i,j,k)*(wedge(i,j,k+1)-wedge(i,j,k))/hz
c    &                    - dth*s(i,j,k)*(vedge(i,j+1,k)-vedge(i,j,k))/hy

                  else

                     if (uedge(i,j,k)*uedge(i+1,j,k).le.zero) then
                        ubar = half*(uedge(i,j,k)+uedge(i+1,j,k))
                        if (ubar.lt.zero) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr1 = ubar*(s(i+inc,j,k)-s(i+inc-1,j,k))/hx
                     else
                        tr1 = half*(uedge(i+1,j,k) + uedge(i,j,k)) *
     &                               (xlo(i+1,j,k) -   xlo(i,j,k)  ) / hx
                     endif
                     if (wedge(i,j,k)*wedge(i,j,k+1).lt.zero) then
                        wbar = half*(wedge(i,j,k)+wedge(i,j,k+1))
                        if (wbar.lt.zero) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr2 = wbar*(s(i,j,k+inc)-s(i,j,k+inc-1))/hz
                     else
                        tr2 = half*(wedge(i,j,k+1) + wedge(i,j,k)) *
     &                               (zlo(i,j,k+1) -   zlo(i,j,k)  ) / hz
                     endif

                     st = -dth*(tr1 + tr2) + dth*tforces(i,j,k)

                  endif

                  stylo(j+1)= spy - dthy*vedge(i,j+1,k)*sy(i,j,k) + st
                  styhi(j  )= smy - dthy*vedge(i,j  ,k)*sy(i,j,k) + st
               end do

               if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin,k).ge.zero) then
                  styhi(jmin) = s(i,jmin-1,k)
                  stylo(jmin) = s(i,jmin-1,k)
               else if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin,k).lt.zero) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP
     &                 .or.bc(2,1).eq.REFLECT_EVEN) then
                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.REFLECT_ODD) then
                  styhi(jmin) = zero
                  stylo(jmin) = zero
               end if
               
               if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
                  stylo(jmax+1) = s(i,jmax+1,k)
                  styhi(jmax+1) = s(i,jmax+1,k)
               else if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP
     &                 .or.bc(2,2).eq.REFLECT_EVEN) then
                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.REFLECT_ODD) then
                  stylo(jmax+1) = zero
                  styhi(jmax+1) = zero
               end if

               do j=jmin,jmax+1
                  ystate(i,j,k) = cvmgp(stylo(j),styhi(j),vedge(i,j,k))
                  ystate(i,j,k) = cvmgt(half*(stylo(j)+styhi(j)),ystate(i,j,k),
     &                 abs(vedge(i,j,k)).lt.eps)
               end do
            end do
      end do
c
c     compute the zedge states
c
      do j = jmin,jmax
         do i = imin,imax
            do k = kmin-1,kmax+1
               
                  spz = s(i,j,k) + half*sz(i,j,k)
                  smz = s(i,j,k) - half*sz(i,j,k)

                  if (iconserv.eq.1) then
                     tr =
     &                    (uedge(i+1,j,k)*xlo(i+1,j,k) - uedge(i,j,k)*xlo(i,j,k))/hx +   
     &                    (vedge(i,j+1,k)*ylo(i,j+1,k) - vedge(i,j,k)*ylo(i,j,k))/hy   
                     
                     st = -dth*tr + dth*(tforces(i,j,k) - s(i,j,k)*divu(i,j,k))
     &                    + dth*s(i,j,k)*(uedge(i+1,j,k)-uedge(i,j,k))/hx
     &                    + dth*s(i,j,k)*(vedge(i,j+1,k)-vedge(i,j,k))/hy
c    &                    - dth*s(i,j,k)*(wedge(i,j,k+1)-wedge(i,j,k))/hz
                  else
                     if (uedge(i,j,k)*uedge(i+1,j,k).le.zero) then
                        ubar = half*(uedge(i,j,k)+uedge(i+1,j,k))
                        if (ubar.lt.zero) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr1 = ubar*(s(i+inc,j,k)-s(i+inc-1,j,k))/hx
                     else
                        tr1 = half*(uedge(i+1,j,k) + uedge(i,j,k)) *
     &                       (xlo(i+1,j,k) - xlo(i,j,k)  ) / hx
                     endif
                     if (vedge(i,j,k)*vedge(i,j+1,k).lt.zero) then
                        vbar = half*(vedge(i,j,k)+vedge(i,j+1,k))
                        if (vbar.lt.zero) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr2 = vbar*(s(i,j+inc,k)-s(i,j+inc-1,k))/hy
                     else
                        tr2 = half*(vedge(i,j+1,k) + vedge(i,j,k)) *
     &                       (ylo(i,j+1,k) - ylo(i,j,k)  ) / hy
                     endif

                     st = -dth*(tr1 + tr2) + dth*tforces(i,j,k)
                  endif

                  stzlo(k+1)= spz - dthz*wedge(i,j,k+1)*sz(i,j,k) + st
                  stzhi(k  )= smz - dthz*wedge(i,j,k  )*sz(i,j,k) + st

               end do

               if (bc(3,1).eq.EXT_DIR .and. wedge(i,j,kmin).ge.zero) then
                  stzlo(kmin) = s(i,j,kmin-1)
                  stzhi(kmin) = s(i,j,kmin-1)
               else if (bc(3,1).eq.EXT_DIR .and. wedge(i,j,kmin).lt.zero) then
                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP
     &                 .or.bc(3,1).eq.REFLECT_EVEN) then
                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1).eq.REFLECT_ODD) then
                  stzlo(kmin) = zero
                  stzhi(kmin) = zero
               end if
               if (bc(3,2).eq.EXT_DIR .and. wedge(i,j,kmax+1).le.zero) then
                  stzlo(kmax+1) = s(i,j,kmax+1)
                  stzhi(kmax+1) = s(i,j,kmax+1)
               else if (bc(3,2).eq.EXT_DIR .and. wedge(i,j,kmax+1).gt.zero) then
                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP
     &                 .or.bc(3,2).eq.REFLECT_EVEN) then
                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2).eq.REFLECT_ODD) then
                  stzlo(kmax+1) = zero
                  stzhi(kmax+1) = zero
               end if

               do k = kmin,kmax+1
                  zstate(i,j,k) = cvmgp(stzlo(k),stzhi(k),wedge(i,j,k))
                  zstate(i,j,k) = cvmgt(half*(stzlo(k)+stzhi(k)),zstate(i,j,k),
     &                 abs(wedge(i,j,k)).lt.eps)
               end do
            end do
      end do

      end

      subroutine FORT_ESTATE_LIN (s, DIMS(s), sn, DIMS(sn), 
     &     tf, rphi, DIMS(rphi),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), 
     &     sz, slzscr, 
     &     wedge, DIMS(wedge), zstate, DIMS(zstate),
     &     DIMS(work),
     &     bc,lo,hi,dt,dx,nscal)
      
c
c ::: This subroutine computes edges states and fluxes. 
c ::: ncomps must be equal to nscal.
c
      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(rphi)
      integer DIMDEC(work)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(xstate),DIMDEC(ystate),DIMDEC(zstate)
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer nscal
      REAL_T   s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6),slzscr(DIM2(s),6)
      REAL_T  uedge(DIMV(uedge)), xstate(DIMV(xstate),nscal)
      REAL_T  vedge(DIMV(vedge)), ystate(DIMV(ystate),nscal)  
      REAL_T  wedge(DIMV(wedge)), zstate(DIMV(zstate),nscal)
      REAL_T  sx(DIMV(work),nscal),  sy(DIMV(work),nscal), sz(DIMV(work),nscal)
      REAL_T  tf(DIMV(work),nscal)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  dx(SDIM), dt

c ::: internal variables
      REAL_T  gravx, gravy, gravz
      REAL_T, ALLOCATABLE :: fluxx(:,:,:,:),fluxy(:,:,:,:),fluxz(:,:,:,:)
      REAL_T, ALLOCATABLE :: xlo(:,:,:,:), xhi(:,:,:,:)
      REAL_T, ALLOCATABLE :: ylo(:,:,:,:), yhi(:,:,:,:)
      REAL_T, ALLOCATABLE :: zlo(:,:,:,:), zhi(:,:,:,:)
      integer nc
      
      gravx = zero
      gravy = zero
      gravz = -gravity

      ALLOCATE(fluxx(DIMV(uedge),nscal))
      ALLOCATE(fluxy(DIMV(vedge),nscal))
      ALLOCATE(fluxz(DIMV(wedge),nscal))
      ALLOCATE(xlo(DIMV(uedge),nscal))
      ALLOCATE(ylo(DIMV(vedge),nscal))
      ALLOCATE(zlo(DIMV(wedge),nscal))
      ALLOCATE(xhi(DIMV(uedge),nscal))
      ALLOCATE(yhi(DIMV(vedge),nscal)) 
      ALLOCATE(zhi(DIMV(wedge),nscal)) 

c ::: compute the slopes
      do nc = 1,nscal
         call FORT_SLOPES(ALL,s(ARG_L1(s),ARG_L2(s),ARG_L3(s),nc),DIMS(s),
     &        sx(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),
     &        sy(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),
     &        sz(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),DIMS(work),
     &        lo,hi,slxscr,slyscr,slzscr,bc)
      end do     

c ::: initial trace
      call trace_lin(s,DIMS(s),sn,DIMS(sn),
     &     rphi,DIMS(rphi),
     &     sx,sy,sz,tf,DIMS(work),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     wedge,zlo,zhi,DIMS(wedge),
     &     dt,dx,lo,hi,bc,nscal)
      
c ::: predictor step: determine final edge states 
      call riemann_solver_lin(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,nscal)
      call riemann_solver_lin(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,nscal)
      call riemann_solver_lin(fluxz,wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,nscal)

      call trace_p_lin(s,DIMS(s),sn,DIMS(sn),
     &     fluxx,uedge,DIMS(uedge),
     &     fluxy,vedge,DIMS(vedge),
     &     fluxz,wedge,DIMS(wedge),
     &     rphi,DIMS(rphi),xlo,xhi,ylo,yhi,zlo,zhi,
     &     dt,dx,lo,hi,bc,nscal)

c ::: compute fluxes based on final edge states
      call riemann_solver_lin(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,nscal)
      call riemann_solver_lin(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,nscal)
      call riemann_solver_lin(fluxz,wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,nscal)


      call flux_copy(xstate,DIMS(xstate),fluxx,DIMS(uedge),
     &     lo,hi,1,0,0,nscal)
      call flux_copy(ystate,DIMS(ystate),fluxy,DIMS(vedge),
     &     lo,hi,0,1,0,nscal)   
      call flux_copy(zstate,DIMS(zstate),fluxz,DIMS(wedge),
     &     lo,hi,0,0,1,nscal) 

      DEALLOCATE(fluxx,fluxy,fluxz)
      DEALLOCATE(xlo,xhi,ylo,yhi,zlo,zhi)

      end

      subroutine FORT_ESTATE_RMN (s,  DIMS(s), sn, DIMS(sn), 
     &     tf, divu, lbd_cc, dlbd_cc,
     &     rphi, DIMS(rphi), kappa, DIMS(kappa),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate),kappax, DIMS(kpx),
     &     eigvx, eiglx, eigrx, 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate),kappay, DIMS(kpy),
     &     eigvy, eigly, eigry, 
     &     sz, slzscr, 
     &     wedge, DIMS(wedge), zstate, DIMS(zstate),kappaz, DIMS(kpz),
     &     eigvz, eiglz, eigrz, 
     &     DIMS(work),
     &     krcoef, DIMS(krcoef), nkrcoef,
     &     bc,lo,hi,dt,dx,use_minion,iconserv,nscal)
c
c ::: This subroutine computes edges states and fluxes for 
c ::: single phase model and two-phase two-component model.
c ::: ncomps must be equal to nscal.
c
      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(work)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(xstate),DIMDEC(ystate),DIMDEC(zstate)
      integer DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer use_minion, iconserv, nscal
      REAL_T  s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6),slzscr(DIM3(s),6)
      REAL_T  uedge(DIMV(uedge)),xstate(DIMV(xstate),nscal)
      REAL_T  vedge(DIMV(vedge)),ystate(DIMV(ystate),nscal)
      REAL_T  wedge(DIMV(wedge)),zstate(DIMV(zstate),nscal)
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy)), kappaz(DIMV(kpz))
      REAL_T  sx(DIMV(work),nscal)
      REAL_T  sy(DIMV(work),nscal)
      REAL_T  sz(DIMV(work),nscal)
      REAL_T  tf(DIMV(work))
      REAL_T  divu(DIMV(work))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  lbd_cc(DIMV(work),ncomps)
      REAL_T  dlbd_cc(DIMV(work),3)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  dt,dx(SDIM)
      REAL_T  eigvx(DIMV(work),2*nscal)
      REAL_T  eigvy(DIMV(work),2*nscal)
      REAL_T  eigvz(DIMV(work),2*nscal)
      REAL_T  eiglx(DIMV(work),nscal,nscal),eigrx(DIMV(work),nscal,nscal)
      REAL_T  eigly(DIMV(work),nscal,nscal),eigry(DIMV(work),nscal,nscal)
      REAL_T  eiglz(DIMV(work),nscal,nscal),eigrz(DIMV(work),nscal,nscal)

c ::: internal variables
      integer nc, corrector
      REAL_T  gravx, gravy, gravz
      REAL_T, ALLOCATABLE :: fluxx(:,:,:,:),fluxy(:,:,:,:),fluxz(:,:,:,:)
      REAL_T, ALLOCATABLE :: xlo(:,:,:,:), xhi(:,:,:,:)
      REAL_T, ALLOCATABLE :: ylo(:,:,:,:), yhi(:,:,:,:)
      REAL_T, ALLOCATABLE :: zlo(:,:,:,:), zhi(:,:,:,:)
      

      gravx = zero
      gravy = zero
      gravz = -gravity

      ALLOCATE(fluxx(DIMV(uedge),nscal))
      ALLOCATE(fluxy(DIMV(vedge),nscal))
      ALLOCATE(fluxz(DIMV(wedge),nscal))
      ALLOCATE(xlo(DIMV(uedge),nscal))
      ALLOCATE(ylo(DIMV(vedge),nscal))
      ALLOCATE(zlo(DIMV(wedge),nscal))
      ALLOCATE(xhi(DIMV(uedge),nscal))
      ALLOCATE(yhi(DIMV(vedge),nscal)) 
      ALLOCATE(zhi(DIMV(wedge),nscal))

c ::: compute eigen solutions
      call eigsol_mp(eigvx,eiglx,eigrx,lbd_cc,dlbd_cc,DIMS(work),
     &     uedge,DIMS(uedge),rphi,DIMS(rphi),kappax,DIMS(kpx),
     &     krcoef,DIMS(krcoef),nkrcoef,gravx,muval,
     &     1,0,0,lo,hi,rhoval,ncomps,nphases,model,bc)
      call eigsol_mp(eigvy,eigly,eigry,lbd_cc,dlbd_cc,DIMS(work),
     &     vedge,DIMS(vedge),rphi,DIMS(rphi),kappay,DIMS(kpy),
     &     krcoef,DIMS(krcoef),nkrcoef,gravy,muval,
     &     0,1,0,lo,hi,rhoval,ncomps,nphases,model,bc)
      call eigsol_mp(eigvz,eiglz,eigrz,lbd_cc,dlbd_cc,DIMS(work),
     &     wedge,DIMS(wedge),rphi,DIMS(rphi),kappaz,DIMS(kpz),
     &     krcoef,DIMS(krcoef),nkrcoef,gravz,muval,
     &     0,0,1,lo,hi,rhoval,ncomps,nphases,model,bc)

c ::: compute the slopes
      do nc = 1,nscal
         call FORT_SLOPES(ALL,s(ARG_L1(s),ARG_L2(s),ARG_L3(s),nc),DIMS(s),
     &        sx(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),
     &        sy(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),
     &        sz(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),DIMS(work),
     &        lo,hi,slxscr,slyscr,slzscr,bc)
      end do     
 
      call trace(s,DIMS(s),sn,DIMS(sn),
     &     eigvx,eiglx,eigrx,eigvy,eigly,eigry,eigvz,eiglz,eigrz,
     &     sx,sy,sz,divu,tf,DIMS(work),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     wedge,zlo,zhi,DIMS(wedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,nscal)

c ::: predictor step: determine final edge states 
      corrector = 0
      call riemann_solver(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),eiglx,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,0,lo,hi,gravx,muval,rhoval,nscal,model,corrector)
      call riemann_solver(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),eigly,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,0,lo,hi,gravy,muval,rhoval,nscal,model,corrector)
      call riemann_solver(fluxz,wedge,zlo,zhi,DIMS(wedge),
     &     kappaz,DIMS(kpz),eiglz,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,0,1,lo,hi,gravz,muval,rhoval,nscal,model,corrector)

      call trace_p(s,DIMS(s),sn,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),wedge,DIMS(wedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),kappaz,DIMS(kpz),
     &     fluxx,DIMS(uedge),fluxy,DIMS(vedge),fluxz,DIMS(wedge),
     &     xlo,xhi,ylo,yhi,zlo,zhi,divu,tf,DIMS(work),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     dt,dx,lo,hi,bc,muval,rhoval,iconserv,
     &     nscal,model,gravz)

c ::: compute fluxes based on final edge states
      corrector = 1
      call riemann_solver(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),eiglx,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,0,lo,hi,gravx,muval,rhoval,nscal,model,corrector)
      call riemann_solver(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),eigly,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,0,lo,hi,gravy,muval,rhoval,nscal,model,corrector)
      call riemann_solver(fluxz,wedge,zlo,zhi,DIMS(wedge),
     &     kappaz,DIMS(kpz),eiglz,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,0,1,lo,hi,gravz,muval,rhoval,nscal,model,corrector)

      call flux_copy(xstate,DIMS(xstate),fluxx,DIMS(uedge),
     &     lo,hi,1,0,0,nscal)
      call flux_copy(ystate,DIMS(ystate),fluxy,DIMS(vedge),
     &     lo,hi,0,1,0,nscal)   
      call flux_copy(zstate,DIMS(zstate),fluxz,DIMS(wedge),
     &     lo,hi,0,0,1,nscal) 

      DEALLOCATE(fluxx,fluxy,fluxz)
      DEALLOCATE(xlo,xhi,ylo,yhi,zlo,zhi)

      end

      subroutine FORT_ESTATE_CPL (s,  DIMS(s), sn, DIMS(sn),
     &     pc,DIMS(pc),
     &     tf, divu, lbd_cc, dlbd_cc,
     &     rphi, DIMS(rphi), kappa, DIMS(kappa),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate),kappax, DIMS(kpx),
     &     eigvx, eiglx, eigrx, 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate),kappay, DIMS(kpy),
     &     eigvy, eigly, eigry, 
     &     sz, slzscr, 
     &     wedge, DIMS(wedge), zstate, DIMS(zstate),kappaz, DIMS(kpz),
     &     eigvz, eiglz, eigrz, 
     &     DIMS(work),
     &     krcoef, DIMS(krcoef), nkrcoef,
     &     bc,lo,hi,dt,dx,use_minion,iconserv,nscal)
c
c ::: This subroutine computes edges states and fluxes for 
c ::: single phase model and two-phase two-component model.
c ::: ncomps must be equal to nscal.
c
      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(work)
      integer DIMDEC(pc)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(xstate),DIMDEC(ystate),DIMDEC(zstate)
      integer DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer use_minion, iconserv, nscal
      REAL_T  s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  pc(DIMV(pc))
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6),slzscr(DIM3(s),6)
      REAL_T  uedge(DIMV(uedge)),xstate(DIMV(xstate),nscal)
      REAL_T  vedge(DIMV(vedge)),ystate(DIMV(ystate),nscal)
      REAL_T  wedge(DIMV(wedge)),zstate(DIMV(zstate),nscal)
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy)), kappaz(DIMV(kpz))
      REAL_T  sx(DIMV(work),nscal)
      REAL_T  sy(DIMV(work),nscal)
      REAL_T  sz(DIMV(work),nscal)
      REAL_T  tf(DIMV(work))
      REAL_T  divu(DIMV(work))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  lbd_cc(DIMV(work),ncomps)
      REAL_T  dlbd_cc(DIMV(work),3)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  dt,dx(SDIM)
      REAL_T  eigvx(DIMV(work),2*nscal)
      REAL_T  eigvy(DIMV(work),2*nscal)
      REAL_T  eigvz(DIMV(work),2*nscal)
      REAL_T  eiglx(DIMV(work),nscal,nscal),eigrx(DIMV(work),nscal,nscal)
      REAL_T  eigly(DIMV(work),nscal,nscal),eigry(DIMV(work),nscal,nscal)
      REAL_T  eiglz(DIMV(work),nscal,nscal),eigrz(DIMV(work),nscal,nscal)

      call FORT_ESTATE_RMN (s,  DIMS(s), sn, DIMS(sn), 
     &     tf, divu, lbd_cc, dlbd_cc,
     &     rphi, DIMS(rphi), kappa, DIMS(kappa),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate),kappax, DIMS(kpx),
     &     eigvx, eiglx, eigrx, 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate),kappay, DIMS(kpy),
     &     eigvy, eigly, eigry, 
     &     sz, slzscr, 
     &     wedge, DIMS(wedge), zstate, DIMS(zstate),kappaz, DIMS(kpz),
     &     eigvz, eiglz, eigrz, 
     &     DIMS(work),
     &     krcoef, DIMS(krcoef), nkrcoef,
     &     bc,lo,hi,dt,dx,use_minion,iconserv,nscal)

      end

      subroutine FORT_ESTATE_PMR (s, u, DIMS(s), 
     &     sn, un, DIMS(sn), tf, divu,
     &     rphi, DIMS(rphi), kappa, DIMS(kappa),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), kappax, DIMS(kpx),
     &     strucx, eigvx, eiglx, eigrx, 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), kappay, DIMS(kpy),
     &     strucy, eigvy, eigly, eigry, 
     &     sz, slzscr, 
     &     wedge, DIMS(wedge), zstate, DIMS(zstate),kappaz, DIMS(kpz),
     &     strucz, eigvz, eiglz, eigrz, 
     &     DIMS(work),
     &     bc,lo,hi,dt,dx,
     &     use_minion,iconserv,grav,eigmax)
c
c ::: This subroutine computes edges states and fluxes for the polymer problem
c     The algorithm goes as follows:
c     1. transform problem into the auxiliary variables
c     2. do everything in the auxiliary variables
c     3. transform back the flux only at the very last step.

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(work)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(xstate),DIMDEC(ystate),DIMDEC(zstate)
      integer DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer use_minion, iconserv
      REAL_T  s (DIMV( s),ncomps)
      REAL_T  sn(DIMV(sn),ncomps)
      REAL_T  u (DIMV(s),ncomps-1)
      REAL_T  un(DIMV(sn),ncomps-1)
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6),slzscr(DIM3(s),6)
      REAL_T  uedge(DIMV(uedge)), xstate(DIMV(xstate),ncomps)
      REAL_T  vedge(DIMV(vedge)), ystate(DIMV(ystate),ncomps)    
      REAL_T  wedge(DIMV(wedge)), zstate(DIMV(zstate),ncomps)
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy)), kappaz(DIMV(kpz))
      REAL_T  sx(DIMV(work),ncomps-1)
      REAL_T  sy(DIMV(work),ncomps-1)
      REAL_T  sz(DIMV(work),ncomps-1)
      REAL_T  tf(DIMV(work))
      REAL_T  divu(DIMV(work))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  dx(SDIM), eigmax(SDIM),  grav
      REAL_T  eigvx(DIMV(work),ncomps-1)
      REAL_T  eigvy(DIMV(work),ncomps-1)
      REAL_T  eigvz(DIMV(work),ncomps-1)
      REAL_T  eiglx(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  eigrx(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  eigly(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  eigry(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  eiglz(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  eigrz(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  strucx(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  strucy(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  strucz(DIMV(work),ncomps-1,ncomps-1)

c ::: internal variables
      integer i,j,k
      integer compcount(nphases)
      integer compidx(nphases,ncomps)
      logical, ALLOCATABLE ::  ldefx(:,:,:), ldefy(:,:,:), ldefz(:,:,:)
      REAL_T, ALLOCATABLE :: fluxx(:,:,:,:),fluxy(:,:,:,:),fluxz(:,:,:,:)
      REAL_T, ALLOCATABLE :: xlo(:,:,:,:), xhi(:,:,:,:)
      REAL_T, ALLOCATABLE :: ylo(:,:,:,:), yhi(:,:,:,:)
      REAL_T, ALLOCATABLE :: zlo(:,:,:,:), zhi(:,:,:,:)
      REAL_T  hx, hy, hz, dt, dth, dthx, dthy, dthz
      REAL_T  gravx,gravy,gravz
      REAL_T  athick
      REAL_T  mu(nphases), rho(nphases)
      
      ALLOCATE(ldefx(DIMV(work)))
      ALLOCATE(ldefy(DIMV(work)))
      ALLOCATE(ldefz(DIMV(work)))
      ALLOCATE(fluxx(DIMV(uedge),ncomps-1))
      ALLOCATE(fluxy(DIMV(vedge),ncomps-1))
      ALLOCATE(fluxz(DIMV(wedge),ncomps-1))
      ALLOCATE(xlo(DIMV(uedge),ncomps-1))
      ALLOCATE(ylo(DIMV(vedge),ncomps-1))
      ALLOCATE(zlo(DIMV(wedge),ncomps-1))
      ALLOCATE(xhi(DIMV(uedge),ncomps-1))
      ALLOCATE(yhi(DIMV(vedge),ncomps-1)) 
      ALLOCATE(zhi(DIMV(wedge),ncomps-1))

      dth   = half*dt
      dthx  = half*dt / dx(1)
      dthy  = half*dt / dx(2)
      dthz  = half*dt / dx(3)
      hx    = dx(1)
      hy    = dx(2)
      hz    = dx(3)
      gravx = zero
      gravy = zero
      gravz = -grav

c ::: identify comp-phase relation (in GODUNOV_F.F)
      call get_mu_rho_phases(rho, mu, athick, rhoval, muval, 
     &     compcount, compidx, phaseidx, ncomps, nphases)
c ::: s --> u
      
      call get_aux(s,u,DIMS(s),compcount,compidx,rhoval,ncomps,nphases)
      call get_aux_n(sn,un,DIMS(sn),compcount,compidx,rhoval,bc,ncomps,nphases)

c ::: compute the equivalent eigensystem
c ::: Jacobian entries computed on the fly

      call eigsol2(u,DIMS(s),uedge,DIMS(uedge),rphi,DIMS(rphi),
     &     kappa,DIMS(work),eigvx,eiglx,eigrx,strucx,ldefx,DIMS(work),
     &     lo,hi,1,0,0,mu,athick,rho,gravx,ncomps-1,nphases,eigmax(1))
      call eigsol2(u,DIMS(s),vedge,DIMS(vedge),rphi,DIMS(rphi),
     &     kappa,DIMS(work),eigvy,eigly,eigry,strucy,ldefy,DIMS(work),
     &     lo,hi,0,1,0,mu,athick,rho,gravy,ncomps-1,nphases,eigmax(2))
      call eigsol2(u,DIMS(s),wedge,DIMS(wedge),rphi,DIMS(rphi),
     &     kappa,DIMS(work),eigvz,eiglz,eigrz,strucz,ldefz,DIMS(work),
     &     lo,hi,0,0,1,mu,athick,rho,gravz,ncomps-1,nphases,eigmax(3))

c ::: compute the slopes
      call slope_alpha(u,DIMS(s),sx,sy,sz,DIMS(work),
     &     eigvx,eiglx,eigrx,strucx,ldefx,
     &     eigvy,eigly,eigry,strucy,ldefy,
     &     eigvz,eiglz,eigrz,strucz,ldefz,DIMS(work),
     &     slxscr,slyscr,slzscr,lo,hi,bc,ncomps-1,nphases)

c ::: initial step
      call trace2(u,DIMS(s),un,DIMS(sn),
     &     eigvx,eigrx,eigvy,eigry,eigvz,eigrz,DIMS(work),
     &     sx,sy,sz,divu,tf,DIMS(work),uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),wedge,zlo,zhi,DIMS(wedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,ncomps-1,nphases)
      
c ::: predictor step: determine final edge states 
      call riemann_solver2(fluxx,DIMS(uedge),u,DIMS(s),xlo,xhi,DIMS(uedge),
     &     eigvx,DIMS(work),uedge,DIMS(uedge),kappax,DIMS(kpx),
     &     1,0,0,lo,hi,athick,gravx,mu,rho,ncomps-1,nphases,dx,dt,bc)
      call riemann_solver2(fluxy,DIMS(vedge),u,DIMS(s),ylo,yhi,DIMS(vedge),
     &     eigvy,DIMS(work),vedge,DIMS(vedge),kappay,DIMS(kpy),
     &     0,1,0,lo,hi,athick,gravy,mu,rho,ncomps-1,nphases,dx,dt,bc)
      call riemann_solver2(fluxz,DIMS(wedge),u,DIMS(s),zlo,zhi,DIMS(wedge),
     &     eigvz,DIMS(work),wedge,DIMS(wedge),kappaz,DIMS(kpz),
     &     0,0,1,lo,hi,athick,gravz,mu,rho,ncomps-1,nphases,dx,dt,bc)

      call trace_p2(u,DIMS(s),un,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),wedge,DIMS(wedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),kappaz,DIMS(kpz),
     &     fluxx,DIMS(uedge),fluxy,DIMS(vedge),fluxz,DIMS(wedge),
     &     xlo,xhi,ylo,yhi,zlo,zhi,divu,tf,DIMS(work),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     dt,dx,lo,hi,bc,mu,rho,athick,iconserv,ncomps-1,nphases,gravz)

c ::: compute fluxes based on final edge states
      call riemann_solver2(fluxx,DIMS(uedge),u,DIMS(s),xlo,xhi,DIMS(uedge),
     &     eigvx,DIMS(work),uedge,DIMS(uedge),kappax,DIMS(kpx),
     &     1,0,0,lo,hi,athick,gravx,mu,rho,ncomps-1,nphases,dx,dt,bc)
      call riemann_solver2(fluxy,DIMS(vedge),u,DIMS(s),ylo,yhi,DIMS(vedge),
     &     eigvy,DIMS(work),vedge,DIMS(vedge),kappay,DIMS(kpy),
     &     0,1,0,lo,hi,athick,gravy,mu,rho,ncomps-1,nphases,dx,dt,bc) 
      call riemann_solver2(fluxz,DIMS(wedge),u,DIMS(s),zlo,zhi,DIMS(wedge),
     &     eigvz,DIMS(work),wedge,DIMS(wedge),kappaz,DIMS(kpz),
     &     0,0,1,lo,hi,athick,gravz,mu,rho,ncomps-1,nphases,dx,dt,bc)

c ::: f(u) -> f(s)
      call flux_convert(xstate,DIMS(xstate),fluxx,DIMS(uedge),uedge,DIMS(uedge),
     &     lo,hi,1,0,0,rho,compcount,compidx,ncomps,nphases,bc)
      call flux_convert(ystate,DIMS(ystate),fluxy,DIMS(vedge),vedge,DIMS(vedge),
     &     lo,hi,0,1,0,rho,compcount,compidx,ncomps,nphases,bc)   
      call flux_convert(zstate,DIMS(zstate),fluxz,DIMS(wedge),wedge,DIMS(wedge),
     &     lo,hi,0,0,1,rho,compcount,compidx,ncomps,nphases,bc)  

c ::: convert the second auxiliary variable to c instead of c*s
      do k = ARG_L3(s), ARG_H3(s)
         do j =  ARG_L2(s), ARG_H2(s)
            do i = ARG_L1(s), ARG_H1(s)
               if (u(i,j,k,1) .gt. 1.d-10) then
                  u(i,j,k,2) = u(i,j,k,2)/u(i,j,k,1)
               end if
            end do
         end do
      end do
      
      DEALLOCATE(ldefx,ldefy,ldefz)
      DEALLOCATE(fluxx,fluxy,fluxz)
      DEALLOCATE(xlo,xhi,ylo,yhi,zlo,zhi)

      end

      subroutine FORT_ESTATE_TRACER (s, DIMS(s),
     &     sn, DIMS(sn), 
     &     st, DIMS(st), stn, DIMS(stn),
     &     sx, slxscr, xlo, xhi, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), 
     &     sy, slyscr, ylo, yhi,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), 
     &     sz, slzscr, zlo, zhi,
     &     wedge, DIMS(wedge), zstate, DIMS(zstate), 
     &     rphi,  DIMS(rphi), DIMS(work),
     &     bc,lo,hi,dt,dx,nscal)
c
c ::: This subroutine computes tracers' fluxes at edges
c
      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(st)
      integer DIMDEC(rphi)
      integer DIMDEC(stn)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(wedge)
      integer DIMDEC(xstate)
      integer DIMDEC(ystate)
      integer DIMDEC(zstate)
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer nscal
      REAL_T   s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  st(DIMV(st))
      REAL_T  stn(DIMV(stn))
      REAL_T  slxscr(DIM1(s),6)
      REAL_T  slyscr(DIM2(s),6)
      REAL_T  slzscr(DIM3(s),6)
      REAL_T  uedge(DIMV(uedge))
      REAL_T  vedge(DIMV(vedge))
      REAL_T  wedge(DIMV(wedge))
      REAL_T  xstate(DIMV(xstate),nscal)      
      REAL_T  ystate(DIMV(ystate),nscal)
      REAL_T  zstate(DIMV(zstate),nscal)
      REAL_T  sx(DIMV(work),nscal)
      REAL_T  sy(DIMV(work),nscal)
      REAL_T  sz(DIMV(work),nscal)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  dx(SDIM)
      REAL_T  xlo(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal)
      REAL_T  zlo(DIMV(wedge),nscal)
      REAL_T  xhi(DIMV(uedge),nscal)
      REAL_T  yhi(DIMV(vedge),nscal)
      REAL_T  zhi(DIMV(wedge),nscal)

      integer nc
      REAL_T hx, hy, hz, dt, dth, dthx, dthy, dthz
      
      dth   = half*dt
      dthx  = half*dt / dx(1)
      dthy  = half*dt / dx(2)
      dthz  = half*dt / dx(3)
      hx    = dx(1)
      hy    = dx(2)
      hz    = dx(3)

c ::: compute the slopes
      do nc = 1,nscal
         call FORT_SLOPES(ALL,s(ARG_L1(s),ARG_L2(s),ARG_L3(s),nc),DIMS(s),
     &        sx(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),
     &        sy(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),
     &        sz(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),DIMS(work),
     &        lo,hi,slxscr,slyscr,slzscr,bc)
      end do 

      call trace_tracers(s,DIMS(s),sn,DIMS(sn),
     &     st,DIMS(st),stn,DIMS(stn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),wedge,DIMS(wedge),
     &     rphi,DIMS(rphi),sx,sy,sz,DIMS(work),
     &     xlo,xhi,ylo,yhi,zlo,zhi,
     &     dt,dx,lo,hi,bc,nscal) 

c ::: predictor step: determine final edge states 
      call advect_tracers(xstate,DIMS(xstate),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,nscal)
      call advect_tracers(ystate,DIMS(ystate),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,nscal)
      call advect_tracers(zstate,DIMS(zstate),wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,nscal)

      end

      subroutine FORT_SYNC_RMN (s,  DIMS(s), sn, DIMS(sn), 
     &     tf, divu, lbd_cc, dlbd_cc,
     &     rphi, DIMS(rphi), kappa, DIMS(kappa),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate),kappax, DIMS(kpx),
     &     eigvx, eiglx, eigrx, 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate),kappay, DIMS(kpy),
     &     eigvy, eigly, eigry, 
     &     sz, slzscr, 
     &     wedge, DIMS(wedge), zstate, DIMS(zstate),kappaz, DIMS(kpz),
     &     eigvz, eiglz, eigrz, 
     &     DIMS(work),
     &     krcoef, DIMS(krcoef), nkrcoef,
     &     bc,lo,hi,dt,dx,use_minion,iconserv,nscal)
c
c ::: This subroutine computes edges states and fluxes
c

      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(work)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(xstate),DIMDEC(ystate),DIMDEC(zstate)
      integer DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer use_minion, iconserv, nscal
      REAL_T  s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6),slzscr(DIM3(s),6)
      REAL_T  uedge(DIMV(uedge)),xstate(DIMV(xstate),nscal)
      REAL_T  vedge(DIMV(vedge)),ystate(DIMV(ystate),nscal)
      REAL_T  wedge(DIMV(wedge)),zstate(DIMV(zstate),nscal)
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy)), kappaz(DIMV(kpz))
      REAL_T  sx(DIMV(work),nscal)
      REAL_T  sy(DIMV(work),nscal)
      REAL_T  sz(DIMV(work),nscal)
      REAL_T  tf(DIMV(work))
      REAL_T  divu(DIMV(work))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  lbd_cc(DIMV(work),ncomps)
      REAL_T  dlbd_cc(DIMV(work),3)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  dx(SDIM)
      REAL_T  eigvx(DIMV(work),2*nscal)
      REAL_T  eigvy(DIMV(work),2*nscal)
      REAL_T  eigvz(DIMV(work),2*nscal)
      REAL_T  eiglx(DIMV(work),nscal,nscal),eigrx(DIMV(work),nscal,nscal)
      REAL_T  eigly(DIMV(work),nscal,nscal),eigry(DIMV(work),nscal,nscal)
      REAL_T  eiglz(DIMV(work),nscal,nscal),eigrz(DIMV(work),nscal,nscal)

c ::: internal variables
      integer nc, corrector
      REAL_T  hx, hy, hz, dt, dth, dthx, dthy, dthz
      REAL_T  gravx, gravy, gravz
      REAL_T, ALLOCATABLE :: fluxx(:,:,:,:),fluxy(:,:,:,:),fluxz(:,:,:,:)
      REAL_T, ALLOCATABLE :: xlo(:,:,:,:), xhi(:,:,:,:)
      REAL_T, ALLOCATABLE :: ylo(:,:,:,:), yhi(:,:,:,:)
      REAL_T, ALLOCATABLE :: zlo(:,:,:,:), zhi(:,:,:,:)
      
      dth   = half*dt
      dthx  = half*dt / dx(1)
      dthy  = half*dt / dx(2)
      dthz  = half*dt / dx(3)
      hx    = dx(1)
      hy    = dx(2)
      hz    = dx(3)
      gravx = zero
      gravy = zero
      gravz = -gravity

      ALLOCATE(fluxx(DIMV(uedge),nscal))
      ALLOCATE(fluxy(DIMV(vedge),nscal))
      ALLOCATE(fluxz(DIMV(wedge),nscal))
      ALLOCATE(xlo(DIMV(uedge),nscal))
      ALLOCATE(ylo(DIMV(vedge),nscal))
      ALLOCATE(zlo(DIMV(wedge),nscal))
      ALLOCATE(xhi(DIMV(uedge),nscal))
      ALLOCATE(yhi(DIMV(vedge),nscal)) 
      ALLOCATE(zhi(DIMV(wedge),nscal))

c ::: compute eigen solutions
      call eigsol_mp(eigvx,eiglx,eigrx,lbd_cc,dlbd_cc,DIMS(work),
     &     uedge,DIMS(uedge),rphi,DIMS(rphi),kappax,DIMS(kpx),    
     &     krcoef,DIMS(krcoef),nkrcoef,gravx,muval,
     &     1,0,0,lo,hi,rhoval,ncomps,nphases,model,bc)
      call eigsol_mp(eigvy,eigly,eigry,lbd_cc,dlbd_cc,DIMS(work),
     &     vedge,DIMS(vedge),rphi,DIMS(rphi),kappay,DIMS(kpy),    
     &     krcoef,DIMS(krcoef),nkrcoef,gravy,muval,
     &     0,1,0,lo,hi,rhoval,ncomps,nphases,model,bc)
      call eigsol_mp(eigvz,eiglz,eigrz,lbd_cc,dlbd_cc,DIMS(work),
     &     wedge,DIMS(wedge),rphi,DIMS(rphi),kappaz,DIMS(kpz),    
     &     krcoef,DIMS(krcoef),nkrcoef,gravz,muval,
     &     0,0,1,lo,hi,rhoval,ncomps,nphases,model,bc)

c ::: compute the slopes
      do nc = 1,nscal
         call FORT_SLOPES(ALL,s(ARG_L1(s),ARG_L2(s),ARG_L3(s),nc),DIMS(s),
     &        sx(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),
     &        sy(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),
     &        sz(ARG_L1(work),ARG_L2(work),ARG_L3(work),nc),DIMS(work),
     &        lo,hi,slxscr,slyscr,slzscr,bc)
      end do
     
      call trace(s,DIMS(s),sn,DIMS(sn),
     &     eigvx,eiglx,eigrx,eigvy,eigly,eigry,eigvz,eiglz,eigrz,
     &     sx,sy,sz,divu,tf,DIMS(work),
     &     uedge,xlo,xhi,DIMS(uedge),vedge,ylo,yhi,DIMS(vedge),
     &     wedge,zlo,zhi,DIMS(wedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,nscal)

c ::: predictor step: determine final edge states 
      corrector = 0
      call riemann_solver(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),eiglx,DIMS(work),    
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,0,lo,hi,gravx,muval,rhoval,nscal,model,corrector)
      call riemann_solver(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),eigly,DIMS(work),    
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,0,lo,hi,gravy,muval,rhoval,nscal,model,corrector)
      call riemann_solver(fluxz,wedge,zlo,zhi,DIMS(wedge),
     &     kappaz,DIMS(kpz),eiglz,DIMS(work),    
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,0,1,lo,hi,gravz,muval,rhoval,nscal,model,corrector)

      call trace_p(s,DIMS(s),sn,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),wedge,DIMS(wedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),kappaz,DIMS(kpz),
     &     fluxx,DIMS(uedge),fluxy,DIMS(vedge),fluxz,DIMS(wedge),
     &     xlo,xhi,ylo,yhi,zlo,zhi,divu,tf,DIMS(work),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),    
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     dt,dx,lo,hi,bc,muval,rhoval,iconserv,
     &     nscal,model,gravz)

c ::: compute fluxes based on final edge states
      call riemann_solver_arg(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),eiglx,DIMS(work),    
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,0,lo,hi,gravx,muval,rhoval,nscal,model)
      call riemann_solver_arg(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),eigly,DIMS(work),    
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,0,lo,hi,gravy,muval,rhoval,nscal,model)
      call riemann_solver_arg(fluxz,wedge,zlo,zhi,DIMS(wedge),
     &     kappaz,DIMS(kpz),eiglz,DIMS(work),    
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,0,1,lo,hi,gravz,muval,rhoval,nscal,model)

      call flux_copy(xstate,DIMS(xstate),fluxx,DIMS(uedge),
     &     lo,hi,1,0,0,nscal)
      call flux_copy(ystate,DIMS(ystate),fluxy,DIMS(vedge),
     &     lo,hi,0,1,0,nscal)   
      call flux_copy(zstate,DIMS(zstate),fluxz,DIMS(wedge),
     &     lo,hi,0,0,1,nscal) 

      DEALLOCATE(fluxx,fluxy,fluxz)
      DEALLOCATE(xlo,xhi,ylo,yhi,zlo,zhi)

      end


      subroutine trace_lin(s,DIMS(s),sn,DIMS(sn),
     &     rphi,DIMS(rphi),
     &     sx,sy,sz,tf,DIMS(w),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     wedge,zlo,zhi,DIMS(wedge),
     &     dt,dx,lo,hi,bc,nscal)

c ::: Trace the cell centered states to edge for linear advection.
c     Valid for single-phase flow.

      implicit none

c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(rphi)
      integer DIMDEC(w)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      integer nscal
      REAL_T   s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  tf(DIMV(w),nscal)
      REAL_T  sx(DIMV(w),nscal) , sy(DIMV(w),nscal), sz(DIMV(w),nscal)
      REAL_T  uedge(DIMV(uedge)), xlo(DIMV(uedge),nscal), xhi(DIMV(uedge),nscal)
      REAL_T  vedge(DIMV(vedge)), ylo(DIMV(vedge),nscal), yhi(DIMV(vedge),nscal)
      REAL_T  wedge(DIMV(wedge)), zlo(DIMV(wedge),nscal), zhi(DIMV(wedge),nscal)
      REAL_T  dt, dx(SDIM)

c ::: local variable
      integer old_time
      integer nc,i,j,k
      REAL_T dth,dthx,dthy,dthz

      old_time = 0
      
      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)

c ::: x-direction
      do nc = 1,nscal
         do k = lo(3)-1,hi(3)+1
            do j = lo(2)-1,hi(2)+1
               do i = lo(1),hi(1)+1
                  xlo(i,j,k,nc) = s(i-1,j,k,nc) + 
     &                 (half-dthx*uedge(i,j,k)/rphi(i-1,j,k))*sx(i-1,j,k,nc) +
     &                 dth*tf(i-1,j,k,nc)
                  xhi(i,j,k,nc) = s(i,j,k,nc) + 
     &                 (-half-dthx*uedge(i,j,k)/rphi(i,j,k))*sx(i,j,k,nc) +
     &                 dth*tf(i,j,k,nc)
               end do
            end do
         end do
      end do

c ::: y-direction
      do nc = 1,nscal
         do k = lo(3)-1,hi(3)+1
            do j = lo(2),hi(2)+1
               do i = lo(1)-1,hi(1)+1
                  ylo(i,j,k,nc) = s(i,j-1,k,nc) + 
     &                 (half-dthy*vedge(i,j,k)/rphi(i,j-1,k))*sy(i,j-1,k,nc) +
     &                 dth*tf(i,j-1,k,nc)
                  yhi(i,j,k,nc) = s(i,j,k,nc) + 
     &                 (-half-dthy*vedge(i,j,k)/rphi(i,j,k))*sy(i,j,k,nc) +
     &                 dth*tf(i,j,k,nc)
               end do
            end do
         end do
      end do


c ::: z-direction
      do nc = 1,nscal
         do k = lo(3),hi(3)+1
            do j = lo(2)-1,hi(2)+1
               do i = lo(1)-1,hi(1)+1
                  zlo(i,j,k,nc) = s(i,j,k-1,nc) + 
     &                 (half-dthz*wedge(i,j,k)/rphi(i,j,k-1))*sz(i,j,k-1,nc) +
     &                 dth*tf(i,j,k-1,nc)
                  zhi(i,j,k,nc) = s(i,j,k,nc) + 
     &                 (-half-dthz*wedge(i,j,k)/rphi(i,j,k))*sz(i,j,k,nc) +
     &                 dth*tf(i,j,k,nc)
               end do
            end do
         end do
      end do

c ::: fix edges at the boundaries.
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,bc,1,nscal,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,bc,1,nscal,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,bc,1,nscal,old_time)

      end

      subroutine trace(s,DIMS(s),sn,DIMS(sn),
     &     eigvx,eiglx,eigrx,
     &     eigvy,eigly,eigry,
     &     eigvz,eiglz,eigrz,
     &     sx,sy,sz,divu,tf,DIMS(w),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     wedge,zlo,zhi,DIMS(wedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,nscal)

c ::: trace the cell centered states to edges.

      implicit none

c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(w)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer use_minion,iconserv
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      integer nscal
      REAL_T   s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  eigvx(DIMV(w),2*nscal)
      REAL_T  eigvy(DIMV(w),2*nscal)
      REAL_T  eigvz(DIMV(w),2*nscal)
      REAL_T  eiglx(DIMV(w),nscal,nscal),eigrx(DIMV(w),nscal,nscal)
      REAL_T  eigly(DIMV(w),nscal,nscal),eigry(DIMV(w),nscal,nscal)
      REAL_T  eiglz(DIMV(w),nscal,nscal),eigrz(DIMV(w),nscal,nscal)
      REAL_T  divu(DIMV(w)),tf(DIMV(w))
      REAL_T  sx(DIMV(w),nscal) , sy(DIMV(w),nscal), sz(DIMV(w),nscal)
      REAL_T  uedge(DIMV(uedge)), xlo(DIMV(uedge),nscal), xhi(DIMV(uedge),nscal)
      REAL_T  vedge(DIMV(vedge)), ylo(DIMV(vedge),nscal), yhi(DIMV(vedge),nscal)
      REAL_T  wedge(DIMV(wedge)), zlo(DIMV(wedge),nscal), zhi(DIMV(wedge),nscal)
      REAL_T  dt, dx(SDIM)

c ::: local variable
      integer old_time
      REAL_T dth,dthx,dthy,dthz

      old_time = 0
      
      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)

c ::: trace to edges in each directions
      call trace_dir(s,DIMS(s),eigvx,eiglx,eigrx,DIMS(w),
     &     xlo,xhi,DIMS(uedge),sx,tf,divu,DIMS(w),lo,hi,
     &     1,0,0,dthx,dth,use_minion,iconserv,nscal)
      call trace_dir(s,DIMS(s),eigvy,eigly,eigry,DIMS(w),
     &     ylo,yhi,DIMS(vedge),sy,tf,divu,DIMS(w),lo,hi,
     &     0,1,0,dthy,dth,use_minion,iconserv,nscal)
      call trace_dir(s,DIMS(s),eigvz,eiglz,eigrz,DIMS(w),
     &     zlo,zhi,DIMS(wedge),sz,tf,divu,DIMS(w),lo,hi,
     &     0,0,1,dthz,dth,use_minion,iconserv,nscal)

c ::: fix edges at the boundaries.
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,bc,1,nscal,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,bc,1,nscal,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,bc,1,nscal,old_time)
      end

      subroutine trace2(s,DIMS(s),sn,DIMS(sn),
     &     eigvx,eigrx,eigvy,eigry,eigvz,eigrz,DIMS(eig),
     &     sx,sy,sz,divu,tf,DIMS(w),uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),wedge,zlo,zhi,DIMS(wedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,ncomps,nphases) 
c ::: trace the cell centered states to edges.
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(w)
      integer DIMDEC(eig)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer ncomps,nphases
      integer use_minion,iconserv
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  sn(DIMV(sn),ncomps)
      REAL_T  divu(DIMV(w)),tf(DIMV(w))
      REAL_T  eigvx(DIMV(eig),ncomps)
      REAL_T  eigvy(DIMV(eig),ncomps)
      REAL_T  eigvz(DIMV(eig),ncomps)
      REAL_T  eigrx(DIMV(eig),ncomps,ncomps)
      REAL_T  eigry(DIMV(eig),ncomps,ncomps)
      REAL_T  eigrz(DIMV(eig),ncomps,ncomps)
      REAL_T  sx(DIMV(w),ncomps)
      REAL_T  sy(DIMV(w),ncomps)
      REAL_T  sz(DIMV(w),ncomps)
      REAL_T  uedge(DIMV(uedge))
      REAL_T  vedge(DIMV(vedge))
      REAL_T  wedge(DIMV(wedge))
      REAL_T  xlo(DIMV(uedge),ncomps), xhi(DIMV(uedge),ncomps)
      REAL_T  ylo(DIMV(vedge),ncomps), yhi(DIMV(vedge),ncomps)
      REAL_T  zlo(DIMV(wedge),ncomps), zhi(DIMV(wedge),ncomps)
      REAL_T  dt, dx(SDIM)

c ::: local variable
      integer old_time
      REAL_T dth,dthx,dthy,dthz

      old_time  = 0
      
      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)

c ::: trace to edges in each directions
      call trace_dir2(s,DIMS(s),eigvx,eigrx,DIMS(eig),
     &     xlo,xhi,DIMS(uedge),sx,tf,divu,DIMS(w),lo,hi,
     &     1,0,0,dthx,dth,use_minion,iconserv,ncomps,nphases)
      call trace_dir2(s,DIMS(s),eigvy,eigry,DIMS(eig),
     &     ylo,yhi,DIMS(vedge),sy,tf,divu,DIMS(w),lo,hi,
     &     0,1,0,dthy,dth,use_minion,iconserv,ncomps,nphases)
      call trace_dir2(s,DIMS(s),eigvz,eigrz,DIMS(eig),
     &     zlo,zhi,DIMS(wedge),sz,tf,divu,DIMS(w),lo,hi,
     &     0,0,1,dthz,dth,use_minion,iconserv,ncomps,nphases)

c ::: fix edges at the boundaries.
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,bc,1,ncomps,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,bc,1,ncomps,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,bc,1,ncomps,old_time)

      end
      
      subroutine trace_p_lin(s,DIMS(s),sn,DIMS(sn),
     &     fx,uedge,DIMS(uedge),fy,vedge,DIMS(vedge),fz,wedge,DIMS(wedge),
     &     rphi,DIMS(rphi),xlo,xhi,ylo,yhi,zlo,zhi,
     &     dt,dx,lo,hi,bc,nscal)

      implicit none
c
c ::: correct the edge values during the predictor stage
c
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(rphi)
      integer bc(SDIM,2), lo(SDIM), hi(SDIM)
      integer nscal  
      REAL_T  s(DIMV(s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  uedge(DIMV(uedge)),vedge(DIMV(vedge)),wedge(DIMV(wedge))
      REAL_T  fx(DIMV(uedge),nscal),fy(DIMV(vedge),nscal),fz(DIMV(wedge),nscal)  
      REAL_T  xlo(DIMV(uedge),nscal),xhi(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal),yhi(DIMV(vedge),nscal)
      REAL_T  zlo(DIMV(wedge),nscal),zhi(DIMV(wedge),nscal)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  dt,dx(SDIM)

c ::: internal variables
      integer i,j,k,nc
      integer half_time
      REAL_T  dth,xhx,xhy,xhz
      REAL_T  st,drmndv

      half_time = 1

      dth  = half*dt
      xhx  = one/dx(1)
      xhy  = one/dx(2)
      xhz  = one/dx(3)

c ::: x-edges
!$omp parallel private(i,j,k,nc,st)
      do nc = 1,nscal
!$omp do
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1)-1,hi(1)+1
                  drmndv = s(i,j,k,nc)/rphi(i,j,k)
                  st = - (fy(i,j+1,k,nc)-fy(i,j,k,nc))*xhy/rphi(i,j,k)
     &                 - (fz(i,j,k+1,nc)-fz(i,j,k,nc))*xhz/rphi(i,j,k) 
     &                 + drmndv*(vedge(i,j+1,k)-vedge(i,j,k))*xhy
     &                 + drmndv*(wedge(i,j,k+1)-wedge(i,j,k))*xhz

                  xlo(i+1,j,k,nc) = xlo(i+1,j,k,nc) + dth*st
                  xhi(i,  j,k,nc) = xhi(i,  j,k,nc) + dth*st
               end do
            end do         
         end do
!$omp end do nowait
      end do
!$omp end parallel

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,bc,1,nscal,half_time)

c ::: y-edges
!$omp parallel private(i,j,k,nc,st)
      do nc = 1,nscal
!$omp do
         do k = lo(3),hi(3)
            do j = lo(2)-1,hi(2)+1
               do i = lo(1),hi(1)
                  drmndv = s(i,j,k,nc)/rphi(i,j,k)
                  st = - (fx(i+1,j,k,nc)-fx(i,j,k,nc))*xhx/rphi(i,j,k)
     &                 - (fz(i,j,k+1,nc)-fz(i,j,k,nc))*xhz/rphi(i,j,k)
     &                 + drmndv*(uedge(i+1,j,k)-uedge(i,j,k))*xhx
     &                 + drmndv*(wedge(i,j,k+1)-wedge(i,j,k))*xhz

                  ylo(i,j+1,k,nc) = ylo(i,j+1,k,nc) + dth*st
                  yhi(i,j,  k,nc) = yhi(i,j,  k,nc) + dth*st 
               end do
            end do            
         end do  
!$omp end do nowait
      end do
!$omp end parallel

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,bc,1,nscal,half_time)

c ::: z-edges
      do nc = 1,nscal
         do k = lo(3)-1,hi(3)+1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  drmndv = s(i,j,k,nc)/rphi(i,j,k)
                  st = - (fx(i+1,j,k,nc)-fx(i,j,k,nc))*xhx/rphi(i,j,k)
     &                 - (fy(i,j+1,k,nc)-fy(i,j,k,nc))*xhy/rphi(i,j,k)
     &                 + drmndv*(uedge(i+1,j,k)-uedge(i,j,k))*xhx
     &                 + drmndv*(vedge(i,j+1,k)-vedge(i,j,k))*xhy

                  zlo(i,j,k+1,nc) = zlo(i,j,k+1,nc) + dth*st
                  zhi(i,j,k,  nc) = zhi(i,j,k,  nc) + dth*st 
               end do
            end do            
         end do  
      end do

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,bc,1,nscal,half_time)

      end

      subroutine trace_p(s,DIMS(s),sn,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),wedge,DIMS(wedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),kappaz,DIMS(kpz),
     &     fx,DIMS(fx),fy,DIMS(fy),fz,DIMS(fz),
     &     xlo,xhi,ylo,yhi,zlo,zhi,divu,tf,DIMS(w),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     krcoef, DIMS(krcoef), nkrcoef,
     &     dt,dx,lo,hi,bc,mu,rhoval,iconserv,nscal,model,grav)

      implicit none

c
c ::: correct the edge values during the predictor stage
c
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(fx), DIMDEC(fy), DIMDEC(fz)
      integer DIMDEC(w)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa),DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer bc(SDIM,2), lo(SDIM), hi(SDIM)
      integer nscal,model
      integer iconserv  
      REAL_T  s(DIMV(s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  uedge(DIMV(uedge)),vedge(DIMV(vedge)),wedge(DIMV(wedge))
      REAL_T  fx(DIMV(fx),nscal),fy(DIMV(fy),nscal),fz(DIMV(fz),nscal)  
      REAL_T  kappax(DIMV(kpx)),kappay(DIMV(kpy)),kappaz(DIMV(kpz))
      REAL_T  xlo(DIMV(uedge),nscal),xhi(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal),yhi(DIMV(vedge),nscal)
      REAL_T  zlo(DIMV(wedge),nscal),zhi(DIMV(wedge),nscal)
      REAL_T  divu(DIMV(w)), tf(DIMV(w))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa)) 
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  mu(nscal),rhoval(nscal)
      REAL_T  dt,dx(SDIM)
      REAL_T  grav

c ::: internal variables
      integer i,j,k,nc
      integer imin,imax,jmin,jmax,kmin,kmax
      integer half_time
      REAL_T  dth,xhx,xhy,xhz
      REAL_T  st

      REAL_T, ALLOCATABLE::drmndv(:,:,:,:)
      REAL_T, ALLOCATABLE::drmndk(:,:,:,:)
      
      allocate(drmndv(DIMV(w),nscal),drmndk(DIMV(w),nscal))

      half_time = 1

      dth  = half*dt
      xhx  = one/dx(1)
      xhy  = one/dx(2)
      xhz  = one/dx(3)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)     
      kmax = hi(3)

      if (iconserv.ne.1) then
         call bl_pd_abort("trans_p: iconserv not equal to 1")
      end if

      call dfluxdv(drmndv,drmndk,DIMS(w),s,DIMS(s),rphi,DIMS(rphi),
     &     kappa,DIMS(kappa),krcoef,DIMS(krcoef),nkrcoef,
     &     mu,rhoval,nscal,model,grav)

c ::: x-edges
!$omp parallel private(i,j,k,nc,st)
      do nc = 1,nscal
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin-1,imax+1
                  st = - (fy(i,j+1,k,nc)-fy(i,j,k,nc))*xhy/rphi(i,j,k)
     &                 - (fz(i,j,k+1,nc)-fz(i,j,k,nc))*xhz/rphi(i,j,k) 
     &                 + drmndv(i,j,k,nc)*(vedge(i,j+1,k)-vedge(i,j,k))*xhy
     &                 + drmndv(i,j,k,nc)*(wedge(i,j,k+1)-wedge(i,j,k))*xhz

                  xlo(i+1,j,k,nc) = xlo(i+1,j,k,nc) + dth*st
                  xhi(i,  j,k,nc) = xhi(i,  j,k,nc) + dth*st
               end do
            end do         
         end do
!$omp end do nowait
      end do
!$omp end parallel

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,bc,1,nscal,half_time)

c ::: y-edges
!$omp parallel private(i,j,k,nc,st)
      do nc = 1,nscal
!$omp do
         do k = kmin,kmax
            do j = jmin-1,jmax+1
               do i = imin,imax
                  st = - (fx(i+1,j,k,nc)-fx(i,j,k,nc))*xhx/rphi(i,j,k)
     &                 - (fz(i,j,k+1,nc)-fz(i,j,k,nc))*xhz/rphi(i,j,k)
     &                 + drmndv(i,j,k,nc)*(uedge(i+1,j,k)-uedge(i,j,k))*xhx
     &                 + drmndv(i,j,k,nc)*(wedge(i,j,k+1)-wedge(i,j,k))*xhz

                  ylo(i,j+1,k,nc) = ylo(i,j+1,k,nc) + dth*st
                  yhi(i,j,  k,nc) = yhi(i,j,  k,nc) + dth*st 
               end do
            end do            
         end do  
!$omp end do nowait
      end do
!$omp end parallel

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,bc,1,nscal,half_time)

c ::: z-edges
      do nc = 1,nscal
         do k = kmin-1,kmax+1
            do j = jmin,jmax
               do i = imin,imax
                  st = - (fx(i+1,j,k,nc)-fx(i,j,k,nc))*xhx/rphi(i,j,k)
     &                 - (fy(i,j+1,k,nc)-fy(i,j,k,nc))*xhy/rphi(i,j,k)
     &                 - drmndk(i,j,k,nc)*(kappaz(i,j,k+1)-kappaz(i,j,k))*xhz
     &                 + drmndv(i,j,k,nc)*(uedge(i+1,j,k)-uedge(i,j,k))*xhx
     &                 + drmndv(i,j,k,nc)*(vedge(i,j+1,k)-vedge(i,j,k))*xhy

                  zlo(i,j,k+1,nc) = zlo(i,j,k+1,nc) + dth*st
                  zhi(i,j,k,  nc) = zhi(i,j,k,  nc) + dth*st 
               end do
            end do            
         end do  
      end do

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,bc,1,nscal,half_time)

      deallocate(drmndv,drmndk)

      end

      subroutine trace_p2(s,DIMS(s),sn,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),wedge,DIMS(wedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),kappaz,DIMS(kpz),
     &     fx,DIMS(fx),fy,DIMS(fy),fz,DIMS(fz),
     &     xlo,xhi,ylo,yhi,zlo,zhi,divu,tf,DIMS(w),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     dt,dx,lo,hi,bc,
     &     mu,rho,athick,iconserv,ncomps,nphases,grav)
c
c ::: correct the edge values during the predictor stage
c
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(uedge),DIMDEC(vedge),DIMDEC(wedge)
      integer DIMDEC(fx),DIMDEC(fy),DIMDEC(fz)
      integer DIMDEC(w)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa),DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer bc(SDIM,2),lo(SDIM),hi(SDIM)
      integer ncomps,nphases,iconserv
      REAL_T   s(DIMV(s),ncomps)
      REAL_T  sn(DIMV(sn),ncomps)
      REAL_T  uedge(DIMV(uedge)),vedge(DIMV(vedge)),wedge(DIMV(wedge))
      REAL_T  kappax(DIMV(kpx)),kappay(DIMV(kpy)),kappaz(DIMV(kpz))
      REAL_T  fx(DIMV(fx),ncomps),fy(DIMV(fy),ncomps),fz(DIMV(fz),ncomps)
      REAL_T  xlo(DIMV(uedge),ncomps),xhi(DIMV(uedge),ncomps)
      REAL_T  ylo(DIMV(vedge),ncomps),yhi(DIMV(vedge),ncomps)
      REAL_T  zlo(DIMV(wedge),ncomps),zhi(DIMV(wedge),ncomps)
      REAL_T  divu(DIMV(w)), tf(DIMV(w))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  mu(nphases),rho(nphases),athick
      REAL_T  dt,dx(SDIM)
      REAL_T  grav
      
c ::: internal variables
      integer i,j,k,nc
      integer imin,imax,jmin,jmax,kmin,kmax
      integer half_time
      REAL_T  dth,hx,hy,hz
      REAL_T  st,u(2),drmndv(2),drmndk(2)

      half_time = 1
      
      dth  = half*dt
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)-1
      jmin = lo(2)-1
      kmin = lo(3)-1
      imax = hi(1)+1
      jmax = hi(2)+1
      kmax = hi(3)+1

      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               do nc = 1,ncomps
                  u(nc) = s(i,j,k,nc)
               end do

               call dfluxdv2(drmndv,drmndk,u,rphi(i,j,k),kappa(i,j,k),
     &              mu(1),mu(2),athick,rho(1),rho(2),grav)
               do nc = 1,ncomps
                  if (iconserv.eq.1) then
                     st = - (fy(i,j+1,k,nc)-fy(i,j,k,nc))/hy/rphi(i,j,k)
     &                    - (fz(i,j,k+1,nc)-fz(i,j,k,nc))/hz/rphi(i,j,k)
     &                    - drmndv(nc)*(uedge(i+1,j,k)-uedge(i,j,k))/hx
                     xlo(i+1,j,k,nc) = xlo(i+1,j,k,nc) + dth*st
                     xhi(i,j,k,nc)   = xhi(i,j,k,nc)   + dth*st

                     st = - (fx(i+1,j,k,nc)-fx(i,j,k,nc))/hx/rphi(i,j,k)
     &                    - (fz(i,j,k+1,nc)-fz(i,j,k,nc))/hz/rphi(i,j,k)
     &                    - drmndv(nc)*(vedge(i,j+1,k)-vedge(i,j,k))/hy
                     ylo(i,j+1,k,nc) = ylo(i,j+1,k,nc) + dth*st
                     yhi(i,j,k,nc)   = yhi(i,j,k,nc)   + dth*st


                     st = - (fx(i+1,j,k,nc)-fx(i,j,k,nc))/hx/rphi(i,j,k)
     &                    - (fy(i,j+1,k,nc)-fy(i,j,k,nc))/hy/rphi(i,j,k)
     &                    - drmndk(nc)*(kappaz(i,j,k+1)-kappaz(i,j,k))/hz
     &                    - drmndv(nc)*(wedge(i,j,k+1)-wedge(i,j,k))/hz 
                     zlo(i,j,k+1,nc) = zlo(i,j,k+1,nc) + dth*st
                     zhi(i,j,k,nc)   = zhi(i,j,k,nc)   + dth*st

                  endif 
               end do
            end do         
         end do
      end do

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,bc,1,ncomps,half_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,bc,1,ncomps,half_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,bc,1,ncomps,half_time)

      end

      subroutine trace_tracers(s,DIMS(s),sn,DIMS(sn),
     &     st,DIMS(st),stn,DIMS(stn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),wedge,DIMS(wedge),
     &     phi,DIMS(phi),sx,sy,sz,DIMS(w),
     &     xlo,xhi,ylo,yhi,zlo,zhi,
     &     dt,dx,lo,hi,bc,nscal) 

c ::: Trace the cell centered states to edges for tracers.
c ::: Assume linear advection based on the phase velocity.

      implicit none

c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(st)
      integer DIMDEC(stn)
      integer DIMDEC(w)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(wedge)
      integer DIMDEC(phi)
      integer nscal
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      REAL_T  s(DIMV(s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  st(DIMV(st)),stn(DIMV(stn))
      REAL_T  uedge(DIMV(uedge))
      REAL_T  vedge(DIMV(vedge))
      REAL_T  wedge(DIMV(wedge))
      REAL_T  phi(DIMV(phi))
      REAL_T  sx(DIMV(w),nscal)
      REAL_T  sy(DIMV(w),nscal)
      REAL_T  sz(DIMV(w),nscal)
      REAL_T  xlo(DIMV(uedge),nscal), xhi(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal), yhi(DIMV(vedge),nscal)
      REAL_T  zlo(DIMV(wedge),nscal), zhi(DIMV(wedge),nscal)
      REAL_T  dt, dx(SDIM)

c ::: local variable
      integer old_time

      REAL_T dth,dthx,dthy,dthz

      old_time = 0
      
      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)

c ::: trace to edges in each directions
      call trace_dir_tracer(s,DIMS(s),st,DIMS(st),stn,DIMS(stn),
     &     uedge,DIMS(uedge),phi,DIMS(phi),
     &     xlo,xhi,DIMS(uedge),sx,DIMS(w),lo,hi,
     &     1,0,0,dthx,nscal)
      call trace_dir_tracer(s,DIMS(s),st,DIMS(st),stn,DIMS(stn),
     &     vedge,DIMS(vedge),phi,DIMS(phi),
     &     ylo,yhi,DIMS(vedge),sy,DIMS(w),lo,hi,
     &     0,1,0,dthy,nscal)
      call trace_dir_tracer(s,DIMS(s),st,DIMS(st),stn,DIMS(stn),
     &     wedge,DIMS(wedge),phi,DIMS(phi),
     &     zlo,zhi,DIMS(wedge),sz,DIMS(w),lo,hi,
     &     0,0,1,dthz,nscal)

c ::: fix edges at the boundaries.
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,0,lo,hi,bc,1,nscal,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,0,lo,hi,bc,1,nscal,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),wedge,zlo,zhi,DIMS(wedge),
     &     0,0,1,lo,hi,bc,1,nscal,old_time)

      end

      subroutine trace_dir(s,DIMS(s),eigv,eigl,eigr,DIMS(eig),
     &     vlo,vhi,DIMS(v),dv,tf,divu,DIMS(w),lo,hi,
     &     ix,iy,iz,dthx,dth,use_minion,iconserv,nscal)

c ::: arguments
      integer DIMDEC(s), DIMDEC(eig), DIMDEC(w), DIMDEC(v)
      integer lo(SDIM),hi(SDIM)
      integer ix,iy,iz,nscal,use_minion,iconserv
      REAL_T  s(DIMV(s),nscal)
      REAL_T  eigv(DIMV(eig),2*nscal)
      REAL_T  eigl(DIMV(eig),nscal,nscal)
      REAL_T  eigr(DIMV(eig),nscal,nscal)
      REAL_T  vlo(DIMV(v),nscal),vhi(DIMV(v),nscal)
      REAL_T  dv(DIMV(w),nscal)
      REAL_T  tf(DIMV(w)),divu(DIMV(w))
      REAL_T  dthx, dth

c ::: internal variables
      integer i,j,k,nc,nd
      integer imin,jmin,kmin,imax,jmax,kmax
      REAL_T  alphalo, alphahi

      imin = lo(1) - 1
      imax = hi(1) + 1
      jmin = lo(2) - 1
      jmax = hi(2) + 1
      kmin = lo(3) - 1
      kmax = hi(3) + 1
      if (ix .eq. 1) then
         imin = lo(1)
      elseif (iy .eq. 1) then
         jmin = lo(2) 
      elseif (iz .eq. 1) then
         kmin = lo(3) 
      end if

!$omp parallel do private(i,j,k,nc,nd,alphalo,alphahi)
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               do nc = 1,nscal
                  vlo(i,j,k,nc) = s(i-ix,j-iy,k-iz,nc)
                  vhi(i,j,k,nc) = s(i,j,k,nc)
               end do
               do nd = 1,nscal
                  alphalo = zero
                  alphahi = zero
                  do nc = 1,nscal
                     alphalo = alphalo + eigl(i-ix,j-iy,k-iz,nd,nc)*
     &                    dv(i-ix,j-iy,k-iz,nc)
                     alphahi = alphahi + eigl(i,j,k,nd,nc)*dv(i,j,k,nc)
                  end do               
                  do nc = 1,nscal                
                     vlo(i,j,k,nc) = vlo(i,j,k,nc) 
     &                 + ( half - dthx*max(zero,eigv(i-ix,j-iy,k-iz,nd)))
     &                 *alphalo*eigr(i-ix,j-iy,k-iz,nc,nd)
                     vhi(i,j,k,nc) = vhi(i,j,k,nc) + (-half
     &                 - dthx*min(zero,eigv(i,j,k,nd)))
     &                 *alphahi*eigr(i,j,k,nc,nd)
                  end do
               end do
            end do
         end do
      end do
!$omp end parallel do

      use_minion = 1

      if (use_minion.eq.1)then
!$omp parallel private(i,j,k,nc)
         do nc = 1, nscal
!$omp do
            do k = kmin,kmax
               do j = jmin,jmax
                  do i = imin,imax
                     vlo(i,j,k,nc) = vlo(i,j,k,nc) + dth*tf(i-ix,j-iy,k-iz)
                     vhi(i,j,k,nc) = vhi(i,j,k,nc) + dth*tf(i,j,k)
                  end do
               end do
            end do
!$omp end do nowait
         end do
!$omp end parallel

         if (iconserv .eq. 1) then
!$omp parallel private(i,j,k,nc)
            do nc = 1, nscal
!$omp do
               do k = kmin,kmax
                  do j = jmin,jmax
                     do i = imin,imax
                        vlo(i,j,k,nc) = vlo(i,j,k,nc) - dth*s(i,j,k,nc)*divu(i-ix,j-iy,k-iz)
                        vhi(i,j,k,nc) = vhi(i,j,k,nc) - dth*s(i,j,k,nc)*divu(i,j,k)
                     end do
                  end do
               end do
!$omp end do nowait
            end do
!$omp end parallel
         end if
      end if

      end

      subroutine trace_dir_tracer(s,DIMS(s),
     &     st,DIMS(st),stn,DIMS(stn),
     &     uedge,DIMS(uedge),phi,DIMS(phi),
     &     vlo,vhi,DIMS(v),dv,DIMS(w),lo,hi,
     &     ix,iy,iz,dthx,nscal)

      implicit none

c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(st)
      integer DIMDEC(stn)
      integer DIMDEC(uedge)
      integer DIMDEC(phi)
      integer DIMDEC(w)
      integer DIMDEC(v)
      integer lo(SDIM),hi(SDIM)
      integer ix,iy,iz,nscal
      REAL_T  s(DIMV(s),nscal)
      REAL_T  st(DIMV(st)),stn(DIMV(stn))
      REAL_T  uedge(DIMV(uedge))
      REAL_T  phi(DIMV(phi))
      REAL_T  vlo(DIMV(v),nscal)
      REAL_T  vhi(DIMV(v),nscal)
      REAL_T  dv(DIMV(w),nscal)
      REAL_T  dthx

c ::: internal variables
      integer i,j,k,nc
      integer imin,jmin,kmin,imax,jmax,kmax
      REAL_T  alphalo,alphahi,sthalf

      imin = lo(1) - 1
      imax = hi(1) + 1
      jmin = lo(2) - 1
      jmax = hi(2) + 1
      kmin = lo(3) - 1
      kmax = hi(3) + 1

!$omp parallel private(i,j,k,nc,sthalf,alphalo,alphahi)
      do nc = 1,nscal
!$omp do
         do k = kmin+iz,kmax
            do j = jmin+iy,jmax
               do i = imin+ix,imax

                  vlo(i,j,k,nc) = s(i-ix,j-iy,k-iz,nc)
                  vhi(i,j,k,nc) = s(i,j,k,nc)

                  sthalf  = half*(st(i-ix,j-iy,k-iz)+stn(i-ix,j-iy,k-iz))
                  alphalo = uedge(i,j,k)/(phi(i-ix,j-iy,k-iz)*sthalf)
                  sthalf  = half*(st(i,j,k)+stn(i,j,k))
                  alphahi = uedge(i,j,k)/(phi(i,j,k)*sthalf)

                  vlo(i,j,k,nc) = vlo(i,j,k,nc) + (half
     &                 - dthx*max(zero,alphalo))*dv(i-ix,j-iy,k-iz,nc)
                  vhi(i,j,k,nc) = vhi(i,j,k,nc) + (-half
     &                 - dthx*min(zero,alphahi))*dv(i,j,k,nc)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end
 
      subroutine trace_dir2(s,DIMS(s),eigv,eigr,DIMS(eig),
     &     vlo,vhi,DIMS(v),dv,tf,divu,DIMS(w),lo,hi,
     &     ix,iy,iz,dthx,dth,use_minion,iconserv,ncomps,nphases)

c ::: arguments
      integer DIMDEC(s), DIMDEC(eig), DIMDEC(w), DIMDEC(v)
      integer lo(SDIM),hi(SDIM)
      integer ix,iy,iz,ncomps,nphases,use_minion,iconserv
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  eigv(DIMV(eig),ncomps)
      REAL_T  eigr(DIMV(eig),ncomps,ncomps) 
      REAL_T  vlo(DIMV(v),ncomps),vhi(DIMV(v),ncomps)
      REAL_T  dv(DIMV(w),ncomps)
      REAL_T  tf(DIMV(w)),divu(DIMV(w))
      REAL_T  dthx, dth

c ::: internal variables
      integer i,j,k,nc,nd
      integer imin,jmin,kmin,imax,jmax,kmax
      REAL_T  alphalo, alphahi

      if (ix .eq. 1) then
         imin = lo(1)
         imax = hi(1) + 1
         jmin = lo(2) - 1
         jmax = hi(2) + 1
         kmin = lo(3) - 1
         kmax = hi(3) + 1
      elseif (iy .eq. 1) then
         imin = lo(1) - 1
         imax = hi(1) + 1
         jmin = lo(2) 
         jmax = hi(2) + 1
         kmin = lo(3) - 1
         kmax = hi(3) + 1
      elseif (iz .eq. 1) then
         imin = lo(1) - 1
         imax = hi(1) + 1
         jmin = lo(2) - 1
         jmax = hi(2) + 1
         kmin = lo(3) 
         kmax = hi(3) + 1
      end if
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax  
               do nc = 1,ncomps
                  vlo(i,j,k,nc) = s(i-ix,j-iy,k-iz,nc)
                  vhi(i,j,k,nc) = s(i,j,k,nc)
               end do
 
               do nd = 1,ncomps     
                  alphalo = cvmgt(zero, one, eigv(i-ix,j-iy,k-iz,nd) .lt. zero)
                  alphahi = cvmgt(zero, one, eigv(i,j,k,nd) .gt. zero)
                  do nc = 1,ncomps   
                     vlo(i,j,k,nc) = vlo(i,j,k,nc) + half*
     &                    (alphalo*(one - 2*dthx*eigv(i-ix,j-iy,k-iz,nd)))
     &                    *dv(i-ix,j-iy,k-iz,nd)*eigr(i-ix,j-iy,k-iz,nc,nd)
                     vhi(i,j,k,nc) = vhi(i,j,k,nc) - half*
     &                    (alphahi*(one + 2*dthx*eigv(i,j,k,nd)))
     &                    *dv(i,j,k,nd)*eigr(i,j,k,nc,nd)      
                  end do
               end do
            end do
         end do
      end do

      end

      subroutine trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vel,vlo,vhi,DIMS(v),
     &     ix,iy,iz,lo,hi,bc,idx_first,nscal,which_time)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the specified direction.
c     (ix,iy,iz) is the unit vector used to specify direction.  
c
      implicit none

#include "probdata.H"

c ::: arguments

      integer DIMDEC(s), DIMDEC(sn), DIMDEC(v)
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      integer ix,iy,iz,idx_first,nscal, which_time
      REAL_T  s(DIMV( s),nscal)
      REAL_T sn(DIMV(sn),nscal)
      REAL_T vel(DIMV(v))
      REAL_T vlo(DIMV(v),nscal)
      REAL_T vhi(DIMV(v),nscal)
      
c ::: local variables

      integer i,j,k,nc
      integer ibc,imin,jmin,kmin,imax,jmax,kmax
      integer old_time, half_time

      old_time  = 0
      half_time = 1
      
      ibc  = ix + 2*iy + 3*iz
      imin = (lo(1)-1)*(1-ix) + lo(1)*ix
      imax = (hi(1)+1)*(1-ix) + lo(1)*ix
      jmin = (lo(2)-1)*(1-iy) + lo(2)*iy
      jmax = (hi(2)+1)*(1-iy) + lo(2)*iy
      kmin = (lo(3)-1)*(1-iz) + lo(3)*iz
      kmax = (hi(3)+1)*(1-iz) + lo(3)*iz 
         
      if (bc(ibc,1).eq.EXT_DIR) then
         if (model .eq. 0 .or. model .eq. 1) then
            if (which_time .eq. old_time) then
               do nc = idx_first,nscal
                  do k = kmin,kmax
                     do j = jmin,jmax
                        do i = imin,imax
                           if (vel(i,j,k) .gt. zero) then
                              vlo(i,j,k,nc) = s(i-ix,j-iy,k-iz,nc)
                           else
                              vlo(i,j,k,nc) = vhi(i,j,k,nc)
                           end if
                        end do
                     end do
                  end do
               end do
            else if (which_time .eq. half_time) then
               do nc = idx_first,nscal
                  do k = kmin,kmax
                     do j = jmin,jmax
                        do i = imin,imax
                           if (vel(i,j,k) .gt. zero) then
                              vlo(i,j,k,nc) = half * ( s(i-ix,j-iy,k-iz,nc) + 
     &                             sn(i-ix,j-iy,k-iz,nc) ) 
                           else
                              vlo(i,j,k,nc) = vhi(i,j,k,nc)
                           end if
                        end do
                     end do
                  end do
               end do
            else
               print *,'BAD WHICH_TIME IN TRANS_BC_RMN ',which_time
               stop
            end if

         elseif (model .eq. 2 .or. model .eq. 3) then
            if (which_time .eq. old_time) then
               do nc = idx_first,nscal
                  do k = kmin,kmax
                     do j = jmin,jmax
                        do i = imin,imax
                           if (vel(i,j,k) .gt. zero) then
                              vlo(i,j,k,nc) = s(i-ix,j-iy,k-iz,nc)
                           else
                              vlo(i,j,k,nc) = vhi(i,j,k,nc)
                           end if
                        end do
                     end do
                  end do
               end do
            else if (which_time .eq. half_time) then
               do nc = idx_first,nscal
                  do k = kmin,kmax
                     do j = jmin,jmax
                        do i = imin,imax
                           if (vel(i,j,k) .gt. zero) then
                              vlo(i,j,k,nc) = half * ( s(i-ix,j-iy,k-iz,nc) + 
     &                             sn(i-ix,j-iy,k-iz,nc) ) 
                           else
                              vlo(i,j,k,nc) = vhi(i,j,k,nc)
                           end if
                        end do
                     end do
                  end do
               end do
            else
               print *,'BAD WHICH_TIME IN TRANS_BC_RMN ',which_time
               stop
            end if
         end if

      else if (bc(ibc,1).eq.SEEPAGE .or. bc(ibc,1).eq.FOEXTRAP.or.
     &        bc(ibc,1).eq.HOEXTRAP .or. bc(ibc,1).eq.REFLECT_EVEN) then
         do nc = idx_first,nscal
            do k = kmin,kmax
               do j = jmin,jmax
                  do i = imin,imax
                     vlo(i,j,k,nc) = vhi(i,j,k,nc)
                  end do
               end do
            end do
         end do
      else if (bc(ibc,1).eq.REFLECT_ODD) then
         do nc = idx_first,nscal
            do k = kmin,kmax
               do j = jmin,jmax
                  do i = imin,imax
                     vhi(i,j,k,nc) = zero
                     vlo(i,j,k,nc) = zero
                  end do
               end do
            end do
         end do
      end if

      imin = (lo(1)-1)*(1-ix) + (hi(1)+1)*ix
      imax = (hi(1)+1)*(1-ix) + (hi(1)+1)*ix
      jmin = (lo(2)-1)*(1-iy) + (hi(2)+1)*iy
      jmax = (hi(2)+1)*(1-iy) + (hi(2)+1)*iy
      kmin = (lo(3)-1)*(1-iz) + (hi(3)+1)*iz
      kmax = (hi(3)+1)*(1-iz) + (hi(3)+1)*iz 

      if (bc(ibc,2).eq.EXT_DIR) then
         if (model .eq. 0 .or. model .eq. 1) then
            if (which_time .eq. old_time) then
               do nc = idx_first,nscal
                  do k = kmin,kmax
                     do j = jmin,jmax
                        do i = imin,imax
                           if (vel(i,j,k) .gt. zero) then
                              vhi(i,j,k,nc) = s(i,j,k,nc)
                           else
                              vhi(i,j,k,nc) = vlo(i,j,k,nc)
                           end if
                        end do
                     end do
                  end do
               end do
            else if (which_time .eq. half_time) then
               do nc = idx_first,nscal
                  do k = kmin,kmax
                     do j = jmin,jmax
                        do i = imin,imax
                           if (vel(i,j,k) .gt. zero) then
                              vhi(i,j,k,nc) = half * ( s(i,j,k,nc) + 
     &                             sn(i,j,k,nc) ) 
                           else
                              vhi(i,j,k,nc) = vlo(i,j,k,nc)
                           end if
                        end do
                     end do
                  end do
               end do
            else
               print *,'BAD WHICH_TIME IN TRANS_BC_RMN ',which_time
               stop
            end if

         elseif (model .eq. 2 .or. model .eq. 3) then
            if (which_time .eq. old_time) then
               do nc = idx_first,nscal
                  do k = kmin,kmax
                     do j = jmin,jmax
                        do i = imin,imax
                           if (vel(i,j,k) .gt. zero) then
                              vhi(i,j,k,nc) = s(i,j,k,nc)
                           else
                              vhi(i,j,k,nc) = vhi(i,j,k,nc)
                           end if
                        end do
                     end do
                  end do
               end do
            else if (which_time .eq. half_time) then
               do nc = idx_first,nscal
                  do k = kmin,kmax
                     do j = jmin,jmax
                        do i = imin,imax
                           if (vel(i,j,k) .gt. zero) then
                              vhi(i,j,k,nc) = half * ( s(i,j,k,nc) + 
     &                             sn(i,j,k,nc) ) 
                           else
                              vhi(i,j,k,nc) = vlo(i,j,k,nc)
                           end if
                        end do
                     end do
                  end do
               end do
            else
               print *,'BAD WHICH_TIME IN TRANS_BC_RMN ',which_time
               stop
            end if
         end if

      else if (bc(ibc,2).eq.SEEPAGE.or.bc(ibc,2).eq.FOEXTRAP.or.
     &         bc(ibc,2).eq.HOEXTRAP.or.bc(ibc,2).eq.REFLECT_EVEN) then
         do nc = idx_first,nscal
            do k = kmin,kmax
               do j = jmin,jmax
                  do i = imin,imax
                     vhi(i,j,k,nc) = vlo(i,j,k,nc)
                  end do
               end do
            end do
         end do

      else if (bc(ibc,2).eq.REFLECT_ODD) then
         do nc = idx_first,nscal
            do k = kmin,kmax
               do j = jmin,jmax
                  do i = imin,imax
                     vhi(i,j,k,nc) = zero
                     vlo(i,j,k,nc) = zero
                  end do
               end do
            end do
         end do
      end if

      end

      subroutine riemann_solver_lin(flux,vel,vlo,vhi,DIMS(v),
     &     ix,iy,iz,lo,hi,nscal)

c ::: Riemann solver for linear advection.  
c     valid for single-phase flow.

      implicit none

      integer DIMDEC(v)
      integer ix,iy,iz
      integer nscal
      integer lo(SDIM),hi(SDIM)
      REAL_T  flux(DIMV(v),nscal)
      REAL_T  vlo(DIMV(v),nscal), vhi(DIMV(v),nscal)
      REAL_T  vel(DIMV(v))

c ::: internal variables
      integer i,j,k,nc    
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1) - 1
      imax = hi(1) + 1
      jmin = lo(2) - 1
      jmax = hi(2) + 1
      kmin = lo(3) - 1
      kmax = hi(3) + 1
      if (ix .eq. 1) then
         imin = lo(1)
      elseif (iy .eq. 1) then
         jmin = lo(2) 
      elseif (iz .eq. 1) then
         kmin = lo(3) 
      end if

!$omp parallel private(i,j,k,nc)
      do nc = 1,nscal
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  flux(i,j,k,nc) = cvmgp(vlo(i,j,k,nc),vhi(i,j,k,nc),
     &                 vel(i,j,k))
                  flux(i,j,k,nc) = cvmgt(half*(vlo(i,j,k,nc)+vhi(i,j,k,nc)),
     &                 flux(i,j,k,nc),abs(vel(i,j,k)).lt.1.d-16)
                  flux(i,j,k,nc) = flux(i,j,k,nc) * vel(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end

      subroutine riemann_solver(flux,vel,vlo,vhi,DIMS(v),
     &     kpedge,DIMS(kp),eigl,DIMS(eig),krcoef,DIMS(krcoef),nkrcoef,
     &     ix,iy,iz,lo,hi,
     &     grav,muval,rhoval,nscal,model,corrector)

c ::: Riemann Solver: wrapper for cmpflx
c ::: arguments

      implicit none

      integer DIMDEC(v)
      integer DIMDEC(kp)
      integer DIMDEC(eig)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ix,iy,iz
      integer nscal,model
      integer lo(SDIM),hi(SDIM)
      integer corrector
      REAL_T  flux(DIMV(v),nscal)
      REAL_T  vlo(DIMV(v),nscal), vhi(DIMV(v),nscal)
      REAL_T  vel(DIMV(v))
      REAL_T  kpedge(DIMV(kp))
      REAL_T  eigl(DIMV(eig),nscal,nscal)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  rhoval(nscal),muval(nscal),grav

c ::: internal variables
      integer i,j,k,nc
      integer imin,jmin,kmin,imax,jmax,kmax
      integer l_kr_type, r_kr_type
      REAL_T  l_kr_coef, l_sat_res, r_kr_coef, r_sat_res
      REAL_T aright, aleft, stmp(2), rtmp(2), st(2), gstar
      logical no_change

      imin = lo(1) - 1 + ix
      imax = hi(1) + 1
      jmin = lo(2) - 1 + iy
      jmax = hi(2) + 1
      kmin = lo(3) - 1 + iz
      kmax = hi(3) + 1

      if (model .eq. 2) then

         rtmp = one

!$omp parallel do private(i,j,k,aright,aleft,nc,st,stmp)
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  r_kr_type = nint(krcoef(i,j,k,1))
                  l_kr_type = nint(krcoef(i-ix,j-iy,k-iz,1))  
                  r_kr_coef = krcoef(i,j,k,2)
                  l_kr_coef = krcoef(i-ix,j-iy,k-iz,2)
                  r_sat_res = krcoef(i,j,k,3)
                  l_sat_res = krcoef(i-ix,j-iy,k-iz,3)
                  no_change = .true.
                  if (r_kr_type .ne. l_kr_type .or. r_kr_coef .ne. l_kr_coef .or.
     &                 r_sat_res .ne. l_sat_res) then
                     no_change = .false.
                  end if

                  aright = zero
                  aleft  = zero

                  do nc = 1,nscal
                     aright = aright + eigl(i,j,k,2,nc)*vhi(i,j,k,nc)
                     aleft  = aleft  + eigl(i-ix,j-iy,k-iz,2,nc)*vlo(i,j,k,nc)
                  end do 

                  stmp(1) = aright
                  stmp(2) = one-stmp(1)
                  call get_st(st,stmp,rtmp,r_sat_res,nscal)
                  aright = st(1)

                  stmp(1) = aleft
                  stmp(2) = one-stmp(1)
                  call get_st(st,stmp,rtmp,l_sat_res,nscal)
                  aleft = st(1)

                  flux(i,j,k,:) = zero
                  gstar = kpedge(i,j,k)*grav*(rhoval(1)-rhoval(2)) 
                  if (no_change) then
                     if (l_kr_type .eq. 2) then
                        call cmpflx(flux(i,j,k,1),aleft,aright,
     &                       vel(i,j,k),muval(1),muval(2),gstar) 
                     else if (l_kr_type .eq. 3) then
                        call cmpflx_vang(flux(i,j,k,1),aleft,aright,
     &                       vel(i,j,k),muval(1),muval(2),gstar,
     &                       l_kr_coef) 
                     end if
                  else
                     call cmpflx_upwind(flux(i,j,k,1),aleft,aright,vel(i,j,k),
     &                    muval,gstar,
     &                    l_kr_type,r_kr_type,l_kr_coef,r_kr_coef)
                  end if

                  flux(i,j,k,2) = vel(i,j,k) - flux(i,j,k,1)

                  if (corrector .eq. 1) then
                     vel(i,j,k) = flux(i,j,k,1)
                  end if

                  flux(i,j,k,1) = rhoval(1) * flux(i,j,k,1)
                  flux(i,j,k,2) = rhoval(2) * flux(i,j,k,2)
               end do
            end do
         end do
!$omp end parallel do

      else
         call bl_pd_abort('model not defined in GODUNOV_3D:riemann_solver');
      end if

      end

      subroutine riemann_solver_arg(flux,vel,vlo,vhi,DIMS(v),
     &     kpedge,DIMS(kp),eigl,DIMS(eig),krcoef,DIMS(krcoef),nkrcoef,
     &     ix,iy,iz,lo,hi,
     &     grav,muval,rhoval,nscal,model)

c ::: Riemann Solver: wrapper for cmpflx arguments

      implicit none

      integer DIMDEC(v)
      integer DIMDEC(kp)
      integer DIMDEC(eig)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ix,iy,iz
      integer nscal,model
      integer lo(SDIM),hi(SDIM)
      REAL_T  flux(DIMV(v),nscal)
      REAL_T  vlo(DIMV(v),nscal), vhi(DIMV(v),nscal)
      REAL_T  vel(DIMV(v))
      REAL_T  kpedge(DIMV(kp))
      REAL_T  eigl(DIMV(eig),nscal,nscal)
      REAL_T krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  rhoval(nscal),muval(nscal),grav

c ::: internal variables
      integer i,j,k,nc    
      integer imin,jmin,kmin,imax,jmax,kmax
      integer l_kr_type, r_kr_type
      REAL_T  l_kr_coef, l_sat_res, r_kr_coef, r_sat_res
      REAL_T aright, aleft, stmp(2), rtmp(2), st(2), gstar, ftmp
      logical no_change

      imin = lo(1) - 1 + ix
      imax = hi(1) + 1
      jmin = lo(2) - 1 + iy
      jmax = hi(2) + 1
      kmin = lo(3) - 1 + iz
      kmax = hi(3) + 1

      if (model .eq. 2) then

         rtmp = one

!$omp parallel do private(i,j,k,aright,aleft,nc,st,stmp)
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  r_kr_type = nint(krcoef(i,j,k,1))
                  l_kr_type = nint(krcoef(i-ix,j-iy,k-iz,1))  
                  r_kr_coef = krcoef(i,j,k,2)
                  l_kr_coef = krcoef(i-ix,j-iy,k-iz,2)
                  r_sat_res = krcoef(i,j,k,3)
                  l_sat_res = krcoef(i-ix,j-iy,k-iz,3)
                  no_change = .true.
                  if (r_kr_type .ne. l_kr_type .or. r_kr_coef .ne. l_kr_coef .or.
     &                 r_sat_res .ne. l_sat_res) then
                     no_change = .false.
                  end if

                  aright = zero
                  aleft  = zero
                  do nc = 1,nscal
                     aright = aright + eigl(i,j,k,2,nc)*vhi(i,j,k,nc)
                     aleft  = aleft  + eigl(i-ix,j-iy,k-iz,2,nc)*vlo(i,j,k,nc)
                  end do 

                  stmp(1) = aright
                  stmp(2) = one-stmp(1)
                  call get_st(st,stmp,rtmp,r_sat_res,nscal)
                  aright = st(1)

                  stmp(1) = aleft
                  stmp(2) = one-stmp(1)
                  call get_st(st,stmp,rtmp,l_sat_res,nscal)
                  aleft = st(1)
                  gstar = kpedge(i,j,k)*grav*(rhoval(1)-rhoval(2)) 
                  if (no_change) then  
                     if (l_kr_type .eq. 2) then
                        call cmpflx_arg(flux(i,j,k,1),aleft,aright,
     &                       vel(i,j,k),muval(1),muval(2),gstar) 
                     else if (l_kr_type .eq. 3) then
                        call cmpflx_vang_arg(flux(i,j,k,1),aleft,aright,
     &                       vel(i,j,k),muval(1),muval(2),gstar,
     &                       l_kr_coef) 
                     end if
                     flux(i,j,k,2) = one  - flux(i,j,k,1)
                  else
                     call cmpflx_upwind(ftmp,flux(i,j,k,:),aleft,aright,
     &                 vel(i,j,k),muval,gstar,
     &                 l_kr_type,r_kr_type,l_kr_coef,r_kr_coef)
                  end if
               end do
            end do
         end do
!$omp end parallel do

      else
         call bl_pd_abort('model not defined in riemann_solver_arg');
      end if

      end

      subroutine riemann_est_eig(s,DIMS(s),lbd,dlbd,d2lbd,DIMS(lbd),
     &     vel,DIMS(vel),kpedge,DIMS(kp),rphi,DIMS(rphi),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     ix,iy,iz,lo,hi,ncomps,grav,rhoval,
     &     muval,dx,bc,eigmax)

c        
c ::: Estimate maximum eigenvalues
c     There are possibly two extrema for the eigenvalues.  In some cases,
c     we must first determine the minimum of the flux function, and find
c     the extrema in both directions.
c
      implicit none

c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(lbd)
      integer DIMDEC(vel)
      integer DIMDEC(kp)
      integer DIMDEC(rphi)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ix,iy,iz
      integer ncomps
      integer lo(SDIM),hi(SDIM)
      integer bc(SDIM,2)
      REAL_T s(DIMV(s),ncomps)
      REAL_T lbd(DIMV(lbd),2)
      REAL_T dlbd(DIMV(lbd),3)
      REAL_T d2lbd(DIMV(lbd),3)
      REAL_T krcoef(DIMV(krcoef),nkrcoef)
      REAL_T vel(DIMV(vel))
      REAL_T kpedge(DIMV(kp))
      REAL_T rphi(DIMV(rphi))
      REAL_T rhoval(ncomps)
      REAL_T muval(ncomps)
      REAL_T grav
      REAL_T eigmax
      REAL_T dx(SDIM)

c ::: internal variables
      integer i,j,k
      integer l_kr_type, r_kr_type
      REAL_T  l_kr_coef, l_sat_res, r_kr_coef, r_sat_res
      REAL_T eigtmp
      REAL_T krhog
      logical no_change

      eigmax = zero

!$omp parallel do private(i,j,k,m,sat,lam,dlam,d2lam,slr,lt,krhog)
!$omp& private(f1,f2,stmp,dtmp,ncycle,nc,df,d2f,dftmp,eigtmp,tmps)
!$omp& reduction(max:eigmax)
      do k = lo(3),hi(3)+iz
         do j = lo(2),hi(2)+iy
            do i = lo(1),hi(1)+ix 
               r_kr_type = nint(krcoef(i,j,k,1))
               l_kr_type = nint(krcoef(i-ix,j-iy,k-iz,1))  
               r_kr_coef = krcoef(i,j,k,2)
               l_kr_coef = krcoef(i-ix,j-iy,k-iy,2)
               r_sat_res = krcoef(i,j,k,3)
               l_sat_res = krcoef(i-ix,j-iy,k-iz,3)
               no_change = .true.
               if (r_kr_type .ne. l_kr_type .or. r_kr_coef .ne. l_kr_coef .or.
     &              r_sat_res .ne. l_sat_res) then
                  no_change = .false.
               end if
               krhog = kpedge(i,j,k)*((rhoval(1)-rhoval(2))*grav) 
               if (no_change) then
                  call find_max_eig(s(i-ix,j-iy,k-iz,:),s(i,j,k,:),ncomps,
     &                 lbd(i-ix,j-iy,k-iz,:),dlbd(i-ix,j-iy,k-iz,:),d2lbd(i-ix,j-iy,k-iz,:),
     &                 lbd(i,j,k,:),dlbd(i,j,k,:),d2lbd(i,j,k,:),
     &                 vel(i,j,k),krcoef(i,j,k,:),nkrcoef,
     &                 krhog,rhoval,muval,eigtmp)
               else     
                  call find_upwind_eig(s(i-ix,j-iy,k-iz,:),s(i,j,k,:),ncomps,
     &                 lbd(i-ix,j-iy,k-iz,:),dlbd(i-ix,j-iy,k-iz,:),
     &                 lbd(i,j,k,:),dlbd(i,j,k,:),
     &                 vel(i,j,k),krcoef(i,j,k,:),nkrcoef,krhog,eigtmp)
               end if
               eigtmp = eigtmp/(half*(rphi(i,j,k)+rphi(i-ix,j-iy,k-iz)))
               eigmax = max(eigmax,dabs(eigtmp))
            end do
         end do
      end do
!$omp end parallel do

      end 

      subroutine get_dlbd_cc(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     ncomps,rhoval,muval,model)
c        
c ::: Determine first and second derivative of flux at cell center 
c
      implicit none

c ::: Arguments
      integer DIMDEC(s)
      integer DIMDEC(lam)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ncomps,model
      REAL_T s(DIMV(s),ncomps)
      REAL_T lam(DIMV(lam),2)
      REAL_T dlam(DIMV(lam),3)
      REAL_T d2lam(DIMV(lam),3)
      REAL_T krcoef(DIMV(krcoef),nkrcoef)
      REAL_T rhoval(ncomps)
      REAL_T muval(ncomps)

c ::: Internal variables
      integer i,j,k
      integer l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T sat(2),rtmp(2), tlam(2), tdlam(3), tmps(2)

      if (model .eq. 2) then

         rtmp = one

!$omp parallel do private(i,j,k,sat,tlam,tdlam,tmps)
         do k = ARG_L3(s),ARG_H3(s)
            do j = ARG_L2(s),ARG_H2(s)
               do i = ARG_L1(s),ARG_H1(s)

                  tmps = s(i,j,k,:)
                  l_kr_type = nint(krcoef(i,j,k,1))
                  l_kr_coef = krcoef(i,j,k,2)
                  l_sat_res = krcoef(i,j,k,3)
                  call get_st(sat,tmps,rhoval,l_sat_res,ncomps)

                  if (l_kr_type .eq. 3) then
                     sat(1) = min(one-1.d-6,sat(1))
                     sat(1) = max(1.d-6,sat(1))
                     sat(2) = one-sat(1)
                  end if

                  call get_ld(tlam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
                  lam(i,j,k,:) = tlam
                  call get_dld(tdlam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
                  dlam(i,j,k,:) = tdlam
                  call get_d2ld(tdlam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
                  d2lam(i,j,k,:) = tdlam
               end do
            end do
         end do
!$omp end parallel do

      end if

      end 

      subroutine riemann_solver2(flux,DIMS(flux),s,DIMS(s),vlo,vhi,DIMS(v),
     &     eigv,DIMS(eig),vel,DIMS(vel),kpedge,DIMS(kp),ix,iy,iz,lo,hi,
     &     athick,grav,muval,rhoval,ncomps,nphases,dx,dt,bc)
c ::: Riemann Solver: wrapper for cmpflx_polymer
c        ncomps is the number of auxiliary variables.
c ::: arguments
      integer DIMDEC(flux)
      integer DIMDEC(s)
      integer DIMDEC(v)
      integer DIMDEC(eig)
      integer DIMDEC(vel)
      integer DIMDEC(kp)
      integer ix,iy,iz,ncomps,nphases
      integer lo(SDIM),hi(SDIM)
      integer bc(SDIM,2)
      REAL_T flux(DIMV(flux),ncomps)
      REAL_T s(DIMV(s),ncomps)
      REAL_T vlo(DIMV(v),ncomps), vhi(DIMV(v),ncomps)
      REAL_T eigv(DIMV(eig),ncomps)
      REAL_T vel(DIMV(vel))
      REAL_T kpedge(DIMV(kp))
      REAL_T tum(2),tup(2),tulft(2),turgt(2),teigr(2),teigl(2),tflx(2)
      REAL_T rhoval(nphases),muval(nphases),grav,athick
      REAL_T dx, dt

c ::: internal variables
      integer i,j,k,nc    
      integer imin,jmin,kmin,imax,jmax,kmax

      if (ix .eq. 1) then
         imin = lo(1)
         imax = hi(1) + 1
         jmin = lo(2) - 1
         jmax = hi(2) + 1
         kmin = lo(3) - 1
         kmax = hi(3) + 1
      elseif (iy .eq. 1) then
         imin = lo(1) - 1
         imax = hi(1) + 1
         jmin = lo(2) 
         jmax = hi(2) + 1
         kmin = lo(3) - 1
         kmax = hi(3) + 1
      elseif (iz .eq. 1) then
         imin = lo(1) - 1
         imax = hi(1) + 1
         jmin = lo(2) - 1
         jmax = hi(2) + 1
         kmin = lo(3) 
         kmax = hi(3) + 1
      end if

      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               do nc = 1,ncomps
                  tum(nc)   = s(i-ix,j-iy,k-iz,nc)
                  tup(nc)   = s(i,j,k,nc)
                  tulft(nc) = vlo(i,j,k,nc)
                  turgt(nc) = vhi(i,j,k,nc)
                  teigl(nc) = eigv(i-ix,j-iy,k-iz,nc)
                  teigr(nc) = eigv(i,j,k,nc)
               end do
               call cmpflx_polymer(tum,tup,tulft,turgt,teigl,teigr,tflx,vel(i,j,k),
     &              grav, muval(1), muval(2), athick, rhoval(1),rhoval(2),
     &              kpedge(i,j,k),dx, dt)
               flux(i,j,k,1) = tflx(1)
               flux(i,j,k,2) = tflx(2)
            end do
         end do
      end do

      if (iy .eq. 1) then
         if (bc(2,2).eq.FOEXTRAP) then
            do k = kmin,kmax
               do i = imin,imax
                  flux(i,jmax+1,k,1) = 0
                  flux(i,jmax+1,k,2) = 0
               end do
            end do
         end if
         if (bc(2,1).eq.FOEXTRAP) then
            do k = kmin,kmax
               do i = imin,imax
                  flux(i,jmin,k,1) = 0
                  flux(i,jmin,k,2) = 0
               end do
            end do
         end if
      end if

      if (iz .eq. 1) then
         if (bc(3,2).eq.FOEXTRAP) then
            do j = jmin,jmax
               do i = imin,imax
                  flux(i,j,kmax+1,1) = 0
                  flux(i,j,kmax+1,2) = 0
               end do
            end do
         end if
         if (bc(3,1).eq.FOEXTRAP) then
            do j = jmin,jmax
               do i = imin,imax
                  flux(i,j,kmin,1) = 0
                  flux(i,j,kmin,2) = 0
               end do
            end do
         end if
      end if

      end

      subroutine advect_tracers(flux,DIMS(flux),uedge,vlo,vhi,DIMS(uedge),
     &     ix,iy,iz,lo,hi,nscal)

c ::: Calculate tracers' flux

c ::: arguments
      integer DIMDEC(uedge)
      integer DIMDEC(flux)
      integer ix,iy,iz
      integer nscal
      integer lo(SDIM),hi(SDIM)
      REAL_T  flux(DIMV(flux),nscal)
      REAL_T  vlo(DIMV(uedge),nscal)
      REAL_T  vhi(DIMV(uedge),nscal)
      REAL_T  uedge(DIMV(uedge))

c ::: internal variables
      integer i,j,k,nc 

!$omp parallel private(i,j,k,nc)
      do nc = 1,nscal
!$omp do
         do k = lo(3),hi(3)+iz
            do j = lo(2),hi(2)+iy
               do i = lo(1),hi(1)+ix
                  flux(i,j,k,nc) = cvmgp(vlo(i,j,k,nc),vhi(i,j,k,nc),uedge(i,j,k))
                  flux(i,j,k,nc) = cvmgt(half*(vlo(i,j,k,nc)+vhi(i,j,k,nc)),
     &                 flux(i,j,k,nc),abs(uedge(i,j,k)).lt. 1.d-16)
                  flux(i,j,k,nc) = flux(i,j,k,nc) * uedge(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end

      subroutine flux_copy(fluxs,DIMS(fluxs),fluxu,DIMS(fluxu),
     &     lo,hi,ix,iy,iz,ncomps)
c ::: arguments
      integer DIMDEC(fluxs)
      integer DIMDEC(fluxu)
      integer lo(SDIM),hi(SDIM),ix,iy,iz,ncomps
      REAL_T  fluxs(DIMV(fluxs),ncomps)
      REAL_T  fluxu(DIMV(fluxu),ncomps)

c ::: internal variables
      integer i,j,k,nc    
      integer imin,jmin,imax,jmax,kmin,kmax
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)

      do nc = 1,ncomps
         do k = kmin,kmax+iz
            do j = jmin,jmax+iy
               do i = imin,imax+ix
                  fluxs(i,j,k,nc) = fluxu(i,j,k,nc)
               end do    
            end do
         end do
      end do
      end

      subroutine flux_convert(fluxs,DIMS(fluxs),fluxu,DIMS(fluxu),vel,DIMS(v),
     &     lo,hi,ix,iy,iz,rho,compcount,compidx,ncomps,nphases,bc)
c ::: arguments
      integer DIMDEC(fluxs)
      integer DIMDEC(fluxu)
      integer DIMDEC(v)
      integer lo(SDIM),hi(SDIM),ix,iy,iz,ncomps,nphases
      integer bc(SDIM,2)
      integer compcount(nphases)
      integer compidx(nphases,ncomps)
      REAL_T  fluxs(DIMV(fluxs),ncomps)
      REAL_T  fluxu(DIMV(fluxu),ncomps-1)
      REAL_T  vel(DIMV(v))
      REAL_T  rho(nphases)


c ::: internal variables
      integer i,j,k,nc    
      integer imin,jmin,imax,jmax,kmin,kmax
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)

      do nc = 1,nphases
         do k = kmin,kmax+iz
            do j = jmin,jmax+iy
               do i = imin,imax+ix
                  if (compcount(nc) .eq. 2) then
                     fluxs(i,j,k,compidx(nc,1)) = rho(1)*
     &                 (fluxu(i,j,k,1)-fluxu(i,j,k,2))
                     fluxs(i,j,k,compidx(nc,2)) = rho(1)*
     &                    fluxu(i,j,k,2)
                  else if (compcount(nc) .eq. 1) then
                     fluxs(i,j,k,compidx(nc,1)) = rho(2)*
     &                 (vel(i,j,k)-fluxu(i,j,k,1))
                  end if
               end do
            end do
         end do
      end do  

      if (iy .eq. 1) then
         if (bc(2,2).eq.FOEXTRAP) then

            do nc = 1,nphases
               do k = kmin,kmax
                  do i = imin,imax
                     if (compcount(nc) == 1) then
                        fluxs(i,jmax+1,k,compidx(nc,1)) = 0
                     end if
                  end do
               end do
            end do
         end if
         if (bc(2,1).eq.FOEXTRAP) then

            do nc = 1,nphases
               do k = kmin,kmax
                  do i = imin,imax
                     if (compcount(nc) == 1) then
                        fluxs(i,jmin,k,compidx(nc,1)) = 0
                     end if
                  end do
               end do
            end do
         end if
      end if

      if (iz .eq. 1) then
         if (bc(3,2).eq.FOEXTRAP) then

            do nc = 1,nphases
               do j = jmin,jmax
                  do i = imin,imax
                     if (compcount(nc) == 1) then
                        fluxs(i,j,kmax+1,compidx(nc,1)) = 0
                     end if
                  end do
               end do
            end do
         end if
         if (bc(3,1).eq.FOEXTRAP) then

            do nc = 1,nphases
               do j = jmin,jmax
                  do i = imin,imax
                     if (compcount(nc) == 1) then
                        fluxs(i,j,kmin,compidx(nc,1)) = 0
                     end if
                  end do
               end do
            end do
         end if
      end if
      
      end

      subroutine eigsol_mp(eigv,eigl,eigr,lbd,dlbd,DIMS(w),
     &     v,DIMS(v),rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,grav,mu,
     &     ix,iy,iz,lo,hi,rhoval,ncomps,nphases,model,bc)

c    
c     eigensystem for multiphase models.
c    

      implicit none

      integer DIMDEC(w)
      integer DIMDEC(v)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ncomps,nphases,model
      integer lo(SDIM),hi(SDIM)
      integer ix,iy,iz
      integer bc(SDIM,2)
      REAL_T  eigv(DIMV(w),2*ncomps)
      REAL_T  eigl(DIMV(w),ncomps,ncomps)
      REAL_T  eigr(DIMV(w),ncomps,ncomps)
      REAL_T  lbd(DIMV(w),ncomps)
      REAL_T  dlbd(DIMV(w),3)
      REAL_T  v(DIMV(v))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  mu(nphases),grav
      REAL_T  rhoval(ncomps)

c ::: internal variables
      integer i,j,k,n
      integer imin,jmin,kmin,imax,jmax,kmax,ibc
      integer l_kr_type
      REAL_T  l_sat_res
      REAL_T  lam(nphases), dlam(3), ltot
      REAL_T  gs1,gs2,vt1,vt2

      ibc  = ix + 2*iy + 3*iz
      imin = lo(1)-1
      imax = hi(1)+1
      jmin = lo(2)-1
      jmax = hi(2)+1
      kmin = lo(3)-1
      kmax = hi(3)+1
      
      if (model == 2) then 

!$omp parallel do private(i,j,k,lam,dlam,ltot,gs1,gs2,vt1,vt2)
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax

                  l_kr_type = nint(krcoef(i,j,k,1))
                  l_sat_res = krcoef(i,j,k,3)
c     eigenvalues for the left (1 to ncomps) and 
c     right (ncomps to 2*ncomps) state.

                  lam  = lbd(i,j,k,:)
                  dlam = dlbd(i,j,k,:)

                  ltot = lam(1)+lam(2)

                  gs1  = kappa(i,j,k)*(rhoval(1)-rhoval(2))*grav
                  gs2  = kappa(i+ix,j+iy,k+iz)*(rhoval(1)-rhoval(2))*grav

                  vt1 = v(i,j,k)
                  vt2 = v(i+ix,j+iy,k+iz)

                  if (ix .eq. 1) then
                     if (bc(ibc,1) .eq. EXT_DIR .and. i .eq. imin) then
                        vt1 = v(i+1,j,k)
                     else if  (bc(ibc,2) .eq. EXT_DIR .and. i .eq. imax) then
                        vt2 = v(i,j,k)
                     end if
                  elseif (iy .eq. 1) then
                     if (bc(ibc,1) .eq. EXT_DIR .and. j .eq. jmin) then
                        vt1 = v(i,j+1,k)
                     else if  (bc(ibc,2) .eq. EXT_DIR .and. j .eq. jmax) then
                        vt2 = v(i,j,k)
                     end if
                  else
                     if (bc(ibc,1) .eq. EXT_DIR .and. k .eq. kmin) then
                        vt1 = v(i,j,k+1)
                     else if  (bc(ibc,2) .eq. EXT_DIR .and. k .eq. kmax) then
                        vt2 = v(i,j,k)
                     end if
                  end if

                  eigv(i,j,k,1) = zero
                  eigv(i,j,k,3) = zero

                  if (l_kr_type .eq. 2) then
                     eigv(i,j,k,2) = vt1*(dlam(1)/ltot+lam(1)*dlam(3))+
     &                    gs1*(dlam(2)/ltot+lam(1)*lam(2)*dlam(3))
                     eigv(i,j,k,2) = eigv(i,j,k,2)/rphi(i,j,k)/(one-l_sat_res)
                     eigv(i,j,k,4) = vt2*(dlam(1)/ltot+lam(1)*dlam(3))+
     &                    gs2*(dlam(2)/ltot+lam(1)*lam(2)*dlam(3))
                     eigv(i,j,k,4) = eigv(i,j,k,4)/rphi(i,j,k)/(one-l_sat_res)
                  elseif (l_kr_type .eq. 3) then
                     eigv(i,j,k,2) = (dlam(1)*lam(2)*(vt1+gs1*lam(2))
     &                    + lam(1)*(-vt1 + gs1*lam(1))*dlam(2))/(ltot**2)
                     eigv(i,j,k,2) = eigv(i,j,k,2)/rphi(i,j,k)/(one-l_sat_res)
                     eigv(i,j,k,4) = (dlam(1)*lam(2)*(vt2+gs2*lam(2))
     &                    + lam(1)*(-vt2 + gs2*lam(1))*dlam(2))/(ltot**2)
                     eigv(i,j,k,4) = eigv(i,j,k,4)/rphi(i,j,k)/(one-l_sat_res)
                  else
                     call bl_pd_abort('eigensolver for the kr_type specified
     &does not exist')
                  end if

c     ::: left eigenvectors
                  eigl(i,j,k,1,1) = one/rhoval(1)
                  eigl(i,j,k,1,2) = one/rhoval(2)
                  eigl(i,j,k,2,1) = one/rhoval(1)
                  eigl(i,j,k,2,2) = zero
                  
c     ::: right eigenvectors
                  eigr(i,j,k,1,1) =  zero
                  eigr(i,j,k,1,2) =  rhoval(1)
                  eigr(i,j,k,2,1) =  rhoval(2)
                  eigr(i,j,k,2,2) = -rhoval(2)
               end do
            end do
         end do
!$omp end parallel do

      elseif (model == 0 .or. model == 1) then
        eigl = 0
        eigr = 0
        eigv = 0
        do k = kmin,kmax
           do j = jmin,jmax
              do i = imin,imax

                  vt1 = v(i,j,k)
                  vt2 = v(i+ix,j+iy,k+iz)

                  if (ix .eq. 1) then
                     if (bc(ibc,1) .eq. EXT_DIR .and. i .eq. imin) then
                        vt1 = v(i+1,j,k)
                     else if  (bc(ibc,2) .eq. EXT_DIR .and. i .eq. imax) then
                        vt2 = v(i,j,k)
                     end if
                  elseif (iy .eq. 1) then
                     if (bc(ibc,1) .eq. EXT_DIR .and. j .eq. jmin) then
                        vt1 = v(i,j+1,k)
                     else if  (bc(ibc,2) .eq. EXT_DIR .and. j .eq. jmax) then
                        vt2 = v(i,j,k)
                     end if
                  else
                     if (bc(ibc,1) .eq. EXT_DIR .and. k .eq. kmin) then
                        vt1 = v(i,j,k+1)
                     else if  (bc(ibc,2) .eq. EXT_DIR .and. k .eq. kmax) then
                        vt2 = v(i,j,k)
                     end if
                  end if

                  do n = 1,ncomps
c     ::: eigenvalues
                     eigv(i,j,k,n) = vt1/rphi(i,j,k)
                     eigv(i,j,k,n+ncomps) = vt2/rphi(i,j,k)
c     ::: left eigenvectors
                     eigl(i,j,k,n,n) = one

c     ::: right eigenvectors
                     eigr(i,j,k,n,n) = one
                  end do
               end do
            end do
         end do
      else
         call bl_pd_abort('model is not defined in GODUNOV_3D:eigsol')
      end if
      
      end

      subroutine eigsol2(s,DIMS(s),v,DIMS(v),rphi,DIMS(rphi),
     &     kappa,DIMS(kappa),eigv,eigl,eigr,struc,ldef,DIMS(w),
     &     lo,hi,ix,iy,iz,mup,athick,rhop,grav,ncomps,nphases,eigmax)
c ::  eigensystem for the problem
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(v)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(w)
      integer ncomps,nphases,lo(SDIM),hi(SDIM),ix,iy,iz
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  v(DIMV(v))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  eigv(DIMV(w),ncomps)
      REAL_T  eigl(DIMV(w),ncomps,ncomps)
      REAL_T  eigr(DIMV(w),ncomps,ncomps)
      REAL_T  struc(DIMV(w),ncomps,ncomps)
      logical ldef(DIMV(w))
      REAL_T  mup(nphases)
      REAL_T  rhop(nphases)
      REAL_T  athick
      REAL_T  grav
      REAL_T  eigmax

c ::: internal variables
      integer i,j,k,n,m
      integer imin,jmin,kmin,imax,jmax,kmax
      REAL_T  vt, gstar
      REAL_T  u(2)
      REAL_T  teigv(2), teigl(2,2), teigr(2,2), dl(2,2)
      logical ldeflocal
      imin = lo(1)-1
      jmin = lo(2)-1
      kmin = lo(3)-1
      imax = hi(1)+1
      jmax = hi(2)+1
      kmax = hi(3)+1
      eigmax = 0

      do k = kmin,kmax      
         do j = jmin,jmax
            do i = imin,imax
               gstar = kappa(i,j,k)*(rhop(1)-rhop(2))*grav
               u(1) = s(i,j,k,1)
               u(2) = s(i,j,k,2)
               vt = half*(v(i+ix,j+iy,k+iz) + v(i,j,k))

               call eigen_polymer(vt,gstar,mup(1),mup(2),
     &              athick,u,teigl,teigr,teigv,dl,ldeflocal)
               ldef(i,j,k) = ldeflocal
               do n = 1,2
                  eigv(i,j,k,n) = teigv(n)/rphi(i,j,k)
                  do m = 1,2
                     eigl(i,j,k,n,m)  = teigl(n,m)
                     eigr(i,j,k,n,m)  = teigr(n,m)
                     struc(i,j,k,n,m) = dl(n,1)*teigr(1,m) + dl(n,2)*teigr(2,m)
                  end do
               end do
            end do
         end do
      end do
      
c     ::: find maximum eigenvalues

      do n = 1,2
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  eigmax = max(eigmax, abs(eigv(i,j,k,n)))
               end do 
            end do
         end do
      end do

      end

      subroutine dfluxdv(drmndv,drmndk,DIMS(drmn),s,DIMS(s),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),krcoef,DIMS(krcoef),nkrcoef,
     &     mu,rhoval,ncomps,model,grav)

c ::: Derivative Flux w.r.t. v_T and kappa.  These are analytical functions.

      implicit none

c ::: arguments      
      integer DIMDEC(drmn)
      integer DIMDEC(s)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ncomps,model
      REAL_T  drmndv(DIMV(drmn),ncomps)
      REAL_T  drmndk(DIMV(drmn),ncomps) 
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))    
      REAL_T krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  mu(ncomps),rhoval(ncomps)
      REAL_T  grav

c ::: internal variables
      integer i,j,k,nc
      integer l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T  kr(ncomps)
      REAL_T  ltot,gstar,xmu1,xmu2

      if (model .eq. 2) then

         gstar = (rhoval(1)-rhoval(2))*grav
         xmu1  = one / mu(1)
         xmu2  = one / mu(2)

!$omp parallel do private(i,j,k,kr,ltot)
         do k = ARG_L3(drmn),ARG_H3(drmn)  
            do j = ARG_L2(drmn),ARG_H2(drmn)
               do i = ARG_L1(drmn),ARG_H1(drmn)
                  l_kr_type = nint(krcoef(i,j,k,1))
                  l_kr_coef = krcoef(i,j,k,2)
                  l_sat_res = krcoef(i,j,k,3)
                  call get_kr(kr,s(i,j,k,1),rhoval,ncomps,l_kr_type,l_kr_coef,l_sat_res)
                  ltot = (kr(1) + kr(2))*xmu1
                  drmndv(i,j,k,1) = kr(1)/ltot/rphi(i,j,k)
                  ltot = kr(1)*mu(2)*xmu1 + kr(2)
                  drmndv(i,j,k,2) = kr(2)/ltot/rphi(i,j,k)
                  drmndk(i,j,k,1) = kr(2)*drmndv(i,j,k,1)*gstar*xmu2
                  drmndk(i,j,k,2) = -kr(1)*drmndv(i,j,k,2)*gstar*xmu1
                  drmndv(i,j,k,1) = rhoval(1)*drmndv(i,j,k,1)
                  drmndv(i,j,k,2) = rhoval(2)*drmndv(i,j,k,2)
                  drmndk(i,j,k,1) = rhoval(1)*drmndk(i,j,k,1)
                  drmndk(i,j,k,2) = rhoval(2)*drmndk(i,j,k,2)
               end do
            end do
         end do
!$omp end parallel do

      end if         

      end

      subroutine dfluxdv2(drmndv,drmndk,s,rphi,kappa,
     &     mu1,mu2, athick,rho1,rho2,grav)
c ::: Derivative Flux w.r.t. v_T.  This is an analytical function
c ::: arguments      

      REAL_T  drmndv(2)
      REAL_T  drmndk(2)
      REAL_T  s(2)
      REAL_T  rphi
      REAL_T  kappa
      REAL_T  mu1, mu2, rho1, rho2, athick
      REAL_T  grav
      
c ::: internal variables
      REAL_T  coef,lambda(2),ltot,gstar

      gstar = (rho1-rho2)*grav
      coef = zero
      if (dabs(s(1)) .gt. 1.d-10) then 
         coef  = s(2)/s(1)
      end if
      lambda(1) = s(1)**2/(mu1+athick*coef)
      lambda(2) = (one-s(1))**2/mu2
      ltot = lambda(1) + lambda(2)

      drmndv(1) = lambda(1)*lambda(2)/ltot/rphi
      drmndv(2) = coef*drmndv(1)
      drmndk(1) = lambda(1)*lambda(2)*gstar/ltot/rphi
      drmndk(2) = coef*drmndk(1)
        
      end


      subroutine get_aux(s, u, DIMS(s),compcount,compidx, 
     &     rhoval, ncomps, nphases)
c
c     ----------------------------------------------------------
c     Loop through a fab to convert from primary to auxiliary
c     variables for the polymer problem.  
c     Use get_aux_i in GODUNOV_F.F
c     ----------------------------------------------------------
c
      integer DIMDEC(s)
      integer ncomps, nphases
      integer compcount(nphases)
      integer compidx(nphases,ncomps)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  u(DIMV(s),ncomps-1)
      REAL_T  rhoval(ncomps)

      integer i,j,k
      REAL_T  stmp(ncomps), utmp(ncomps-1)

      do k =  ARG_L3(s), ARG_H3(s)      
         do j =  ARG_L2(s), ARG_H2(s)
            do i = ARG_L1(s), ARG_H1(s)
               utmp = zero
               stmp = s(i,j,k,:)
               call get_aux_i(stmp,utmp,compcount,compidx,
     &              rhoval,ncomps,nphases)    
               u(i,j,k,:) = utmp
            end do
         end do
      end do

      end

      subroutine get_aux_n(s, u, DIMS(s), compcount, compidx, 
     &     rhoval, bc, ncomps, nphases)
c
c     ----------------------------------------------------------
c     Loop through a fab to convert from primary to auxiliary
c     variables for the polymer problem.  
c     Use get_aux_i in GODUNOV_F.F
c     ----------------------------------------------------------
c 
      integer DIMDEC(s)
      integer ncomps, nphases
      integer compcount(nphases)
      integer compidx(nphases,ncomps)
      integer bc(SDIM,2)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  u(DIMV(s),ncomps-1)
      REAL_T  rhoval(ncomps)

      integer i,j,k,ibc
      integer imin,jmin,imax,jmax,kmin,kmax 
      integer idx(SDIM,SDIM)
      REAL_T  stmp(ncomps), utmp(ncomps-1)
      
      idx = 0
      do i=1,3
         idx(i,i) = 1
      end do
      do ibc = 1,3
         if (bc(ibc,1) .eq. EXT_DIR) then
            imin = ARG_L1(s)*(1-idx(ibc,1)) + ARG_L1(s)*idx(ibc,1)
            imax = ARG_H1(s)*(1-idx(ibc,1)) + ARG_L1(s)*idx(ibc,1)
            jmin = ARG_L2(s)*(1-idx(ibc,2)) + ARG_L2(s)*idx(ibc,2)
            jmax = ARG_H2(s)*(1-idx(ibc,2)) + ARG_L2(s)*idx(ibc,2)
            kmin = ARG_L3(s)*(1-idx(ibc,3)) + ARG_L3(s)*idx(ibc,3)
            kmax = ARG_H3(s)*(1-idx(ibc,3)) + ARG_L3(s)*idx(ibc,3)

            do k = kmin,kmax
               do j =  jmin,jmax
                  do i = imin,imax
                     utmp = zero
                     stmp = s(i,j,k,:)
                     call get_aux_i(stmp,utmp,compcount,compidx,
     &                    rhoval,ncomps,nphases)    
                     u(i,j,k,:) = utmp
                  end do
               end do
            end do
         end if
      
         if (bc(ibc,2) .eq. EXT_DIR) then
            imin = ARG_L1(s)*(1-idx(ibc,1)) + ARG_H1(s)*idx(ibc,1)
            imax = ARG_H1(s)*(1-idx(ibc,1)) + ARG_H1(s)*idx(ibc,1)
            jmin = ARG_L2(s)*(1-idx(ibc,2)) + ARG_H2(s)*idx(ibc,2)
            jmax = ARG_H2(s)*(1-idx(ibc,2)) + ARG_H2(s)*idx(ibc,2)
            kmin = ARG_L3(s)*(1-idx(ibc,3)) + ARG_H3(s)*idx(ibc,3)
            kmax = ARG_H3(s)*(1-idx(ibc,3)) + ARG_H3(s)*idx(ibc,3)

            do k = kmin,kmax
               do j =  jmin,jmax
                  do i = imin,imax
                     utmp = zero
                     stmp = s(i,j,k,:)
                     call get_aux_i(stmp,utmp,compcount,compidx,
     &                    rhoval,ncomps,nphases)    
                     u(i,j,k,:) = utmp
                  end do
               end do
            end do
         end if
      end do

      end

      subroutine FORT_GET_EIG(s,DIMS(s),eigv,DIMS(w),ncomps,idx) 
c ::: trace the cell centered states to edges.
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(w)
      integer ncomps,idx
      REAL_T  s(DIMV(s))
      REAL_T  eigv(DIMV(w),2*ncomps)
      integer i,j,k

      do k = ARG_L3(w),ARG_H3(w)
         do j = ARG_L2(w),ARG_H2(w)
            do i = ARG_L1(w),ARG_H1(w)
                  s(i,j,k) = eigv(i,j,k,idx)
               end do
            end do
         end do
      
      end

      subroutine FORT_ADV_FORCING(
     &     aofs,DIMS(aofs),
     &     xflux,DIMS(xflux),
     &     uedge,DIMS(uedge),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vedge,DIMS(vedge),
     &     areay,DIMS(ay),
     &     zflux,DIMS(zflux),
     &     wedge,DIMS(wedge),
     &     areaz,DIMS(az),
     &     vol,DIMS(vol),
     &     lo,hi,iconserv )
c
c     This subroutine uses scalar edge states to compute
c     an advective tendency
c
      implicit none
      integer i,j,k
      integer iconserv
      REAL_T divux,divuy,divuz
      integer imin,jmin,kmin,imax,jmax,kmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(aofs)
      integer DIMDEC(vol)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(wedge)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      REAL_T aofs(DIMV(aofs))
      REAL_T vol(DIMV(vol))
      REAL_T uedge(DIMV(uedge))
      REAL_T vedge(DIMV(vedge))
      REAL_T wedge(DIMV(wedge))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T zflux(DIMV(zflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
      REAL_T areaz(DIMV(az))

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     if nonconservative initialize the advective tendency as -U*grad(S)
c
      if ( iconserv .ne. 1 ) then
!$omp parallel do private(i,j,k,divux,divuy,divuz)
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  divux = (
     &                 areax(i+1,j,k)*uedge(i+1,j,k)-
     &                 areax(i,  j,k)*uedge(i,  j,k))/vol(i,j,k)
                  divuy = (
     &                 areay(i,j+1,k)*vedge(i,j+1,k)-
     &                 areay(i,j,  k)*vedge(i,j,  k))/vol(i,j,k)
                  divuz = (
     &                 areaz(i,j,k+1)*wedge(i,j,k+1)-
     &                 areaz(i,j,k  )*wedge(i,j,k  ))/vol(i,j,k)
                  aofs(i,j,k) =
     &                 - divux*half*(xflux(i+1,j,k)+xflux(i,j,k))
     &                 - divuy*half*(yflux(i,j+1,k)+yflux(i,j,k))
     &                 - divuz*half*(zflux(i,j,k+1)+zflux(i,j,k))
               end do
            end do
         end do
!$omp end parallel do
      end if
c
c     convert edge states to fluxes
c
!$omp parallel private(i,j,k)
!$omp do
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax+1
               xflux(i,j,k) = xflux(i,j,k)*uedge(i,j,k)*areax(i,j,k)
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = kmin,kmax
         do j = jmin,jmax+1
            do i = imin,imax
               yflux(i,j,k) = yflux(i,j,k)*vedge(i,j,k)*areay(i,j,k)
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = kmin,kmax+1
         do j = jmin,jmax
            do i = imin,imax
               zflux(i,j,k) = zflux(i,j,k)*wedge(i,j,k)*areaz(i,j,k)
            end do
         end do
      end do
!$omp end do nowait
!$omp end parallel
c
c     compute the part of the advective tendency 
c     that depends on the flux convergence
c
      if ( iconserv .ne. 1 ) then
!$omp parallel do private(i,j,k)
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  aofs(i,j,k) = aofs(i,j,k) + (
     &                 xflux(i+1,j,k) - xflux(i,j,k) +
     &                 yflux(i,j+1,k) - yflux(i,j,k) +
     &                 zflux(i,j,k+1) - zflux(i,j,k))/vol(i,j,k)
               end do
            end do
         end do
!$omp end parallel do
      else
!$omp parallel do private(i,j,k)
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  aofs(i,j,k) = (
     &                 xflux(i+1,j,k) - xflux(i,j,k) +
     &                 yflux(i,j+1,k) - yflux(i,j,k) +
     &                 zflux(i,j,k+1) - zflux(i,j,k))/vol(i,j,k)
               end do
            end do
         end do
!$omp end parallel do
      end if

      end

      subroutine FORT_ADV_RMN_FORCING(
     &     aofs,DIMS(aofs),
     &     xflux,DIMS(xflux),uedge,DIMS(uedge),areax,DIMS(ax),
     &     yflux,DIMS(yflux),vedge,DIMS(vedge),areay,DIMS(ay),
     &     zflux,DIMS(zflux),wedge,DIMS(wedge),areaz,DIMS(az),
     &     vol,DIMS(vol),lo,hi,nscal)
c
c     This subroutine uses scalar edge states to compute
c     an advective tendency
c
      implicit none

#include "probdata.H"

      integer i,j,k,nc
      integer nscal
      integer imin,jmin,kmin,imax,jmax,kmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(aofs)
      integer DIMDEC(vol)
      integer DIMDEC(uedge), DIMDEC(vedge), DIMDEC(wedge)
      integer DIMDEC(xflux), DIMDEC(yflux), DIMDEC(zflux)
      integer DIMDEC(ax),DIMDEC(ay),DIMDEC(az)
      REAL_T  aofs(DIMV(aofs),nscal)
      REAL_T  vol(DIMV(vol))
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge)), wedge(DIMV(wedge))
      REAL_T  xflux(DIMV(xflux),nscal)
      REAL_T  yflux(DIMV(yflux),nscal)
      REAL_T  zflux(DIMV(zflux),nscal)
      REAL_T  areax(DIMV(ax)), areay(DIMV(ay)), areaz(DIMV(az))

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     if nonconservative initialize the advective tendency as -U*grad(S)
c
     
c
c     convert edge states to fluxes
c
!$omp parallel private(i,j,k,nc)
      do nc = 1,nscal
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax+1
                  xflux(i,j,k,nc) = xflux(i,j,k,nc)*areax(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax+1
               do i = imin,imax
                  yflux(i,j,k,nc) = yflux(i,j,k,nc)*areay(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
!$omp do
         do k = kmin,kmax+1
            do j = jmin,jmax
               do i = imin,imax
                  zflux(i,j,k,nc) = zflux(i,j,k,nc)*areaz(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel
c
c     compute the part of the advective tendency 
c     that depends on the flux convergence
c     
!$omp parallel private(i,j,k,nc)
      do nc = 1,nscal
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  aofs(i,j,k,nc) = (
     &                 xflux(i+1,j,k,nc) - xflux(i,j,k,nc) +
     &                 yflux(i,j+1,k,nc) - yflux(i,j,k,nc) +
     &                 zflux(i,j,k+1,nc) - zflux(i,j,k,nc))/vol(i,j,k)   
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end

      subroutine FORT_SYNC_ADV_FORCING(
     &     sync ,DIMS(sync),
     &     xflux,DIMS(xflux),
     &     ucor ,DIMS(ucor),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vcor ,DIMS(vcor),
     &     areay,DIMS(ay),
     &     zflux,DIMS(zflux),
     &     wcor ,DIMS(wcor),
     &     areaz,DIMS(az),
     &     vol ,DIMS(vol),
     &     lo,hi,iconserv )
c
c     This subroutine computes the sync advective tendency
c     for a state variable
c
      implicit none
      integer i,j,k
      integer iconserv
      integer imin,jmin,kmin,imax,jmax,kmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(sync)
      integer DIMDEC(vol)
      integer DIMDEC(ucor)
      integer DIMDEC(vcor)
      integer DIMDEC(wcor)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      REAL_T sync(DIMV(sync))
      REAL_T vol(DIMV(vol))
      REAL_T ucor(DIMV(ucor))
      REAL_T vcor(DIMV(vcor))
      REAL_T wcor(DIMV(wcor))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T zflux(DIMV(zflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
      REAL_T areaz(DIMV(az))

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     compute corrective fluxes from edge states 
c     and perform conservative update
c
!$omp parallel private(i,j,k)
!$omp do
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax+1
               xflux(i,j,k) = xflux(i,j,k)*ucor(i,j,k)*areax(i,j,k)
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = kmin,kmax
         do j = jmin,jmax+1
            do i = imin,imax
               yflux(i,j,k) = yflux(i,j,k)*vcor(i,j,k)*areay(i,j,k)
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = kmin,kmax+1
         do j = jmin,jmax
            do i = imin,imax
               zflux(i,j,k) = zflux(i,j,k)*wcor(i,j,k)*areaz(i,j,k)
            end do
         end do
      end do
!$omp end do
!$omp end parallel

!$omp parallel do private(i,j,k)
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               sync(i,j,k) = sync(i,j,k) + (
     &              xflux(i+1,j,k)-xflux(i,j,k) +
     &              yflux(i,j+1,k)-yflux(i,j,k) +
     &              zflux(i,j,k+1)-zflux(i,j,k) )/vol(i,j,k)
            end do
         end do
      end do
!$omp end parallel do

      end

      subroutine FORT_SYNC_RMN_FORCING(
     &     sync ,DIMS(sync),
     &     xflux,DIMS(xflux),
     &     ucor ,DIMS(ucor),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vcor ,DIMS(vcor),
     &     areay,DIMS(ay),
     &     zflux,DIMS(zflux),
     &     wcor ,DIMS(wcor),
     &     areaz,DIMS(az),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     vol ,DIMS(vol),
     &     lo,hi,nscal)
c
c     This subroutine computes the sync advective tendency
c     for a state variable
c
      implicit none

#include "probdata.H"

      integer nscal
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(sync)
      integer DIMDEC(vol)
      integer DIMDEC(ucor)
      integer DIMDEC(vcor)
      integer DIMDEC(wcor)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      integer DIMDEC(krcoef)
      integer nkrcoef
      REAL_T sync(DIMV(sync),nscal)
      REAL_T vol(DIMV(vol))
      REAL_T ucor(DIMV(ucor))
      REAL_T vcor(DIMV(vcor))
      REAL_T wcor(DIMV(wcor))
      REAL_T xflux(DIMV(xflux),nscal)
      REAL_T yflux(DIMV(yflux),nscal)
      REAL_T zflux(DIMV(zflux),nscal)
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
      REAL_T areaz(DIMV(az))
      REAL_T krcoef(DIMV(krcoef),nkrcoef)

      integer i,j,k,nc
      integer imin,jmin,kmin,imax,jmax,kmax
      integer l_kr_type
      REAL_T  l_kr_coef
      REAL_T  lam(2),rtmp(2)

      rtmp = one
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     compute corrective fluxes from edge states 
c     and perform conservative update
c
      if (model .eq. 2) then
!$omp parallel private(i,j,k,lam)
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax+1
                  l_kr_type = nint(krcoef(i,j,k,1))
                  l_kr_coef = krcoef(i,j,k,2)
                  call get_ld(lam,xflux(i,j,k,:),rtmp,muval,ncomps,
     &                 l_kr_type,l_kr_coef,zero)
                  xflux(i,j,k,1) = areax(i,j,k)*rhoval(1)*lam(1)/(lam(1)+lam(2))*
     &                 ucor(i,j,k)
                  xflux(i,j,k,2) = areax(i,j,k)*rhoval(2)*lam(2)/(lam(1)+lam(2))*
     &                 ucor(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax+1
               do i = imin,imax
                  l_kr_type = nint(krcoef(i,j,k,1))
                  l_kr_coef = krcoef(i,j,k,2)
                  call get_ld(lam,yflux(i,j,k,:),rtmp,muval,ncomps,
     &                 l_kr_type,l_kr_coef,zero)
                  yflux(i,j,k,1) = areay(i,j,k)*rhoval(1)*lam(1)/(lam(1)+lam(2))*
     &                 vcor(i,j,k)
                  yflux(i,j,k,2) = areay(i,j,k)*rhoval(2)*lam(2)/(lam(1)+lam(2))*
     &                 vcor(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
!$omp do
         do k = kmin,kmax+1
            do j = jmin,jmax
               do i = imin,imax
                  l_kr_type = nint(krcoef(i,j,k,1))
                  l_kr_coef = krcoef(i,j,k,2)
                  call get_ld(lam,zflux(i,j,k,:),rtmp,muval,ncomps,
     &                 l_kr_type,l_kr_coef,zero)
                  zflux(i,j,k,1) = areaz(i,j,k)*rhoval(1)*lam(1)/(lam(1)+lam(2))*
     &                 wcor(i,j,k)
                  zflux(i,j,k,2) = areaz(i,j,k)*rhoval(2)*lam(2)/(lam(1)+lam(2))*
     &                 wcor(i,j,k)
               end do
            end do
         end do
!$omp end do
!$omp end parallel
      else
         call bl_pd_abort("SYNC_RMN_FORCING: model not defined.")
      end if

!$omp parallel private(i,j,k,nc)
      do nc = 1,nscal
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  sync(i,j,k,nc) = sync(i,j,k,nc) + (
     &                 xflux(i+1,j,k,nc)-xflux(i,j,k,nc) +
     &                 yflux(i,j+1,k,nc)-yflux(i,j,k,nc) +
     &                 zflux(i,j,k+1,nc)-zflux(i,j,k,nc) )/vol(i,j,k)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end

      subroutine trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(xx),uad,DIMS(uad),
     &     lo,hi,n,xbc,eps)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the x direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      integer DIMDEC(s)
      REAL_T s(DIMV(s))
      integer DIMDEC(xx)
      integer DIMDEC(uad)
      REAL_T xlo(DIMV(xx))
      REAL_T xhi(DIMV(xx))
      REAL_T uad(DIMV(uad))
      REAL_T eps
      integer lo(SDIM), hi(SDIM)
      integer n
      integer xbc(SDIM,2)

      REAL_T stx
      logical ltest
      integer j,k
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     -------------- the lower x boundary
c
      if (xbc(1,1).eq.EXT_DIR) then
         if ( n .eq. XVEL ) then
            do k = kmin-1,kmax+1
               do j = jmin-1,jmax+1
                  if (uad(imin,j,k) .ge. zero) then
                     xlo(imin,j,k) = s(imin-1,j,k)
                     xhi(imin,j,k) = s(imin-1,j,k)
                  else
                     xlo(imin,j,k) = xhi(imin,j,k)
                  endif
               end do
            end do
         else
            do k = kmin-1,kmax+1
               do j = jmin-1,jmax+1
                  ltest = uad(imin,j,k).le.eps
                  stx   = cvmgt(xhi(imin,j,k),s(imin-1,j,k),ltest)
                  xlo(imin,j,k) = stx
                  xhi(imin,j,k) = stx
               end do
            end do
         end if
      else if (xbc(1,1).eq.FOEXTRAP.or.xbc(1,1).eq.HOEXTRAP
     &        .or.xbc(1,1).eq.REFLECT_EVEN) then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               xlo(imin,j,k) = xhi(imin,j,k)
            end do
         end do
      else if (xbc(1,1).eq.REFLECT_ODD) then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               xhi(imin,j,k) = zero
               xlo(imin,j,k) = zero
            end do
         end do
      end if
c
c     -------------- the upper x boundary
c
      if (xbc(1,2).eq.EXT_DIR) then
         if ( n .eq. XVEL ) then
            do k = kmin-1,kmax+1
               do j = jmin-1,jmax+1
                  if (uad(imax+1,j,k) .le. zero) then
                     xlo(imax+1,j,k) = s(imax+1,j,k)
                     xhi(imax+1,j,k) = s(imax+1,j,k)
                  else
                     xhi(imax+1,j,k) = xlo(imax+1,j,k)
                  endif
               end do
            end do
         else
            do k = kmin-1,kmax+1
               do j = jmin-1,jmax+1
                  ltest = uad(imax+1,j,k).ge.-eps
                  stx   = cvmgt(xlo(imax+1,j,k),s(imax+1,j,k),ltest)
                  xlo(imax+1,j,k) = stx
                  xhi(imax+1,j,k) = stx
               end do
            end do
         end if
      else if (xbc(1,2).eq.FOEXTRAP.or.xbc(1,2).eq.HOEXTRAP
     &        .or.xbc(1,2).eq.REFLECT_EVEN) then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               xhi(imax+1,j,k) = xlo(imax+1,j,k)
            end do
         end do
      else if (xbc(1,2).eq.REFLECT_ODD) then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               xhi(imax+1,j,k) = zero
               xlo(imax+1,j,k) = zero
            end do
         end do
      end if

      end

      subroutine trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(yy),vad,DIMS(vad),
     &     lo,hi,n,ybc,eps)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the y direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      integer DIMDEC(s)
      REAL_T s(DIMV(s))
      integer DIMDEC(yy)
      integer DIMDEC(vad)
      REAL_T ylo(DIMV(yy))
      REAL_T yhi(DIMV(yy))
      REAL_T vad(DIMV(vad))
      REAL_T eps
      integer lo(SDIM), hi(SDIM)
      integer n
      integer ybc(SDIM,2)

      REAL_T sty
      logical ltest
      integer i,k
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     -------------- the lower y boundary
c
      if (ybc(2,1).eq.EXT_DIR) then
         if ( n .eq. YVEL ) then
            do k = kmin-1,kmax+1
               do i = imin-1,imax+1
                  if (vad(i,jmin,k) .ge. zero) then
                     ylo(i,jmin,k) = s(i,jmin-1,k)
                     yhi(i,jmin,k) = s(i,jmin-1,k)
                  else
                     ylo(i,jmin,k) = yhi(i,jmin,k)
                  endif
               end do
            end do
         else
            do k = kmin-1,kmax+1
               do i = imin-1,imax+1
                  ltest = vad(i,jmin,k).le.eps
                  sty   = cvmgt(yhi(i,jmin,k),s(i,jmin-1,k),ltest)
                  ylo(i,jmin,k) = sty
                  yhi(i,jmin,k) = sty
               end do
            end do
         end if
      else if (ybc(2,1).eq.FOEXTRAP.or.ybc(2,1).eq.HOEXTRAP
     &        .or.ybc(2,1).eq.REFLECT_EVEN) then
         do k = kmin-1,kmax+1
            do i = imin-1,imax+1
               ylo(i,jmin,k) = yhi(i,jmin,k)
            end do
         end do
      else if (ybc(2,1).eq.REFLECT_ODD) then
         do k = kmin-1,kmax+1
            do i = imin-1,imax+1
               yhi(i,jmin,k) = zero
               ylo(i,jmin,k) = zero
            end do
         end do
      end if
c
c     -------------- the upper y boundary
c
      if (ybc(2,2).eq.EXT_DIR) then
         if ( n .eq. YVEL ) then
            do k = kmin-1,kmax+1
               do i = imin-1,imax+1
                  if (vad(i,jmax+1,k) .le. zero) then
                     ylo(i,jmax+1,k) = s(i,jmax+1,k)
                     yhi(i,jmax+1,k) = s(i,jmax+1,k)
                  else
                     yhi(i,jmax+1,k) = ylo(i,jmax+1,k)
                  endif
               end do
            end do
         else
            do k = kmin-1,kmax+1
               do i = imin-1,imax+1
                  ltest = vad(i,jmax+1,k).ge.-eps
                  sty   = cvmgt(ylo(i,jmax+1,k),s(i,jmax+1,k),ltest)
                  ylo(i,jmax+1,k) = sty
                  yhi(i,jmax+1,k) = sty
               end do
            end do
         end if
      else if (ybc(2,2).eq.FOEXTRAP.or.ybc(2,2).eq.HOEXTRAP
     &        .or.ybc(2,2).eq.REFLECT_EVEN) then
         do k = kmin-1,kmax+1
            do i = imin-1,imax+1
               yhi(i,jmax+1,k) = ylo(i,jmax+1,k)
            end do
         end do
      else if (ybc(2,2).eq.REFLECT_ODD) then
         do k = kmin-1,kmax+1
            do i = imin-1,imax+1
               ylo(i,jmax+1,k) = zero
               yhi(i,jmax+1,k) = zero
            end do
         end do
      end if

      end

      subroutine trans_zbc(
     &     s,DIMS(s),
     &     zlo,zhi,DIMS(zz),wad,DIMS(wad),
     &     lo,hi,n,zbc,eps)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the z direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      integer DIMDEC(s)
      REAL_T s(DIMV(s))
      integer DIMDEC(zz)
      integer DIMDEC(wad)
      REAL_T zlo(DIMV(zz))
      REAL_T zhi(DIMV(zz))
      REAL_T wad(DIMV(wad))
      REAL_T eps
      integer lo(SDIM), hi(SDIM)
      integer n
      integer zbc(SDIM,2)

      REAL_T stz
      logical ltest
      integer i,j
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     -------------- the lower z boundary
c
      if (zbc(3,1).eq.EXT_DIR) then
         if ( n .eq. ZVEL ) then
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
               if (wad(i,j,kmin) .ge. zero) then
                  zhi(i,j,kmin) = s(i,j,kmin-1)
                  zlo(i,j,kmin) = s(i,j,kmin-1)
               else
                  zlo(i,j,kmin) = zhi(i,j,kmin)
               endif
             end do
            end do
         else
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  ltest = wad(i,j,kmin).le.eps
                  stz   = cvmgt(zhi(i,j,kmin),s(i,j,kmin-1),ltest)
                  zhi(i,j,kmin) = stz
                  zlo(i,j,kmin) = stz
               end do
            end do
         end if
      else if (zbc(3,1).eq.FOEXTRAP.or.zbc(3,1).eq.HOEXTRAP
     &        .or.zbc(3,1).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               zlo(i,j,kmin) = zhi(i,j,kmin)
            end do
         end do
      else if (zbc(3,1).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               zhi(i,j,kmin) = zero
               zlo(i,j,kmin) = zero
            end do
         end do
      end if
c
c     -------------- the upper z boundary
c
      if (zbc(3,2).eq.EXT_DIR) then
         if ( n .eq. ZVEL ) then
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
               if (wad(i,j,kmax+1) .le. zero) then
                  zlo(i,j,kmax+1) = s(i,j,kmax+1)
                  zhi(i,j,kmax+1) = s(i,j,kmax+1)
               else
                  zhi(i,j,kmax+1) = zlo(i,j,kmax+1)
               endif
             end do
            end do
         else
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  ltest = wad(i,j,kmax+1).ge.-eps
                  stz   = cvmgt(zlo(i,j,kmax+1),s(i,j,kmax+1),ltest)
                  zhi(i,j,kmax+1) = stz
                  zlo(i,j,kmax+1) = stz
               end do
            end do
         end if
      else if (zbc(3,2).eq.FOEXTRAP.or.zbc(3,2).eq.HOEXTRAP
     &        .or.zbc(3,2).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               zhi(i,j,kmax+1) = zlo(i,j,kmax+1)
            end do
         end do
      else if (zbc(3,2).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               zlo(i,j,kmax+1) = zero
               zhi(i,j,kmax+1) = zero
            end do
         end do
      end if

      end

      subroutine FORT_SLOPES( dir,
     &     s,DIMS(s),
     &     slx,sly,slz,DIMS(sl),
     &     lo,hi,slxscr,slyscr,slzscr,bc)
c 
c     this subroutine computes first or forth order slopes of
c     a 3D scalar field.
c
c     (dir) is used to eliminate calculating extra slopes in transvel
c
c     Boundary conditions on interior slopes are handled automatically
c     by the ghost cells
c
c     Boundary conditions on EXT_DIR and HOEXTRAP slopes are implemented
c     by setting them to zero outside of the domain and using a
c     one-sided derivative from the interior
c
      implicit none

#include "GODCOMM_F.H"

      integer dir
      integer DIMDEC(s)
      REAL_T     s(DIMV(s))
      integer DIMDEC(sl)
      REAL_T   slx(DIMV(sl))
      REAL_T   sly(DIMV(sl))
      REAL_T   slz(DIMV(sl))
      integer lo(SDIM), hi(SDIM)
      REAL_T slxscr(DIM1(s), 4)
      REAL_T slyscr(DIM2(s), 4)
      REAL_T slzscr(DIM3(s), 4)
      integer bc(SDIM,2)

      integer imin,jmin,kmin,imax,jmax,kmax,i,j,k
      integer ng
      REAL_T dpls,dmin,ds
      REAL_T del,slim,sflg
      integer cen,lim,flag,fromm
      parameter( cen = 1   )
      parameter( lim = 2   )
      parameter( flag = 3  )
      parameter( fromm = 4 )

      ng = lo(1) - ARG_L1(s)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     Added to prevent underflow for small s values.
c
!$omp parallel do private(i,j,k)
      do k = lo(3)-ng, hi(3)+ng
         do j = lo(2)-ng, hi(2)+ng
            do i = lo(1)-ng, hi(1)+ng
               s(i,j,k) = cvmgt(s(i,j,k), zero, abs(s(i,j,k)).gt.1.0D-20)
            end do
         end do
      end do
!$omp end parallel do
c
c     COMPUTE 0TH order slopes
c
      if (slope_order.eq.1) then
        if (ng .lt. 1) then
	   call bl_pd_abort("FORT_SLOPES: not enough bndry cells for 1st order")
        end if
        do k = kmin-1, kmax+1
           do j = jmin-1, jmax+1 
              do i = imin-1, imax+1
                 slx(i,j,k) = zero
                 sly(i,j,k) = zero
                 slz(i,j,k) = zero
              end do
           end do
        end do
        return
      end if
c
c     COMPUTE 2ND order slopes
c
      if (slope_order.eq.2) then
        if (ng .lt. 2) then
	   call bl_pd_abort("FORT_SLOPES: not enough bndry cells for 2nd order")
        end if
c
c     ------------------------ x slopes
c
        if ( (dir.eq.XVEL) .or. (dir.eq.ALL) ) then
         if (use_unlimited_slopes) then
            do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  slx(i,j,k)= half*(s(i+1,j,k) - s(i-1,j,k))
               end do
            end do
            end do
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do j = jmin-1, jmax+1
                  slx(imin-1,j,k) = zero
                  slx(imin  ,j,k) = (s(imin+1,j,k)+three*s(imin,j,k)-four*s(imin-1,j,k))/three
               end do
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do j = jmin-1, jmax+1
                  slx(imax+1,j,k) = zero
                  slx(imax  ,j,k) = -(s(imax-1,j,k)+three*s(imax,j,k)-four*s(imax+1,j,k))/three
               end do
               end do
            end if
         else
            do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i+1,j,k) - s(i-1,j,k))
                  dpls =  2*(s(i+1,j,k) - s(i  ,j,k))
                  dmin =  2*(s(i  ,j,k) - s(i-1,j,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slx(i,j,k)= sflg*min(slim,abs(del))
               end do
            end do
            end do
            
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do j = jmin-1, jmax+1
                  slx(imin-1,j,k) = zero
                  del  = (s(imin+1,j,k)+three*s(imin,j,k)-four*s(imin-1,j,k))/three
                  dpls = 2*(s(imin+1,j,k) - s(imin  ,j,k))
                  dmin = 2*(s(imin  ,j,k) - s(imin-1,j,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slx(imin,j,k)= sflg*min(slim,abs(del))
               end do
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do j = jmin-1, jmax+1
                  slx(imax+1,j,k) = zero
                  del  = -(s(imax-1,j,k)+three*s(imax,j,k)-four*s(imax+1,j,k))/three
                  dpls = 2*(s(imax+1,j,k) - s(imax  ,j,k))
                  dmin = 2*(s(imax  ,j,k) - s(imax-1,j,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slx(imax,j,k)= sflg*min(slim,abs(del))
               end do
               end do
            end if
          end if
         end if
c
c     ------------------------ y slopes
c
        if ( (dir.eq.YVEL) .or. (dir.eq.ALL) ) then
         if (use_unlimited_slopes) then
            do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  sly(i,j,k) = half*(s(i,j+1,k)-s(i,j-1,k))
               end do
            end do
            end do
            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do i = imin-1, imax+1
                  sly(i,jmin-1,k) = zero
                  sly(i,jmin  ,k) = (s(i,jmin+1,k)+three*s(i,jmin,k)-four*s(i,jmin-1,k))/three
               end do
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do i = imin-1, imax+1
                  sly(i,jmax+1,k) = zero
                  sly(i,jmax  ,k) = -(s(i,jmax-1,k)+three*s(i,jmax,k)-four*s(i,jmax+1,k))/three
               end do
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do k = kmin-1,kmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i,j+1,k) - s(i,j-1,k))
                  dpls =  2*(s(i,j+1,k) - s(i,j  ,k))
                  dmin =  2*(s(i,j  ,k) - s(i,j-1,k))
                  slim = min(abs(dpls),abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  sly(i,j,k)= sflg*min(slim,abs(del))
               end do
               end do
            end do

            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do i = imin-1, imax+1
                  sly(i,jmin-1,k) = zero
                  del  = (s(i,jmin+1,k)+three*s(i,jmin,k)-four*s(i,jmin-1,k))/three
                  dpls = 2*(s(i,jmin+1,k) - s(i,jmin  ,k))
                  dmin = 2*(s(i,jmin  ,k) - s(i,jmin-1,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  sly(i,jmin,k)= sflg*min(slim,abs(del))
               end do
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do i = imin-1, imax+1
                  sly(i,jmax+1,k) = zero
                  del  = -(s(i,jmax-1,k)+three*s(i,jmax,k)-four*s(i,jmax+1,k))/three
                  dpls = 2*(s(i,jmax+1,k) - s(i,jmax  ,k))
                  dmin = 2*(s(i,jmax  ,k) - s(i,jmax-1,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  sly(i,jmax,k)= sflg*min(slim,abs(del))
               end do
               end do
            end if
         end if
        end if

c
c     ------------------------ z slopes
c
        if ( (dir.eq.ZVEL) .or. (dir.eq.ALL) ) then
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               do k = kmin-1,kmax+1
                  slz(i,j,k) = half*(s(i,j,k+1)-s(i,j,k-1))
               end do
            end do
            end do
            if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
               do i = imin-1, imax+1
                  slz(i,j,kmin-1) = zero
                  slz(i,j,kmin  ) = (s(i,j,kmin+1)+three*s(i,j,kmin)-four*s(i,j,kmin-1))/three
               end do
               end do
            end if
            if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
               do i = imin-1, imax+1
                  slz(i,j,kmax+1) = zero
                  slz(i,j,kmax  ) = -(s(i,j,kmax-1)+three*s(i,j,kmax)-four*s(i,j,kmax+1))/three
               end do
               end do
            end if
         else
            do k = kmin-1,kmax+1
               do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i,j,k+1) - s(i,j,k-1))
                  dpls =  2*(s(i,j,k+1) - s(i,j,k  ))
                  dmin =  2*(s(i,j,k  ) - s(i,j,k-1))
                  slim = min(abs(dpls),abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slz(i,j,k)= sflg*min(slim,abs(del))
               end do
               end do
            end do

            if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
               do i = imin-1, imax+1
                  slz(i,j,kmin-1) = zero
                  del  = (s(i,j,kmin+1)+three*s(i,j,kmin)-four*s(i,j,kmin-1))/three
                  dpls = 2*(s(i,j,kmin+1) - s(i,j,kmin  ))
                  dmin = 2*(s(i,j,kmin  ) - s(i,j,kmin-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slz(i,j,kmin)= sflg*min(slim,abs(del))
               end do
               end do
            end if
            if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
               do i = imin-1, imax+1
                  slz(i,j,kmax+1) = zero
                  del  = -(s(i,j,kmax-1)+three*s(i,j,kmax)-four*s(i,j,kmax+1))/three
                  dpls = 2*(s(i,j,kmax+1) - s(i,j,kmax  ))
                  dmin = 2*(s(i,j,kmax  ) - s(i,j,kmax-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slz(i,j,kmax)= sflg*min(slim,abs(del))
               end do
               end do
            end if
         end if
        end if
c
c ... end, if slope_order .eq. 2
c
      end if
c
c     COMPUTE 4TH order slopes
c
      if (slope_order.eq.4)then
        if (ng .lt. 3) then
	   call bl_pd_abort("SLOPE_3D: not enough bndry cells for 4th order")
        end if
c
c     ------------------------ x slopes
c
        if ( (dir.eq.XVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do k = kmin-1,kmax+1
                 do j = jmin-1,jmax+1
                    do i = imin-2,imax+2
                       slxscr(i,cen)  = half*(s(i+1,j,k)-s(i-1,j,k))
                    end do
                    do i = imin-1,imax+1
                       slx(i,j,k) = 2 * two3rd * slxscr(i,cen) -
     &                      sixth * (slxscr(i+1,cen) + slxscr(i-1,cen))
                    end do
                 end do
              end do
              
              if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                 do k = kmin-1,kmax+1
                    do j = jmin-1, jmax+1
                       slx(imin,j,k) = -sixteen/fifteen*s(imin-1,j,k) + half*s(imin,j,k) + 
     &                      two3rd*s(imin+1,j,k) - tenth*s(imin+2,j,k)
                       slx(imin-1,j,k) = zero
                    end do
                 end do
              end if
              if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                 do k = kmin-1,kmax+1
                    do j = jmin-1, jmax+1
                       slx(imax,j,k) = -( -sixteen/fifteen*s(imax+1,j,k) + half*s(imax,j,k) + 
     &                      two3rd*s(imax-1,j,k) - tenth*s(imax-2,j,k) )
                       slx(imax+1,j,k) = zero
                    end do
                 end do
              end if
           else
!$omp parallel do private(i,j,k,slxscr,dmin,dpls,ds,del,slim,sflg)
              do k = kmin-1,kmax+1
                 do j = jmin-1,jmax+1 
                    do i = imin-2,imax+2
                       dmin           =  2*(s(i,  j,k)-s(i-1,j,k))
                       dpls           =  2*(s(i+1,j,k)-s(i,  j,k))
                       slxscr(i,cen)  = half*(s(i+1,j,k)-s(i-1,j,k))
                       slxscr(i,lim)  = min(abs(dmin),abs(dpls))
                       slxscr(i,lim)  = cvmgp(slxscr(i,lim),zero,dpls*dmin)
                       slxscr(i,flag) = sign(one,slxscr(i,cen))
                       slxscr(i,fromm)= slxscr(i,flag)*
     &                      min(slxscr(i,lim),abs(slxscr(i,cen)))
                    end do
                    do i = imin-1,imax+1
                       ds = 2 * two3rd * slxscr(i,cen) - 
     &                      sixth * (slxscr(i+1,fromm) + slxscr(i-1,fromm))
                       slx(i,j,k) = slxscr(i,flag)*min(abs(ds),slxscr(i,lim))
                    end do

                    if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                       del  = -sixteen/fifteen*s(imin-1,j,k) + half*s(imin,j,k) + 
     &                      two3rd*s(imin+1,j,k) -  tenth*s(imin+2,j,k)
                       dmin = 2*(s(imin  ,j,k)-s(imin-1,j,k))
                       dpls = 2*(s(imin+1,j,k)-s(imin  ,j,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       slx(imin-1,j,k) = zero
                       slx(imin,  j,k) = sflg*min(slim,abs(del))

c                      Recalculate the slope at imin+1 using the revised slxscr(imin,fromm)
                       slxscr(imin,fromm) = slx(imin,j,k)
                       ds = 2 * two3rd * slxscr(imin+1,cen) -
     $                    sixth * (slxscr(imin+2,fromm) + slxscr(imin,fromm))
                       slx(imin+1,j,k) = slxscr(imin+1,flag)*min(abs(ds),slxscr(imin+1,lim))
                    end if

                    if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                       del  = -( -sixteen/fifteen*s(imax+1,j,k) + half*s(imax,j,k) + 
     &                      two3rd*s(imax-1,j,k) - tenth*s(imax-2,j,k) )
                       dmin = 2*(s(imax  ,j,k)-s(imax-1,j,k))
                       dpls = 2*(s(imax+1,j,k)-s(imax  ,j,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       slx(imax,  j,k) = sflg*min(slim,abs(del))
                       slx(imax+1,j,k) = zero

c                      Recalculate the slope at imax-1 using the revised slxscr(imax,fromm)
                       slxscr(imax,fromm) = slx(imax,j,k)
                       ds = 2 * two3rd * slxscr(imax-1,cen) -
     $                    sixth * (slxscr(imax-2,fromm) + slxscr(imax,fromm))
                       slx(imax-1,j,k) = slxscr(imax-1,flag)*min(abs(ds),slxscr(imax-1,lim))
                    end if
                 end do
              end do
!$omp end parallel do
           end if
        end if
c
c     ------------------------ y slopes
c
        if ( (dir.eq.YVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do k = kmin-1,kmax+1
                 do i = imin-1,imax+1
                    do j = jmin-2,jmax+2
                       slyscr(j,cen)  = half*(s(i,j+1,k)-s(i,j-1,k))
                    end do
                    do j = jmin-1,jmax+1
                       sly(i,j,k) = 2 * two3rd * slyscr(j,cen) -
     &                      sixth * (slyscr(j+1,cen) + slyscr(j-1,cen))
                    end do
                 end do
              end do
              
              if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                 do k = kmin-1,kmax+1
                    do i = imin-1, imax+1
                       sly(i,jmin-1,k) = zero
                       sly(i,jmin,k) = -sixteen/fifteen*s(i,jmin-1,k) + half*s(i,jmin,k) + 
     &                      two3rd*s(i,jmin+1,k) - tenth*s(i,jmin+2,k)
                    end do
                 end do
              end if
              if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                 do k = kmin-1,kmax+1
                    do i = imin-1, imax+1
                       sly(i,jmax,k) = -( -sixteen/fifteen*s(i,jmax+1,k) + half*s(i,jmax,k) + 
     &                      two3rd*s(i,jmax-1,k) - tenth*s(i,jmax-2,k) )
                       sly(i,jmax+1,k) = zero
                    end do
                 end do
              end if
           else
!$omp parallel do private(i,j,k,slyscr,dmin,dpls,ds,del,slim,sflg)
              do k = kmin-1,kmax+1
                 do i = imin-1,imax+1 
                    do j = jmin-2,jmax+2
                       dmin           =  2*(s(i,j,  k)-s(i,j-1,k))
                       dpls           =  2*(s(i,j+1,k)-s(i,j,  k))
                       slyscr(j,cen)  = half*(s(i,j+1,k)-s(i,j-1,k))
                       slyscr(j,lim)  = min(abs(dmin),abs(dpls))
                       slyscr(j,lim)  = cvmgp(slyscr(j,lim),zero,dpls*dmin)
                       slyscr(j,flag) = sign(one,slyscr(j,cen))
                       slyscr(j,fromm)= slyscr(j,flag)*
     &                      min(slyscr(j,lim),abs(slyscr(j,cen)))
                    end do
                    do j = jmin-1,jmax+1
                       ds = 2 * two3rd * slyscr(j,cen) - 
     &                      sixth * (slyscr(j+1,fromm) + slyscr(j-1,fromm))
                       sly(i,j,k) = slyscr(j,flag)*min(abs(ds),slyscr(j,lim))
                    end do
c
                    if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                       del  = -sixteen/fifteen*s(i,jmin-1,k) + half*s(i,jmin,k) + 
     &                      two3rd*s(i,jmin+1,k) - tenth*s(i,jmin+2,k)
                       dmin = 2*(s(i,jmin  ,k)-s(i,jmin-1,k))
                       dpls = 2*(s(i,jmin+1,k)-s(i,jmin  ,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       sly(i,jmin-1,k) = zero
                       sly(i,jmin,  k) = sflg*min(slim,abs(del))

c                      Recalculate the slope at jmin+1 using the revised slyscr(jmin,fromm)
                       slyscr(jmin,fromm) = sly(i,jmin,k)
                       ds = 2 * two3rd * slyscr(jmin+1,cen) -
     $                    sixth * (slyscr(jmin+2,fromm) + slyscr(jmin,fromm))
                       sly(i,jmin+1,k) = slyscr(jmin+1,flag)*min(abs(ds),slyscr(jmin+1,lim))
                    end if
                    if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                       del  = -( -sixteen/fifteen*s(i,jmax+1,k) + half*s(i,jmax,k) +
     &                      two3rd*s(i,jmax-1,k) - tenth*s(i,jmax-2,k) )
                       dmin = 2*(s(i,jmax  ,k)-s(i,jmax-1,k))
                       dpls = 2*(s(i,jmax+1,k)-s(i,jmax  ,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       sly(i,jmax, k)  = sflg*min(slim,abs(del))
                       sly(i,jmax+1,k) = zero

c                      Recalculate the slope at jmax-1 using the revised slyscr(jmax,fromm)
                       slyscr(jmax,fromm) = sly(i,jmax,k)
                       ds = 2 * two3rd * slyscr(jmax-1,cen) -
     $                    sixth * (slyscr(jmax-2,fromm) + slyscr(jmax,fromm))
                       sly(i,jmax-1,k) = slyscr(jmax-1,flag)*min(abs(ds),slyscr(jmax-1,lim))
                    end if
                 end do
              end do
!$omp end parallel do
           end if
        end if
c
c     ------------------------ z slopes
c
        if ( (dir.eq.ZVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do j = jmin-1,jmax+1
                 do i = imin-1,imax+1
                    do k = kmin-2,kmax+2
                       slzscr(k,cen)  = half*(s(i,j,k+1)-s(i,j,k-1))
                    end do
                    do k = kmin-1,kmax+1
                       slz(i,j,k) = 2 * two3rd * slzscr(k,cen) -
     &                      sixth * (slzscr(k+1,cen) + slzscr(k-1,cen))
                    end do
                 end do
              end do
              
              if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
                 do j = jmin-1, jmax+1
                    do i = imin-1, imax+1
                       slz(i,j,kmin-1) = zero
                       slz(i,j,kmin) = -sixteen/fifteen*s(i,j,kmin-1) + half*s(i,j,kmin) + 
     &                      two3rd*s(i,j,kmin+1) - tenth*s(i,j,kmin+2)
                    end do
                 end do
              end if
              if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
                 do j = jmin-1, jmax+1
                    do i = imin-1, imax+1
                       slz(i,j,kmax) = -( -sixteen/fifteen*s(i,j,kmax+1) + half*s(i,j,kmax) + 
     &                      two3rd*s(i,j,kmax-1) - tenth*s(i,j,kmax-2) )
                       slz(i,j,kmax+1) = zero
                    end do
                 end do
              end if
           else
!$omp parallel do private(i,j,k,slzscr,dmin,dpls,ds,del,slim,sflg)
              do j = jmin-1,jmax+1
                 do i = imin-1,imax+1
                    do k = kmin-2,kmax+2
                       dmin           =  2*(s(i,j,k  )-s(i,j,k-1))
                       dpls           =  2*(s(i,j,k+1)-s(i,j,k  ))
                       slzscr(k,cen)  = half*(s(i,j,k+1)-s(i,j,k-1))
                       slzscr(k,lim)  = min(abs(dmin),abs(dpls))
                       slzscr(k,lim)  = cvmgp(slzscr(k,lim),zero,dpls*dmin)
                       slzscr(k,flag) = sign(one,slzscr(k,cen))
                       slzscr(k,fromm)= slzscr(k,flag)*
     &                      min(slzscr(k,lim),abs(slzscr(k,cen)))
                    end do
                    do k = kmin-1,kmax+1
                       ds = 2 * two3rd * slzscr(k,cen) -
     &                      sixth * (slzscr(k+1,fromm) + slzscr(k-1,fromm))
                       slz(i,j,k) = slzscr(k,flag)*min(abs(ds),slzscr(k,lim))
                    end do
c
                    if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
                       del  = -sixteen/fifteen*s(i,j,kmin-1) + half*s(i,j,kmin) +
     &                      two3rd*s(i,j,kmin+1) - tenth*s(i,j,kmin+2)
                       dmin = 2*(s(i,j,kmin  )-s(i,j,kmin-1))
                       dpls = 2*(s(i,j,kmin+1)-s(i,j,kmin  ))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       slz(i,j,kmin-1) = zero
                       slz(i,j,kmin  ) = sflg*min(slim,abs(del))

c                      Recalculate the slope at jmin+1 using the revised slzscr(kmin,fromm)
                       slzscr(kmin,fromm) = slz(i,j,kmin)
                       ds = 2 * two3rd * slzscr(kmin+1,cen) -
     $                    sixth * (slzscr(kmin+2,fromm) + slzscr(kmin,fromm))
                       slz(i,j,kmin+1) = slzscr(kmin+1,flag)*min(abs(ds),slzscr(kmin+1,lim))
                    end if
                    if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
                       del  = sixteen/fifteen*s(i,j,kmax+1) - half*s(i,j,kmax)
     &                      - two3rd*s(i,j,kmax-1) + tenth*s(i,j,kmax-2)
                       dmin = 2*(s(i,j,kmax  )-s(i,j,kmax-1))
                       dpls = 2*(s(i,j,kmax+1)-s(i,j,kmax  ))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       slz(i,j,kmax  ) = sflg*min(slim,abs(del))
                       slz(i,j,kmax+1) = zero

c                      Recalculate the slope at jmax-1 using the revised slzscr(kmax,fromm)
                       slzscr(kmax,fromm) = slz(i,j,kmax)
                       ds = 2 * two3rd * slzscr(kmax-1,cen) -
     $                    sixth * (slzscr(kmax-2,fromm) + slzscr(kmax,fromm))
                       slz(i,j,kmax-1) = slzscr(kmax-1,flag)*min(abs(ds),slzscr(kmax-1,lim))
                    end if
                 end do
              end do
!$omp end parallel do
           end if
        end if
c
c ... end, if slope_order .eq. 4
c
      end if

      end

      subroutine FORT_CONVSCALMINMAX(s,DIMS(s),sn,DIMS(sn),
     &                               smin,smax,DIMS(smin),lo,hi,bc)
c
c     correct an advected field for under/over shoots
c
      integer  i, j, k, imin, imax, jmin, jmax, kmin, kmax
      integer  DIMDEC(s)
      integer  DIMDEC(sn)
      integer  DIMDEC(smin)
      integer  lo(SDIM), hi(SDIM)
      integer  bc(SDIM,2)
      REAL_T   s(DIMV(s))
      REAL_T   sn(DIMV(sn))
      integer  km, kk, kp
      REAL_T   smn, smx
      REAL_T   smin(DIM12(smin),0:2)
      REAL_T   smax(DIM12(smin),0:2)

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)
c
c     correct the 8 corners
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imin-1,jmin-1,kmin-1) = s(imin,jmin,kmin)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imax+1,jmin-1,kmin-1) = s(imax,jmin,kmin)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imin-1,jmax+1,kmin-1) = s(imin,jmax,kmin)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imax+1,jmax+1,kmin-1) = s(imax,jmax,kmin)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imin-1,jmin-1,kmax+1) = s(imin,jmin,kmax)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imax+1,jmin-1,kmax+1) = s(imax,jmin,kmax)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imin-1,jmax+1,kmax+1) = s(imin,jmax,kmax)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imax+1,jmax+1,kmax+1) = s(imax,jmax,kmax)
      end if
c
c     correct the 12 edges
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
          s(imin-1,jmin-1,k) = s(imin,jmin,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
          s(imax+1,jmin-1,k) = s(imax,jmin,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
          s(imin-1,jmax+1,k) = s(imin,jmax,k)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
          s(imax+1,jmax+1,k) = s(imax,jmax,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
          s(imin-1,j,kmin-1) = s(imin,j,kmin)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
          s(imax+1,j,kmin-1) = s(imax,j,kmin)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
          s(imin-1,j,kmax+1) = s(imin,j,kmax)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
          s(imax+1,j,kmax+1) = s(imax,j,kmax)
        end do
      end if
      if (bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmin-1,kmin-1) = s(i,jmin,kmin)
        end do
      end if
      if (bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmax+1,kmin-1) = s(i,jmax,kmin)
        end do
      end if
      if (bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmin-1,kmax+1) = s(i,jmin,kmax)
        end do
      end if
      if (bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmax+1,kmax+1) = s(i,jmax,kmax)
        end do
      end if
c
c     ::::: compute min/max a slab at a time
c     ::::: compute min and max of neighbors on kmin-1 slab
c
      km = 0
      kk = 1
      kp = 2

      k = kmin-1
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,km) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,km) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
         end do         
      end do
c
c     ::::: compute min and max of neighbors on kmin slab
c
      k = kmin
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,kk) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,kk) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
         end do         
      end do

      do k = kmin, kmax
c
c        ::::: compute min and max of neighbors on k+1 slab
c
         do j = jmin, jmax     
            do i = imin, imax   
               smin(i,j,kp) = min(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &              s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &              s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
               smax(i,j,kp) = max(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &              s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &              s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
c
c        ::::: compute min/max of cell
c
               smn = min(smin(i,j,km),smin(i,j,kk),smin(i,j,kp))
               smx = max(smax(i,j,km),smax(i,j,kk),smax(i,j,kp))
               sn(i,j,k) = max(sn(i,j,k),smn)
               sn(i,j,k) = min(sn(i,j,k),smx)
               
            end do
         end do
c
c        ::::: roll indices for next slab
c
         km = mod(km+1,3)
         kk = mod(kk+1,3)
         kp = mod(kp+1,3)
      end do

      end

      subroutine FORT_CONSSCALMINMAX(s,sn,DIMS(s),
     &                               smin,smax,DIMS(smin),lo,hi,bc)
c
c     correct an conservatively-advected field for under/over shoots
c
      integer  i, j, k, imin, imax, jmin, jmax, kmin, kmax
      integer  DIMDEC(s)
      integer  DIMDEC(smin)
      integer  lo(SDIM), hi(SDIM)
      integer  bc(SDIM,2)
      REAL_T      s(DIMV(s))
      REAL_T     sn(DIMV(s))
      integer  km, kk, kp
      REAL_T   smn, smx
      REAL_T   smin(DIM12(smin),0:2)
      REAL_T   smax(DIM12(smin),0:2)

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)
c
c     correct the 8 corners
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
           s(imin-1,jmin-1,kmin-1) =   s(imin,jmin,kmin)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
           s(imax+1,jmin-1,kmin-1) =   s(imax,jmin,kmin)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
           s(imin-1,jmax+1,kmin-1) =   s(imin,jmax,kmin)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
           s(imax+1,jmax+1,kmin-1) =   s(imax,jmax,kmin)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
           s(imin-1,jmin-1,kmax+1) =   s(imin,jmin,kmax)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
           s(imax+1,jmin-1,kmax+1) =   s(imax,jmin,kmax)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
           s(imin-1,jmax+1,kmax+1) =   s(imin,jmax,kmax)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
           s(imax+1,jmax+1,kmax+1) =   s(imax,jmax,kmax)
      end if
c
c     correct the 12 edges
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
            s(imin-1,jmin-1,k) =   s(imin,jmin,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
            s(imax+1,jmin-1,k) =   s(imax,jmin,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
            s(imin-1,jmax+1,k) =   s(imin,jmax,k)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
            s(imax+1,jmax+1,k) =   s(imax,jmax,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
            s(imin-1,j,kmin-1) =   s(imin,j,kmin)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
            s(imax+1,j,kmin-1) =   s(imax,j,kmin)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
            s(imin-1,j,kmax+1) =   s(imin,j,kmax)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
            s(imax+1,j,kmax+1) =   s(imax,j,kmax)
        end do
      end if
      if (bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
            s(i,jmin-1,kmin-1) =   s(i,jmin,kmin)
        end do
      end if
      if (bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
            s(i,jmax+1,kmin-1) =   s(i,jmax,kmin)
        end do
      end if
      if (bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
            s(i,jmin-1,kmax+1) =   s(i,jmin,kmax)
        end do
      end if
      if (bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
            s(i,jmax+1,kmax+1) =   s(i,jmax,kmax)
        end do
      end if
c
c     ::::: compute min/max a slab at a time
c     ::::: compute min and max of neighbors on kmin-1 slab
c
      km = 0
      kk = 1
      kp = 2

      k = kmin-1
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,km) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,km) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k)) 
         end do         
      end do
c
c     ::::: compute min and max of neighbors on kmin slab
c
      k = kmin
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,kk) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,kk) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
         end do         
      end do

      do k = kmin, kmax
c
c        ::::: compute min and max of neighbors on k+1 slab
c
         do j = jmin, jmax     
            do i = imin, imax   
               smin(i,j,kp) = min(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &              s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &              s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
               smax(i,j,kp) = max(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &              s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &              s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
c
c        ::::: compute min/max of cell
c
               smn = min(smin(i,j,km),smin(i,j,kk),smin(i,j,kp))
               smx = max(smax(i,j,km),smax(i,j,kk),smax(i,j,kp))
               sn(i,j,k) = max(sn(i,j,k),smn)
               sn(i,j,k) = min(sn(i,j,k),smx)
               
            end do
         end do
c
c        ::::: roll indices for next slab
c
         km = mod(km+1,3)
         kk = mod(kk+1,3)
         kp = mod(kp+1,3)
      end do

      end

      subroutine FORT_SUM_TF_DIVU(
     &     S,DIMS(S),
     &     tforces,DIMS(tf),
     &     divu,DIMS(divu),
     &     lo,hi,nvar,iconserv )
c
c     sum tforces, viscous forcing and divU*S into tforces
c     depending on the value of iconserv
c
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, k, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)

      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T divu(DIMV(divu))

      if ( iconserv .eq. 1 ) then
!$omp parallel private(i,j,k,n)
         do n = 1, nvar
!$omp do
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = tforces(i,j,k,n) - S(i,j,k,n)*divu(i,j,k)
                  end do
               end do
            end do
!$omp end do nowait
         end do
!$omp end parallel
      end if

      end

      subroutine FORT_SUM_TF_DIVU_VISC(
     &     S,DIMS(S),
     &     tforces,DIMS(tf),
     &     divu,DIMS(divu),
     &     visc,DIMS(visc),
     &     lo,hi,nvar,iconserv )
c
c     sum tforces, viscous forcing and divU*S into tforces
c     depending on the value of iconserv
c
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, k, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(visc)

      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T divu(DIMV(divu))
      REAL_T visc(DIMV(visc),nvar)

!$omp parallel private(i,j,k,n)
      if ( iconserv .eq. 1 ) then
         do n = 1, nvar
!$omp do
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = tforces(i,j,k,n)
     &                    +  visc(i,j,k,n) -  S(i,j,k,n)*divu(i,j,k)
                  end do
               end do
            end do
!$omp end do nowait
         end do
      else
         do n = 1, nvar
!$omp do
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = tforces(i,j,k,n) +  visc(i,j,k,n) 
                  end do
               end do
            end do
!$omp end do nowait
         end do
      end if
!$omp end parallel

      end

      subroutine FORT_UPDATE_AOFS_TF(
     &     s,       DIMS(s),
     &     sn,      DIMS(sn),
     &     aofs,    DIMS(aofs),
     &     tforces, DIMS(tf),
     &     rockphi, DIMS(rockphi),
     &     lo,hi,dt,nvar)
c
c     update a field with an advective tendency
c     and a forcing term
c
      integer i, j, k, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(aofs)
      integer DIMDEC(tf)
      integer DIMDEC(rockphi)
      integer lo(SDIM), hi(SDIM)
      REAL_T dt
      REAL_T s(DIMV(s),nvar)
      REAL_T sn(DIMV(sn),nvar)
      REAL_T aofs(DIMV(aofs),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T rockphi(DIMV(rockphi))

!$omp parallel private(i,j,k,n)
      do n = 1,nvar
!$omp do
         do k = lo(3), hi(3)
            do j = lo(2), hi(2) 
               do i = lo(1), hi(1)  
                  sn(i,j,k,n) = s(i,j,k,n)
     &                 - dt*aofs(i,j,k,n)/rockphi(i,j,k)
     &                 + dt*tforces(i,j,k,n)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end

      subroutine FORT_UPDATE_AOFS_TRACER(
     &     s,       DIMS(s),
     &     sn,      DIMS(sn),
     &     aofs,    DIMS(aofs),
     &     tforces, DIMS(tf),
     &     rockphi, DIMS(rockphi),
     &     lo,hi, idx_total, idx_total_n, dt, nvar)
c
c     update tracer with an advective tendency
c     and a forcing term
c
      implicit none
      
#include "probdata.H"

      integer nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(aofs)
      integer DIMDEC(tf)
      integer DIMDEC(rockphi)
      integer lo(SDIM), hi(SDIM) 
      integer idx_total_n
      integer idx_total(idx_total_n)
      REAL_T dt
      REAL_T s(DIMV(s),nvar)
      REAL_T sn(DIMV(sn),nvar)
      REAL_T aofs(DIMV(aofs),nvar)
      REAL_T tforces(DIMV(tf),nvar-ncomps)
      REAL_T rockphi(DIMV(rockphi))

      integer i, j, k, n, idx
      REAL_T  sat, satn

!$omp parallel do private(i,j,k,n,idx,sat,satn)
      do k = lo(3), hi(3)
         do j = lo(2), hi(2) 
            do i = lo(1), hi(1) 
               do n = ncomps+1,nvar
                  sn(i,j,k,n) = s(i,j,k,n)
               end do

               sat = s(i,j,k,1)/rhoval(1)
               satn = sn(i,j,k,1)/rhoval(1)
               do n = 1,idx_total_n
                  idx = idx_total(n)
                  if (satn .gt. zero) then
                     sn(i,j,k,idx) = (sat*s(i,j,k,idx)
     &                    - dt*aofs(i,j,k,idx)/rockphi(i,j,k))/satn
     &                    + dt*tforces(i,j,k,idx-ncomps)
                  end if
               end do
            end do
         end do
      end do
!$omp end parallel do

      end

      subroutine slope_alpha(s,DIMS(s),slx,sly,slz,DIMS(sl),
     &     eigvx,eiglx,eigrx,strucx,ldefx,
     &     eigvy,eigly,eigry,strucy,ldefy,
     &     eigvz,eiglz,eigrz,strucz,ldefz,
     &     DIMS(eig),
     &     slxscr,slyscr,slzscr,lo,hi,bc,ncomps,nphases)
c ::: slope limiting for alpha.
c ::: ncomps is adjusted to be ncomps-1
#include "GODCOMM_F.H"

      integer DIMDEC(s)
      integer DIMDEC(sl)
      integer DIMDEC(eig)
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer ncomps, nphases
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  slx(DIMV(sl),ncomps)
      REAL_T  sly(DIMV(sl),ncomps)
      REAL_T  slz(DIMV(sl),ncomps)
      REAL_T  slxscr(DIM1(s),3,2)
      REAL_T  slyscr(DIM2(s),3,2)
      REAL_T  slzscr(DIM3(s),3,2)
      REAL_T  eigvx(DIMV(eig),ncomps)
      REAL_T  eigvy(DIMV(eig),ncomps)
      REAL_T  eigvz(DIMV(eig),ncomps)
      REAL_T  eiglx(DIMV(eig),ncomps,ncomps)
      REAL_T  eigly(DIMV(eig),ncomps,ncomps)
      REAL_T  eiglz(DIMV(eig),ncomps,ncomps)
      REAL_T  eigrx(DIMV(eig),ncomps,ncomps)
      REAL_T  eigry(DIMV(eig),ncomps,ncomps)
      REAL_T  eigrz(DIMV(eig),ncomps,ncomps)
      REAL_T  strucx(DIMV(eig),ncomps,ncomps)
      REAL_T  strucy(DIMV(eig),ncomps,ncomps)
      REAL_T  strucz(DIMV(eig),ncomps,ncomps)
      logical ldefx(DIMV(eig)),ldefy(DIMV(eig)),ldefz(DIMV(eig))
      integer ng,n,m,i,j,k
      integer imin,jmin,kmin,imax,jmax,kmax
      REAL_T alpl(2), alpc(2), alpr(2)
      REAL_T anslp, den, faclim, cdef, esnew, slrmin, def

      ng = lo(1) - ARG_L1(s)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      cdef = 1.0d-1

      if (ng .lt. 3) then
         call bl_pd_abort("SLOPE_3D: not enough bndry cells for 4th order")
      end if

c ::: added to prevent underflow for small s values
      do n = 1,ncomps
         do k = kmin-ng,kmax+ng
            do j = jmin-ng, jmax+ng
               do i = imin-ng, imax+ng
                  s(i,j,k,n) = cvmgt(s(i,j,k,n), zero, abs(s(i,j,k,n)).gt.1.0D-20)
               end do
            end do
         end do
      end do


      if (use_unlimited_slopes) then
         call bl_pd_abort("cannot use unlimited slopes")
      else

c :::::: x-slopes
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1 
                  do n = 1,ncomps
                     slxscr(i,1,n)  = s(i,j,k,n)   - s(i-1,j,k,n)
                     slxscr(i,3,n)  = s(i+1,j,k,n) - s(i,j,k,n)
                     slxscr(i,2,n)  = 2*(s(i+1,j,k,n)-s(i-1,j,k,n))/3
     &                 -(s(i+2,j,k,n)-s(i-2,j,k,n))/12
                  end do
               end do
            

               if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slxscr(imin,2,n) = -sixteen/fifteen*s(imin-1,j,k,n) + 
     &                    half*s(imin,j,k,n) + two3rd*s(imin+1,j,k,n) - tenth*s(imin+2,j,k,n)
                     slxscr(imin-1,2,n) = zero
                     slx(imin-1,j,k,n) = zero
                  end do
               end if
               if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slxscr(imax,2,n) = -( -sixteen/fifteen*s(imax+1,j,k,n) +  
     &                    half*s(imax,j,k,n) + two3rd*s(imax-1,j,k,n) - tenth*s(imax-2,j,k,n) )
                     slxscr(imax+1,2,n) = zero
                     slx(imax+1,j,k,n) = zero
                  end do
               end if

               do i = imin-1,imax+1
                  do n = 1,ncomps
                     alpl(n) = zero
                     alpc(n) = zero
                     alpr(n) = zero                  
                     do m = 1,ncomps
                        alpl(n) = alpl(n) + eiglx(i,j,k,n,m)*slxscr(i,1,m)
                        alpc(n) = alpc(n) + eiglx(i,j,k,n,m)*slxscr(i,2,m)
                        alpr(n) = alpr(n) + eiglx(i,j,k,n,m)*slxscr(i,3,m)
                     end do
                  end do

                  do n = 1,ncomps
                     faclim = two
                     if (i .ge. imin .and. i .le. imax) then
                        if(strucx(i,j,k,n,n)*strucx(i+1,j,k,n,n).lt.-1.e-10)faclim=1.
                        if(strucx(i,j,k,n,n)*strucx(i-1,j,k,n,n).lt.-1.e-10)faclim=1.
                     end if
                     anslp = dmin1(faclim*dabs(alpl(n)),
     &                    dabs(alpc(n)),faclim*dabs(alpr(n)))
     &                    *dmax1(zero,dsign(one,alpl(n)*alpr(n)))
                     den=cvmgt(dabs(alpc(n)),one,dabs(alpc(n)).ne.zero)
                     slx(i,j,k,n)=anslp/den
                  end do
               
                  def=dabs(alpc(1))*cdef*dabs(strucx(i,j,k,1,1)-strucx(i,j,k,2,1))
     &                 +dabs(alpc(2))*cdef*dabs(strucx(i,j,k,1,2)-strucx(i,j,k,2,2))
     &                 -dabs(eigvx(i,j,k,1)-eigvx(i,j,k,2))

                  if(def .gt. 1.0d-10 .or. ldefx(i,j,k))then
                     if(eigvx(i,j,k,1)*eigvx(i,j,k,2).lt.0)then
                        slx(i,j,k,1)=zero
                        slx(i,j,k,2)=zero
                     else
                        esnew=half*(eigvx(i,j,k,1)+eigvx(i,j,k,2))
                        slrmin=dmin1(slx(i,j,k,1),slx(i,j,k,2))
                        slx(i,j,k,1)=slrmin
                        slx(i,j,k,2)=slrmin
                        eigvx(i,j,k,1)=esnew
                        eigvx(i,j,k,2)=esnew
                     endif
                  endif
                  slx(i,j,k,1) = alpc(1)*slx(i,j,k,1)
                  slx(i,j,k,2) = alpc(2)*slx(i,j,k,2)
               end do

            end do
         end do

c :::::: y-slopes
         do k = kmin-1,kmax+1
            do i = imin-1,imax+1 
               do j = jmin-1,jmax+1      
                  do n = 1,ncomps
                     slyscr(j,1,n)  = s(i,j,k,n)   - s(i,j-1,k,n)
                     slyscr(j,3,n)  = s(i,j+1,k,n) - s(i,j,k,n)
                     slyscr(j,2,n)  = 2*(s(i,j+1,k,n)-s(i,j-1,k,n))/3
     &                    -(s(i,j+2,k,n)-s(i,j-2,k,n))/12
                  end do
               end do

               if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slyscr(jmin,2,n) = -sixteen/fifteen*s(i,jmin-1,k,n) + half*s(i,jmin,k,n) + 
     &                    two3rd*s(i,jmin+1,k,n) - tenth*s(i,jmin+2,k,n)
                     slyscr(jmin-1,2,n) = zero
                     sly(i,jmin-1,k,n) = zero
                  end do
               end if
               if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slyscr(jmax,2,n) = -( -sixteen/fifteen*s(i,jmax+1,k,n) + half*s(i,jmax,k,n) + 
     &                    two3rd*s(i,jmax-1,k,n) - tenth*s(i,jmax-2,k,n) )
                     slyscr(jmax+1,2,n) = zero
                     sly(i,jmax+1,k,n) = zero
                  end do
               end if
               do j = jmin-1,jmax+1 
                  do n = 1,ncomps
                     alpl(n) = zero
                     alpc(n) = zero
                     alpr(n) = zero                  
                     do m = 1,2
                        alpl(n) = alpl(n) + eigly(i,j,k,n,m)*slyscr(j,1,m)
                        alpc(n) = alpc(n) + eigly(i,j,k,n,m)*slyscr(j,2,m)
                        alpr(n) = alpr(n) + eigly(i,j,k,n,m)*slyscr(j,3,m)
                     end do
                  end do

                  do n = 1,ncomps
                     faclim = two
                     if (j .ge. jmin .and. j .le. jmax) then
                        if(strucy(i,j,k,n,n)*strucy(i,j+1,k,n,n).lt.-1.e-10)faclim=1.
                        if(strucy(i,j,k,n,n)*strucy(i,j-1,k,n,n).lt.-1.e-10)faclim=1.
                     end if
                     anslp = dmin1(faclim*dabs(alpl(n)),
     &                    dabs(alpc(n)),faclim*dabs(alpr(n)))
     &                    *dmax1(zero,dsign(one,alpl(n)*alpr(n)))
                     den=cvmgt(dabs(alpc(n)),one,dabs(alpc(n)).ne.zero)
                     sly(i,j,k,n)=anslp/den
                  end do
               
                  def= dabs(alpc(1))*cdef*dabs(strucy(i,j,k,1,1)-strucy(i,j,k,2,1))
     &                 +dabs(alpc(2))*cdef*dabs(strucy(i,j,k,1,2)-strucy(i,j,k,2,2))
     &                 -dabs(eigvy(i,j,k,1)-eigvy(i,j,k,2))

                  if(def.gt. 1.0d-10 .or. ldefy(i,j,k) )then
                     if(eigvy(i,j,k,1)*eigvy(i,j,k,2).lt.0)then
                        sly(i,j,k,1)=zero
                        sly(i,j,k,2)=zero
                     else
                        esnew=half*(eigvy(i,j,k,1)+eigvy(i,j,k,2))
                        slrmin=dmin1(sly(i,j,k,1),sly(i,j,k,2))
                        sly(i,j,k,1)=slrmin
                        sly(i,j,k,2)=slrmin
                        eigvy(i,j,k,1)=esnew
                        eigvy(i,j,k,2)=esnew
                     endif
                  endif
                  sly(i,j,k,1) = alpc(1)*sly(i,j,k,1)
                  sly(i,j,k,2) = alpc(2)*sly(i,j,k,2)
               end do

            end do
          end do

c :::::: z-slopes

          do i = imin-1,imax+1 
             do j = jmin-1,jmax+1 
                do k = kmin-1,kmax+1
                  do n = 1,ncomps
                     slzscr(k,1,n)  = s(i,j,k,n)   - s(i,j,k-1,n)
                     slzscr(k,3,n)  = s(i,j,k+1,n) - s(i,j,k,n)
                     slzscr(k,2,n)  = 2*(s(i,j,k+1,n)-s(i,j,k-1,n))/3
     &                    -(s(i,j,k+2,n)-s(i,j,k-2,n))/12
                  end do
               end do

               if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slzscr(kmin,2,n) = -sixteen/fifteen*s(i,j,kmin-1,n) + half*s(i,j,kmin,n) + 
     &                    two3rd*s(i,j,kmin+1,n) - tenth*s(i,j,kmin+2,n)
                     slzscr(kmin-1,2,n) = zero
                     slz(i,j,kmin-1,n) = zero
                  end do
               end if
               if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slzscr(kmax,2,n) = -( -sixteen/fifteen*s(i,j,kmax+1,n) + 
     &                    half*s(i,j,kmax,n) + two3rd*s(i,j,kmax-1,n) - tenth*s(i,j,kmax-2,n) )
                     slzscr(kmax+1,2,n) = zero
                     slz(i,j,kmax+1,n) = zero
                  end do
               end if
               do k = kmin-1,kmax+1 
                  do n = 1,ncomps
                     alpl(n) = zero
                     alpc(n) = zero
                     alpr(n) = zero                  
                     do m = 1,2
                        alpl(n) = alpl(n) + eiglz(i,j,k,n,m)*slzscr(k,1,m)
                        alpc(n) = alpc(n) + eiglz(i,j,k,n,m)*slzscr(k,2,m)
                        alpr(n) = alpr(n) + eiglz(i,j,k,n,m)*slzscr(k,3,m)
                     end do
                  end do

                  do n = 1,ncomps
                     faclim = two
                     if (k .ge. kmin .and. k .le. kmax) then
                        if(strucz(i,j,k,n,n)*strucz(i,j,k+1,n,n).lt.-1.e-10)faclim=1.
                        if(strucz(i,j,k,n,n)*strucz(i,j,k-1,n,n).lt.-1.e-10)faclim=1.
                     end if
                     anslp = dmin1(faclim*dabs(alpl(n)),
     &                    dabs(alpc(n)),faclim*dabs(alpr(n)))
     &                    *dmax1(zero,dsign(one,alpl(n)*alpr(n)))
                     den=cvmgt(dabs(alpc(n)),one,dabs(alpc(n)).ne.zero)
                     slz(i,j,k,n)=anslp/den
                  end do
               
                  def= dabs(alpc(1))*cdef*dabs(strucz(i,j,k,1,1)-strucz(i,j,k,2,1))
     &                 +dabs(alpc(2))*cdef*dabs(strucz(i,j,k,1,2)-strucz(i,j,k,2,2))
     &                 -dabs(eigvz(i,j,k,1)-eigvz(i,j,k,2))

                  if(def.gt. 1.0d-10 .or. ldefz(i,j,k) )then
                     if(eigvz(i,j,k,1)*eigvz(i,j,k,2).lt.zero)then
                        slz(i,j,k,1)=zero
                        slz(i,j,k,2)=zero
                     else
                        esnew=half*(eigvz(i,j,k,1)+eigvz(i,j,k,2))
                        slrmin=dmin1(slz(i,j,k,1),slz(i,j,k,2))
                        slz(i,j,k,1)=slrmin
                        slz(i,j,k,2)=slrmin
                        eigvz(i,j,k,1)=esnew
                        eigvz(i,j,k,2)=esnew
                     endif 
                  end if
                     
                  slz(i,j,k,1) = alpc(1)*slz(i,j,k,1)
                  slz(i,j,k,2) = alpc(2)*slz(i,j,k,2)
                     
               end do

            end do
          end do

       end if

       end

