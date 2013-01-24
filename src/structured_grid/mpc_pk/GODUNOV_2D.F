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

#define SDIM 2
#define XVEL 1
#define YVEL 2

      subroutine FORT_EST_EIG (s, DIMS(s), lam, dlam, d2lam, DIMS(lam),
     &     rphi,  DIMS(rphi),  krcoef, DIMS(krcoef), nkrcoef,
     &     uedge, DIMS(uedge), kappax, DIMS(kpx),
     &     vedge, DIMS(vedge), kappay, DIMS(kpy),
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
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(kpx),DIMDEC(kpy)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM)
      integer bc(SDIM,2)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  lam(DIMV(lam),2)
      REAL_T  dlam(DIMV(lam),3)
      REAL_T  d2lam(DIMV(lam),3)
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge))
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  dx(SDIM)
      REAL_T  eigmax(SDIM)

c ::: internal variables
      REAL_T  gravx, gravy
      
      gravx = zero
      gravy = -gravity

      call get_dlbd_cc(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     ncomps,rhoval,muval,model)
      call riemann_est_eig(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     uedge,DIMS(uedge),kappax,DIMS(kpx),rphi,DIMS(rphi),     
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,lo,hi,ncomps,gravx,rhoval,muval,dx,bc,eigmax(1))
      call riemann_est_eig(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     vedge,DIMS(vedge),kappay,DIMS(kpy),rphi,DIMS(rphi),     
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,lo,hi,ncomps,gravy,rhoval,muval,dx,bc,eigmax(2))
      end      

      subroutine FORT_EST_EIG_LIN (rphi,  DIMS(rphi), 
     &     uedge, DIMS(uedge),
     &     vedge, DIMS(vedge), 
     &     lo,hi,eigmax)
c
c ::: This subroutine estimates the eigenvalues for timestep prediction
c
      implicit none
      
c ::: arguments
      integer DIMDEC(rphi)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer lo(SDIM),hi(SDIM)
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  eigmax(SDIM)

c ::: internal variables
      integer i,j
      REAL_T  eigtmp
      
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)+1
            eigtmp = uedge(i,j)/(half*(rphi(i,j)+rphi(i-1,j)))
            eigmax(1) = max(eigmax(1),dabs(eigtmp))
         end do
      end do
      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)
            eigtmp = vedge(i,j)/(half*(rphi(i,j)+rphi(i,j-1)))
            eigmax(2) = max(eigmax(2),dabs(eigtmp))
         end do
      end do

      end 

      subroutine FORT_EST_EIG_TRC (s, DIMS(s), nc,
     &     rphi,  DIMS(rphi), 
     &     uedge, DIMS(uedge),
     &     vedge, DIMS(vedge), 
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
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer lo(SDIM),hi(SDIM)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  eigmax(SDIM)

c ::: internal variables
      integer i,j
      REAL_T  eigtmp,sat
      
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)+1
            sat = half*(s(i,j,nc)+s(i-1,j,nc))/rhoval(nc)
            eigtmp = uedge(i,j)/(half*(rphi(i,j)+rphi(i-1,j)))/sat
            eigmax(1) = max(eigmax(1),dabs(eigtmp))
         end do
      end do
      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)
            sat = half*(s(i,j,nc)+s(i,j-1,nc))/rhoval(nc)
            eigtmp = vedge(i,j)/(half*(rphi(i,j)+rphi(i,j-1)))/sat
            eigmax(2) = max(eigmax(2),dabs(eigtmp))
         end do
      end do

      end 

      subroutine FORT_EST_EIG_CPL (s, DIMS(s), lam, dlam, d2lam, DIMS(lam),
     &     rphi,  DIMS(rphi),  krcoef, DIMS(krcoef), nkrcoef,
     &     pc, DIMS(pc), uedge, DIMS(uedge), kappax, DIMS(kpx),
     &     vedge, DIMS(vedge), kappay, DIMS(kpy),
     &     lo,hi,dx,bc,eigmax)
c
c ::: This subroutine estimates the eigenvalues for timestep prediction
c
      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(pc)
      integer DIMDEC(lam)
      integer DIMDEC(rphi)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(kpx),DIMDEC(kpy)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM)
      integer bc(SDIM,2)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  lam(DIMV(lam),2)
      REAL_T  dlam(DIMV(lam),3)
      REAL_T  d2lam(DIMV(lam),3)
      REAL_T  pc(DIMV(pc))
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge))
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  dx(SDIM)
      REAL_T  eigmax(SDIM)

c ::: internal variables
      REAL_T  gravx, gravy
      
      gravx = zero
      gravy = -gravity

      call get_dlbd_cc(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     ncomps,rhoval,muval,model)
      call riemann_est_eig_cpl(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     pc,DIMS(pc),uedge,DIMS(uedge),kappax,DIMS(kpx),rphi,DIMS(rphi),     
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,lo,hi,ncomps,gravx,rhoval,muval,dx,bc,eigmax(1))
      call riemann_est_eig_cpl(s,DIMS(s),lam,dlam,d2lam,DIMS(lam),
     &     pc,DIMS(pc),vedge,DIMS(vedge),kappay,DIMS(kpy),rphi,DIMS(rphi),     
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,lo,hi,ncomps,gravy,rhoval,muval,dx,bc,eigmax(2))
      end

      subroutine FORT_ESTATE_FPU(s, tforces, divu, DIMS(s),
     &     xlo, xhi, sx, slxscr, stxlo, stxhi,
     &     uedge, DIMS(uedge), xstate, DIMS(xstate),

     &     ylo, yhi, sy, slyscr, stylo, styhi,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate),

     &     DIMS(work),
     &     bc,lo,hi,dt,dx,n,use_minion,iconserv)
      implicit none
c
c     This subroutine computes edges states, right now it uses
c     a lot of memory, but there becomes a trade off between
c     simplicity-efficiency in the new way of computing states
c     and complexity in the old way.  By eliminating loops over
c     state components though, the new way uses much less memory.
c
      integer i,j,n
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer imin,jmin,imax,jmax
      REAL_T hx, hy, dt, dth, dthx, dthy
      REAL_T tr,ubar,vbar,stx,sty,fu,fv,dx(SDIM)
      REAL_T eps,eps_for_bc
      parameter( eps        = zero )
      parameter( eps_for_bc = zero )

      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(xstate)
      integer DIMDEC(vedge)
      integer DIMDEC(ystate)

      REAL_T s(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)

      REAL_T uedge(DIMV(uedge)), xstate(DIMV(xstate))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(ystate))

      REAL_T xlo(DIMV(work)), xhi(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work))
      REAL_T  sx(DIMV(work))
      REAL_T  sy(DIMV(work))
      REAL_T tforces(DIMV(work))
      REAL_T divu(DIMV(work))

     
      REAL_T spx,smx,spy,smy,st

      integer use_minion, iconserv
      integer inc

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
c     compute the slopes
c
      call FORT_SLOPES(
     &     s,DIMS(s),
     &     sx,sy,DIMS(work),
     &     lo,hi,slxscr,slyscr,bc)
c
c     trace the state to the cell edges
c
      do j = jmin-1,jmax+1
         do i = imin,  imax+1
            xlo(i,j) = s(i-1,j) + (half  - dthx*uedge(i,j))*sx(i-1,j)
            xhi(i,j) = s(i,  j) + (-half - dthx*uedge(i,j))*sx(i,  j)
         end do
      end do

      if(use_minion.eq.1)then
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = xlo(i,j) + dth*tforces(i-1,j)
               xhi(i,j) = xhi(i,j) + dth*tforces(i,  j)
            end do
         end do
         if (iconserv .eq. 1) then
           do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = xlo(i,j) - dth*s(i-1,j)*divu(i-1,j)
               xhi(i,j) = xhi(i,j) - dth*s(i  ,j)*divu(i,  j)
            end do
           end do
         end if
      end if

      call trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(work),uedge,DIMS(uedge),
     &     lo,hi,n,bc,eps_for_bc)

      do j = jmin-1,jmax+1
         do i = imin,  imax+1
            fu  = merge(zero,one,abs(uedge(i,j)).lt.eps)
            stx = merge(xlo(i,j),xhi(i,j),uedge(i,j)>=zero)
            xlo(i,j) = fu*stx + (one - fu)*half*(xhi(i,j)+xlo(i,j))
         end do
      end do

      do j = jmin,  jmax+1
         do i = imin-1,imax+1
            ylo(i,j) = s(i,j-1) + (half  - dthy*vedge(i,j))*sy(i,j-1)
            yhi(i,j) = s(i,j)   + (-half - dthy*vedge(i,j))*sy(i,j)
         end do
      end do

      if (use_minion.eq.1)then
         do j = jmin, jmax+1
            do i = imin-1,  imax+1
               ylo(i,j) = ylo(i,j) + dth*tforces(i,j-1)
               yhi(i,j) = yhi(i,j) + dth*tforces(i,j  )
            end do
         end do
         if (iconserv .eq. 1) then
           do j = jmin-1,jmax+1
            do i = imin,  imax+1
               ylo(i,j) = ylo(i,j) - dth*s(i,j-1)*divu(i,j-1)
               yhi(i,j) = yhi(i,j) - dth*s(i,j  )*divu(i,j  )
            end do
           end do
         end if
      end if

      call trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(work),vedge,DIMS(vedge),
     &     lo,hi,n,bc,eps_for_bc)

      do j = jmin,  jmax+1
         do i = imin-1,imax+1
            fv  = merge(zero,one,abs(vedge(i,j)).lt.eps)
            sty = merge(ylo(i,j),yhi(i,j),vedge(i,j)>=zero)
            ylo(i,j) = fv*sty + (one - fv)*half*(yhi(i,j)+ylo(i,j))
         end do
      end do

c
c     compute the xedge states
c
      do j = jmin,jmax

            do i = imin-1,imax+1
               spx = s(i,j) + half * sx(i,j)
               smx = s(i,j) - half * sx(i,j)

               if (iconserv.eq.1) then

                  st = -(vedge(i,j+1)*ylo(i,j+1) - vedge(i,j)*ylo(i,j))/hy
     &                 + s(i,j)*(vedge(i,j+1)-vedge(i,j))/hy
     &                 - s(i,j)*divu(i,j)

               else

#if 1
                  if (vedge(i,j)*vedge(i,j+1).le.0.d0) then
                     vbar = 0.5d0*(vedge(i,j)+vedge(i,j+1))
                     if (vbar.lt.0.d0) then
                        inc = 1
                     else
                        inc = 0
                     endif
                     tr = vbar*(s(i,j+inc)-s(i,j+inc-1))/hy
                  else
                     tr = half*(vedge(i,j+1) + vedge(i,j)) *
     &                    (  ylo(i,j+1) - ylo(i,j)  ) / hy
                  endif

#else
                  tr = half*(vedge(i,j+1) + vedge(i,j)) * extra_work(i,j,2)
#endif
                  st =  -tr

               endif

               stxlo(i+1)= spx - dthx*uedge(i+1,j)*sx(i,j) + dth*(st + tforces(i,j))
               stxhi(i  )= smx - dthx*uedge(i  ,j)*sx(i,j) + dth*(st + tforces(i,j))

            end do

            if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j).ge.zero) then
               stxhi(imin) = s(imin-1,j)
               stxlo(imin) = s(imin-1,j)
            else if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j).lt.zero) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP
     &              .or.bc(1,1).eq.REFLECT_EVEN) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1).eq.REFLECT_ODD) then
               stxhi(imin) = zero
               stxlo(imin) = zero
            end if
            if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j).le.zero) then
               stxlo(imax+1) = s(imax+1,j)
               stxhi(imax+1) = s(imax+1,j)
            else if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j).gt.zero) then
               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP
     &              .or.bc(1,2).eq.REFLECT_EVEN) then
               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2).eq.REFLECT_ODD) then
               stxlo(imax+1) = zero
               stxhi(imax+1) = zero
            end if
            
            do i = imin, imax+1
               xstate(i,j) = merge(stxlo(i),stxhi(i),uedge(i,j)>=zero)
               xstate(i,j) = merge(half*(stxlo(i)+stxhi(i)),xstate(i,j)
     &              ,abs(uedge(i,j)).lt.eps)
            end do
      end do
c
c     compute the yedge states
c
      do i = imin,imax
            
            do j = jmin-1,jmax+1
               spy = s(i,j) + half * sy(i,j)
               smy = s(i,j) - half * sy(i,j)

               if (iconserv.eq.1) then

                  st = -(uedge(i+1,j)*xlo(i+1,j) - uedge(i,j)*xlo(i,j))/hx
     &                 + s(i,j)*(uedge(i+1,j)-uedge(i,j))/hx
     &                 - s(i,j)*divu(i,j)

               else

#if 1
                  if (uedge(i,j)*uedge(i+1,j).le.0.d0) then
                     ubar = 0.5d0*(uedge(i,j)+uedge(i+1,j))
                     if (ubar.lt.0.d0) then
                        inc = 1
                     else
                        inc = 0
                     endif
                     tr = ubar*(s(i+inc,j)-s(i+inc-1,j))/hx
                  else

                     tr = half*(uedge(i+1,j) + uedge(i,j)) *
     &                           (xlo(i+1,j) -   xlo(i,j)  ) / hx
                  endif
#else
                  tr = half*(uedge(i+1,j) + uedge(i,j)) * extra_work(i,j,1)
#endif
                  st = -tr

               endif

               stylo(j+1)= spy - dthy*vedge(i,j+1)*sy(i,j) + dth*(st + tforces(i,j))
               styhi(j  )= smy - dthy*vedge(i,j  )*sy(i,j) + dth*(st + tforces(i,j))
            end do
            
            if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin).ge.zero) then
               styhi(jmin) = s(i,jmin-1)
               stylo(jmin) = s(i,jmin-1)
            else if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin).lt.zero) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP
     &              .or.bc(2,1).eq.REFLECT_EVEN) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1).eq.REFLECT_ODD) then
               styhi(jmin) = zero
               stylo(jmin) = zero
            end if
            
            if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1).le.zero) then
               stylo(jmax+1) = s(i,jmax+1)
               styhi(jmax+1) = s(i,jmax+1)
            else if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1).gt.zero) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP
     &              .or.bc(2,2).eq.REFLECT_EVEN) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2).eq.REFLECT_ODD) then
               stylo(jmax+1) = zero
               styhi(jmax+1) = zero
            end if
            do j=jmin,jmax+1
               ystate(i,j) = merge(stylo(j),styhi(j),vedge(i,j)>=zero)
               ystate(i,j) = merge(half*(stylo(j)+styhi(j)),ystate(i,j),
     &              abs(vedge(i,j)).lt.eps)
            end do
      end do
      
      end

      subroutine FORT_ESTATE_LIN (s, DIMS(s), sn, DIMS(sn), 
     &     tf, rphi, DIMS(rphi),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), 
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
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(xstate),DIMDEC(ystate)
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer nscal
      REAL_T   s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6)
      REAL_T  uedge(DIMV(uedge)), xstate(DIMV(xstate),nscal)
      REAL_T  vedge(DIMV(vedge)), ystate(DIMV(ystate),nscal)  
      REAL_T  sx(DIMV(work),nscal),  sy(DIMV(work),nscal)
      REAL_T  tf(DIMV(work),nscal)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  dx(SDIM), dt

c ::: internal variables
      REAL_T, ALLOCATABLE :: fluxx(:,:,:), fluxy(:,:,:)
      REAL_T, ALLOCATABLE :: xlo(:,:,:), xhi(:,:,:)
      REAL_T, ALLOCATABLE :: ylo(:,:,:), yhi(:,:,:)
      REAL_T  gravx, gravy
      integer nc
      
      gravx = zero
      gravy = -gravity

      ALLOCATE(fluxx(DIMV(uedge),nscal))
      ALLOCATE(fluxy(DIMV(vedge),nscal))
      ALLOCATE(xlo(DIMV(uedge),nscal))
      ALLOCATE(ylo(DIMV(vedge),nscal))
      ALLOCATE(xhi(DIMV(uedge),nscal))
      ALLOCATE(yhi(DIMV(vedge),nscal)) 


c ::: compute the slopes
      do nc = 1,nscal
         call FORT_SLOPES(s(ARG_L1(s),ARG_L2(s),nc),DIMS(s),
     &        sx(ARG_L1(work),ARG_L2(work),nc),
     &        sy(ARG_L1(work),ARG_L2(work),nc),DIMS(work),
     &        lo,hi,slxscr,slyscr,bc)
      end do 

c ::: initial trace
      call trace_lin(s,DIMS(s),sn,DIMS(sn),
     &     rphi,DIMS(rphi),
     &     sx,sy,tf,DIMS(work),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     dt,dx,lo,hi,bc,nscal)
      
c ::: predictor step: determine final edge states 
      call riemann_solver_lin(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,nscal)
      call riemann_solver_lin(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,nscal)

      call trace_p_lin(s,DIMS(s),sn,DIMS(sn),
     &     fluxx,uedge,DIMS(uedge),fluxy,vedge,DIMS(vedge),
     &     rphi,DIMS(rphi),xlo,xhi,ylo,yhi,
     &     dt,dx,lo,hi,bc,nscal)

c ::: compute fluxes based on final edge states
      call riemann_solver_lin(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,nscal)
      call riemann_solver_lin(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,nscal)

      call flux_copy(xstate,DIMS(xstate),fluxx,DIMS(uedge),
     &     lo,hi,1,0,nscal)
      call flux_copy(ystate,DIMS(ystate),fluxy,DIMS(vedge),
     &     lo,hi,0,1,nscal)

      DEALLOCATE(fluxx,fluxy)
      DEALLOCATE(xlo,xhi,ylo,yhi)

      end

      subroutine FORT_ESTATE_RMN (s, DIMS(s), sn, DIMS(sn), 
     &     tf, divu, lbd_cc, dlbd_cc,
     &     rphi, DIMS(rphi), kappa, DIMS(kappa),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), kappax, DIMS(kpx),
     &     eigvx, eiglx, eigrx, 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), kappay, DIMS(kpy),
     &     eigvy, eigly, eigry, 
     &     DIMS(work),
     &     krcoef, DIMS(krcoef), nkrcoef,
     &     bc,lo,hi,dt,dx,
     &     use_minion,iconserv,nscal)
      
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
      integer DIMDEC(kappa)
      integer DIMDEC(work)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(xstate),DIMDEC(ystate)
      integer DIMDEC(kpx),DIMDEC(kpy)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer use_minion, iconserv
      integer nscal
      REAL_T   s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6)
      REAL_T  uedge(DIMV(uedge)), xstate(DIMV(xstate),nscal)
      REAL_T  vedge(DIMV(vedge)), ystate(DIMV(ystate),nscal)    
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy))
      REAL_T  sx(DIMV(work),nscal),  sy(DIMV(work),nscal)
      REAL_T  tf(DIMV(work),nscal)
      REAL_T  divu(DIMV(work))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  lbd_cc(DIMV(work),ncomps)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  dlbd_cc(DIMV(work),3)
      REAL_T  dx(SDIM), dt

c ::: internal variables
      REAL_T  eigvx(DIMV(work),2*nscal),eigvy(DIMV(work),2*nscal)
      REAL_T  eiglx(DIMV(work),nscal,nscal),eigly(DIMV(work),nscal,nscal)
      REAL_T  eigrx(DIMV(work),nscal,nscal),eigry(DIMV(work),nscal,nscal)
      REAL_T, ALLOCATABLE :: fluxx(:,:,:), fluxy(:,:,:)
      REAL_T, ALLOCATABLE :: xlo(:,:,:), xhi(:,:,:)
      REAL_T, ALLOCATABLE :: ylo(:,:,:), yhi(:,:,:)
      REAL_T  gravx, gravy
      integer nc, corrector
      
      gravx = zero
      gravy = -gravity

      ALLOCATE(fluxx(DIMV(uedge),nscal))
      ALLOCATE(fluxy(DIMV(vedge),nscal))
      ALLOCATE(xlo(DIMV(uedge),nscal))
      ALLOCATE(ylo(DIMV(vedge),nscal))
      ALLOCATE(xhi(DIMV(uedge),nscal))
      ALLOCATE(yhi(DIMV(vedge),nscal)) 

c ::: compute eigen solutions
      call eigsol_mp(eigvx,eiglx,eigrx,lbd_cc,dlbd_cc,DIMS(work),
     &     uedge,DIMS(uedge),rphi,DIMS(rphi),kappax,DIMS(kpx),
     &     krcoef,DIMS(krcoef),nkrcoef,gravx,muval,
     &     1,0,lo,hi,rhoval,ncomps,nphases,model,bc)
      call eigsol_mp(eigvy,eigly,eigry,lbd_cc,dlbd_cc,DIMS(work),
     &     vedge,DIMS(vedge),rphi,DIMS(rphi),kappay,DIMS(kpy),
     &     krcoef,DIMS(krcoef),nkrcoef,gravy,muval,
     &     0,1,lo,hi,rhoval,ncomps,nphases,model,bc)
c ::: compute the slopes
      do nc = 1,nscal
         call FORT_SLOPES(s(ARG_L1(s),ARG_L2(s),nc),DIMS(s),
     &        sx(ARG_L1(work),ARG_L2(work),nc),
     &        sy(ARG_L1(work),ARG_L2(work),nc),DIMS(work),
     &        lo,hi,slxscr,slyscr,bc)
      end do 

c ::: initial step
      call trace(s,DIMS(s),sn,DIMS(sn),
     &     eigvx,eiglx,eigrx,eigvy,eigly,eigry,
     &     sx,sy,divu,tf,DIMS(work),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,nscal)
      
c ::: predictor step: determine final edge states 
      corrector = 0
      call riemann_solver(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),eiglx,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,lo,hi,gravx,muval,rhoval,nscal,model,corrector)
      call riemann_solver(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),eigly,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,lo,hi,gravy,muval,rhoval,nscal,model,corrector)

      call trace_p(s,DIMS(s),sn,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),
     &     fluxx,DIMS(uedge),fluxy,DIMS(vedge),
     &     xlo,xhi,ylo,yhi,divu,tf,DIMS(work),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     dt,dx,lo,hi,bc,muval,rhoval,iconserv,
     &     nscal,model,gravy)

c ::: compute fluxes based on final edge states
      corrector = 1
      call riemann_solver(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),eiglx,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,lo,hi,gravx,muval,rhoval,nscal,model,corrector)
      call riemann_solver(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),eigly,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,lo,hi,gravy,muval,rhoval,nscal,model,corrector) 

      call flux_copy(xstate,DIMS(xstate),fluxx,DIMS(uedge),
     &     lo,hi,1,0,nscal)
      call flux_copy(ystate,DIMS(ystate),fluxy,DIMS(vedge),
     &     lo,hi,0,1,nscal)

      DEALLOCATE(fluxx,fluxy)
      DEALLOCATE(xlo,xhi,ylo,yhi)

      end

      subroutine FORT_ESTATE_CPL (s, DIMS(s), sn, DIMS(sn), 
     &     pc, DIMS(pc),
     &     tf, divu, lbd_cc, dlbd_cc,
     &     rphi, DIMS(rphi), kappa, DIMS(kappa),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), kappax, DIMS(kpx),
     &     lambdax, DIMS(lambdax), eigvx, eiglx, eigrx, 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), kappay, DIMS(kpy),
     &     lambday, DIMS(lambday), eigvy, eigly, eigry, 
     &     DIMS(work),
     &     krcoef, DIMS(krcoef), nkrcoef,
     &     bc,lo,hi,dt,dx,
     &     use_minion,iconserv,nscal)
      
c
c ::: This subroutine computes edges states and fluxes. 
c ::: ncomps must be equal to nscal.
c
      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(pc)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(work)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(xstate),DIMDEC(ystate)
      integer DIMDEC(kpx),DIMDEC(kpy)
      integer DIMDEC(lambdax),DIMDEC(lambday)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer use_minion, iconserv
      integer nscal
      REAL_T   s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6)
      REAL_T  uedge(DIMV(uedge)), xstate(DIMV(xstate),nscal)
      REAL_T  vedge(DIMV(vedge)), ystate(DIMV(ystate),nscal)    
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy))
      REAL_T  lambdax(DIMV(lambdax)) , lambday(DIMV(lambday))
      REAL_T  sx(DIMV(work),nscal),  sy(DIMV(work),nscal)
      REAL_T  tf(DIMV(work),nscal)
      REAL_T  divu(DIMV(work))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  pc(DIMV(pc))
      REAL_T  lbd_cc(DIMV(work),ncomps)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  dlbd_cc(DIMV(work),3)
      REAL_T  dx(SDIM), dt

c ::: internal variables
      REAL_T eigvx(DIMV(work),2*nscal),eigvy(DIMV(work),2*nscal)
      REAL_T eiglx(DIMV(work),nscal,nscal),eigly(DIMV(work),nscal,nscal)
      REAL_T eigrx(DIMV(work),nscal,nscal),eigry(DIMV(work),nscal,nscal)
      REAL_T, ALLOCATABLE :: fluxx(:,:,:), fluxy(:,:,:)
      REAL_T, ALLOCATABLE :: xlo(:,:,:), xhi(:,:,:)
      REAL_T, ALLOCATABLE :: ylo(:,:,:), yhi(:,:,:)
      REAL_T  gravx, gravy
      integer nc, corrector
      
      gravx = zero
      gravy = -gravity

      ALLOCATE(fluxx(DIMV(uedge),nscal))
      ALLOCATE(fluxy(DIMV(vedge),nscal))
      ALLOCATE(xlo(DIMV(uedge),nscal))
      ALLOCATE(ylo(DIMV(vedge),nscal))
      ALLOCATE(xhi(DIMV(uedge),nscal))
      ALLOCATE(yhi(DIMV(vedge),nscal)) 

c ::: compute eigen solutions
      call eigsol_mp_cpl(eigvx,eiglx,eigrx,lbd_cc,dlbd_cc,DIMS(work),
     &     pc,DIMS(pc),uedge,DIMS(uedge),rphi,DIMS(rphi),kappax,DIMS(kpx),
     &     krcoef,DIMS(krcoef),nkrcoef,gravx,muval,
     &     1,0,lo,hi,rhoval,ncomps,nphases,model,dx,bc)
      call eigsol_mp_cpl(eigvy,eigly,eigry,lbd_cc,dlbd_cc,DIMS(work),
     &     pc,DIMS(pc),vedge,DIMS(vedge),rphi,DIMS(rphi),kappay,DIMS(kpy),
     &     krcoef,DIMS(krcoef),nkrcoef,gravy,muval,
     &     0,1,lo,hi,rhoval,ncomps,nphases,model,dx,bc)
c ::: compute the slopes
      do nc = 1,nscal
         call FORT_SLOPES(s(ARG_L1(s),ARG_L2(s),nc),DIMS(s),
     &        sx(ARG_L1(work),ARG_L2(work),nc),
     &        sy(ARG_L1(work),ARG_L2(work),nc),DIMS(work),
     &        lo,hi,slxscr,slyscr,bc)
      end do 

c ::: initial step
      call trace(s,DIMS(s),sn,DIMS(sn),
     &     eigvx,eiglx,eigrx,eigvy,eigly,eigry,
     &     sx,sy,divu,tf,DIMS(work),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,nscal)
      
c ::: predictor step: determine final edge states 
      corrector = 0
      call riemann_solver_cpl(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),lambdax,DIMS(lambdax),pc,DIMS(pc),eiglx,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,lo,hi,dx,gravx,muval,rhoval,nscal,model,bc,corrector)
      call riemann_solver_cpl(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),lambday,DIMS(lambday),pc,DIMS(pc),eigly,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,lo,hi,dx,gravy,muval,rhoval,nscal,model,bc,corrector)

      call trace_p(s,DIMS(s),sn,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),
     &     fluxx,DIMS(uedge),fluxy,DIMS(vedge),
     &     xlo,xhi,ylo,yhi,divu,tf,DIMS(work),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     dt,dx,lo,hi,bc,muval,rhoval,iconserv,
     &     nscal,model,gravy)

c ::: compute fluxes based on final edge states
      corrector = 1
      call riemann_solver_cpl(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),lambdax,DIMS(lambdax),pc,DIMS(pc),eiglx,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,lo,hi,dx,gravx,muval,rhoval,nscal,model,bc,corrector)
      call riemann_solver_cpl(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),lambday,DIMS(lambday),pc,DIMS(pc),eigly,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,lo,hi,dx,gravy,muval,rhoval,nscal,model,bc,corrector) 

      call flux_copy(xstate,DIMS(xstate),fluxx,DIMS(uedge),
     &     lo,hi,1,0,nscal)
      call flux_copy(ystate,DIMS(ystate),fluxy,DIMS(vedge),
     &     lo,hi,0,1,nscal)

      DEALLOCATE(fluxx,fluxy)
      DEALLOCATE(xlo,xhi,ylo,yhi)

      end

      subroutine FORT_ESTATE_PMR (s, u, DIMS(s),  
     &     sn, un, DIMS(sn),  tf, divu,
     &     rphi, DIMS(rphi), kappa, DIMS(kappa),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), kappax, DIMS(kpx),
     &     strucx, eigvx, eiglx, eigrx, 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), kappay, DIMS(kpy),
     &     strucy, eigvy, eigly, eigry, 
     &     DIMS(work),
     &     bc,lo,hi,dt,dx,
     &     use_minion,iconserv,grav,eigmax)
      implicit none
c
c     This subroutine computes edges states and fluxes for the polymer problem
c     The algorithm goes as follows:
c     1. transform problem into the auxiliary variables
c     2. do everything in the auxiliary variables
c     3. transform back the flux only at the very last step.
c     Here index of rho and mu are based on follows:
c       1: phase with 2 components
c       2: phase with 1 component
c     The above (confusing) indexing is due to the use of auxiliary 
c     variables that are not as flexible as components.
c
#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(work)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(xstate),DIMDEC(ystate)
      integer DIMDEC(kpx),DIMDEC(kpy)
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer use_minion, iconserv
      REAL_T  s (DIMV( s),ncomps)
      REAL_T  sn(DIMV(sn),ncomps)
      REAL_T  u(DIMV(s),ncomps-1)
      REAL_T  un(DIMV(sn),ncomps-1)
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6)
      REAL_T  uedge(DIMV(uedge)), xstate(DIMV(xstate),ncomps)
      REAL_T  vedge(DIMV(vedge)), ystate(DIMV(ystate),ncomps)    
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy))
      REAL_T  sx(DIMV(work),ncomps-1)  , sy(DIMV(work),ncomps-1)
      REAL_T  tf(DIMV(work))
      REAL_T  divu(DIMV(work))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  dx(SDIM), eigmax(SDIM),  grav
      REAL_T  eigvx(DIMV(work),ncomps),eigvy(DIMV(work),ncomps)
      REAL_T  eiglx(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  eigly(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  eigrx(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  eigry(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  strucx(DIMV(work),ncomps-1,ncomps-1)
      REAL_T  strucy(DIMV(work),ncomps-1,ncomps-1)

c ::: internal variables
      integer i,j
      integer compcount(nphases)
      integer compidx(nphases,ncomps)
      logical, ALLOCATABLE ::  ldefx(:,:), ldefy(:,:)
      REAL_T, ALLOCATABLE :: fluxx(:,:,:), fluxy(:,:,:)
      REAL_T, ALLOCATABLE :: xlo(:,:,:), xhi(:,:,:)
      REAL_T, ALLOCATABLE :: ylo(:,:,:), yhi(:,:,:)
      REAL_T  hx, hy, dt, dth, dthx, dthy, gravx, gravy
      REAL_T  athick
      REAL_T  mu(2), rho(2)
      
      ALLOCATE(ldefx(DIMV(work)))
      ALLOCATE(ldefy(DIMV(work)))
      ALLOCATE(fluxx(DIMV(uedge),ncomps-1))
      ALLOCATE(fluxy(DIMV(vedge),ncomps-1))
      ALLOCATE(xlo(DIMV(uedge),ncomps-1))
      ALLOCATE(ylo(DIMV(vedge),ncomps-1))
      ALLOCATE(xhi(DIMV(uedge),ncomps-1))
      ALLOCATE(yhi(DIMV(vedge),ncomps-1)) 

      dth   = half*dt
      dthx  = half*dt / dx(1)
      dthy  = half*dt / dx(2)
      hx    = dx(1)
      hy    = dx(2)
      gravx = zero
      gravy = -grav
c ::: identify comp-phase relation (in GODUNOV_F.F)
      call get_mu_rho_phases(rho, mu, athick, rhoval, muval, 
     &     compcount, compidx, phaseidx, ncomps, nphases)

c ::: s --> u
      call get_aux  (s,  u,  DIMS(s),  compcount, compidx, 
     &               rhoval, ncomps, nphases)
      call get_aux_n(sn, un, DIMS(sn), compcount, compidx, 
     &               rhoval, bc, ncomps, nphases)

c ::: compute the equivalent eigensystem
c ::: Jacobian entries computed on the fly
      call eigsol2(u,DIMS(s),uedge,DIMS(uedge),rphi,DIMS(rphi),
     &     kappa,DIMS(kappa),eigvx,eiglx,eigrx,strucx,ldefx,DIMS(work),
     &     lo,hi,1,0,mu,athick,rho,gravx,ncomps-1,nphases,bc,eigmax(1))
      call eigsol2(u,DIMS(s),vedge,DIMS(vedge),rphi,DIMS(rphi),
     &     kappa,DIMS(kappa),eigvy,eigly,eigry,strucy,ldefy,DIMS(work),
     &     lo,hi,0,1,mu,athick,rho,gravy,ncomps-1,nphases,bc,eigmax(2)) 

c ::: compute the slopes
      call slope_alpha(u,DIMS(s),sx,sy,DIMS(work),
     &     eigvx,eiglx,eigrx,strucx,ldefx,
     &     eigvy,eigly,eigry,strucy,ldefy,DIMS(work),
     &     slxscr,slyscr,lo,hi,bc,ncomps-1,nphases) 
 
c ::: initial step
      call trace2(u,DIMS(s),un,DIMS(sn),
     &     eigvx,eiglx,eigrx,eigvy,eigly,eigry,DIMS(work),
     &     sx,sy,divu,tf,DIMS(work),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,ncomps-1,nphases)

c ::: predictor step: determine final edge states 
      call riemann_solver2(fluxx,DIMS(uedge),u,DIMS(s),xlo,xhi,DIMS(uedge),
     &     eigvx,DIMS(work),uedge,DIMS(uedge),kappax,DIMS(kpx),
     &     1,0,lo,hi,athick,gravx,mu,rho,ncomps-1,nphases,dx,dt,bc)
      call riemann_solver2(fluxy,DIMS(vedge),u,DIMS(s),ylo,yhi,DIMS(vedge),
     &     eigvy,DIMS(work),vedge,DIMS(vedge),kappay,DIMS(kpy),
     &     0,1,lo,hi,athick,gravy,mu,rho,ncomps-1,nphases,dx,dt,bc)

      call trace_p2(u,DIMS(s),un,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),
     &     fluxx,DIMS(uedge),fluxy,DIMS(vedge),
     &     xlo,xhi,ylo,yhi,divu,tf,DIMS(work),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     dt,dx,lo,hi,bc,mu,rho,athick,iconserv,ncomps-1,nphases,gravy)

c ::: compute fluxes based on final edge states
      call riemann_solver2(fluxx,DIMS(uedge),u,DIMS(s),xlo,xhi,DIMS(uedge),
     &     eigvx,DIMS(work),uedge,DIMS(uedge),kappax,DIMS(kpx),
     &     1,0,lo,hi,athick,gravx,mu,rho,ncomps-1,nphases,dx,dt,bc)
      call riemann_solver2(fluxy,DIMS(vedge),u,DIMS(s),ylo,yhi,DIMS(vedge),
     &     eigvy,DIMS(work),vedge,DIMS(vedge),kappay,DIMS(kpy),
     &     0,1,lo,hi,athick,gravy,mu,rho,ncomps-1,nphases,dx,dt,bc) 

c ::: f(u) -> f(s)
      call flux_convert(xstate,DIMS(xstate),fluxx,DIMS(uedge),
     &     uedge,DIMS(uedge),lo,hi,1,0,rho,
     &     compcount,compidx,ncomps,nphases,bc)
      call flux_convert(ystate,DIMS(ystate),fluxy,DIMS(vedge),
     &     vedge,DIMS(vedge),lo,hi,0,1,rho,
     &     compcount,compidx,ncomps,nphases,bc) 
        
c ::: convert the second auxiliary variable to c instead of c*s
      do i = ARG_L1(s), ARG_H1(s)
         do j =  ARG_L2(s), ARG_H2(s)
            if (u(i,j,1) .gt. 1.d-10) then
               u(i,j,2) = u(i,j,2)/u(i,j,1)
            end if
         end do
      end do

      DEALLOCATE(ldefx,ldefy)
      DEALLOCATE(fluxx,fluxy)
      DEALLOCATE(xlo,xhi,ylo,yhi)

      end

      subroutine FORT_ESTATE_TRACER (s, DIMS(s),
     &     sn, DIMS(sn), 
     &     st, DIMS(st), stn, DIMS(stn),
     &     sx, slxscr, xlo, xhi,
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), 
     &     sy, slyscr, ylo, yhi,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), 
     &     rphi,  DIMS(rphi),DIMS(work),
     &     bc,lo,hi,dt,dx,nscal)
      implicit none
c
c ::: This subroutine computes tracers' fluxes at edges
c

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(st)
      integer DIMDEC(stn)
      integer DIMDEC(rphi)
      integer DIMDEC(work)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(xstate),DIMDEC(ystate)
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer nscal
      REAL_T   s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  st(DIMV(st))
      REAL_T  stn(DIMV(stn))
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6)
      REAL_T  uedge(DIMV(uedge))
      REAL_T  xstate(DIMV(xstate),nscal)
      REAL_T  vedge(DIMV(vedge))
      REAL_T  ystate(DIMV(ystate),nscal)   
      REAL_T  sx(DIMV(work),nscal),  sy(DIMV(work),nscal)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  dx(SDIM)

c ::: internal variables
      REAL_T  xlo(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal)
      REAL_T  xhi(DIMV(uedge),nscal)
      REAL_T  yhi(DIMV(vedge),nscal) 
      REAL_T hx, hy, dt, dth, dthx, dthy
      integer nc
      
      dth   = half*dt
      dthx  = half*dt / dx(1)
      dthy  = half*dt / dx(2)
      hx    = dx(1)
      hy    = dx(2)

c ::: compute the slopes
      do nc = 1,nscal
         call FORT_SLOPES(s(ARG_L1(s),ARG_L2(s),nc),DIMS(s),
     &        sx(ARG_L1(work),ARG_L2(work),nc),
     &        sy(ARG_L1(work),ARG_L2(work),nc),DIMS(work),
     &        lo,hi,slxscr,slyscr,bc)
      end do 

      call trace_tracers(s,DIMS(s),sn,DIMS(sn),
     &     st,DIMS(st),stn,DIMS(stn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),
     &     rphi,DIMS(rphi),sx,sy,DIMS(work),xlo,xhi,ylo,yhi,
     &     dt,dx,lo,hi,bc,nscal) 

c ::: predictor step: determine final edge states 
      call advect_tracers(xstate,DIMS(xstate),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,nscal)
      call advect_tracers(ystate,DIMS(ystate),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,nscal)

      end


      subroutine FORT_SYNC_RMN (s, DIMS(s), sn, DIMS(sn), 
     &     tf, divu, lbd_cc, dlbd_cc,
     &     rphi, DIMS(rphi), kappa, DIMS(kappa),
     &     sx, slxscr, 
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), kappax, DIMS(kpx),
     &     eigvx, eiglx, eigrx, 
     &     sy, slyscr,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), kappay, DIMS(kpy),
     &     eigvy, eigly, eigry, 
     &     DIMS(work),
     &     krcoef, DIMS(krcoef), nkrcoef,
     &     bc,lo,hi,dt,dx,
     &     use_minion,iconserv,nscal)
      
c
c ::: This subroutine computes edges states and fluxes.  It only makes sense 
c ::: if ncomps = nscal.
c
      implicit none

#include "probdata.H" 
      
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(work)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(xstate),DIMDEC(ystate)
      integer DIMDEC(kpx),DIMDEC(kpy)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer use_minion, iconserv
      integer nscal
      REAL_T   s(DIMV( s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  slxscr(DIM1(s),6),slyscr(DIM2(s),6)
      REAL_T  uedge(DIMV(uedge)), xstate(DIMV(xstate),nscal)
      REAL_T  vedge(DIMV(vedge)), ystate(DIMV(ystate),nscal)    
      REAL_T  kappax(DIMV(kpx)) , kappay(DIMV(kpy))
      REAL_T  sx(DIMV(work),nscal),  sy(DIMV(work),nscal)
      REAL_T  tf(DIMV(work))
      REAL_T  divu(DIMV(work))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  lbd_cc(DIMV(work),ncomps)
      REAL_T  dlbd_cc(DIMV(work),3)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  dx(SDIM)

c ::: internal variables
      REAL_T eigvx(DIMV(work),2*nscal),eigvy(DIMV(work),2*nscal)
      REAL_T eiglx(DIMV(work),nscal,nscal),eigly(DIMV(work),nscal,nscal)
      REAL_T eigrx(DIMV(work),nscal,nscal),eigry(DIMV(work),nscal,nscal)
      REAL_T, ALLOCATABLE :: fluxx(:,:,:), fluxy(:,:,:)
      REAL_T, ALLOCATABLE :: xlo(:,:,:), xhi(:,:,:)
      REAL_T, ALLOCATABLE :: ylo(:,:,:), yhi(:,:,:)
      REAL_T hx, hy, dt, dth, dthx, dthy, gravx, gravy
      integer nc, corrector
      
      dth   = half*dt
      dthx  = half*dt / dx(1)
      dthy  = half*dt / dx(2)
      hx    = dx(1)
      hy    = dx(2)
      gravx = zero
      gravy = -gravity

      ALLOCATE(fluxx(DIMV(uedge),nscal))
      ALLOCATE(fluxy(DIMV(vedge),nscal))
      ALLOCATE(xlo(DIMV(uedge),nscal))
      ALLOCATE(ylo(DIMV(vedge),nscal))
      ALLOCATE(xhi(DIMV(uedge),nscal))
      ALLOCATE(yhi(DIMV(vedge),nscal)) 

c ::: compute eigen solutions
      call eigsol_mp(eigvx,eiglx,eigrx,lbd_cc,dlbd_cc,DIMS(work),
     &     uedge,DIMS(uedge),rphi,DIMS(rphi),kappax,DIMS(kpx),     
     &     krcoef,DIMS(krcoef),nkrcoef,gravx,muval,
     &     1,0,lo,hi,rhoval,ncomps,nphases,model,bc)
      call eigsol_mp(eigvy,eigly,eigry,lbd_cc,dlbd_cc,DIMS(work),
     &     vedge,DIMS(vedge),rphi,DIMS(rphi),kappay,DIMS(kpy),
     &     krcoef,DIMS(krcoef),nkrcoef,gravy,muval,
     &     0,1,lo,hi,rhoval,ncomps,nphases,model,bc)

c ::: compute the slopes
      do nc = 1,nscal
         call FORT_SLOPES(s(ARG_L1(s),ARG_L2(s),nc),DIMS(s),
     &        sx(ARG_L1(work),ARG_L2(work),nc),
     &        sy(ARG_L1(work),ARG_L2(work),nc),DIMS(work),
     &        lo,hi,slxscr,slyscr,bc)
      end do 

c ::: initial step
      call trace(s,DIMS(s),sn,DIMS(sn),
     &     eigvx,eiglx,eigrx,eigvy,eigly,eigry,
     &     sx,sy,divu,tf,DIMS(work),
     &     uedge,xlo,xhi,DIMS(uedge),vedge,ylo,yhi,DIMS(vedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,nscal)
      
c ::: predictor step: determine final edge states 
      corrector = 0
      call riemann_solver(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),eiglx,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,lo,hi,gravx,muval,rhoval,nscal,model,corrector)
      call riemann_solver(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),eigly,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,lo,hi,gravy,muval,rhoval,nscal,model,corrector)

      call trace_p(s,DIMS(s),sn,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),
     &     fluxx,DIMS(uedge),fluxy,DIMS(vedge),
     &     xlo,xhi,ylo,yhi,divu,tf,DIMS(work),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     dt,dx,lo,hi,bc,muval,rhoval,iconserv,
     &     nscal,model,gravy)

c ::: compute fluxes based on final edge states
      call riemann_solver_arg(fluxx,uedge,xlo,xhi,DIMS(uedge),
     &     kappax,DIMS(kpx),eiglx,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     1,0,lo,hi,gravx,muval,rhoval,nscal,model)
      call riemann_solver_arg(fluxy,vedge,ylo,yhi,DIMS(vedge),
     &     kappay,DIMS(kpy),eigly,DIMS(work),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     0,1,lo,hi,gravy,muval,rhoval,nscal,model) 

      call flux_copy(xstate,DIMS(xstate),fluxx,DIMS(uedge),
     &     lo,hi,1,0,nscal)
      call flux_copy(ystate,DIMS(ystate),fluxy,DIMS(vedge),
     &     lo,hi,0,1,nscal)

      DEALLOCATE(fluxx,fluxy)
      DEALLOCATE(xlo,xhi,ylo,yhi)

      end


      subroutine trace_lin(s,DIMS(s),sn,DIMS(sn),
     &     rphi,DIMS(rphi),
     &     sx,sy,tf,DIMS(w),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     dt,dx,lo,hi,bc,nscal)
      implicit none

c ::: Trace the cell centered states to edge for linear advection.
c     Valid for single-phase flow.

c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(rphi)
      integer DIMDEC(w)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer nscal
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      REAL_T  s(DIMV(s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  tf(DIMV(w),nscal)
      REAL_T  sx(DIMV(w),nscal),  sy(DIMV(w),nscal)
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge))
      REAL_T  xlo(DIMV(uedge),nscal), xhi(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal), yhi(DIMV(vedge),nscal)
      REAL_T  dt, dx(2)

c ::: local variable
      integer old_time, i,j,nc
      REAL_T dth,dthx,dthy

      old_time = 0
      
      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)

c ::: x-direction
      do nc = 1,nscal
         do j = lo(2)-1,hi(2)+1
            do i = lo(1),hi(1)+1
               xlo(i,j,nc) = s(i-1,j,nc) + 
     &              (half-dthx*uedge(i,j)/rphi(i-1,j))*sx(i-1,j,nc) +
     &              dth*tf(i-1,j,nc)
               xhi(i,j,nc) = s(i,j,nc) + 
     &              (-half-dthx*uedge(i,j)/rphi(i,j))*sx(i,j,nc) +
     &              dth*tf(i,j,nc)
            end do
         end do
      end do

c ::: y-direction
      do nc = 1,nscal
         do j = lo(2),hi(2)+1
            do i = lo(1)-1,hi(1)+1
               ylo(i,j,nc) = s(i,j-1,nc) + 
     &              (half-dthy*vedge(i,j)/rphi(i,j-1))*sy(i,j-1,nc) +
     &              dth*tf(i,j-1,nc)
               yhi(i,j,nc) = s(i,j,nc) + 
     &              (-half-dthy*vedge(i,j)/rphi(i,j))*sy(i,j,nc) +
     &              dth*tf(i,j,nc)
            end do
         end do
      end do

c ::: fix edges at the boundaries.
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,bc,1,nscal,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,bc,1,nscal,old_time)
     
      end

      subroutine trace(s,DIMS(s),sn,DIMS(sn),
     &     eigvx,eiglx,eigrx,eigvy,eigly,eigry,
     &     sx,sy,divu,tf,DIMS(w),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,nscal)
      implicit none

c ::: trace the cell centered states to edges.

c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(w)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer nscal
      integer use_minion,iconserv
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      REAL_T  s(DIMV(s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  eigvx(DIMV(w),2*nscal),eigvy(DIMV(w),2*nscal)
      REAL_T  eiglx(DIMV(w),nscal,nscal),eigly(DIMV(w),nscal,nscal)
      REAL_T  eigrx(DIMV(w),nscal,nscal),eigry(DIMV(w),nscal,nscal)
      REAL_T  divu(DIMV(w))
      REAL_T  tf(DIMV(w),nscal)
      REAL_T  sx(DIMV(w),nscal),  sy(DIMV(w),nscal)
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge))
      REAL_T  xlo(DIMV(uedge),nscal), xhi(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal), yhi(DIMV(vedge),nscal)
      REAL_T  dt, dx(2)

c ::: local variable
      integer old_time
      REAL_T dth,dthx,dthy

      old_time = 0
      
      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)

c ::: trace to edges in each directions
      call trace_dir(s,DIMS(s),eigvx,eiglx,eigrx,DIMS(w),
     &     xlo,xhi,DIMS(uedge),sx,tf,divu,DIMS(w),lo,hi,
     &     1,0,dthx,dth,use_minion,iconserv,nscal)
      call trace_dir(s,DIMS(s),eigvy,eigly,eigry,DIMS(w),
     &     ylo,yhi,DIMS(vedge),sy,tf,divu,DIMS(w),lo,hi,
     &     0,1,dthy,dth,use_minion,iconserv,nscal)

c ::: fix edges at the boundaries.
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,bc,1,nscal,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,bc,1,nscal,old_time)
     
      end

      subroutine trace2(s,DIMS(s),sn,DIMS(sn),
     &     eigvx,eiglx,eigrx,eigvy,eigly,eigry,DIMS(eig),
     &     sx,sy,divu,tf,DIMS(w),
     &     uedge,xlo,xhi,DIMS(uedge),
     &     vedge,ylo,yhi,DIMS(vedge),
     &     dt,dx,lo,hi,bc,use_minion,iconserv,ncomps,nphases)
      implicit none
c ::: trace the cell centered states to edges.
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(w)
      integer DIMDEC(eig)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer ncomps,nphases
      integer use_minion,iconserv
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  sn(DIMV(sn),ncomps)
      REAL_T  divu(DIMV(w)),tf(DIMV(w))
      REAL_T  eigvx(DIMV(eig),ncomps),eigvy(DIMV(eig),ncomps)
      REAL_T  eigrx(DIMV(eig),ncomps,ncomps),eigry(DIMV(eig),ncomps,ncomps)
      REAL_T  eiglx(DIMV(eig),ncomps,ncomps),eigly(DIMV(eig),ncomps,ncomps)
      REAL_T  sx(DIMV(w),ncomps),  sy(DIMV(w),ncomps)
      REAL_T  uedge(DIMV(uedge))
      REAL_T  vedge(DIMV(vedge))
      REAL_T  xlo(DIMV(uedge),ncomps), xhi(DIMV(uedge),ncomps)
      REAL_T  ylo(DIMV(vedge),ncomps), yhi(DIMV(vedge),ncomps)
      REAL_T  dt, dx(SDIM)

c ::: local variable
      integer old_time
      REAL_T dth,dthx,dthy

      old_time  = 0
      
      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)

c ::: trace to edges in each directions
      call trace_dir2(s,DIMS(s),eigvx,eiglx,eigrx,DIMS(eig),
     &     xlo,xhi,DIMS(uedge),sx,tf,divu,DIMS(w),lo,hi,
     &     1,0,dthx,dth,use_minion,iconserv,ncomps,nphases)
      call trace_dir2(s,DIMS(s),eigvy,eigly,eigry,DIMS(eig),
     &     ylo,yhi,DIMS(vedge),sy,tf,divu,DIMS(w),lo,hi,
     &     0,1,dthy,dth,use_minion,iconserv,ncomps,nphases)
      
c ::: fix edges at the boundaries.
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,bc,1,ncomps,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,bc,1,ncomps,old_time)

      end

      subroutine trace_tracers(s,DIMS(s),sn,DIMS(sn),
     &     st,DIMS(st),stn,DIMS(stn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),
     &     phi,DIMS(phi),sx,sy,DIMS(w),xlo,xhi,ylo,yhi,
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
      integer DIMDEC(phi)
      integer nscal
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      REAL_T  s(DIMV(s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  st(DIMV(st)),stn(DIMV(stn))
      REAL_T  uedge(DIMV(uedge)),vedge(DIMV(vedge))
      REAL_T  phi(DIMV(phi))
      REAL_T  sx(DIMV(w),nscal),  sy(DIMV(w),nscal)
      REAL_T  xlo(DIMV(uedge),nscal), xhi(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal), yhi(DIMV(vedge),nscal)
      REAL_T  dt, dx(2)

c ::: local variable
      integer old_time

      REAL_T dth,dthx,dthy

      old_time = 0
      
      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)

c ::: trace to edges in each directions
      call trace_dir_tracer(s,DIMS(s),st,DIMS(st),stn,DIMS(stn),
     &     uedge,DIMS(uedge),phi,DIMS(phi),
     &     xlo,xhi,DIMS(uedge),sx,DIMS(w),lo,hi,
     &     1,0,dthx,nscal)
      call trace_dir_tracer(s,DIMS(s),st,DIMS(st),stn,DIMS(stn),
     &     vedge,DIMS(vedge),phi,DIMS(phi),
     &     ylo,yhi,DIMS(vedge),sy,DIMS(w),lo,hi,
     &     0,1,dthy,nscal)

c ::: fix edges at the boundaries.
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,bc,1,nscal,old_time)
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,bc,1,nscal,old_time)
     
      end

      subroutine trace_p_lin(s,DIMS(s),sn,DIMS(sn),
     &     fx,uedge,DIMS(uedge),fy,vedge,DIMS(vedge),
     &     rphi,DIMS(rphi),xlo,xhi,ylo,yhi,
     &     dt,dx,lo,hi,bc,nscal)

      implicit none
c
c ::: correct the edge values during based on the predictor stage
c     valid for single-phase flow.
c
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(rphi)
      integer bc(SDIM,2),lo(SDIM),hi(SDIM)
      integer nscal
      REAL_T  s(DIMV(s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  uedge(DIMV(uedge)),vedge(DIMV(vedge))
      REAL_T  fx(DIMV(uedge),nscal),fy(DIMV(vedge),nscal)
      REAL_T  xlo(DIMV(uedge),nscal),xhi(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal),yhi(DIMV(vedge),nscal)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  dt,dx(SDIM)
      
c ::: internal variables
      integer i,j,nc
      integer half_time
      REAL_T  dth,xhx,xhy
      REAL_T  st,drmndv

      half_time = 1
      
      dth  = half*dt
      xhx   = one/dx(1)
      xhy   = one/dx(2)

c ::: x-edges
      do nc = 1,nscal
         do j = lo(2),hi(2)
            do i = lo(1)-1,hi(1)+1
               drmndv = s(i,j,nc)/rphi(i,j)
               st = -(fy(i,j+1,nc)-fy(i,j,nc))/rphi(i,j)*xhy
c                 BEGIN HACK :: THIS HACK ASSUMES INCOMPRESSIBLE
c     &                - drmndv*(uedge(i+1,j)-uedge(i,j))*xhx
     &              + drmndv*(vedge(i,j+1)-vedge(i,j))*xhy
c                  END HACK
               xlo(i+1,j,nc) = xlo(i+1,j,nc) + dth*st
               xhi(i,j,nc)   = xhi(i,j,nc)   + dth*st
            end do  
         end do

         
      end do

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,bc,1,nscal,half_time)

c ::: y-edges
      do nc = 1,nscal
         do j = lo(2)-1,hi(2)+1
            do i = lo(1),hi(1)  
               drmndv = s(i,j,nc)/rphi(i,j)
               st = - (fx(i+1,j,nc)-fx(i,j,nc))/rphi(i,j)*xhx
c     &                - drmndv*(vedge(i,j+1)-vedge(i,j))/hy
     &              + drmndv*(uedge(i+1,j)-uedge(i,j))*xhx
               ylo(i,j+1,nc) = ylo(i,j+1,nc) + dth*st
               yhi(i,j,nc)   = yhi(i,j,nc)   + dth*st
            end do     
         end do        
      end do

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,bc,1,nscal,half_time)
      
      end

      subroutine trace_p(s,DIMS(s),sn,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),
     &     fx,DIMS(fx),fy,DIMS(fy),
     &     xlo,xhi,ylo,yhi,divu,tf,DIMS(w),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     krcoef, DIMS(krcoef), nkrcoef,
     &     dt,dx,lo,hi,bc,
     &     mu,rhoval,iconserv,nscal,model,grav)
c
c ::: correct the edge values during the predictor stage
c
      implicit none

c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(fx)   , DIMDEC(fy)
      integer DIMDEC(w)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa),DIMDEC(kpx),DIMDEC(kpy)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer bc(SDIM,2),lo(SDIM),hi(SDIM)
      integer nscal,model
      integer iconserv
      REAL_T  s(DIMV(s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  uedge(DIMV(uedge)),vedge(DIMV(vedge))
      REAL_T  kappax(DIMV(kpx)),kappay(DIMV(kpy))
      REAL_T  fx(DIMV(fx),nscal),fy(DIMV(fy),nscal)
      REAL_T  xlo(DIMV(uedge),nscal),xhi(DIMV(uedge),nscal)
      REAL_T  ylo(DIMV(vedge),nscal),yhi(DIMV(vedge),nscal)
      REAL_T  divu(DIMV(w)), tf(DIMV(w))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))    
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  mu(nscal),rhoval(nscal)
      REAL_T  dt,dx(SDIM)
      REAL_T  grav
      
c ::: internal variables
      integer i,j,nc
      integer imin,imax,jmin,jmax
      integer half_time
      REAL_T  dth,hx,hy
      REAL_T  st
      
      REAL_T, ALLOCATABLE::drmndv(:,:,:)
      REAL_T, ALLOCATABLE::drmndk(:,:,:)
      
      allocate(drmndv(DIMV(w),nscal),drmndk(DIMV(w),nscal))

      half_time = 1
      
      dth  = half*dt
      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      call dfluxdv(drmndv,drmndk,DIMS(w),s,DIMS(s),rphi,DIMS(rphi),
     &     kappa,DIMS(kappa),krcoef,DIMS(krcoef),nkrcoef,
     &     mu,rhoval,nscal,model,grav) 
c ::: x-edges
      do nc = 1,nscal
         do i = imin-1,imax+1
            do j = jmin,jmax
               if (iconserv.eq.1) then
                  st = -(fy(i,j+1,nc)-fy(i,j,nc))/hy/rphi(i,j)
c                 BEGIN HACK :: THIS HACK ASSUMES INCOMPRESSIBLE
c     &                 - drmndv(i,j,nc)*(uedge(i+1,j)-uedge(i,j))/hx
     &                 + drmndv(i,j,nc)*(vedge(i,j+1)-vedge(i,j))/hy
c                 END HACK
                  xlo(i+1,j,nc) = xlo(i+1,j,nc) + dth*st
                  xhi(i,j,nc)   = xhi(i,j,nc) + dth*st
                  
               endif 
            end do  
         end do

      end do
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,bc,1,nscal,half_time)

c ::: y-edges
      do nc = 1,nscal
         do j = jmin-1,jmax+1
            do i = imin,imax        
               if (iconserv.eq.1) then              
                  st = - (fx(i+1,j,nc)-fx(i,j,nc))/hx/rphi(i,j)
     &                 - drmndk(i,j,nc)*(kappay(i,j+1)-kappay(i,j))/hy
c     &                 - drmndv(i,j,nc)*(vedge(i,j+1)-vedge(i,j))/hy
     &			+ drmndv(i,j,nc)*(uedge(i+1,j)-uedge(i,j))/hx
                  ylo(i,j+1,nc) = ylo(i,j+1,nc) + dth*st
                  yhi(i,j,nc)   = yhi(i,j,nc) + dth*st
               end if
            end do     
         end do        
      end do

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,bc,1,nscal,half_time)

      deallocate(drmndv,drmndk)
      
      end

      subroutine trace_p2(s,DIMS(s),sn,DIMS(sn),
     &     uedge,DIMS(uedge),vedge,DIMS(vedge),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),
     &     fx,DIMS(fx),fy,DIMS(fy),
     &     xlo,xhi,ylo,yhi,divu,tf,DIMS(w),
     &     rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     dt,dx,lo,hi,bc,
     &     mu,rho,athick,iconserv,ncomps,nphases,grav)
      implicit none
c
c ::: correct the edge values during the predictor stage
c
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(uedge),DIMDEC(vedge)
      integer DIMDEC(fx)   , DIMDEC(fy)
      integer DIMDEC(w)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa),DIMDEC(kpx),DIMDEC(kpy)
      integer bc(SDIM,2),lo(SDIM),hi(SDIM)
      integer ncomps,nphases,iconserv
      REAL_T   s(DIMV(s),ncomps)
      REAL_T  sn(DIMV(sn),ncomps)
      REAL_T  uedge(DIMV(uedge)), vedge(DIMV(vedge))
      REAL_T  kappax(DIMV(kpx)),kappay(DIMV(kpy))
      REAL_T  fx(DIMV(fx),ncomps),fy(DIMV(fy),ncomps)
      REAL_T  xlo(DIMV(uedge),ncomps),xhi(DIMV(uedge),ncomps)
      REAL_T  ylo(DIMV(vedge),ncomps),yhi(DIMV(vedge),ncomps)
      REAL_T  divu(DIMV(w)), tf(DIMV(w))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  mu(nphases),rho(nphases),athick
      REAL_T  dt,dx(SDIM)
      REAL_T  grav
      
c ::: internal variables
      integer i,j,nc
      integer imin,imax,jmin,jmax
      integer half_time
      REAL_T  dth,hx,hy
      REAL_T  st,u(2),drmndv(2),drmndk(2)

      half_time = 1
      
      dth  = half*dt
      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)     

      do i = imin,imax+1
         do j = jmin,jmax
            do nc = 1,ncomps
               u(nc) = s(i,j,nc)
            end do

            call dfluxdv2(drmndv,drmndk,u,rphi(i,j),kappa(i,j),
     &           mu(1),mu(2),athick,rho(1),rho(2),grav)
            do nc = 1,ncomps
               if (iconserv.eq.1) then
                  st = -(fy(i,j+1,nc)-fy(i,j,nc))/hy/rphi(i,j)
     &                 - drmndv(nc)*(uedge(i+1,j)-uedge(i,j))/hx

                  xlo(i+1,j,nc) = xlo(i+1,j,nc) + dth*st
                  xhi(i,j,nc)   = xhi(i,j,nc) + dth*st

               end if
            end do
         end do
      end do

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(uedge),
     &     1,0,lo,hi,bc,1,ncomps,half_time)
      
      do j = jmin,jmax+1
         do i = imin,imax
            do nc = 1,ncomps
               u(nc) = s(i,j,nc)
            end do
            
            call dfluxdv2(drmndv,drmndk,u,rphi(i,j),kappa(i,j),
     &           mu(1),mu(2),athick,rho(1),rho(2),grav)
            
            do nc = 1,ncomps
               if (iconserv.eq.1) then
                  st = - (fx(i+1,j,nc)-fx(i,j,nc))/hx/rphi(i,j)
     &                 - drmndk(nc)*(kappay(i,j+1)-kappay(i,j))/hy
     &                 - drmndv(nc)*(vedge(i,j+1)-vedge(i,j))/hy
                  ylo(i,j+1,nc) = ylo(i,j+1,nc) + dth*st
                  yhi(i,j,nc)   = yhi(i,j,nc)   + dth*st
               endif
            end do         
         end do
      end do
      
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(vedge),
     &     0,1,lo,hi,bc,1,ncomps,half_time)

      end


      subroutine trace_p_tracer(s,DIMS(s),sn,DIMS(sn),
     &     st,DIMS(st),stn,DIMS(stn),
     &     fx,uedge,xlo,xhi,DIMS(fx),fy,vedge,ylo,yhi,DIMS(fy),
     &     rphi,DIMS(rphi),dt,dx,lo,hi,bc,nscal)
      
c
c ::: correct the edge values during the predictor stage
c

      implicit none

c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(st)
      integer DIMDEC(stn)
      integer DIMDEC(fx), DIMDEC(fy)
      integer DIMDEC(rphi)
      integer bc(SDIM,2),lo(SDIM),hi(SDIM)
      integer nscal
      REAL_T  s(DIMV(s),nscal)
      REAL_T  sn(DIMV(sn),nscal)
      REAL_T  st(DIMV(st))
      REAL_T  stn(DIMV(stn))
      REAL_T  uedge(DIMV(fx)),vedge(DIMV(fy))
      REAL_T  fx(DIMV(fx),nscal),fy(DIMV(fy),nscal)
      REAL_T  xlo(DIMV(fx),nscal),xhi(DIMV(fx),nscal)
      REAL_T  ylo(DIMV(fy),nscal),yhi(DIMV(fy),nscal)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  dt,dx(SDIM)
      
c ::: internal variables
      integer i,j,nc
      integer imin,imax,jmin,jmax
      integer half_time
      REAL_T  dth,hx,hy
      REAL_T  stp

      half_time = 1
      
      dth  = half*dt
      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

c ::: x-edges
      do nc = 1,nscal
         do i = imin-1,imax+1
            do j = jmin,jmax
               stp = -(fy(i,j+1,nc)-fy(i,j,nc))/hy/rphi(i,j)
               xlo(i+1,j,nc) = xlo(i+1,j,nc) + dth*stp
               xhi(i,j,nc)   = xhi(i,j,nc) + dth*stp
            end do  
         end do

      end do
      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),uedge,xlo,xhi,DIMS(fx),
     &     1,0,lo,hi,bc,1,nscal,half_time)

c ::: y-edges
      do nc = 1,nscal
         do j = jmin-1,jmax+1
            do i = imin,imax            
               stp = - (fx(i+1,j,nc)-fx(i,j,nc))/hx/rphi(i,j)
               ylo(i,j+1,nc) = ylo(i,j+1,nc) + dth*stp
               yhi(i,j,nc)   = yhi(i,j,nc) + dth*stp
            end do     
         end do        
      end do

      call trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vedge,ylo,yhi,DIMS(fy),
     &     0,1,lo,hi,bc,1,nscal,half_time)
      
      end



      subroutine trace_dir(s,DIMS(s),eigv,eigl,eigr,DIMS(eig),
     &     vlo,vhi,DIMS(v),dv,tf,divu,DIMS(w),lo,hi,
     &     ix,iy,dthx,dth,use_minion,iconserv,nscal)
      implicit none
c ::: arguments
      integer DIMDEC(s), DIMDEC(eig), DIMDEC(w)
      integer DIMDEC(v)
      integer lo(SDIM),hi(SDIM)
      integer ix,iy,nscal,use_minion,iconserv
      REAL_T  s(DIMV(s),nscal)
      REAL_T  eigv(DIMV(eig),2*nscal)
      REAL_T  eigl(DIMV(eig),nscal,nscal)
      REAL_T  eigr(DIMV(eig),nscal,nscal) 
      REAL_T  vlo(DIMV(v),nscal)
      REAL_T  vhi(DIMV(v),nscal)
      REAL_T  dv(DIMV(w) ,nscal)
      REAL_T  tf(DIMV(w),nscal)
      REAL_T  divu(DIMV(w))
      REAL_T  dthx, dth

c ::: internal variables
      integer i,j,nc,nd
      integer imin,jmin,imax,jmax
      REAL_T  alphalo,alphahi
c     REAL_T  tlo,thi

      if (ix .eq. 1) then
         imin = lo(1)
         imax = hi(1) + 1
         jmin = lo(2) - 1
         jmax = hi(2) + 1
      elseif (iy .eq. 1) then
         imin = lo(1) - 1
         imax = hi(1) + 1
         jmin = lo(2) 
         jmax = hi(2) + 1
      end if         

      do j = jmin,jmax
         do i = imin,imax
            do nc = 1,nscal
               vlo(i,j,nc) = s(i-ix,j-iy,nc)
               vhi(i,j,nc) = s(i,j,nc)
            end do

            do nd = 1,nscal
               alphalo=0
               alphahi=0
               do nc = 1,nscal
                  alphalo = alphalo + eigl(i-ix,j-iy,nd,nc)*dv(i-ix,j-iy,nc)
                  alphahi = alphahi + eigl(i,j,nd,nc)*dv(i,j,nc)
               end do               
               do nc = 1,nscal                  
                  vlo(i,j,nc) = vlo(i,j,nc) + ( half
     &                 - dthx*max(zero,eigv(i-ix,j-iy,nscal+nd))) 
     &                 *alphalo*eigr(i-ix,j-iy,nc,nd)
                  vhi(i,j,nc) = vhi(i,j,nc) + (-half
     &                 - dthx*min(zero,eigv(i,j,nd)))
     &                 *alphahi*eigr(i,j,nc,nd)
               end do
            end do 
         end do
      end do

      do nc = 1, nscal
         do j = jmin,jmax
            do i = imin,imax
               vlo(i,j,nc) = vlo(i,j,nc) + dth*tf(i-ix,j-iy,nc)
               vhi(i,j,nc) = vhi(i,j,nc) + dth*tf(i,j,nc)
            end do
         end do

      end do

      if (iconserv .eq. 1) then
         do nc = 1, nscal
            do j = jmin,jmax
               do i = imin,imax
                  vlo(i,j,nc)=vlo(i,j,nc) 
     &                 - dth*s(i-ix,j-iy,nc)*divu(i-ix,j-iy)
                  vhi(i,j,nc)=vhi(i,j,nc) 
     &                 - dth*s(i,j,nc)*divu(i,j)
               end do
            end do               
         end do
      end if
      
      end

      subroutine trace_dir_tracer(s,DIMS(s),
     &     st,DIMS(st),stn,DIMS(stn),
     &     uedge,DIMS(uedge),phi,DIMS(phi),
     &     vlo,vhi,DIMS(v),dv,DIMS(w),lo,hi,
     &     ix,iy,dthx,nscal)

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
      integer ix,iy,nscal
      REAL_T  s(DIMV(s),nscal)
      REAL_T  st(DIMV(st)),stn(DIMV(stn))
      REAL_T  uedge(DIMV(uedge))
      REAL_T  phi(DIMV(phi))
      REAL_T  vlo(DIMV(v),nscal)
      REAL_T  vhi(DIMV(v),nscal)
      REAL_T  dv(DIMV(w),nscal)
      REAL_T  dthx

c ::: internal variables
      integer i,j,nc
      REAL_T  alphalo,alphahi,sthalf

      do j = lo(2)-ix,hi(2)+1
         do i = lo(1)-iy,hi(1)+1
            do nc = 1,nscal
               vlo(i,j,nc) = s(i-ix,j-iy,nc)
               vhi(i,j,nc) = s(i,j,nc)

               sthalf  = half*(st(i-ix,j-iy)+stn(i-ix,j-iy))
               alphalo = uedge(i,j)/(phi(i-ix,j-iy)*sthalf)
               sthalf  = half*(st(i,j)+stn(i,j))
               alphahi = uedge(i,j)/(phi(i,j)*sthalf)

               vlo(i,j,nc) = vlo(i,j,nc) + (half
     &                 - dthx*max(zero,alphalo))*dv(i-ix,j-iy,nc)
               vhi(i,j,nc) = vhi(i,j,nc) + (-half
     &                 - dthx*min(zero,alphahi))*dv(i,j,nc)
            end do
         end do
      end do

      end
 
      subroutine trace_dir2(s,DIMS(s),eigv,eigl,eigr,DIMS(eig),
     &     vlo,vhi,DIMS(v),dv,tf,divu,DIMS(w),lo,hi,
     &     ix,iy,dthx,dth,use_minion,iconserv,ncomps,nphases)
      implicit none
c ::: arguments
      integer DIMDEC(s), DIMDEC(eig), DIMDEC(w), DIMDEC(v)
      integer lo(SDIM),hi(SDIM)
      integer ix,iy,ncomps,nphases,use_minion,iconserv
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  eigv(DIMV(eig),ncomps)
      REAL_T  eigl(DIMV(eig),ncomps,ncomps) 
      REAL_T  eigr(DIMV(eig),ncomps,ncomps) 
      REAL_T  vlo(DIMV(v),ncomps),vhi(DIMV(v),ncomps)
      REAL_T  dv(DIMV(w),ncomps)
      REAL_T  tf(DIMV(w)),divu(DIMV(w))
      REAL_T  dthx, dth

c ::: internal variables
      integer i,j,nc,nd
      integer imin,jmin,imax,jmax
      REAL_T  alphalo, alphahi
 
      if (ix .eq. 1) then
         imin = lo(1)
         imax = hi(1) + 1
         jmin = lo(2) - 1
         jmax = hi(2) + 1
      elseif (iy .eq. 1) then
         imin = lo(1) - 1
         imax = hi(1) + 1
         jmin = lo(2) 
         jmax = hi(2) + 1
      end if 

      do j = jmin,jmax
         do i = imin,imax
            do nc = 1,ncomps
               vlo(i,j,nc) = s(i-ix,j-iy,nc)
               vhi(i,j,nc) = s(i,j,nc)
            end do
 
            do nd = 1,ncomps  
               alphalo = merge(zero, one, eigv(i-ix,j-iy,nd) .lt. zero)
               alphahi = merge(zero, one, eigv(i,j,nd) .gt. zero)
         
              do nc = 1,ncomps   
                  vlo(i,j,nc) = vlo(i,j,nc) + half*
     &                (alphalo*(one - two*dthx*eigv(i-ix,j-iy,nd)))
     &                 *dv(i-ix,j-iy,nd)*eigr(i-ix,j-iy,nc,nd)
                  vhi(i,j,nc) = vhi(i,j,nc) - half*
     &                (alphahi*(one + two*dthx*eigv(i,j,nd)))
     &                 *dv(i,j,nd)*eigr(i,j,nc,nd)
               end do
            end do
         end do
      end do

      if(use_minion.eq.1)then
         do nc = 1, ncomps
            do j = jmin,jmax
               do i = imin,imax
                  vlo(i+ix,j+iy,nc) = vlo(i+ix,j+iy,nc) + dth*tf(i,j)
                  vhi(i,j,nc) = vhi(i,j,nc) + dth*tf(i,j)
               end do
            end do
         end do

         if (iconserv .eq. 1) then
            do nc = 1, ncomps
            do j = jmin,jmax
               do i = imin,imax
                     vlo(i+ix,j+iy,nc)=vlo(i+ix,j+iy,nc) 
     &                    - dth*s(i,j,nc)*divu(i,j)
                     vhi(i,j,nc)=vhi(i,j,nc) 
     &                    - dth*s(i,j,nc)*divu(i,j)
                  end do
               end do               
            end do
         end if
      end if

      end

      subroutine trans_bc_rmn(s,DIMS(s),sn,DIMS(sn),vel,vlo,vhi,DIMS(v),
     &                        ix,iy,lo,hi,bc,idx_first,nscal,which_time)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the specified direction based on ix and iy.  
c

      implicit none

#include "probdata.H"

c ::: arguments
      integer DIMDEC(s), DIMDEC(sn), DIMDEC(v)
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      integer ix,iy,idx_first,nscal
      integer which_time
      REAL_T  s(DIMV( s),nscal)
      REAL_T sn(DIMV(sn),nscal)
      REAL_T vel(DIMV(v))
      REAL_T vlo(DIMV(v),nscal)
      REAL_T vhi(DIMV(v),nscal)
      
c ::: local variables
      integer i,j,nc
      integer ibc,imin,jmin,imax,jmax
      integer old_time, half_time
      REAL_T  sat_threshold

      sat_threshold = one - 1.d-2

       old_time = 0
      half_time = 1

      ibc  = ix + 2*iy 
      imin = (lo(1)-1)*(1-ix) + lo(1)*ix
      imax = (hi(1)+1)*(1-ix) + lo(1)*ix
      jmin = (lo(2)-1)*(1-iy) + lo(2)*iy
      jmax = (hi(2)+1)*(1-iy) + lo(2)*iy
         
      if (bc(ibc,1).eq.EXT_DIR) then
         if (model .eq. 0 .or. model .eq. 1) then
            if (which_time .eq. old_time) then
               do nc = idx_first,nscal
                  do i = imin,imax
                     do j = jmin,jmax
c                        if (vel(i,j) .gt. zero) then
c                           vlo(i,j,nc) = s(i-ix,j-iy,nc)
c                        else
c                           vlo(i,j,nc) = vhi(i,j,nc)
c                        end if
                        vlo(i,j,nc) = s(i-ix,j-iy,nc)
                     end do
                  end do
               end do
            else if (which_time .eq. half_time) then
               do nc = idx_first,nscal
                  do i = imin,imax
                     do j = jmin,jmax
c                        if (vel(i,j)  .gt. zero) then
                           vlo(i,j,nc) = 0.5d0 * (s(i-ix,j-iy,nc) + sn(i-ix,j-iy,nc))
c                        else
c                           vlo(i,j,nc) = vhi(i,j,nc)
c                        end if
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
                  do i = imin,imax
                     do j = jmin,jmax
c                        if (vel(i,j)  .gt. zero) then
                           vlo(i,j,nc) = s(i-ix,j-iy,nc)
c                        else
c                           vlo(i,j,nc) = vhi(i,j,nc)
c                        end if
                     end do
                  end do
               end do
            else if (which_time .eq. half_time) then
               do nc = idx_first,nscal
                   do i = imin,imax
                     do j = jmin,jmax
c                        if (vel(i,j)  .gt. zero) then
                           vlo(i,j,nc) = 0.5d0 * (s(i-ix,j-iy,nc) + sn(i-ix,j-iy,nc))
c                        else
c                           vlo(i,j,nc) = vhi(i,j,nc)
c                        end if
                     end do
                  end do
               end do
            else 
               print *,'BAD WHICH_TIME IN TRANS_BC_RMN ',which_time
               stop
            end if
         end if

c$$$      else if (bc(ibc,1).eq.SEEPAGE) then
c$$$         if (model .eq. 0 .or. model .eq. 1) then
c$$$            do i = imin,imax
c$$$               do j = jmin,jmax
c$$$                  do nc = idx_first,nscal
c$$$                     vlo(i,j,nc) = vhi(i,j,nc)
c$$$                  end do
c$$$               end do
c$$$            end do
c$$$
c$$$         elseif (model .eq. 2 .or. model .eq. 3) then
c$$$            do i = imin,imax
c$$$               do j = jmin,jmax
c$$$                  if (which_time .eq. old_time) then
c$$$                     call get_st(sat,s(i,j,:),rhoval,sat_residual,ncomps)
c$$$                  else 
c$$$                     call get_st(sat,sn(i,j,:),rhoval,sat_residual,ncomps)
c$$$                  end if
c$$$                  vlo(i,j,1) = merge(zero, rhoval(1), sat(1).lt.sat_threshold)
c$$$                  vlo(i,j,2) = merge(rhoval(2), zero, sat(1).lt.sat_threshold)
c$$$                  vhi(i,j,1:2) = vlo(i,j,1:2)
c$$$               end do
c$$$            end do
c$$$         end if

      else if (bc(ibc,1).eq.SEEPAGE  .or. bc(ibc,1).eq.FOEXTRAP .or. 
     &         bc(ibc,1).eq.HOEXTRAP .or. bc(ibc,1).eq.REFLECT_EVEN) then
         do i = imin,imax
            do j = jmin,jmax
               do nc = idx_first,nscal
                  vlo(i,j,nc) = vhi(i,j,nc)
               end do
            end do
         end do

      else if (bc(ibc,1).eq.REFLECT_ODD) then
         do i = imin,imax
            do j = jmin,jmax
               do nc = idx_first,nscal
                  vhi(i,j,nc) = zero
                  vlo(i,j,nc) = zero
               end do
            end do
         end do
      end if

      imin = (lo(1)-1)*(1-ix) + (hi(1)+1)*ix
      imax = (hi(1)+1)*(1-ix) + (hi(1)+1)*ix
      jmin = (lo(2)-1)*(1-iy) + (hi(2)+1)*iy
      jmax = (hi(2)+1)*(1-iy) + (hi(2)+1)*iy

      if (bc(ibc,2).eq.EXT_DIR) then
         if (model .eq. 0 .or. model .eq. 1) then
            if (which_time .eq. old_time) then
               do nc = idx_first,nscal
                  do i = imin,imax
                     do j = jmin,jmax
c                        if (vel(i,j)  .lt. zero) then
                           vhi(i,j,nc) = s(i,j,nc) 
c                        else
c                           vhi(i,j,nc) = vlo(i,j,nc)
c                        end if
                     end do
                  end do
               end do
            else if (which_time .eq. half_time) then
               do nc = idx_first,nscal
                  do i = imin,imax
                     do j = jmin,jmax
c                       if (vel(i,j)  .lt. zero) then
                           vhi(i,j,nc) = 0.5d0 * (s(i,j,nc)+sn(i,j,nc))
c                       else
c                           vhi(i,j,nc) = vlo(i,j,nc)
c                       end if
                     end do
                  end do
               end do
            else 
               print *,'BAD WHICH_TIME IN TRANS_BC_RMN ',which_time
               stop
            end if
         elseif (model  .eq. 2 .or. model .eq. 3) then
            if (which_time .eq. old_time) then
               do nc = idx_first,nscal
                  do i = imin,imax
                     do j = jmin,jmax
c                        if (vel(i,j)  .lt. zero) then
                           vhi(i,j,nc) = s(i,j,nc) 
c                        else
c                           vhi(i,j,nc) = vlo(i,j,nc)
c                        end if
                     end do
                  end do
               end do
            else if (which_time .eq. half_time) then
               do nc = idx_first,nscal
                  do i = imin,imax
                     do j = jmin,jmax
c                        if (vel(i,j)  .lt. zero) then
                           vhi(i,j,nc) = 0.5d0 * (s(i,j,nc)+sn(i,j,nc))
c                        else
c                           vhi(i,j,nc) = vlo(i,j,nc)
c                        end if
                     end do
                  end do
               end do
            else 
               print *,'BAD WHICH_TIME IN TRANS_BC_RMN ',which_time
               stop
            end if
         end if


c      else if (bc(ibc,2).eq.SEEPAGE) then
c         if (model .eq. 0 .or. model .eq. 1) then
c            do i = imin,imax
c               do j = jmin,jmax
c                  do nc = idx_first,nscal
c                     vhi(i,j,nc) = vlo(i,j,nc)
c                  end do
c               end do
c            end do

c         elseif (model .eq. 2 .or. model .eq. 3) then
c            do i = imin,imax
c               do j = jmin,jmax
c                  if (which_time .eq. old_time) then
c                     call get_st(sat,s(i-ix,j-iy,:),rhoval,sat_residual,ncomps)
c                  else
c                     call get_st(sat,sn(i-ix,j-iy,:),rhoval,sat_residual,ncomps)
c                  end if
c                  vlo(i,j,1) = merge(zero, rhoval(1), sat(1).lt.sat_threshold)
c                  vlo(i,j,2) = merge(rhoval(2), zero, sat(1).lt.sat_threshold)
c                  vhi(i,j,1:2) = vlo(i,j,1:2)
c                  vhi(i,j,1) = zero
c                  if (nscal .gt. 1) then
c                     vhi(i,j,2) = rhoval(2)
c                  end if
c               end do
c            end do
c         end if

      else if (bc(ibc,2).eq.SEEPAGE.or.bc(ibc,2).eq.FOEXTRAP.or.
     &         bc(ibc,2).eq.HOEXTRAP.or.bc(ibc,2).eq.REFLECT_EVEN) then
         do i = imin,imax
            do j = jmin,jmax
               do nc = idx_first,nscal
                  vhi(i,j,nc) = vlo(i,j,nc)
               end do
            end do
         end do
      else if (bc(ibc,2).eq.REFLECT_ODD) then
         do i = imin,imax
            do j = jmin,jmax
               do nc = idx_first,nscal
                  vhi(i,j,nc) = zero
                  vlo(i,j,nc) = zero
               end do
            end do
         end do
      end if

      end

      subroutine vel_bc_rmn(vel,DIMS(vel),ix,iy,lo,hi,bc)
      implicit none
c
c     This subroutine replaces the velocity at lower edge of ghost shell to that of 
c     higher edge but only for those on domain boundary.  
c
c ::: arguments
      integer DIMDEC(vel)
      integer lo(SDIM), hi(SDIM), bc(SDIM,2)
      integer ix,iy
      REAL_T vel(DIMV(vel))
      
c ::: local variables
      integer i,j
      integer ibc,imin,jmin,imax,jmax

      ibc  = ix + 2*iy 
         
      if (bc(ibc,1) .eq. EXT_DIR .or. bc(ibc,1) .eq. FOEXTRAP  .or.
     &    bc(ibc,1) .eq. FOEXTRAP .or. bc(ibc,1) .eq. REFLECT_EVEN) then

         imin = (lo(1)-1)*(1-ix) + lo(1)*ix
         imax = (hi(1)+1)*(1-ix) + lo(1)*ix
         jmin = (lo(2)-1)*(1-iy) + lo(2)*iy
         jmax = (hi(2)+1)*(1-iy) + lo(2)*iy

         do i = imin,imax
            do j = jmin,jmax
               vel(i-ix,j-iy) = vel(i,j)
            end do
         end do

      elseif (bc(ibc,2) .eq. EXT_DIR .or. bc(ibc,2) .eq. FOEXTRAP .or.
     &    bc(ibc,2) .eq. FOEXTRAP .or. bc(ibc,2) .eq. REFLECT_EVEN) then

         imin = (lo(1)-1)*(1-ix) + (hi(1)+1)*ix
         imax = (hi(1)+1)*(1-ix) + (hi(1)+1)*ix
         jmin = (lo(2)-1)*(1-iy) + (hi(2)+1)*iy
         jmax = (hi(2)+1)*(1-iy) + (hi(2)+1)*iy

         do i = imin,imax
            do j = jmin,jmax
               vel(i+ix,j+iy) = vel(i,j)
            end do
         end do
      endif

      end

      subroutine FORT_ADV_FORCING(
     &     aofs,DIMS(aofs),
     &     xflux,DIMS(xflux),
     &     uedge,DIMS(uedge),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vedge,DIMS(vedge),
     &     areay,DIMS(ay),
     &     vol,DIMS(vol),
     &     lo,hi,iconserv )
c
c     Uses scalar edge states to compute an advective tendency
c     This is only used with the older version.
c
      implicit none
      integer i,j
      integer iconserv
      REAL_T divux,divuy
      integer imin,jmin,imax,jmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(aofs)
      integer DIMDEC(vol)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T aofs(DIMV(aofs))
      REAL_T vol(DIMV(vol))
      REAL_T uedge(DIMV(uedge))
      REAL_T vedge(DIMV(vedge))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
 
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
c     if nonconservative initialize the advective tendency as -U*grad(S)
c
      if ( iconserv .ne. 1 ) then
         do j = jmin,jmax
            do i = imin,imax
               divux = (
     &              areax(i+1,j)*uedge(i+1,j) -
     &              areax(i,  j)*uedge(i,  j) )/vol(i,j)
               divuy = (
     &              areay(i,j+1)*vedge(i,j+1) -
     &              areay(i,j  )*vedge(i,j  ) )/vol(i,j)
               aofs(i,j) =
     &              - divux*half*(xflux(i+1,j) + xflux(i,j))
     &              - divuy*half*(yflux(i,j+1) + yflux(i,j))              
            end do
         end do
      end if
c
c     convert edge states to fluxes
c  
      do j = jmin,jmax
         do i = imin,imax+1
            xflux(i,j) = xflux(i,j)*uedge(i,j)*areax(i,j)
         end do
      end do
      do j = jmin,jmax+1
         do i = imin,imax
            yflux(i,j) = yflux(i,j)*vedge(i,j)*areay(i,j)
         end do
      end do
c     
c     compute part of advective tendency that depends on the flux convergence
c     
      if ( iconserv .ne. 1 ) then
         do j = jmin,jmax
            do i = imin,imax
               aofs(i,j) = aofs(i,j) + (
     &              xflux(i+1,j) - xflux(i,j) +
     &              yflux(i,j+1) - yflux(i,j))/vol(i,j)
            end do
         end do
      else
         do j = jmin,jmax
            do i = imin,imax
               aofs(i,j) = (
     &              xflux(i+1,j) - xflux(i,j) +
     &              yflux(i,j+1) - yflux(i,j))/vol(i,j)
            end do
         end do
      end if

      end

      subroutine FORT_ADV_RMN_FORCING(
     &     aofs,DIMS(aofs),
     &     xflux,DIMS(xflux),
     &     uedge,DIMS(uedge),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vedge,DIMS(vedge),
     &     areay,DIMS(ay),
     &     vol,DIMS(vol),
     &     lo,hi,nscal)
      implicit none
c
c     This subroutine uses scalar edge states to compute
c     an advective tendency
c
#include "probdata.H"

      integer i,j,nc
      integer imin,jmin,imax,jmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(aofs)
      integer DIMDEC(vol)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer nscal
      REAL_T aofs(DIMV(aofs),nscal)
      REAL_T vol(DIMV(vol))
      REAL_T uedge(DIMV(uedge))
      REAL_T vedge(DIMV(vedge))
      REAL_T xflux(DIMV(xflux),nscal)
      REAL_T yflux(DIMV(yflux),nscal)
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

c ::: fluxes times area
      do nc = 1,nscal
         do j = jmin,jmax
            do i = imin,imax+1
               xflux(i,j,nc) = xflux(i,j,nc)*areax(i,j)
            end do
         end do
         do j = jmin,jmax+1
            do i = imin,imax
               yflux(i,j,nc) = yflux(i,j,nc)*areay(i,j)
            end do
         end do
      end do

c ::: compute part of advective tendency that depends on the flux convergence
      do j = jmin,jmax
         do i = imin,imax
            do nc = 1,nscal
               aofs(i,j,nc) = (
     &              xflux(i+1,j,nc) - xflux(i,j,nc) +
     &              yflux(i,j+1,nc) - yflux(i,j,nc))/vol(i,j)
            end do
         end do
      end do
      end

      subroutine riemann_solver_lin(flux,vel,vlo,vhi,DIMS(v),
     &     ix,iy,lo,hi,nscal)

c ::: Riemann solver for linear advection.  
c     valid for single-phase flow.

      implicit none

c ::: arguments
      integer DIMDEC(v)
      integer ix,iy
      integer nscal
      integer lo(SDIM),hi(SDIM)
      REAL_T  flux(DIMV(v),nscal)
      REAL_T  vel(DIMV(v))
      REAL_T  vlo(DIMV(v),nscal), vhi(DIMV(v),nscal)

c ::: internal variables
      integer i,j,nc    
      integer imin,jmin,imax,jmax

      if (ix .eq. 1) then
         imin = lo(1)
         imax = hi(1) + 1
         jmin = lo(2) - 1
         jmax = hi(2) + 1
      elseif (iy .eq. 1) then
         imin = lo(1) - 1 
         imax = hi(1) + 1 
         jmin = lo(2) 
         jmax = hi(2) + 1
      end if

      do nc = 1,nscal
         do j = jmin,jmax
            do i = imin,imax 
               flux(i,j,nc) = merge(vlo(i,j,nc),vhi(i,j,nc),vel(i,j)>=zero)
               flux(i,j,nc) = merge(half*(vlo(i,j,nc)+vhi(i,j,nc)),
     &              flux(i,j,nc),abs(vel(i,j)).lt. 1.d-16)
               flux(i,j,nc) = flux(i,j,nc) * vel(i,j)
            end do
         end do
      end do
         
      end

      subroutine riemann_solver(flux,vel,vlo,vhi,DIMS(v),
     &     kpedge,DIMS(kp),
     &     eigl,DIMS(eig),krcoef,DIMS(krcoef),nkrcoef,
     &     ix,iy,lo,hi,
     &     grav,muval,rhoval,nscal,model,corrector)

c ::: Riemann Solver for two phase flow.
c     Determine the flux values along edges.
c     eigv(1:2) is overwrite with the velocity of each phase 
c     for subsequent advection of tracers.

      implicit none

c ::: arguments
      integer DIMDEC(v)
      integer DIMDEC(kp)
      integer DIMDEC(eig)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ix,iy
      integer nscal,model
      integer lo(SDIM),hi(SDIM)
      integer corrector
      REAL_T  flux(DIMV(v),nscal)
      REAL_T  vel(DIMV(v))
      REAL_T  vlo(DIMV(v),nscal), vhi(DIMV(v),nscal)
      REAL_T  kpedge(DIMV(kp))
      REAL_T  eigl(DIMV(eig),nscal,nscal)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  rhoval(nscal),muval(nscal),grav

c ::: internal variables
      integer i,j,nc 
      integer imin,jmin,imax,jmax 
      integer l_kr_type, r_kr_type
      REAL_T  l_kr_coef, l_sat_res, r_kr_coef, r_sat_res
      REAL_T  aright, aleft, stmp(2), rtmp(2), st(2),gstar
      logical no_change

      imin = lo(1) - 1 + ix
      imax = hi(1) + 1
      jmin = lo(2) - 1 + iy
      jmax = hi(2) + 1

      if (model .eq. 2) then
         rtmp(:) = one
         do j = jmin,jmax
            do i = imin,imax
               r_kr_type = nint(krcoef(i,j,1))
               l_kr_type = nint(krcoef(i-ix,j-iy,1))  
               r_kr_coef = krcoef(i,j,2)
               l_kr_coef = krcoef(i-ix,j-iy,2)
               r_sat_res = krcoef(i,j,3)
               l_sat_res = krcoef(i-ix,j-iy,3)
               no_change = .true.

               aright = zero
               aleft  = zero
               do nc = 1,nscal
                  aright = aright + eigl(i,j,2,nc)*vhi(i,j,nc)
                  aleft  = aleft  + eigl(i-ix,j-iy,2,nc)*vlo(i,j,nc)
               end do 

               stmp(1) = aright
               stmp(2) = one-stmp(1)
               call get_st(st,stmp,rtmp,r_sat_res,nscal)
               aright = st(1)

               stmp(1) = aleft
               stmp(2) = one-stmp(1) 
               call get_st(st,stmp,rtmp,l_sat_res,nscal)
               aleft = st(1)               

               if (r_kr_type .ne. l_kr_type .or. r_kr_coef .ne. l_kr_coef .or.
     &              r_sat_res .ne. l_sat_res) then
                  no_change = .false.
               end if

               flux(i,j,:) = zero
               gstar = kpedge(i,j)*grav*(rhoval(1)-rhoval(2)) 
               if (no_change) then                   
                  if (l_kr_type .eq. 2) then
                     call cmpflx(flux(i,j,1),aleft,aright,vel(i,j),
     &                    muval(1),muval(2),gstar) 
                  else if (l_kr_type .eq. 3) then
                     call cmpflx_vang(flux(i,j,1),aleft,aright,vel(i,j),
     &                    muval(1),muval(2),gstar,l_kr_coef) 
                  else 
                     call bl_abort('riemann solver for kr_type specified 
     &                              does not exist')
                  end if
               else
                  call cmpflx_upwind(flux(i,j,1),aleft,aright,vel(i,j),
     &                 muval,gstar,
     &                 l_kr_type,r_kr_type,l_kr_coef,r_kr_coef)
               end if

               flux(i,j,2) = vel(i,j)  - flux(i,j,1)

               if (corrector .eq. 1) then
                  vel(i,j) = flux(i,j,1)
               end if

               flux(i,j,1) = rhoval(1) * flux(i,j,1)
               flux(i,j,2) = rhoval(2) * flux(i,j,2)

            end do
         end do
      else
         call bl_abort('model not defined in GODUNOV_2D:riemann_solver');
         
      end if

      end

      subroutine riemann_solver_cpl(flux,vel,vlo,vhi,DIMS(v),
     &     kpedge,DIMS(kp),lambda,DIMS(lambda),pc,DIMS(pc),
     &     eigl,DIMS(eig),krcoef,DIMS(krcoef),nkrcoef,
     &     ix,iy,lo,hi,dx,
     &     grav,muval,rhoval,nscal,model,bc,corrector)

c ::: Riemann Solver for two phase flow.
c     Determine the flux values along edges.
c     eigv(1:2) is overwrite with the velocity of each phase 
c     for subsequent advection of tracers.

      implicit none

c ::: arguments
      integer DIMDEC(v)
      integer DIMDEC(lambda)
      integer DIMDEC(pc)
      integer DIMDEC(kp)
      integer DIMDEC(eig)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ix,iy
      integer nscal,model
      integer lo(SDIM),hi(SDIM)
      integer bc(SDIM,2)
      integer corrector
      REAL_T  flux(DIMV(v),nscal)
      REAL_T  vel(DIMV(v))
      REAL_T  vlo(DIMV(v),nscal), vhi(DIMV(v),nscal)
      REAL_T  kpedge(DIMV(kp))
      REAL_T  pc(DIMV(pc))
      REAL_T  lambda(DIMV(lambda))
      REAL_T  eigl(DIMV(eig),nscal,nscal)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  rhoval(nscal),muval(nscal),grav
      REAL_T  dx(SDIM)

c ::: internal variables
      integer i,j,nc,ibc
      integer imin,jmin,imax,jmax
      integer l_kr_type, r_kr_type
      REAL_T  l_kr_coef, l_sat_res, r_kr_coef, r_sat_res
      REAL_T  aright, aleft, stmp(2), rtmp(2), st(2), gstar, lam(2),dpc
      logical no_change
      ibc = ix+2*iy
      imin = lo(1) - 1 + ix
      imax = hi(1) + 1
      jmin = lo(2) - 1 + iy
      jmax = hi(2) + 1

      if (model .eq. 2) then
         rtmp(:) = one
         do j = jmin,jmax
            do i = imin,imax

               r_kr_type = nint(krcoef(i,j,1))
               l_kr_type = nint(krcoef(i-ix,j-iy,1))  
               r_kr_coef = krcoef(i,j,2)
               l_kr_coef = krcoef(i-ix,j-iy,2)
               r_sat_res = krcoef(i,j,3)
               l_sat_res = krcoef(i-ix,j-iy,3)
               no_change = .true.
               if (r_kr_type .ne. l_kr_type .or. r_kr_coef .ne. l_kr_coef .or.
     &             r_sat_res .ne. l_sat_res) then
                  no_change = .false.
               end if
               l_kr_type = nint(krcoef(i,j,1))
               l_kr_coef = krcoef(i,j,2)

               aright = zero
               aleft  = zero
               do nc = 1,nscal
                  aright = aright + eigl(i,j,2,nc)*vhi(i,j,nc)
                  aleft  = aleft  + eigl(i-ix,j-iy,2,nc)*vlo(i,j,nc)
               end do 

               stmp(1) = aright
               stmp(2) = one-stmp(1)
               call get_st(st,stmp,rtmp,r_sat_res,nscal)
               aright = st(1)

               stmp(1) = aleft
               stmp(2) = one-stmp(1) 
               call get_st(st,stmp,rtmp,l_sat_res,nscal)
               aleft = st(1)

               flux(i,j,:) = zero
               dpc = (pc(i,j)-pc(i-ix,j-iy))/dx(ibc)

               if (ix .eq. 1) then
                  if (i .eq. imin .and. bc(1,1) .ne. INT_DIR) then
                     dpc = dpc*two
                  else if (i .eq. imax .and. bc(1,2) .ne. INT_DIR) then
                     dpc = dpc*two
                  end if
               else if (iy .eq. 1) then
                  if (j .eq. jmin .and. bc(2,1) .ne. INT_DIR) then
                     dpc = dpc*two
                  elseif (j.eq. jmax .and. bc(2,2) .ne. INT_DIR) then
                     dpc = dpc*two
                  end if
               end if

               gstar = kpedge(i,j)*(grav*(rhoval(1)-rhoval(2)) + dpc)
               if (no_change) then
                  if (l_kr_type .eq. 2) then
                     call cmpflx_arg(flux(i,j,1),stmp(1),aleft,aright,vel(i,j),
     &                    muval(1),muval(2),gstar)
                  else if (l_kr_type .eq. 3) then
                     call cmpflx_vang_arg(flux(i,j,1),stmp(1),aleft,aright,vel(i,j),
     &                    muval(1),muval(2),gstar,l_kr_coef) 
                  else 
                     call bl_abort('riemann solver for kr_type specified 
     &                              does not exist')
                  end if
                  stmp(2) = one - stmp(1)
               else
                  call cmpflx_upwind_arg(flux(i,j,1),stmp,aleft,aright,vel(i,j),muval,gstar,
     &                 l_kr_type,r_kr_type,l_kr_coef,r_kr_coef)
               end if
               flux(i,j,2) = vel(i,j)  - flux(i,j,1)
               if (corrector .eq. 1) then
                  if (i .ge. ARG_L1(lambda) .and. i .le. ARG_H1(lambda) .and.
     &                j .gt. ARG_L2(lambda) .and. j .lt. ARG_H2(lambda)) then
                     call get_ld(lam,stmp,rtmp,muval,nscal,l_kr_type,
     &                 l_kr_coef,zero)
                     lambda(i,j) = kpedge(i,j)*lam(1)*lam(2)/(lam(1)+lam(2))
                  end if
                  vel(i,j) = flux(i,j,1)
               end if

               flux(i,j,1) = rhoval(1) * flux(i,j,1)
               flux(i,j,2) = rhoval(2) * flux(i,j,2)
            end do
         end do
      else
         call bl_abort('model not defined in GODUNOV_2D:riemann_solver');
         
      end if

      end

      subroutine riemann_solver_arg(flux,vel,vlo,vhi,DIMS(v),
     &     kpedge,DIMS(kp),
     &     eigl,DIMS(eig),krcoef,DIMS(krcoef),nkrcoef,ix,iy,lo,hi,
     &     grav,muval,rhoval,nscal,model)
c      
c ::: Riemann Solver: finds the argument for the solution.
c
      implicit none

c ::: arguments
      integer DIMDEC(v)
      integer DIMDEC(kp)
      integer DIMDEC(eig)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ix,iy
      integer nscal,model
      integer lo(SDIM),hi(SDIM)
      REAL_T  flux(DIMV(v),nscal)
      REAL_T  vel(DIMV(v))
      REAL_T  vlo(DIMV(v),nscal), vhi(DIMV(v),nscal)
      REAL_T  kpedge(DIMV(kp))
      REAL_T  eigl(DIMV(eig),nscal,nscal)
      REAL_T krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  rhoval(nscal),muval(nscal),grav

c ::: internal variables
      integer i,j,nc    
      integer imin,jmin,imax,jmax
      integer l_kr_type, r_kr_type
      REAL_T  l_kr_coef, l_sat_res, r_kr_coef, r_sat_res
      REAL_T  aright, aleft, stmp(2), rtmp(2), st(2), gstar, ftmp
      logical no_change

      imin = lo(1) - 1 + ix
      imax = hi(1) + 1
      jmin = lo(2) - 1 + iy
      jmax = hi(2) + 1

      if (model .eq. 2) then
         rtmp(:) = one
         do j = jmin,jmax
            do i = imin,imax
               r_kr_type = nint(krcoef(i,j,1))
               l_kr_type = nint(krcoef(i-ix,j-iy,1))  
               r_kr_coef = krcoef(i,j,2)
               l_kr_coef = krcoef(i-ix,j-iy,2)
               r_sat_res = krcoef(i,j,3)
               l_sat_res = krcoef(i-ix,j-iy,3)
               no_change = .true.
               if (r_kr_type .ne. l_kr_type .or. r_kr_coef .ne. l_kr_coef .or.
     &             r_sat_res .ne. l_sat_res) then
                  no_change = .false.
               end if

               aright = 0
               aleft  = 0
               do nc = 1,nscal
                  aright = aright + eigl(i,j,2,nc)*vhi(i,j,nc)
                  aleft  = aleft  + eigl(i-ix,j-iy,2,nc)*vlo(i,j,nc)
               end do 

               stmp(1) = aright
               stmp(2) = 1-stmp(1)
               call get_st(st,stmp,rtmp,r_sat_res,nscal)
               aright = st(1)

               stmp(1) = aleft
               stmp(2) = 1-stmp(1)
               call get_st(st,stmp,rtmp,l_sat_res,nscal)
               aleft = st(1)

               gstar = kpedge(i,j)*grav*(rhoval(1)-rhoval(2))
               if (no_change) then
                  if (l_kr_type .eq. 2) then
                     call cmpflx_arg(ftmp,flux(i,j,1),aleft,aright,vel(i,j),
     &                    muval(1),muval(2),gstar)  
                  else if (l_kr_type .eq. 3) then
                     call cmpflx_vang_arg(ftmp,flux(i,j,1),aleft,aright,vel(i,j),
     &                    muval(1),muval(2),gstar,l_kr_coef)  
                  else
                     call bl_abort('riemann solver for the kr_type specified
     &                              does not exist')
                  end if
                  flux(i,j,2) = one - flux(i,j,1)
               else
                   call cmpflx_upwind(ftmp,flux(i,j,:),aleft,aright,vel(i,j),muval,gstar,
     &                 l_kr_type,r_kr_type,l_kr_coef,r_kr_coef)
                end if
             end do
         end do
      else
         call bl_abort('model not defined in GODUNOV_2D:riemann_solver');
         
      end if

      end

      subroutine riemann_est_eig(s,DIMS(s),lbd,dlbd,d2lbd,DIMS(lbd),
     &     vel,DIMS(vel),kpedge,DIMS(kp),rphi,DIMS(rphi),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     ix,iy,lo,hi,ncomps,grav,rhoval,
     &     muval,dx,bc,eigmax)
c        
c ::: Estimate maximum eigenvalues
c     Two cases:
c     1. If the kr_type changes, we use upwind values.
c     2. Otherwise, there are possibly two extrema for the eigenvalues.  
c        In those cases, we determine the minimum of the flux function,
c        and find the extrema in both directions.
c     The values are given by (dv1/ds)/\phi
c
      implicit none

c ::: Arguments
      integer DIMDEC(s)
      integer DIMDEC(lbd)
      integer DIMDEC(vel)
      integer DIMDEC(kp)
      integer DIMDEC(rphi)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ix,iy
      integer ncomps
      integer lo(2),hi(2)
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
      REAL_T dx(SDIM)
      REAL_T grav
      REAL_T eigmax

c ::: Internal variables
      integer i,j
      integer imin,jmin,imax,jmax
      integer l_kr_type, r_kr_type
      REAL_T  l_kr_coef, l_sat_res, r_kr_coef, r_sat_res
      REAL_T  eigtmp,krhog
      logical no_change

c ::: The indices are edge-based with no grown cell.
      imin = lo(1) 
      imax = hi(1) + ix
      jmin = lo(2) 
      jmax = hi(2) + iy

c     For two-phase system, the values are given by (dv1/ds)/\phi
         do j = jmin,jmax
            do i = imin,imax
               r_kr_type = nint(krcoef(i,j,1))
               l_kr_type = nint(krcoef(i-ix,j-iy,1))  
               r_kr_coef = krcoef(i,j,2)
               l_kr_coef = krcoef(i-ix,j-iy,2)
               r_sat_res = krcoef(i,j,3)
               l_sat_res = krcoef(i-ix,j-iy,3)
               no_change = .true.
               if (r_kr_type .ne. l_kr_type .or. r_kr_coef .ne. l_kr_coef .or.
     &             r_sat_res .ne. l_sat_res) then
                  no_change = .false.
               end if
               krhog = kpedge(i,j)*((rhoval(1)-rhoval(2))*grav) 
               if (no_change) then
                  call find_max_eig(s(i-ix,j-iy,:),s(i,j,:),ncomps,
     &                 lbd(i-ix,j-iy,:),dlbd(i-ix,j-iy,:),d2lbd(i-ix,j-iy,:),
     &                 lbd(i,j,:),dlbd(i,j,:),d2lbd(i,j,:),
     &                 vel(i,j),krcoef(i,j,:),nkrcoef,
     &                 krhog,rhoval,muval,eigtmp)
               else     
                  call find_upwind_eig(s(i-ix,j-iy,:),s(i,j,:),ncomps,
     &                 lbd(i-ix,j-iy,:),dlbd(i-ix,j-iy,:),
     &                 lbd(i,j,:),dlbd(i,j,:),
     &                 vel(i,j),krcoef(i,j,:),nkrcoef,krhog,eigtmp)
               end if
               eigtmp = eigtmp/(half*(rphi(i,j)+rphi(i-ix,j-iy)))
               eigmax = max(eigmax,dabs(eigtmp))
            end do
         end do
      end 

      subroutine riemann_est_eig_cpl(s,DIMS(s),lbd,dlbd,d2lbd,DIMS(lbd),
     &     pc,DIMS(pc),vel,DIMS(vel),kpedge,DIMS(kp),rphi,DIMS(rphi),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     ix,iy,lo,hi,ncomps,grav,rhoval,
     &     muval,dx,bc,eigmax)
c        
c ::: Estimate maximum eigenvalues
c     There are possibly two extrema for the eigenvalues.  In some cases,
c     we must first determine the minimum of the flux function, and find
c     the extrema in both directions.
c
      implicit none

c ::: Arguments
      integer DIMDEC(s)
      integer DIMDEC(lbd)
      integer DIMDEC(pc)
      integer DIMDEC(vel)
      integer DIMDEC(kp)
      integer DIMDEC(rphi)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ix,iy,ibc
      integer ncomps
      integer lo(2),hi(2)
      integer bc(SDIM,2)
      REAL_T s(DIMV(s),ncomps)
      REAL_T lbd(DIMV(lbd),2)
      REAL_T dlbd(DIMV(lbd),3)
      REAL_T d2lbd(DIMV(lbd),3)
      REAL_T pc(DIMV(pc))
      REAL_T krcoef(DIMV(krcoef),nkrcoef)
      REAL_T vel(DIMV(vel))
      REAL_T kpedge(DIMV(kp))
      REAL_T rphi(DIMV(rphi))
      REAL_T rhoval(ncomps)
      REAL_T muval(ncomps)
      REAL_T dx(SDIM)
      REAL_T grav
      REAL_T eigmax

c ::: Internal variables
      integer i,j
      integer imin,jmin,imax,jmax
      integer l_kr_type, r_kr_type
      REAL_T  l_kr_coef, l_sat_res, r_kr_coef, r_sat_res
      REAL_T eigtmp
      REAL_T krhog
      REAL_T dpc
      logical no_change

c ::: The indices are edge-based with no grown cell.
      ibc  = ix + 2*iy
      imin = lo(1) 
      imax = hi(1) + ix
      jmin = lo(2) 
      jmax = hi(2) + iy

c     For two-phase system, the values are given by (dv1/ds)/\phi
      do j = jmin,jmax
         do i = imin,imax
            r_kr_type = nint(krcoef(i,j,1))
            l_kr_type = nint(krcoef(i-ix,j-iy,1))  
            r_kr_coef = krcoef(i,j,2)
            l_kr_coef = krcoef(i-ix,j-iy,2)
            r_sat_res = krcoef(i,j,3)
            l_sat_res = krcoef(i-ix,j-iy,3)
            no_change = .true.
            if (r_kr_type .ne. l_kr_type .or. r_kr_coef .ne. l_kr_coef .or.
     &           r_sat_res .ne. l_sat_res) then
               no_change = .false.
            end if
            dpc = (pc(i,j)-pc(i-ix,j-iy))/dx(ibc)
            if (ix .eq. 1) then
               if (i .eq. imin .and. bc(1,1) .ne. INT_DIR) then
                  dpc = dpc*two
               else if (i .eq. imax .and. bc(1,2) .ne. INT_DIR) then
                  dpc = dpc*two
               end if
            else if (iy .eq. 1) then
               if (j .eq. jmin .and. bc(2,1) .ne. INT_DIR) then
                  dpc = dpc*two
               elseif (j.eq. jmax .and. bc(2,2) .ne. INT_DIR) then
                  dpc = dpc*two
               end if
            end if
            krhog = kpedge(i,j)*((rhoval(1)-rhoval(2))*grav + dpc)  
            if (no_change) then
               call find_max_eig(s(i-ix,j-iy,:),s(i,j,:),ncomps,
     &              lbd(i-ix,j-iy,:),dlbd(i-ix,j-iy,:),d2lbd(i-ix,j-iy,:),
     &              lbd(i,j,:),dlbd(i,j,:),d2lbd(i,j,:),
     &              vel(i,j),krcoef(i,j,:),nkrcoef,
     &              krhog,rhoval,muval,eigtmp)
            else     
               call find_upwind_eig(s(i-ix,j-iy,:),s(i,j,:),ncomps,
     &              lbd(i-ix,j-iy,:),dlbd(i-ix,j-iy,:),
     &              lbd(i,j,:),dlbd(i,j,:),
     &              vel(i,j),krcoef(i,j,:),nkrcoef,krhog,eigtmp)
            end if
            eigtmp = eigtmp/(half*(rphi(i,j)+rphi(i-ix,j-iy)))
            eigmax = max(eigmax,dabs(eigtmp))
         end do
      end do
      
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
      integer i,j
      integer l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T sat(2),rtmp(2)

      if (model .eq. 2) then

         rtmp   = one
         do j = ARG_L2(s),ARG_H2(s)
            do i = ARG_L1(s),ARG_H1(s)
               l_kr_type = nint(krcoef(i,j,1))
               l_kr_coef = krcoef(i,j,2)
               l_sat_res = krcoef(i,j,3)

               call get_st(sat,s(i,j,:),rhoval,l_sat_res,ncomps)
               if (l_kr_type .eq. 3) then
                  sat(1) = min(one-1.d-6,sat(1))
                  sat(1) = max(1.d-6,sat(1))
                  sat(2) = one-sat(1)
               end if
               
               call get_ld(lam(i,j,:),sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
               call get_dld(dlam(i,j,:),sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
               call get_d2ld(d2lam(i,j,:),sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)

            end do
         end do
      end if

      end 

      subroutine riemann_solver2(flux,DIMS(flux),s,DIMS(s),vlo,vhi,DIMS(v),
     &     eigv,DIMS(eig),vel,DIMS(vel),kpedge,DIMS(kp),ix,iy,lo,hi,
     &     athick,grav,muval,rhoval,ncomps,nphases,dx,dt,bc)
      implicit none
c ::: Riemann Solver: wrapper for cmpflx_polymer
c        ncomps is the number of auxiliary variables.
c ::: arguments
      integer DIMDEC(flux)
      integer DIMDEC(s)
      integer DIMDEC(v)
      integer DIMDEC(eig)
      integer DIMDEC(vel)
      integer DIMDEC(kp)
      integer ix,iy,ncomps,nphases
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
      REAL_T dx(SDIM), dt, dx_local

c ::: internal variables
      integer i,j,nc    
      integer imin,jmin,imax,jmax
      if (ix .eq. 1) then
         imin = lo(1)
         imax = hi(1) + 1
         jmin = lo(2) - 1
         jmax = hi(2) + 1
         dx_local = dx(1)
      elseif (iy .eq. 1) then
         imin = lo(1) - 1
         imax = hi(1) + 1
         jmin = lo(2) 
         jmax = hi(2) + 1
         dx_local = dx(2)
      end if

      do j = jmin,jmax
         do i = imin,imax
            do nc = 1,ncomps
               tum(nc)   = s(i-ix,j-iy,nc)
               tup(nc)   = s(i,j,nc)
               tulft(nc) = vlo(i,j,nc)
               turgt(nc) = vhi(i,j,nc)
               teigl(nc) = eigv(i-ix,j-iy,nc)
               teigr(nc) = eigv(i,j,nc)
            end do
            call cmpflx_polymer(tum,tup,tulft,turgt,teigl,teigr,tflx,vel(i,j),
     &           grav, muval(1), muval(2), athick, rhoval(1),rhoval(2),
     &           kpedge(i,j),dx_local, dt)
            flux(i,j,1) = tflx(1)
            flux(i,j,2) = tflx(2)   
         end do
      end do

      if (iy .eq. 1) then
         if (bc(2,2).eq.FOEXTRAP) then
            do i = imin,imax
               flux(i,jmax+1,1) = 0
               flux(i,jmax+1,2) = 0
            end do
         end if
         if (bc(2,1).eq.FOEXTRAP) then
            do i = imin,imax
               flux(i,jmin,1) = 0
               flux(i,jmin,2) = 0
            end do
         end if

      end if

      end

      subroutine advect_tracers(flux,DIMS(flux),uedge,vlo,vhi,DIMS(uedge),
     &     ix,iy,lo,hi,nscal)
      implicit none

c ::: Calculate tracers' flux

c ::: arguments
      integer DIMDEC(uedge)
      integer DIMDEC(flux)
      integer ix,iy
      integer nscal
      integer lo(SDIM),hi(SDIM)
      REAL_T  flux(DIMV(flux),nscal)
      REAL_T  vlo(DIMV(uedge),nscal)
      REAL_T  vhi(DIMV(uedge),nscal)
      REAL_T  uedge(DIMV(uedge))

c ::: internal variables
      integer i,j,nc 

      do j = lo(2),hi(2)+iy
         do i = lo(1),hi(1)+ix
            do nc = 1,nscal
               flux(i,j,nc) = merge(vlo(i,j,nc),vhi(i,j,nc),uedge(i,j)>=zero)
               flux(i,j,nc) = merge(half*(vlo(i,j,nc)+vhi(i,j,nc)),
     &              flux(i,j,nc),abs(uedge(i,j)).lt. 1.d-16)
               flux(i,j,nc) = flux(i,j,nc) * uedge(i,j)
            end do
         end do
      end do         

      end


      subroutine flux_copy(fluxs,DIMS(fluxs),fluxu,DIMS(fluxu),
     &     lo,hi,ix,iy,ncomps)
      implicit none

c ::: arguments
      integer DIMDEC(fluxs)
      integer DIMDEC(fluxu)
      integer lo(SDIM),hi(SDIM),ix,iy,ncomps
      REAL_T  fluxs(DIMV(fluxs),ncomps)
      REAL_T  fluxu(DIMV(fluxu),ncomps)

c ::: internal variables
      integer i,j,nc    
      integer imin,jmin,imax,jmax
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)

      do j = jmin,jmax+iy
         do i = imin,imax+ix
            do nc = 1,ncomps
               fluxs(i,j,nc) = fluxu(i,j,nc)
            end do    
         end do
      end do

      end

      subroutine flux_convert(fluxs,DIMS(fluxs),fluxu,DIMS(fluxu),vel,DIMS(v),
     &     lo,hi,ix,iy,rho,compcount,compidx,ncomps,nphases,bc)
      implicit none
c ::: arguments
      integer DIMDEC(fluxs)
      integer DIMDEC(fluxu)
      integer DIMDEC(v)
      integer lo(SDIM),hi(SDIM),ix,iy,ncomps,nphases
      integer bc(SDIM,2)
      integer compcount(nphases)
      integer compidx(nphases,ncomps)
      REAL_T  fluxs(DIMV(fluxs),ncomps)
      REAL_T  fluxu(DIMV(fluxu),ncomps-1)
      REAL_T  vel(DIMV(v))
      REAL_T  rho(nphases)


c ::: internal variables
      integer i,j,nc    
      integer imin,jmin,imax,jmax
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)

      do j = jmin,jmax+iy
         do i = imin,imax+ix
            do nc = 1,nphases
               if (compcount(nc) .eq. 2) then
                  fluxs(i,j,compidx(nc,1)) = rho(1)*
     &                 (fluxu(i,j,1)-fluxu(i,j,2))
                  fluxs(i,j,compidx(nc,2)) = rho(1)*
     &                 fluxu(i,j,2)
               else if (compcount(nc) .eq. 1) then
                  fluxs(i,j,compidx(nc,1)) = rho(2)*
     &                 (vel(i,j)-fluxu(i,j,1))
               end if
            end do
         end do
      end do  

      if (iy .eq. 1) then
         if (bc(2,1).eq.FOEXTRAP) then
            do i = imin,imax
               do nc = 1,nphases
                  if (compcount(nc) .eq. 1) then
                     fluxs(i,jmin,compidx(nc,1)) = 0
                  end if
               end do               
            end do
         end if
         if (bc(2,2).eq.FOEXTRAP) then
            do i = imin,imax
               do nc = 1,nphases
                  if (compcount(nc) .eq. 1) then
                     fluxs(i,jmax+1,compidx(nc,1)) = 0
                  end if
               end do               
            end do
         end if
      end if
      
      end

      subroutine eigsol_mp(eigv,eigl,eigr,lbd,dlbd,DIMS(w),
     &     v,DIMS(v),rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,grav,mu,
     &     ix,iy,lo,hi,rhoval,ncomps,nphases,model,bc)
      
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
      integer ix,iy
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
      integer i,j
      integer imin,jmin,imax,jmax,ibc
      integer l_kr_type
      REAL_T  l_sat_res
      REAL_T  lam(nphases), dlam(3), ltot 
      REAL_T  gs1,gs2,vt1,vt2

      ibc  = ix + 2*iy 
      imin = lo(1)-1
      imax = hi(1)+1
      jmin = lo(2)-1
      jmax = hi(2)+1
      if (model .eq. 2) then 

         do i = imin,imax
            do j = jmin,jmax
               l_kr_type = nint(krcoef(i,j,1))
               l_sat_res = krcoef(i,j,3)
c     eigenvalues for the left (1 to ncomps) and 
c     right (ncomps to 2*ncomps) state.
               
               lam = lbd(i,j,:)
               dlam = dlbd(i,j,:)

               ltot = lam(1)+lam(2)

               gs1 = kappa(i,j)*(rhoval(1)-rhoval(2))*grav
               gs2 = kappa(i+ix,j+iy)*(rhoval(1)-rhoval(2))*grav

               vt1 = v(i,j)
               vt2 = v(i+ix,j+iy)

               if (ix .eq. 1) then
                  if (bc(ibc,1) .eq. EXT_DIR .and. i .eq. imin) then
                     vt1 = v(i+1,j)
                  else if  (bc(ibc,2) .eq. EXT_DIR .and. i .eq. imax) then
                     vt2 = v(i,j)
                  end if
               else
                  if (bc(ibc,1) .eq. EXT_DIR .and. j .eq. jmin) then
                    vt1 = v(i,j+1)
                  else if  (bc(ibc,2) .eq. EXT_DIR .and. j .eq. jmax) then
                     vt2 = v(i,j)
                  end if
               end if

               eigv(i,j,1) = zero
               eigv(i,j,3) = zero

               if (l_kr_type .eq. 2) then
                  eigv(i,j,2) = vt1*(dlam(1)/ltot+lam(1)*dlam(3))+
     &                 gs1*(dlam(2)/ltot+lam(1)*lam(2)*dlam(3))
                  eigv(i,j,2) = eigv(i,j,2)/rphi(i,j)/(one-l_sat_res)
                  eigv(i,j,4) = vt2*(dlam(1)/ltot+lam(1)*dlam(3))+
     &                 gs2*(dlam(2)/ltot+lam(1)*lam(2)*dlam(3))
                  eigv(i,j,4) = eigv(i,j,4)/rphi(i,j)/(one-l_sat_res)
               elseif (l_kr_type .eq. 3) then
                  eigv(i,j,2) = (dlam(1)*lam(2)*(vt1+gs1*lam(2))
     .                 + lam(1)*(-vt1 + gs1*lam(1))*dlam(2))/(ltot**2)
                  eigv(i,j,2) = eigv(i,j,2)/rphi(i,j)/(one-l_sat_res)
                  eigv(i,j,4) = (dlam(1)*lam(2)*(vt2+gs2*lam(2))
     .                 + lam(1)*(-vt2 + gs2*lam(1))*dlam(2))/(ltot**2)
                  eigv(i,j,4) = eigv(i,j,4)/rphi(i,j)/(one-l_sat_res)
               else
                  print *,'l_kr_type,krcoef(i,j,1) : ',l_kr_type,krcoef(i,j,1)
                  call bl_abort('unknown kr_type specified')
               end if
c     ::: left eigenvectors
               eigl(i,j,1,1) = 1/rhoval(1)
               eigl(i,j,1,2) = 1/rhoval(2)
               eigl(i,j,2,1) = 1/rhoval(1)
               eigl(i,j,2,2) = 0

c     ::: right eigenvectors
               eigr(i,j,1,1) =  0
               eigr(i,j,1,2) =  rhoval(1)
               eigr(i,j,2,1) =  rhoval(2)
               eigr(i,j,2,2) = -rhoval(2)

            end do
         end do

      else
         call bl_abort('model is not defined in GODUNOV_2D:eigsol')
      end if

      end


      subroutine eigsol_mp_cpl(eigv,eigl,eigr,lbd,dlbd,DIMS(w),
     &     pc,DIMS(pc),v,DIMS(v),rphi,DIMS(rphi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,grav,mu,
     &     ix,iy,lo,hi,rhoval,ncomps,nphases,model,dx,bc)
      
c    
c     eigensystem for multiphase models.
c     
      implicit none

      integer DIMDEC(w)
      integer DIMDEC(pc)
      integer DIMDEC(v)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer ncomps,nphases,model
      integer lo(SDIM),hi(SDIM)
      integer ix,iy
      integer bc(SDIM,2)
      REAL_T  eigv(DIMV(w),2*ncomps)
      REAL_T  eigl(DIMV(w),ncomps,ncomps)
      REAL_T  eigr(DIMV(w),ncomps,ncomps)
      REAL_T  pc(DIMV(pc))
      REAL_T  lbd(DIMV(w),ncomps)
      REAL_T  dlbd(DIMV(w),3)
      REAL_T  v(DIMV(v))
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  mu(nphases),grav
      REAL_T  rhoval(ncomps)
      REAL_T  dx(SDIM)

c ::: internal variables
      integer i,j
      integer imin,jmin,imax,jmax,ibc
      integer l_kr_type
      REAL_T  l_sat_res
      REAL_T  lam(nphases), dlam(3), ltot 
      REAL_T  gs1,gs2,vt1,vt2, dpc

      ibc  = ix + 2*iy 
      imin = lo(1)-1
      imax = hi(1)+1
      jmin = lo(2)-1
      jmax = hi(2)+1
      
      if (model .eq. 2) then 

         do i = imin,imax
            do j = jmin,jmax
               l_kr_type = nint(krcoef(i,j,1))
               l_sat_res = krcoef(i,j,3)
c     eigenvalues for the left (1 to ncomps) and 
c     right (ncomps to 2*ncomps) state.
               
               lam  = lbd(i,j,:)
               dlam = dlbd(i,j,:)

               ltot = lam(1)+lam(2)

               vt1 = v(i,j)
               vt2 = v(i+ix,j+iy)
               dpc = (pc(i+ix,j+iy)-pc(i,j))/dx(ibc)
               if (ix .eq. 1) then
                  if (i .le. lo(1) .and. bc(ibc,1) .ne. INT_DIR) then
                     dpc = dpc*two
                  elseif (i .ge. hi(1) .and. bc(ibc,2) .ne. INT_DIR) then
                     dpc = dpc*two
                  end if
                  if (i .eq. imin .and. bc(ibc,1) .eq. EXT_DIR) then 
                     vt1 = v(i+1,j)
                  elseif  (i .eq. imax .and. bc(ibc,2) .eq. EXT_DIR) then 
                     vt2 = v(i,j)
                  end if
               else
                  if (j .le. lo(2) .and. bc(ibc,1) .ne. INT_DIR) then
                     dpc = dpc*two
                  elseif (j .ge. hi(2) .and. bc(ibc,2) .ne. INT_DIR) then
                     dpc = dpc*two
                  end if
                  if (j .eq. jmin .and. bc(ibc,1) .eq. EXT_DIR) then 
                     vt1 = v(i,j+1)
                  elseif  (j .eq. jmax .and. bc(ibc,2) .eq. EXT_DIR) then 
                     vt2 = v(i,j)
                  end if 
               end if

               gs1 = kappa(i,j)*((rhoval(1)-rhoval(2))*grav + dpc)
               gs2 = kappa(i+ix,j+iy)*((rhoval(1)-rhoval(2))*grav + dpc)

               eigv(i,j,1) = zero
               eigv(i,j,3) = zero

               if (l_kr_type .eq. 2) then
                  eigv(i,j,2) = vt1*(dlam(1)/ltot+lam(1)*dlam(3))+
     &                 gs1*(dlam(2)/ltot+lam(1)*lam(2)*dlam(3))
                  eigv(i,j,2) = eigv(i,j,2)/rphi(i,j)/(one-l_sat_res)
                  eigv(i,j,4) = vt2*(dlam(1)/ltot+lam(1)*dlam(3))+
     &                 gs2*(dlam(2)/ltot+lam(1)*lam(2)*dlam(3))
                  eigv(i,j,4) = eigv(i,j,4)/rphi(i,j)/(one-l_sat_res)
               elseif (l_kr_type .eq. 3) then
                  eigv(i,j,2) = (dlam(1)*lam(2)*(vt1+gs1*lam(2))
     &                 + lam(1)*(-vt1 + gs1*lam(1))*dlam(2))/(ltot**2)
                  eigv(i,j,2) = eigv(i,j,2)/rphi(i,j)/(one-l_sat_res)
                  eigv(i,j,4) = (dlam(1)*lam(2)*(vt2+gs2*lam(2))
     &                 + lam(1)*(-vt2 + gs2*lam(1))*dlam(2))/(ltot**2)
                  eigv(i,j,4) = eigv(i,j,4)/rphi(i,j)/(one-l_sat_res)
               else
                  call bl_abort('eigensolver for the kr_type specified
     &                           does not exist')
               end if

c     ::: left eigenvectors
               eigl(i,j,1,1) = 1/rhoval(1)
               eigl(i,j,1,2) = 1/rhoval(2)
               eigl(i,j,2,1) = 1/rhoval(1)
               eigl(i,j,2,2) = 0

c     ::: right eigenvectors
               eigr(i,j,1,1) =  0
               eigr(i,j,1,2) =  rhoval(1)
               eigr(i,j,2,1) =  rhoval(2)
               eigr(i,j,2,2) = -rhoval(2)

            end do
         end do

      else
         call bl_abort('model is not defined in GODUNOV_2D:eigsol')
      end if

      end

      subroutine eigsol2(s,DIMS(s),v,DIMS(v),rphi,DIMS(rphi),
     &     kappa,DIMS(kappa),eigv,eigl,eigr,struc,ldef,DIMS(w),
     &     lo,hi,ix,iy,mup,athick,rhop,grav,ncomps,nphases,bc,eigmax)
      implicit none
c ::  eigensystem for the problem
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(v)
      integer DIMDEC(rphi)
      integer DIMDEC(kappa)
      integer DIMDEC(w)
      integer ncomps,nphases,lo(SDIM),hi(SDIM),ix,iy
      integer bc(SDIM,2)
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
      integer i,j,n,m
      integer imin,jmin,imax,jmax,ibc
      REAL_T  vt, gstar
      REAL_T  u(2)
      REAL_T  teigv(2), teigl(2,2), teigr(2,2), dl(2,2)

      ibc  = ix + 2*iy 
      imin = lo(1)-1
      jmin = lo(2)-1
      imax = hi(1)+1
      jmax = hi(2)+1
      eigmax = 0

      do i = imin,imax
         do j = jmin,jmax
            gstar = kappa(i,j)*(rhop(1)-rhop(2))*grav
            u(1) = s(i,j,1)
            u(2) = s(i,j,2)

            vt = half*(v(i+ix,j+iy) + v(i,j))
            if (ix .eq. 1) then
               if (bc(ibc,1) .eq. EXT_DIR .and. i .eq. imin) then
                  vt = v(i+ix,j+iy)
               else if  (bc(ibc,2) .eq. EXT_DIR .and. i .eq. imax) then
                  vt = v(i,j)
               end if
            else
               if (bc(ibc,1) .eq. EXT_DIR .and. j .eq. jmin) then
                  vt = v(i+ix,j+iy)
               else if  (bc(ibc,2) .eq. EXT_DIR .and. j .eq. jmax) then
                  vt = v(i,j)
               end if
            end if
            
           
            call eigen_polymer(vt,gstar,mup(1),mup(2),
     &           athick,u,teigl,teigr,teigv,dl,ldef(i,j))
            do n = 1,2
               eigv(i,j,n) = teigv(n)/rphi(i,j)
 
               do m = 1,2
                  eigl(i,j,n,m)  = teigl(n,m)
                  eigr(i,j,n,m)  = teigr(n,m)
                  struc(i,j,n,m) = dl(n,1)*teigr(1,m) + dl(n,2)*teigr(2,m)
               end do
            end do
         end do
      end do

c     ::: find maximum eigenvalues
      do i = lo(1),hi(1)
         do j = lo(2),hi(2)
            do n = 1,2
               eigmax = max(eigmax, abs(eigv(i,j,n)))
            end do 
         end do
      end do
      
      end

      subroutine dfluxdv(drmndv,drmndk,DIMS(drmn),s,DIMS(s),rphi,DIMS(rphi),
     &     kappa,DIMS(kappa),krcoef,DIMS(krcoef),nkrcoef,mu,rhoval,ncomps,model,grav)
      
c ::: Derivative Flux w.r.t. v_T.  This is an analytical function

c ::: arguments 

      implicit none
     
      integer DIMDEC(drmn)
      integer DIMDEC(s)
      integer DIMDEC(rphi)
      integer DIMDEC(krcoef)
      integer nkrcoef
      integer DIMDEC(kappa)
      integer ncomps,model
      REAL_T  drmndv(DIMV(drmn),ncomps)
      REAL_T  drmndk(DIMV(drmn),ncomps)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  rphi(DIMV(rphi))
      REAL_T  kappa(DIMV(kappa))     
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  mu(ncomps),rhoval(ncomps)
      REAL_T  grav

c ::: internal variables
      integer i,j
      integer l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T  kr(ncomps)
      REAL_T  ltot,gstar
      
      if (model .eq. 2) then
c     This depends on the model and not on kr_type
         gstar = (rhoval(1)-rhoval(2))*grav
         do i = ARG_L1(drmn),ARG_H1(drmn)
            do j = ARG_L2(drmn),ARG_H2(drmn)
               l_kr_type = nint(krcoef(i,j,1))
               l_kr_coef = krcoef(i,j,2)
               l_sat_res = krcoef(i,j,3)
               call get_kr(kr,s(i,j,:),rhoval,ncomps,l_kr_type,l_kr_coef,l_sat_res)
               ltot = kr(1)+kr(2)*mu(1)/mu(2)
               drmndv(i,j,1) = kr(1)/ltot/rphi(i,j)
               ltot = kr(1)*mu(2)/mu(1) + kr(2)
               drmndv(i,j,2) = kr(2)/ltot/rphi(i,j)
               drmndk(i,j,1) = kr(2)*drmndv(i,j,1)*gstar/mu(2)
               drmndk(i,j,2) =  -kr(1)*drmndv(i,j,2)*gstar/mu(1)
               drmndv(i,j,1) = rhoval(1)*drmndv(i,j,1)
               drmndv(i,j,2) = rhoval(2)*drmndv(i,j,2)
               drmndk(i,j,1) = rhoval(1)*drmndk(i,j,1)
               drmndk(i,j,2) = rhoval(2)*drmndk(i,j,2)
            end do
         end do
      else
         call bl_abort('model is not defined in GODUNOV_2D:dfluxdv')
      end if         

      end

      subroutine dfluxdv2(drmndv,drmndk,s,rphi,kappa,
     &     mu1,mu2, athick,rho1,rho2,grav)
      implicit none
c     ::: Derivative Flux w.r.t. v_T.  This is an analytical function
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
      lambda(2) = (1-s(1))**2/mu2
      ltot = lambda(1) + lambda(2)

      drmndv(1) = lambda(1)*lambda(2)/ltot/rphi
      drmndv(2) = coef*drmndv(1)
      drmndk(1) = lambda(1)*lambda(2)*gstar/ltot/rphi
      drmndk(2) = coef*drmndk(1)
        
      end

      

      subroutine get_aux(s,u,DIMS(s),compcount,compidx, 
     &     rhoval, ncomps, nphases)
      implicit none
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

      integer i,j
      REAL_T  stmp(ncomps), utmp(ncomps-1)
      
      do i = ARG_L1(s), ARG_H1(s)
         do j =  ARG_L2(s), ARG_H2(s)
            utmp = zero
            stmp = s(i,j,:)
            call get_aux_i(stmp,utmp,compcount,compidx,
     &           rhoval,ncomps,nphases)    
            u(i,j,:) = utmp
         end do
      end do

      end

      subroutine get_aux_n(s,u,DIMS(s),compcount,compidx,
     &     rhoval,bc,ncomps,nphases)
      implicit none
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

      integer i,j
      REAL_T  stmp(ncomps), utmp(ncomps-1)


      if (bc(1,1) .eq. EXT_DIR) then
         i = ARG_L1(s)
         do j =  ARG_L2(s), ARG_H2(s)
            utmp = zero
            stmp = s(i,j,:)
            call get_aux_i(stmp,utmp,compcount,compidx,
     &           rhoval,ncomps,nphases)    
            u(i,j,:) = utmp
         end do
      end if
      
      if (bc(1,2) .eq. EXT_DIR) then
         i = ARG_H1(s)
         do j =  ARG_L2(s), ARG_H2(s)
            utmp = zero
            stmp = s(i,j,:)
            call get_aux_i(stmp,utmp,compcount,compidx,
     &           rhoval,ncomps,nphases)    
            u(i,j,:) = utmp
         end do
      end if

      if (bc(2,1) .eq. EXT_DIR) then
         j = ARG_L2(s)
         do i =  ARG_L1(s), ARG_H1(s)
            utmp = zero
            stmp = s(i,j,:)
            call get_aux_i(stmp,utmp,compcount,compidx,
     &           rhoval,ncomps,nphases)    
            u(i,j,:) = utmp
         end do
      end if

      if (bc(2,2) .eq. EXT_DIR) then
         j = ARG_H2(s)
         do i =  ARG_L1(s), ARG_H1(s)
            utmp = zero
            stmp = s(i,j,:)
            call get_aux_i(stmp,utmp,compcount,compidx,
     &           rhoval,ncomps,nphases)    
            u(i,j,:) = utmp
         end do
      end if
      end

      subroutine FORT_GET_EIG(s,DIMS(s),eigv,DIMS(w),ncomps,idx) 
      implicit none
c ::: trace the cell centered states to edges.
c ::: arguments
      integer DIMDEC(s)
      integer DIMDEC(w)
      integer ncomps,idx
      REAL_T  s(DIMV(s))
      REAL_T  eigv(DIMV(w),2*ncomps)
      integer i,j

      do i = ARG_L1(w),ARG_H1(w)
         do j = ARG_L2(w),ARG_H2(w)
            s(i,j) = eigv(i,j,idx)
         end do
      end do
      
      end 
      
      subroutine FORT_SYNC_ADV_FORCING(
     &     sync,DIMS(sync),
     &     xflux,DIMS(xflux),
     &     ucor,DIMS(ucor),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vcor,DIMS(vcor),
     &     areay,DIMS(ay),
     &     vol,DIMS(vol),
     &     lo,hi,iconserv )
c
c     This subroutine computes the sync advective tendency
c     for a state variable
c     Note that unlike the advection term, we do not divide by rphi
c     because this will be used directly in the diffusion equation.
c     There might be instances where this should be divided by rphi
c     but at present we are fine without the rphi.
c
      implicit none
      integer i,j
      integer iconserv
      integer imin,jmin,imax,jmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(sync)
      integer DIMDEC(vol)
      integer DIMDEC(ucor)
      integer DIMDEC(vcor)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T sync(DIMV(sync))
      REAL_T vol(DIMV(vol))
      REAL_T ucor(DIMV(ucor))
      REAL_T vcor(DIMV(vcor))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
c     compute corrective fluxes from edge states 
c     and perform conservative update
c
      do j = jmin,jmax
         do i = imin,imax+1
            xflux(i,j) = xflux(i,j)*ucor(i,j)*areax(i,j)
         end do
      end do
      do j = jmin,jmax+1
         do i = imin,imax
            yflux(i,j) = yflux(i,j)*vcor(i,j)*areay(i,j)
         end do
      end do

      do j = jmin,jmax
         do i = imin,imax
            sync(i,j) = sync(i,j) + (
     &           xflux(i+1,j) - xflux(i,j) +
     &           yflux(i,j+1) - yflux(i,j) )/vol(i,j)
         end do
      end do

      end

      subroutine FORT_SYNC_RMN_FORCING(
     &     sync,DIMS(sync),
     &     xflux,DIMS(xflux),
     &     ucor,DIMS(ucor),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vcor,DIMS(vcor),
     &     areay,DIMS(ay),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     vol,DIMS(vol),
     &     lo,hi,nscal)

      implicit none

#include "probdata.H"      
    
      integer nscal
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(sync)
      integer DIMDEC(vol)
      integer DIMDEC(ucor)
      integer DIMDEC(vcor)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(krcoef)
      integer nkrcoef
      REAL_T sync(DIMV(sync),nscal)
      REAL_T vol(DIMV(vol))
      REAL_T ucor(DIMV(ucor))
      REAL_T vcor(DIMV(vcor))
      REAL_T xflux(DIMV(xflux),nscal)
      REAL_T yflux(DIMV(yflux),nscal)
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
      REAL_T krcoef(DIMV(krcoef),nkrcoef)

      integer i,j,nc
      integer imin,jmin,imax,jmax
      integer l_kr_type
      REAL_T  l_kr_coef
      REAL_T  lam(2),rtmp(2)

      rtmp = one
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
c     compute corrective fluxes from edge states 
c     and perform conservative update
c
      if (model .eq. 2 ) then
         do j = jmin,jmax
            do i = imin,imax+1
               l_kr_type = nint(krcoef(i,j,1))
               l_kr_coef = krcoef(i,j,2)
               call get_ld(lam,xflux(i,j,:),rtmp,muval,ncomps,
     &              l_kr_type,l_kr_coef,zero)
               xflux(i,j,1) = areax(i,j)*rhoval(1)*lam(1)/(lam(1)+lam(2))*
     &              ucor(i,j)
               xflux(i,j,2) = areax(i,j)*rhoval(2)*lam(2)/(lam(1)+lam(2))*
     &              ucor(i,j)
               
            end do
         end do
         do j = jmin,jmax+1
            do i = imin,imax
               l_kr_type = nint(krcoef(i,j,1))
               l_kr_coef = krcoef(i,j,2)
               call get_ld(lam,yflux(i,j,:),rtmp,muval,ncomps,
     &              l_kr_type,l_kr_coef,zero)
               yflux(i,j,1) = areay(i,j)*rhoval(1)*lam(1)/(lam(1)+lam(2))*
     &              vcor(i,j)
               yflux(i,j,2) = areay(i,j)*rhoval(2)*lam(2)/(lam(1)+lam(2))*
     &              vcor(i,j) 
            end do
         end do
      else
         call bl_abort("SYNC_RMN_FORCING: model not defined.")
      end if

      do nc = 1,nscal
         do j = jmin,jmax
            do i = imin,imax
               sync(i,j,nc) = sync(i,j,nc) + (
     &              xflux(i+1,j,nc) - xflux(i,j,nc) +
     &              yflux(i,j+1,nc) - yflux(i,j,nc) )/vol(i,j)
               
            end do
         end do
      end do

      end


      subroutine trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(xx),uad,DIMS(uad),
     &     lo,hi,n,xbc,eps)
      implicit none
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the x direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      integer DIMDEC(s), DIMDEC(xx), DIMDEC(uad)
      integer n
      integer lo(SDIM), hi(SDIM)
      integer xbc(SDIM,2)
      REAL_T stx
      REAL_T s(DIMV(s))
      REAL_T xlo(DIMV(xx))
      REAL_T xhi(DIMV(xx))
      REAL_T uad(DIMV(uad))
      REAL_T eps
      
      logical ltest
      integer j
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
      
      if (xbc(1,1).eq.EXT_DIR) then
         if (n .eq. XVEL) then
            do j = jmin-1,jmax+1
             if (uad(imin,j) .ge. zero) then
                xlo(imin,j) = s(imin-1,j)
                xlo(imin,j) = s(imin-1,j)
              else 
                 xlo(imin,j) = xhi(imin,j)
              end if
            end do
         else
           do j = jmin-1,jmax+1
               ltest = uad(imin,j).le.eps
               stx   = merge(xhi(imin,j),s(imin-1,j),ltest)
               xhi(imin,j) = stx
               xlo(imin,j) = stx
            end do
         end if

      else if (xbc(1,1).eq.FOEXTRAP.or.xbc(1,1).eq.HOEXTRAP
     &        .or.xbc(1,1).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            xlo(imin,j) = xhi(imin,j)
         end do
      else if (xbc(1,1).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            xhi(imin,j) = zero
            xlo(imin,j) = zero
         end do
      end if

      if (xbc(1,2).eq.EXT_DIR) then
         if (n .eq. XVEL) then
            do j = jmin-1,jmax+1
              if (uad(imax+1,j) .le. zero) then
                 xhi(imax+1,j) = s(imax+1,j)
                 xlo(imax+1,j) = s(imax+1,j)
               else 
                 xhi(imax+1,j) = xlo(imax+1,j)
               end if
             end do
         else
            do j = jmin-1,jmax+1
               ltest = uad(imax+1,j).ge.-eps
               stx   = merge(xlo(imax+1,j),s(imax+1,j),ltest)
               xhi(imax+1,j) = stx
               xlo(imax+1,j) = stx
            end do
         end if
      else if (xbc(1,2).eq.FOEXTRAP.or.xbc(1,2).eq.HOEXTRAP
     &        .or.xbc(1,2).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            xhi(imax+1,j) = xlo(imax+1,j)
         end do
      else if (xbc(1,2).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            xhi(imax+1,j) = zero
            xlo(imax+1,j) = zero
         end do
      end if

      end

      subroutine trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(yy),vad,DIMS(vad),
     &     lo,hi,n,ybc,eps)
      implicit none
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the y direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      integer DIMDEC(s), DIMDEC(yy), DIMDEC(vad)
      integer n
      integer lo(SDIM), hi(SDIM)
      integer ybc(SDIM,2)
      REAL_T sty
      REAL_T s(DIMV(s))
      REAL_T ylo(DIMV(yy))
      REAL_T yhi(DIMV(yy))
      REAL_T vad(DIMV(vad))
      REAL_T eps

      logical ltest
      integer i
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
         
      if (ybc(2,1).eq.EXT_DIR) then
         if (n .eq. YVEL) then
            do i = imin-1,imax+1
              if (vad(i,jmin).ge.zero) then
                 yhi(i,jmin) = s(i,jmin-1)
                 ylo(i,jmin) = s(i,jmin-1)
              else
                 ylo(i,jmin) = yhi(i,jmin)
              end if
            end do
         else
            do i = imin-1,imax+1
               ltest = vad(i,jmin).le.eps
               sty   = merge(yhi(i,jmin),s(i,jmin-1),ltest)
               yhi(i,jmin) = sty
               ylo(i,jmin) = sty
            end do
         end if
      else if (ybc(2,1).eq.FOEXTRAP.or.ybc(2,1).eq.HOEXTRAP
     &        .or.ybc(2,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            ylo(i,jmin) = yhi(i,jmin)
         end do
      else if (ybc(2,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            yhi(i,jmin) = zero
            ylo(i,jmin) = zero
         end do
      end if

      if (ybc(2,2).eq.EXT_DIR) then
         if (n .eq. YVEL) then
            do i = imin-1,imax+1
              if (vad(i,jmax+1).le.zero) then
                 ylo(i,jmax+1) = s(i,jmax+1)
                 yhi(i,jmax+1) = s(i,jmax+1)
              else
                 yhi(i,jmax+1) = ylo(i,jmax+1)
              end if
            end do
         else
            do i = imin-1,imax+1
               ltest = vad(i,jmax+1).ge.-eps
               sty   = merge(ylo(i,jmax+1),s(i,jmax+1),ltest)
               yhi(i,jmax+1) = sty
               ylo(i,jmax+1) = sty
            end do
         end if
      else if (ybc(2,2).eq.FOEXTRAP.or.ybc(2,2).eq.HOEXTRAP
     &        .or.ybc(2,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            yhi(i,jmax+1) = ylo(i,jmax+1)
         end do
      else if (ybc(2,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            ylo(i,jmax+1) = zero
            yhi(i,jmax+1) = zero
         end do
      end if

      end
      
      subroutine FORT_SLOPES (s,DIMS(s),slx,sly,DIMS(sl),
     &     lo,hi,slxscr,slyscr,bc)
c 
c     this subroutine computes first, second or forth order slopes of
c     a 2D scalar field.
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

      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      integer DIMDEC(s)
      integer DIMDEC(sl)
      REAL_T     s(DIMV(s))
      REAL_T   slx(DIMV(sl))
      REAL_T   sly(DIMV(sl))
      REAL_T slxscr(DIM1(s), 4)
      REAL_T slyscr(DIM2(s), 4)

      integer cen,lim,flag,fromm
      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      integer imin,jmin,imax,jmax,i,j
      integer ng
      REAL_T dpls,dmin,ds
      REAL_T del,slim,sflg

      ng = lo(1) - ARG_L1(s)

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
c ::: ::::: added to prevent underflow for small s values
c
      do j = lo(2)-ng, hi(2)+ng
        do i = lo(1)-ng, hi(1)+ng
	  s(i,j) = merge(s(i,j), zero, abs(s(i,j)).gt.1.0D-20)
        end do
      end do
c
c     COMPUTE 0TH order slopes
c
      if (slope_order.eq.1) then
        if (ng .lt. 1) then
	   call bl_abort("FORT_SLOPES: not enough bndry cells for 1st order")
        end if
        do j = jmin-1,jmax+1
           do i = imin-1,imax+1
              slx(i,j) = zero
              sly(i,j) = zero
  	   end do
        end do
      end if
c
c     COMPUTE 2nd order slopes
c
      if (slope_order.eq.2) then
         if (ng .lt. 2) then
            call bl_abort('SLOPE_2D: not enough bndry cells for 
     .                    2nd order')
         end if
c
c     ------------------------ x slopes
c
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  slx(i,j) = half*(s(i+1,j)-s(i-1,j))
               end do
            end do
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imin-1,j) = zero
                  slx(imin,j)   = (s(imin+1,j)+three*s(imin,j)-four*s(imin-1,j))/three
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imax+1,j) = zero
                  slx(imax,j)   = -(s(imax-1,j)+three*s(imax,j)-four*s(imax+1,j))/three
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i+1,j)-s(i-1,j))
                  dpls = two*(s(i+1,j) - s(i ,j))
                  dmin = two*(s(i ,j) - s(i-1,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  slx(i,j)= sflg*min(slim,abs(del))
               end do
            end do
            
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imin-1,j) = zero
                  del  = (s(imin+1,j)+three*s(imin,j)-four*s(imin-1,j))/three
                  dpls = two*(s(imin+1,j) - s(imin  ,j))
                  dmin = two*(s(imin  ,j) - s(imin-1,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  slx(imin,j)= sflg*min(slim,abs(del))
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imax+1,j) = zero
                  del  = -(s(imax-1,j)+three*s(imax,j)-four*s(imax+1,j))/three
                  dpls = two*(s(imax+1,j) - s(imax  ,j))
                  dmin = two*(s(imax  ,j) - s(imax-1,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  slx(imax,j)= sflg*min(slim,abs(del))
               end do
            end if
         end if
c
c     ------------------------ y slopes
c
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  sly(i,j) = half*(s(i,j+1)-s(i,j-1))
               end do
            end do
            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmin-1) = zero
                  sly(i,jmin) = (s(i,jmin+1)+three*s(i,jmin)-four*s(i,jmin-1))/three
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmax+1) = zero
                  sly(i,jmax) = -(s(i,jmax-1)+three*s(i,jmax)-four*s(i,jmax+1))/three
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i,j+1)-s(i,j-1))
                  dpls = two*(s(i,j+1) - s(i,j ))
                  dmin = two*(s(i,j ) - s(i,j-1))
                  slim = min(abs(dpls),abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  sly(i,j)= sflg*min(slim,abs(del))
               end do
            end do

            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmin-1) = zero
                  del  = (s(i,jmin+1)+three*s(i,jmin)-four*s(i,jmin-1))/three
                  dpls = two*(s(i,jmin+1) - s(i,jmin  ))
                  dmin = two*(s(i,jmin  ) - s(i,jmin-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  sly(i,jmin)= sflg*min(slim,abs(del))
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmax+1) = zero
                  del  = -(s(i,jmax-1)+three*s(i,jmax)-four*s(i,jmax+1))/three
                  dpls = two*(s(i,jmax+1) - s(i,jmax  ))
                  dmin = two*(s(i,jmax  ) - s(i,jmax-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  sly(i,jmax)= sflg*min(slim,abs(del))
               end do
            end if
         end if
c
c ... end, if slope_order .eq. 2
c
      end if
c
c     COMPUTE 4TH order slopes
c
      if (slope_order.eq.4) then
         if (ng .lt. 3) then
            call bl_abort("SLOPE_2D: not enough bndry cells for 4th order")
         end if
c         
c     ------------------------ x slopes
c
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
               do i = imin-2,imax+2
                  slxscr(i,cen)  = half*(s(i+1,j)-s(i-1,j))
               end do
               do i = imin-1,imax+1
                  slx(i,j) = two * two3rd * slxscr(i,cen) -
     &                 sixth * (slxscr(i+1,cen) + slxscr(i-1,cen))
               end do
            end do
            
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imin,j) = -sixteen/fifteen*s(imin-1,j) + half*s(imin,j) + 
     &                 two3rd*s(imin+1,j) - tenth*s(imin+2,j)
                  slx(imin-1,j) = zero
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imax,j) = -( -sixteen/fifteen*s(imax+1,j) + half*s(imax,j) + 
     &                 two3rd*s(imax-1,j) - tenth*s(imax-2,j) )
                  slx(imax+1,j) = zero
               end do
            end if
         else        
            do j = jmin-1,jmax+1
               do i = imin-2,imax+2
                  dmin           = two*(s(i,  j)-s(i-1,j))
                  dpls           = two*(s(i+1,j)-s(i  ,j))
                  slxscr(i,cen)  = half*(s(i+1,j)-s(i-1,j))
                  slxscr(i,lim)  = min(abs(dmin),abs(dpls))
                  slxscr(i,lim)  = merge(slxscr(i,lim),zero,dpls*dmin>=zero)
                  slxscr(i,flag) = sign(one,slxscr(i,cen))
                  slxscr(i,fromm)= slxscr(i,flag)*
     &                 min(slxscr(i,lim),abs(slxscr(i,cen)))
               end do

c              At boundaries, we use a less steep gradient (1 instead of 2)

               do i = imin-1,imax+1
                  ds = two * two3rd * slxscr(i,cen) -
     &                 sixth * (slxscr(i+1,fromm) + slxscr(i-1,fromm))
                  slx(i,j) = slxscr(i,flag)*min(abs(ds),slxscr(i,lim))
               end do

               if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                  del  = -sixteen/fifteen*s(imin-1,j) + half*s(imin,j) + 
     &                 two3rd*s(imin+1,j) - tenth*s(imin+2,j)
                  dmin = (s(imin  ,j)-s(imin-1,j))
                  dpls = (s(imin+1,j)-s(imin  ,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  slx(imin-1,j) = zero
                  slx(imin,  j) = sflg*min(slim,abs(del))

c                 Recalculate the slope at imin+1 using the revised slxscr(imin,fromm)
                  slxscr(imin,fromm) = slx(imin,j)
                  ds = two * two3rd * slxscr(imin+1,cen) -
     $               sixth * (slxscr(imin+2,fromm) + slxscr(imin,fromm))
                  slx(imin+1,j) = slxscr(imin+1,flag)*min(abs(ds),slxscr(imin+1,lim))
               end if

               if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                  del  = -( -sixteen/fifteen*s(imax+1,j) + half*s(imax,j) + 
     &                 two3rd*s(imax-1,j) - tenth*s(imax-2,j) )
                  dmin = (s(imax  ,j)-s(imax-1,j))
                  dpls = (s(imax+1,j)-s(imax  ,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  slx(imax,  j) = sflg*min(slim,abs(del))
                  slx(imax+1,j) = zero

c                 Recalculate the slope at imax-1 using the revised slxscr(imax,fromm)
                  slxscr(imax,fromm) = slx(imax,j)
                  ds = two * two3rd * slxscr(imax-1,cen) -
     $               sixth * (slxscr(imax-2,fromm) + slxscr(imax,fromm))
                  slx(imax-1,j) = slxscr(imax-1,flag)*min(abs(ds),slxscr(imax-1,lim))
               end if
            end do
         end if
c
c     ------------------------ y slopes
c
         if (use_unlimited_slopes) then
            do i = imin-1,imax+1
               do j = jmin-2,jmax+2
                  slyscr(j,cen)  = half*(s(i,j+1)-s(i,j-1))
               end do
               do j = jmin-1,jmax+1
                  sly(i,j) = two * two3rd * slyscr(j,cen) -
     &                 sixth * (slyscr(j+1,cen) + slyscr(j-1,cen))
               end do
            end do
            
            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmin-1) = zero
                  sly(i,jmin) = -sixteen/fifteen*s(i,jmin-1) + half*s(i,jmin) + 
     &                 two3rd*s(i,jmin+1) - tenth*s(i,jmin+2)
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmax) = -( -sixteen/fifteen*s(i,jmax+1) + half*s(i,jmax) + 
     &                 two3rd*s(i,jmax-1) - tenth*s(i,jmax-2) )
                  sly(i,jmax+1) = zero
               end do
            end if
         else
            do i = imin-1,imax+1
              do j = jmin-2,jmax+2
                  dmin           = two* (s(i,j  )-s(i,j-1))
                  dpls           = two* (s(i,j+1)-s(i,j  ))
                  slyscr(j,cen)  = half*(s(i,j+1)-s(i,j-1))
                  slyscr(j,lim)  = min(abs(dmin),abs(dpls))
                  slyscr(j,lim)  = merge(slyscr(j,lim),zero,dpls*dmin>=zero)
                  slyscr(j,flag) = sign(one,slyscr(j,cen))
                  slyscr(j,fromm)= slyscr(j,flag)*
     &                 min(slyscr(j,lim),abs(slyscr(j,cen)))
               end do
               do j = jmin-1,jmax+1
                  ds = two * two3rd * slyscr(j,cen) -
     &                 sixth * (slyscr(j+1,fromm) + slyscr(j-1,fromm))
                  sly(i,j) = slyscr(j,flag)*min(abs(ds),slyscr(j,lim))
               end do

c              At boundaries, we use a less steep gradient (1 instead of 2)

               if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                  del  = -sixteen/fifteen*s(i,jmin-1) + half*s(i,jmin) + 
     &                 two3rd*s(i,jmin+1) - tenth*s(i,jmin+2)
                  dmin = (s(i,jmin  )-s(i,jmin-1))
                  dpls = (s(i,jmin+1)-s(i,jmin  ))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  sly(i,jmin-1) = zero
                  sly(i,jmin  ) = sflg*min(slim,abs(del))

c                 Recalculate the slope at jmin+1 using the revised slyscr(jmin,fromm)
                  slyscr(jmin,fromm) = sly(i,jmin)
                  ds = two * two3rd * slyscr(jmin+1,cen) -
     $               sixth * (slyscr(jmin+2,fromm) + slyscr(jmin,fromm))
                  sly(i,jmin+1) = slyscr(jmin+1,flag)*min(abs(ds),slyscr(jmin+1,lim))
               end if

               if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                  del  = -( -sixteen/fifteen*s(i,jmax+1) + half*s(i,jmax) + 
     &                 two3rd*s(i,jmax-1) - tenth*s(i,jmax-2) )
                  dmin = (s(i,jmax  )-s(i,jmax-1))
                  dpls = (s(i,jmax+1)-s(i,jmax  ))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, dpls*dmin>=zero)
                  sflg = sign(one,del)
                  sly(i,jmax  ) = sflg*min(slim,abs(del))
                  sly(i,jmax+1) = zero

c                 Recalculate the slope at jmax-1 using the revised slyscr(jmax,fromm)
                  slyscr(jmax,fromm) = sly(i,jmax)
                  ds = two * two3rd * slyscr(jmax-1,cen) -
     $               sixth * (slyscr(jmax-2,fromm) + slyscr(jmax,fromm))
                  sly(i,jmax-1) = slyscr(jmax-1,flag)*min(abs(ds),slyscr(jmax-1,lim))
               end if
           end do
         end if
c
c ... end, if slope_order .eq. 4
c
      end if

      end

      subroutine FORT_CONVSCALMINMAX(s,DIMS(s),sn,DIMS(sn),
     &                               lo,hi,bc)
      implicit none
c
c     correct an convectively-advected field for under/over shoots
c
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      REAL_T s(DIMV(s))
      REAL_T sn(DIMV(sn))
      integer  i, j, imin, imax, jmin, jmax
      REAL_T   smin, smax

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
c
c     set corner values for s
c
c     :: corner ilo, jlo
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
         s(imin-1,jmin-1) = s(imin,jmin)
      end if
c
c     :: corner ihi, jlo
c
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
         s(imax+1,jmin-1) = s(imax,jmin)
      end if
c
c     :: corner ilo, jhi
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
         s(imin-1,jmax+1) = s(imin,jmax)
      end if
c
c     :: corner ihi, jhi
c
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
         s(imax+1,jmax+1) = s(imax,jmax)
      end if
c
c     compute extrema of s
c
      do j = jmin, jmax         
         do i = imin, imax
            smin = min(
     &           s(i-1,j-1),s(i  ,j-1),s(i+1,j-1),
     &           s(i-1,j  ),s(i  ,j  ),s(i+1,j  ),
     &           s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            smax = max(
     &           s(i-1,j-1),s(i  ,j-1),s(i+1,j-1),
     &           s(i-1,j  ),s(i  ,j  ),s(i+1,j  ),
     &           s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            sn(i,j) = max(sn(i,j),smin)
            sn(i,j) = min(sn(i,j),smax)
         end do
      end do
      end

      subroutine FORT_CONSSCALMINMAX(s,sn,DIMS(s),lo,hi,bc)
      implicit none
c
c     correct an conservatively-advected field for under/over shoots
c
      integer DIMDEC(s)
      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      REAL_T    s(DIMV(s))
      REAL_T   sn(DIMV(s))
      integer  i, j, imin, imax, jmin, jmax
      REAL_T   slo_min,smd_min,shi_min
      REAL_T   slo_max,smd_max,shi_max
      REAL_T   smin, smax

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
c
c     set corner values for s
c
c     :: corner ilo, jlo
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
           s(imin-1,jmin-1) =   s(imin,jmin)
      end if
c
c     :: corner ihi, jlo
c
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
           s(imax+1,jmin-1) =   s(imax,jmin)
      end if
c
c     :: corner ilo, jhi
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
           s(imin-1,jmax+1) =   s(imin,jmax)
      end if
c
c     :: corner ihi, jhi
c
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
           s(imax+1,jmax+1) =   s(imax,jmax)
      end if
c
c     compute extrema of s
c
      do j = jmin, jmax         
         do i = imin, imax
            slo_min = min(s(i-1,j-1),s(i  ,j-1),s(i+1,j-1))
            smd_min = min(s(i-1,j  ),s(i  ,j  ),s(i+1,j  ))
            shi_min = min(s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            smin = min(slo_min,smd_min,shi_min)

            slo_max = max(s(i-1,j-1),s(i  ,j-1),s(i+1,j-1))
            smd_max = max(s(i-1,j  ),s(i  ,j  ),s(i+1,j  ))
            shi_max = max(s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            smax = max(slo_max,smd_max,shi_max)

            sn(i,j) = max(sn(i,j),smin)
            sn(i,j) = min(sn(i,j),smax)
         end do
      end do

      end

      subroutine FORT_SUM_TF_DIVU(
     &     S,DIMS(S),
     &     tforces,DIMS(tf),
     &     divu,DIMS(divu),
     &     lo,hi,nvar,iconserv )

      implicit none

c
c     sum divU*S into tforces 
c     depending on the value of iconserv
c
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)

      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf))
      REAL_T divu(DIMV(divu))

      if ( iconserv .eq. 1 ) then
         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j) = tforces(i,j) - S(i,j,n)*divu(i,j)
               end do
            end do
         end do
      end if

      end

      subroutine FORT_SUM_TF_DIVU_VISC(
     &     S,DIMS(S),
     &     tforces,DIMS(tf),
     &     divu,DIMS(divu),
     &     visc,DIMS(visc),
     &     lo,hi,nvar,iconserv )

      implicit none
c
c     sum tforces, viscous forcing and divU*S into tforces
c     depending on the value of iconserv
c
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(visc)

      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T divu(DIMV(divu))
      REAL_T visc(DIMV(visc),nvar)
 
      if ( iconserv .eq. 1 ) then
         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = tforces(i,j,n)
     &                 + visc(i,j,n) - S(i,j,n)*divu(i,j)
               end do
            end do
         end do
      else
         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = tforces(i,j,n) + visc(i,j,n) 
               end do
            end do
         end do
      end if

      end

      subroutine FORT_UPDATE_AOFS_TF(
     &     s,       DIMS(s),
     &     sn,      DIMS(sn),
     &     aofs,    DIMS(aofs),
     &     tforces, DIMS(tf),
     &     rockphi, DIMS(rockphi),
     &     lo,hi,dt,nvar)
      implicit none
c
c     update a field with an advective tendency
c     and a forcing term
c

#include "probdata.H"

      integer i, j, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(aofs)
      integer DIMDEC(tf)
      integer DIMDEC(rockphi)
      integer lo(SDIM), hi(SDIM)
      REAL_T  dt
      REAL_T  s(DIMV(s),nvar)
      REAL_T  sn(DIMV(sn),nvar)
      REAL_T  aofs(DIMV(aofs),nvar)
      REAL_T  tforces(DIMV(tf),nvar)
      REAL_T  rockphi(DIMV(rockphi))

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 1,nvar
               sn(i,j,n) = s(i,j,n) - dt*aofs(i,j,n)/rockphi(i,j)
     &              + dt*tforces(i,j,n)
            end do
         end do
      end do

      end 

      subroutine FORT_UPDATE_AOFS_TRACER(
     &     co,      DIMS(co),
     &     cn,      DIMS(cn),
     &     so,      DIMS(so),
     &     sn,      DIMS(sn),
     &     aofs,    DIMS(aofs),
     &     tforces, DIMS(tf),
     &     rockphi, DIMS(rockphi),
     &     lo, hi,idx_total, idx_total_n, dt, nc, nt)

      implicit none
c
c     update tracer with an advective tendency
c     and a forcing term
c

#include "probdata.H"

      integer nc, nt
      integer DIMDEC(co)
      integer DIMDEC(cn)
      integer DIMDEC(so)
      integer DIMDEC(sn)
      integer DIMDEC(aofs)
      integer DIMDEC(tf)
      integer DIMDEC(rockphi)
      integer lo(SDIM), hi(SDIM)
      integer idx_total_n
      integer idx_total(idx_total_n)
      REAL_T  dt
      REAL_T  co(DIMV(co),nt)
      REAL_T  cn(DIMV(cn),nt)
      REAL_T  so(DIMV(so),nc)
      REAL_T  sn(DIMV(sn),nc)
      REAL_T  aofs(DIMV(aofs),nt)
      REAL_T  tforces(DIMV(tf),nt)
      REAL_T  rockphi(DIMV(rockphi))

      integer i, j, n, idx
      integer cWater

c     First update all tracers with no change
      cn(:,:,1:nt) = co(:,:,1:nt)

c     Then, advect/force those in the "Total" group
c     FIXME: Assume first component is water, and that all tracers of interest are in water
      cWater = 1

c     idx_total is an array of 0-based indices selecting a subset of the tracers to update

      do n = 1,idx_total_n
         idx = idx_total(n) + 1

         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (sn(i,j,cWater) .gt. zero) then
                  cn(i,j,idx) = (so(i,j,cWater)*co(i,j,idx)
     &                 - dt*(aofs(i,j,idx)+tforces(i,j,idx))/rockphi(i,j))/sn(i,j,cWater)
               end if
            end do
         end do
      end do

      end 

      subroutine slope_alpha(s,DIMS(s),slx,sly,DIMS(sl),
     &     eigvx,eiglx,eigrx,strucx, ldefx,
     &     eigvy,eigly,eigry,strucy, ldefy,
     &     DIMS(eig),
     &     slxscr,slyscr,lo,hi,bc,ncomps,nphases)
      implicit none
c ::: slope limiting in single direction
c ::: replaced sl with the limited alpha instead
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
      REAL_T  slxscr(DIM1(s),3,2)
      REAL_T  slyscr(DIM2(s),3,2)
      REAL_T  eigvx(DIMV(eig),ncomps), eigvy(DIMV(eig),ncomps)
      REAL_T  eiglx(DIMV(eig),ncomps,ncomps),  eigly(DIMV(eig),ncomps,ncomps)
      REAL_T  eigrx(DIMV(eig),ncomps,ncomps),  eigry(DIMV(eig),ncomps,ncomps)
      REAL_T  strucx(DIMV(eig),ncomps,ncomps), strucy(DIMV(eig),ncomps,ncomps)
      logical ldefx(DIMV(eig)), ldefy(DIMV(eig))
      integer ng,n,m,i,j
      integer imin,jmin,imax,jmax
      REAL_T alpl(2), alpc(2), alpr(2)
      REAL_T anslp, den, faclim, cdef, esnew, slrmin, def

      ng = lo(1) - ARG_L1(s)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      cdef = 1.0d-1
c
c     COMPUTE 0TH order slopes
c
      if (slope_order.eq.1) then
        if (ng .lt. 1) then
	   call bl_abort("FORT_SLOPES: not enough bndry cells for 1st order")
        end if
        do j = jmin-1,jmax+1
           do i = imin-1,imax+1
              slx(i,j,:) = zero
              sly(i,j,:) = zero
  	   end do
        end do
      end if
c
c     COMPUTE 0TH order slopes
c
      if (slope_order .eq. 4) then
         if (ng .lt. 3) then
            call bl_abort("SLOPE_2D: not enough bndry cells for 4th order")
         end if
c     ::: added to prevent underflow for small s values
         do n = 1,ncomps
            do j = jmin-ng, jmax+ng
               do i = imin-ng, imax+ng
                  s(i,j,n) = merge(s(i,j,n), zero, abs(s(i,j,n)).gt.1.0D-20)
               end do
            end do
         end do

         if (use_unlimited_slopes) then
            call bl_abort("cannot used unlimited slopes")
         else
c     :::::: x-slopes
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1 
                  do n = 1,ncomps
                     slxscr(i,1,n)  = s(i,j,n)   - s(i-1,j,n)
                     slxscr(i,3,n)  = s(i+1,j,n) - s(i,j,n)
                     slxscr(i,2,n)  = two*(s(i+1,j,n)-s(i-1,j,n))/3.d0
     &                    -one*(s(i+2,j,n)-s(i-2,j,n))/12.d0  
                  end do
               end do
               
               if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slxscr(imin,2,n) = -sixteen/fifteen*s(imin-1,j,n) + half*s(imin,j,n) + 
     &                    two3rd*s(imin+1,j,n) - tenth*s(imin+2,j,n)
                     slxscr(imin-1,2,n) = zero
                     slx(imin-1,j,n) = zero
                  end do
               end if
               if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slxscr(imax,2,n) = -( -sixteen/fifteen*s(imax+1,j,n) + half*s(imax,j,n) + 
     &                 two3rd*s(imax-1,j,n) - tenth*s(imax-2,j,n) )
                     slxscr(imax+1,2,n) = zero
                     slx(imax+1,j,n) = zero
                  end do
               end if
                                
               do i = imin-1,imax+1
                  do n = 1,ncomps
                     alpl(n) = zero
                     alpc(n) = zero
                     alpr(n) = zero                  
                     do m = 1,ncomps
                        alpl(n) = alpl(n) + eiglx(i,j,n,m)*slxscr(i,1,m)
                        alpc(n) = alpc(n) + eiglx(i,j,n,m)*slxscr(i,2,m)
                        alpr(n) = alpr(n) + eiglx(i,j,n,m)*slxscr(i,3,m)
                     end do
                  end do

                     if (dabs(alpl(2)) .lt. 1.d-5) then
                        alpl(2) = zero
                     end if
                     if (dabs(alpc(2)) .lt. 1.d-5) then
                        alpc(2) = zero
                     end if
                     if (dabs(alpr(2)) .lt. 1.d-5) then
                        alpr(2) = zero
                     end if
                  
                  do n = 1,ncomps
                     faclim = 2.
                     if (i .ge. imin .and. i .le. imax) then
                        if(strucx(i,j,n,n)*strucx(i+1,j,n,n).lt.-1.e-10)faclim=one
                        if(strucx(i,j,n,n)*strucx(i-1,j,n,n).lt.-1.e-10)faclim=one
                     end if
                     anslp = dmin1(faclim*dabs(alpl(n)),
     &                    dabs(alpc(n)),faclim*dabs(alpr(n)))
     &                    *dmax1(zero,dsign(one,alpl(n)*alpr(n)))
                     den=merge(dabs(alpc(n)),one,dabs(alpc(n)).gt.1.d-15)
                     slx(i,j,n)=anslp/den

                  end do

                  def=dabs(alpc(1))*cdef*dabs(strucx(i,j,1,1)-strucx(i,j,2,1))
     &                 +dabs(alpc(2))*cdef*dabs(strucx(i,j,1,2)-strucx(i,j,2,2))
     &                 -dabs(eigvx(i,j,1)-eigvx(i,j,2))

                  if(def.gt. zero .or. ldefx(i,j))then
                     if(eigvx(i,j,1)*eigvx(i,j,2).lt. zero)then
                        slx(i,j,1)=0.
                        slx(i,j,2)=0.
                     else
                        esnew=half*(eigvx(i,j,1)+eigvx(i,j,2))
                        slrmin=dmin1(slx(i,j,1),slx(i,j,2))
                        slx(i,j,1)=slrmin
                        slx(i,j,2)=slrmin
                        eigvx(i,j,1)=esnew
                        eigvx(i,j,2)=esnew
                     end if
                  end if
                 
                  slx(i,j,1) = alpc(1)*slx(i,j,1)
                  slx(i,j,2) = alpc(2)*slx(i,j,2)
               end do
            end do

c :::::: y-slopes
        
            do i = imin-1,imax+1 
               do j = jmin-1,jmax+1      
                  do n = 1,ncomps
                     
                     slyscr(j,1,n)  = s(i,j,n)   - s(i,j-1,n)
                     slyscr(j,3,n)  = s(i,j+1,n) - s(i,j,n)
                     slyscr(j,2,n)  = two*(s(i,j+1,n)-s(i,j-1,n))/3.d0
     &                    -one*(s(i,j+2,n)-s(i,j-2,n))/12.d0
                  end do
               end do
               
               if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slyscr(jmin,2,n) = -sixteen/fifteen*s(i,jmin-1,n) + half*s(i,jmin,n) + 
     &                    two3rd*s(i,jmin+1,n) - tenth*s(i,jmin+2,n)
                     slyscr(jmin-1,2,n) = zero
                     sly(i,jmin-1,n) = zero
                  end do
               end if
               if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                  do n = 1,ncomps
                     slyscr(jmax,2,n)   = -( -sixteen/fifteen*s(i,jmax+1,n) + half*s(i,jmax,n) + 
     &                    two3rd*s(i,jmax-1,n) - tenth*s(i,jmax-2,n) )
                     slyscr(jmax+1,2,n) = zero
                     sly(i,jmax+1,n)    = zero
                  end do
               end if
               
               do j = jmin-1,jmax+1 
                  do n = 1,ncomps
                     alpl(n) = zero
                     alpc(n) = zero
                     alpr(n) = zero                  
                     do m = 1,ncomps
                        alpl(n) = alpl(n) + eigly(i,j,n,m)*slyscr(j,1,m)
                        alpc(n) = alpc(n) + eigly(i,j,n,m)*slyscr(j,2,m)
                        alpr(n) = alpr(n) + eigly(i,j,n,m)*slyscr(j,3,m)
                     end do
                  end do
                     if (dabs(alpl(2)) .lt. 1.d-5) then
                        alpl(2) = zero
                     end if
                     if (dabs(alpc(2)) .lt. 1.d-5) then
                        alpc(2) = zero
                     end if
                     if (dabs(alpr(2)) .lt. 1.d-5) then
                        alpr(2) = zero
                     end if
                  
                  do n = 1,ncomps
                     faclim = 2.
                     if (j .ge. jmin .and. j .le. jmax) then
                        if(strucy(i,j,n,n)*strucy(i,j+1,n,n).lt.-1.d-10)faclim=1.
                        if(strucy(i,j,n,n)*strucy(i,j-1,n,n).lt.-1.d-10)faclim=1.
                     end if
                     anslp = dmin1(faclim*dabs(alpl(n)),
     &                    dabs(alpc(n)),faclim*dabs(alpr(n)))
     &                    *dmax1(0.d0,dsign(1.d0,alpl(n)*alpr(n)))
                     den=merge(dabs(alpc(n)),1.d0,dabs(alpc(n)).ne.0.)
                     sly(i,j,n)=anslp/den
                  end do
                  
                  def= dabs(alpc(1)) *cdef*dabs(strucy(i,j,1,1)-strucy(i,j,2,1))
     &                 +dabs(alpc(2))*cdef*dabs(strucy(i,j,1,2)-strucy(i,j,2,2))
     &                 -dabs(eigvy(i,j,1)-eigvy(i,j,2))
                  
                  if(def.gt. 1.0d-10 .or. ldefy(i,j) )then
                     if(eigvy(i,j,1)*eigvy(i,j,2).lt.0)then
                        sly(i,j,1)=0.
                        sly(i,j,2)=0.
                     else
                        esnew=0.5*(eigvy(i,j,1)+eigvy(i,j,2))
                        slrmin=dmin1(sly(i,j,1),sly(i,j,2))
                        sly(i,j,1)=slrmin
                        sly(i,j,2)=slrmin
                        eigvy(i,j,1)=esnew
                        eigvy(i,j,2)=esnew
                     endif
                  endif
                  sly(i,j,1) = alpc(1)*sly(i,j,1)
                  sly(i,j,2) = alpc(2)*sly(i,j,2)
                  
               end do
            end do
         end if
      end if
      
      end
      
