c
c $Id: GODUNOV_F.F,v 1.43 2011-06-06 06:09:30 gpau Exp $
c
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


      subroutine FORT_SET_PARAMS(slope_order_in,use_unlim_in)
c
c     ----------------------------------------------------------
c     Set the slope order for interpolation to edges.
c     ----------------------------------------------------------
c 
      implicit none
      integer slope_order_in,use_unlim_in

#include "GODCOMM_F.H"

      slope_order = slope_order_in

      use_unlimited_slopes = (use_unlim_in .eq. 1)

      if (slope_order.ne.1.and.
#if (BL_SPACEDIM==2)
     &   slope_order.ne.2.and.
#endif
     &   slope_order.ne.4)then
        write(6,*)'FORT_SET_PARAMS : illegal value of slope_order = ',
     &            slope_order

        call bl_abort(" ")
      end if

      end

      subroutine get_mu_rho_phases(rho, mu, athick, rhoval, muval, 
     &     compcount, compidx, phaseidx, ncomps, nphases) 
      implicit none
c
c     ----------------------------------------------------------
c     Index of rho and mu given by:
c       1: phase with 2 components
c       2: phase with 1 component
c     ----------------------------------------------------------
c
      integer ncomps, nphases
      integer compcount(nphases)
      integer compidx(nphases,ncomps)
      integer phaseidx(ncomps)

      REAL_T  athick
      REAL_T  rho(nphases)
      REAL_T  rhoval(ncomps)
      REAL_T  mu(nphases)
      REAL_T  muval(ncomps)

      integer k

      compcount = 0
      compidx   = 0
      do k = 1,ncomps
         compcount(phaseidx(k)+1) = compcount(phaseidx(k)+1) + 1
         compidx(phaseidx(k)+1,compcount(phaseidx(k)+1)) = k
      end do
      do k = 1,nphases
         if (compcount(k) .eq. 1) then
            mu(2) = muval(compidx(k,1))
            rho(2) = rhoval(compidx(k,1))
         elseif (compcount(k) .eq. 2) then
c           take for now rho in a phase to be the same
c           take mu of the first component in the phase as base mu
            mu(1) = muval(compidx(k,1))
            rho(1) = rhoval(compidx(k,1))
            athick = muval(compidx(k,2))
         else
            print *, "ERROR: compcount(k) > 2)"
            stop
         end if

      end do

      end

c     =========================================================
c     Two-component Two-phase Model
c     =========================================================

      subroutine cmpflx(flux,ulp,urp,vtot,amua,amul,gstar)

c     compute flux based on quadratic relative permeability

      implicit none

      REAL_T flux, ulp, urp, vtot, amua, amul
      REAL_T gstar 

      if (vtot .ge. zero ) then
         call riemann(flux,ulp,urp,amua,amul,vtot,gstar)
      else
         call riemann(flux,urp,ulp,amua,amul,-vtot,-gstar)
         flux = -flux
      endif

      end

      subroutine cmpflx_arg(flux,u_arg,ulp,urp,vtot,amua,amul,gstar)
  
c     compute argument for flux based on quadratic relative permeability

      implicit none

      REAL_T flux, u_arg, ulp, urp, vtot, amua, amul
      REAL_T gstar, fg

      if (vtot .ge. zero) then
         call riemann_arg(u_arg,ulp,urp,amua,amul, vtot,gstar)
         flux = fg(u_arg,amua,amul,vtot,gstar)
      else
         call riemann_arg(u_arg,urp,ulp,amua,amul,-vtot,-gstar)
         flux = fg(u_arg,amua,amul,vtot,gstar)
      endif

      end

      subroutine cmpflx_vang(flux,ulp,urp,vtot,amua,amul,gstar,kr_coef)

c     compute flux based on van Genuchten rel. perm. function

      implicit none

      REAL_T flux, ulp, urp, vtot, amua, amul
      REAL_T gstar, kr_coef

      if (vtot .ge. zero ) then
         call riemann_vang(flux,ulp,urp,amua,amul,vtot,gstar,kr_coef)
      else
         call riemann_vang(flux,urp,ulp,amua,amul,-vtot,-gstar,kr_coef)
         flux = -flux
      endif
      end

      subroutine cmpflx_vang_arg(flux,u_arg,ulp,urp,vtot,
     &     amua,amul,gstar,kr_coef)

c ::: compute argument for flux based on van Genuchten rel. perm. function

      implicit none

      REAL_T flux, u_arg, ulp, urp, vtot, amua, amul
      REAL_T gstar, kr_coef, fg_vang
      
      if (vtot .ge. zero) then
         call riemann_vang_arg(u_arg,ulp,urp,amua,amul, vtot,gstar,kr_coef)
         flux = fg_vang(u_arg,amua,amul,vtot,gstar,kr_coef)
      else
         call riemann_vang_arg(u_arg,urp,ulp,amua,amul,-vtot,-gstar,kr_coef)
         flux = fg_vang(u_arg,amua,amul,vtot,gstar,kr_coef)
      endif
      end

      subroutine cmpflx_upwind(flux,ulp,urp,vtot,mu,gstar,
     &     l_kr_type,r_kr_type,l_kr_coef,r_kr_coef)

c     compute flux based on upwind scheme

      implicit none

      integer l_kr_type, r_kr_type
      REAL_T flux, ulp, urp, vtot, mu(2)
      REAL_T gstar, l_kr_coef, r_kr_coef
      REAL_T ld(2),le(2),rtmp(2),vlr(2),lam(2),vtmp

      rtmp = one
      call get_ld(le,urp,rtmp,mu,2,r_kr_type,r_kr_coef,zero)
      call get_ld(ld,ulp,rtmp,mu,2,l_kr_type,l_kr_coef,zero)
      vlr(1) = ld(1)/(ld(1)+ld(2))*(vtot+ld(2)*gstar)
      vlr(2) = le(1)/(le(1)+le(2))*(vtot+le(2)*gstar)
      vtmp   = half*(vlr(1)+vlr(2))
      vtmp   = cvmgt(zero,vtmp,dabs(vtmp) .lt. 1.d-15)
      lam(1) = cvmgp(ld(1),le(1),vtmp)
      vlr(1) = ld(2)/(ld(1)+ld(2))*(vtot-ld(1)*gstar)
      vlr(2) = le(2)/(le(1)+le(2))*(vtot-le(1)*gstar)
      vtmp   = half*(vlr(1)+vlr(2))
      vtmp   = cvmgt(zero,vtmp,dabs(vtmp) .lt. 1.d-15)
      lam(2) = cvmgp(ld(2),le(2),vtmp)
      flux = lam(1)/(lam(1)+lam(2))*(vtot+lam(2)*gstar)
      
      end

      subroutine cmpflx_upwind_arg(flux,u_arg,ulp,urp,vtot,
     &     mu,gstar,l_kr_type,r_kr_type,l_kr_coef,r_kr_coef)

c ::: compute argument for upwind flux

      implicit none

      integer l_kr_type, r_kr_type
      REAL_T  flux, u_arg(2), ulp, urp, vtot, mu(2)
      REAL_T  gstar,l_kr_coef, r_kr_coef
      REAL_T  ld(2),le(2),rtmp(2),vlr(2),lam(2),vtmp

      rtmp = one
      call get_ld(le,urp,rtmp,mu,2,r_kr_type,r_kr_coef,zero)
      call get_ld(ld,ulp,rtmp,mu,2,l_kr_type,l_kr_coef,zero)
      vlr(1) = ld(1)/(ld(1)+ld(2))*(vtot+ld(2)*gstar)
      vlr(2) = le(1)/(le(1)+le(2))*(vtot+le(2)*gstar)
      vtmp   = half*(vlr(1)+vlr(2))
      vtmp   = cvmgt(zero,vtmp,dabs(vtmp) .lt. 1.d-15)
      lam(1) = cvmgp(ld(1),le(1),vtmp)
      u_arg(1) = cvmgp(urp(1),ulp(1),vtmp)
      vlr(1) = ld(2)/(ld(1)+ld(2))*(vtot-ld(1)*gstar)
      vlr(2) = le(2)/(le(1)+le(2))*(vtot-le(1)*gstar)
      vtmp   = half*(vlr(1)+vlr(2))
      vtmp   = cvmgt(zero,vtmp,dabs(vtmp) .lt. 1.d-15)
      lam(2) = cvmgp(ld(2),le(2),vtmp)
      u_arg(2) = cvmgp(urp(2),ulp(2),vtmp)
      flux   = lam(1)/(lam(1)+lam(2))*(vtot+lam(2)*gstar)

      end

      subroutine riemann(flux,ul,ur,amua,amul,vtot,gstar)

c ::: main riemann solver for quadratic relative permeability

      implicit none
      REAL_T flux, ul, ur, amua, amul, vtot, gstar
      REAL_T fg, smin

      if ((gstar.eq.zero) .and. (vtot.eq.zero)) then
         flux = zero
      else 
         if (gstar .gt. zero ) then
            if (vtot-gstar/amua .ge. zero) then
               flux = fg(ul,amua,amul,vtot,gstar)
            else
               if (ul .lt. ur) then
                  flux = min(fg(ul,amua,amul,vtot,gstar),
     &                 fg(ur,amua,amul,vtot,gstar))
               else
                  call findmin(smin,vtot,gstar,amua,amul)
                  flux = max(fg(ul,amua,amul,vtot,gstar),
     &                 fg(ur,amua,amul,vtot,gstar))
                  if ((ul-smin)*(ur-smin) .lt. zero) then
                     flux = fg(smin,amua,amul,vtot,gstar)
                  end if
               end if
            end if
         else
            if (vtot+gstar/amul .ge. zero ) then
               flux = fg(ul,amua,amul,vtot,gstar)
            else
               if (ul.gt.ur) then
                  flux = max(fg(ul,amua,amul,vtot,gstar),
     &                 fg(ur,amua,amul,vtot,gstar))
               else
                  call findmin(smin,vtot,gstar,amua,amul)
                  flux = min(fg(ul,amua,amul,vtot,gstar),
     &                 fg(ur,amua,amul,vtot,gstar))
                  if ((ul-smin)*(ur-smin).lt. zero ) then
                     flux = fg(smin,amua,amul,vtot,gstar)
                  end if
               end if
            end if
         end if
      endif
      end

      subroutine riemann_vang(flux,ul,ur,amua,amul,
     &     vtot,gstar,kr_coef)

c ::: main riemann solver for van Genuchten permeability functions

      implicit none
      REAL_T flux, ul, ur, amua, amul, vtot, gstar
      REAL_T fg_vang, smin, kr_coef

      if ((gstar.eq.zero) .and. (vtot.eq.zero)) then
         flux = zero
      else 
         if (gstar .gt. zero ) then
            if (vtot-gstar/amua .ge. zero) then
               flux = fg_vang(ul,amua,amul,vtot,gstar,kr_coef)
            else
               if (ul .le. ur) then
                  flux = min(fg_vang(ul,amua,amul,vtot,gstar,kr_coef),
     &                 fg_vang(ur,amua,amul,vtot,gstar,kr_coef))
               else
                  call findext_vang2(smin,vtot,gstar,amua,amul,kr_coef)
                  flux = max(fg_vang(ul,amua,amul,vtot,gstar,kr_coef),
     &                 fg_vang(ur,amua,amul,vtot,gstar,kr_coef))
                  if ((ul-smin)*(ur-smin) .lt. zero) then
                     flux = fg_vang(smin,amua,amul,vtot,gstar,kr_coef)
                  end if
               end if
            end if
         else
            if (vtot+gstar/amul .ge. zero ) then
               flux = fg_vang(ul,amua,amul,vtot,gstar,kr_coef)
            else
               if (ul.ge.ur) then
                  flux = max(fg_vang(ul,amua,amul,vtot,gstar,kr_coef),
     &                 fg_vang(ur,amua,amul,vtot,gstar,kr_coef))
               else
                  call findext_vang2(smin,vtot,gstar,amua,amul,kr_coef)
                  flux = min(fg_vang(ul,amua,amul,vtot,gstar,kr_coef),
     &                 fg_vang(ur,amua,amul,vtot,gstar,kr_coef))
                  if ((ul-smin)*(ur-smin).lt. zero ) then
                     flux = fg_vang(smin,amua,amul,vtot,gstar,kr_coef)
                  end if
               end if
            end if
         end if
      endif
      end
      
      subroutine riemann_arg(u_arg,ul,ur,amua,amul,vtot,gstar)
     
c     ::: compute the argument for riemann solver for quadratic permeability

      implicit none
      REAL_T u_arg, ul, ur, amua, amul, vtot, gstar
      REAL_T fg, smin, flux_l, flux_r

      if (gstar .gt. zero ) then
         if (vtot-gstar/amua .gt. zero) then
            u_arg = ul
         else
            flux_l = fg(ul,amua,amul,vtot,gstar)
            flux_r = fg(ur,amua,amul,vtot,gstar)
            
            if (ul .lt. ur) then
               u_arg  = cvmgt(ul,ur,flux_l .lt. flux_r)
            else
               u_arg  = cvmgt(ul,ur,flux_l .gt. flux_r)
               call findmin(smin,vtot,gstar,amua,amul)
               if ((ul-smin)*(ur-smin) .lt. zero) then
                  u_arg = smin
               end if
            end if
         end if
      else
         if (vtot+gstar/amul .gt. zero ) then
            u_arg = ul
         else
            flux_l = fg(ul,amua,amul,vtot,gstar)
            flux_r = fg(ur,amua,amul,vtot,gstar)
            if (ul .gt. ur) then
               u_arg  = cvmgt(ul,ur,flux_l .gt. flux_r)
            else
               u_arg  = cvmgt(ul,ur,flux_l .lt. flux_r)
               call findmin(smin,vtot,gstar,amua,amul)
               if ((ul-smin)*(ur-smin).lt. zero ) then
                  u_arg = smin
               end if
            end if
         end if
      end if
      end

      subroutine riemann_vang_arg(u_arg,ul,ur,amua,amul,vtot,gstar,kr_coef)
      
c     ::: compute the argument for riemann solver

      implicit none
      REAL_T u_arg, ul, ur, amua, amul, vtot, gstar
      REAL_T fg_vang, smin, flux_l, flux_r
      REAL_T kr_coef

      if (gstar .gt. zero ) then
         if (vtot-gstar/amua .ge. zero) then
            u_arg = ul
         else
            flux_l = fg_vang(ul,amua,amul,vtot,gstar,kr_coef)
            flux_r = fg_vang(ur,amua,amul,vtot,gstar,kr_coef)
            
            if (ul .lt. ur) then
               u_arg  = cvmgt(ul,ur,flux_l .lt. flux_r)
            else
               u_arg  = cvmgt(ul,ur,flux_l .gt. flux_r)
               call findext_vang2(smin,vtot,gstar,amua,amul,kr_coef)
               if ((ul-smin)*(ur-smin) .lt. zero) then
                  u_arg = smin
               end if
            end if
         end if
      else
         if (vtot+gstar/amul .ge. zero ) then
            u_arg = ul
         else
            flux_l = fg_vang(ul,amua,amul,vtot,gstar,kr_coef)
            flux_r = fg_vang(ur,amua,amul,vtot,gstar,kr_coef)
            if (ul .gt. ur) then
               u_arg  = cvmgt(ul,ur,flux_l .gt. flux_r)
            else
               u_arg  = cvmgt(ul,ur,flux_l .lt. flux_r)
               call findext_vang2(smin,vtot,gstar,amua,amul,kr_coef)
               if ((ul-smin)*(ur-smin).lt. zero ) then
                  u_arg = smin
               end if
            end if
         end if
      end if
      end
      
      subroutine findmin(smin,vtot,gstar,amua,amul)

      implicit none
      
      REAL_T smin, vtot, gstar, amua, amul
      REAL_T dels, tol
      integer it
      parameter (tol = 1.D-14)
      
      it   = 0
      smin = amua**third
      smin = smin/(smin+amul**third)
      dels = one

 100  if ((abs(dels) .gt. tol) .and. (it .lt. 20)) then             
         dels = -fgpmin(smin,amua,amul,vtot,gstar)/
     &        fgminpp(smin,amua,amul,vtot,gstar)      
         smin = smin+dels
         it = it + 1
         goto 100
      end if

      contains

        function fgpmin(x,amua,amul,vtot,gstar) result(r)
        REAL_T x, amua, amul, vtot, gstar, r
        r = vtot-gstar*(x**3/amua-(one-x)**3/amul)
        end function fgpmin
        !
        ! ::: another one
        ! ::: George's way of numerically taking the second derivative
        !
        function fgminpp(x,amua,amul,vtot,gstar) result (r)
        REAL_T x, amua, amul, vtot, gstar, eps, r
        parameter (eps = 1.d-14)

        r = (fgpmin(x+eps,amua,amul,vtot,gstar)
     &       - fgpmin(x-eps,amua,amul,vtot,gstar))/(2*eps)
      
        end function fgminpp

      end
      
      subroutine findext_vang2(smin,vtot,gstar,amua,amul,kr_coef)
c     Try bisection using the whole interval... Newton iteration kept failing 

      implicit none

      integer nc
      REAL_T s1,s2,df1,df2,fgpmin_vang
      REAL_T dftmp,smin,vtot,gstar,amua,amul
      REAL_T kr_coef
      REAL_T, parameter :: eps = 1.0d-14

      s1  = eps
      s2  = 0.99999d0
      df1 = fgpmin_vang(s1,amua,amul,vtot,gstar,kr_coef)
      df2 = fgpmin_vang(s2,amua,amul,vtot,gstar,kr_coef)

      smin  = zero
      dftmp = zero
      if (sign(one,df1)-sign(one,df2) .ne. zero) then

c     ## This assumes only one sign change in the interval

         do nc = 1,100
            smin = half*(s1+s2)

            dftmp = fgpmin_vang(smin,amua,amul,vtot,gstar,kr_coef)

            if (dftmp*df1 .lt. zero) then
               df2 = dftmp
               s2 = smin
            elseif (dftmp*df2 .lt. zero) then
               df1 = dftmp
               s1 = smin
            else 
c     # landed on the place where dftmp==zero
               goto 70
            endif

            if (dabs(s1-s2) .lt. eps) then
               goto 70
            endif
         enddo
      else
         smin = one
      endif

 70   end
      
      REAL_T function fg(x,amua,amul,vtot,gstar)
c     ::: determine the flux
      REAL_T x, amua, amul, vtot, gstar
      REAL_T alpha, f
      
      alpha = amua/amul
      f=x*x/(x*x+alpha*(one-x)**2)
      fg=f*(vtot+gstar*(one-x)**2/amul)

      end

      REAL_T function fg_vang(x,amua,amul,vtot,gstar,kr_coef)
c     
c     Van Genuchten flux function
c     x is the current guess for the saturation where f
c     reaches a minimum
c
      implicit none

      REAL_T x, amua, amul, vtot, gstar
      REAL_T rtmp(2),mu(2),ld(2),st(2)
      REAL_T lam1,lam2,ltot
      REAL_T kr_coef

      st(1) = x
      st(2) = one-x
      rtmp  = one
      mu(1) = amua
      mu(2) = amul

      call get_ld(ld,st,rtmp,mu,2,3,kr_coef,zero)
      lam1 = ld(1)
      lam2 = ld(2)
      ltot = lam1+lam2
      fg_vang = lam1/ltot * (vtot+gstar*lam2)

      end

      REAL_T function fgpmin_vang(x,amua,amul,vtot,gstar,kr_coef)
c ::: The first derivative of the flux function with van Genuchten
c ::: permeability functions

      implicit none

      REAL_T x, amua, amul, vtot, gstar
      REAL_T mu(2),rtmp(2),ld(2),dld(3),st(2)
      REAL_T ltot,lam1,lam2,dlam1,dlam2,kr_coef

      st(1) = x
      st(2) = one-x
      rtmp  = one
      mu(1) = amua
      mu(2) = amul

c ::: helper functions from EOS_F.F
      call get_ld(ld,st,rtmp,mu,2,3,kr_coef,zero)
      call get_dld(dld,st,rtmp,mu,2,3,kr_coef,zero)

      lam1  = ld(1)
      lam2  = ld(2)
      ltot  = lam1+lam2
      dlam1 = dld(1)
      dlam2 = dld(2)

      fgpmin_vang = (dlam1*lam2*(vtot+gstar*lam2) 
     .     + lam1*(-vtot + gstar*lam1)*dlam2)/ltot**2

      end

      subroutine find_max_eig(s_l,s_r,ncomps,
     &     lbd_l,dlbd_l,d2lbd_l,lbd_r,dlbd_r,d2lbd_r,
     &     vel,krcoef,nkrcoef,
     &     krhog,rhoval,muval,eigmax)
c        
c ::: Estimate maximum eigenvalues
c        There are possibly two extrema for the eigenvalues.  
c        In those cases, we determine the minimum of the flux function,
c        and find the extrema in both directions.
c        The value is given by (dv1/ds)
c
      implicit none

c ::: Arguments
      integer nkrcoef
      integer ncomps
      REAL_T s_l(ncomps), s_r(ncomps)
      REAL_T lbd_l(2), lbd_r(2)
      REAL_T dlbd_l(3),  dlbd_r(3)
      REAL_T d2lbd_l(3), d2lbd_r(3)
      REAL_T krcoef(nkrcoef)
      REAL_T vel
      REAL_T rhoval(ncomps)
      REAL_T muval(ncomps)
      REAL_T krhog
      REAL_T eigmax

c ::: Internal variables
      integer m,nc,ncycle   
      integer l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T eigtmp
      REAL_T lt,f1,f2
      REAL_T slr(2,2),df(2,2),d2f(2,2),dftmp, stmp(2),dtmp(2)
      REAL_T sat(2),lam(2),dlam(3),d2lam(3),rtmp(2)
c     REAL_T fg,ipoint,fgp,d2tmp,fgtmp

      eigmax = zero
      rtmp   = one
      l_kr_type = nint(krcoef(1))
      l_kr_coef = krcoef(2)
      l_sat_res = krcoef(3)                                 

c     :: right hand side          
      call get_st(sat,s_r,rhoval,l_sat_res,ncomps)
      if (l_kr_type .eq. 3) then
         sat(1) = min(one-1.d-6,sat(1))
         sat(1) = max(1.d-6,sat(1))
         sat(2) = one-sat(1)
      end if
      lam      = lbd_r
      dlam     = dlbd_r
      d2lam    = d2lbd_r
      slr(1,1) = sat(1)
      lt       = lam(1)+lam(2)
      if (l_kr_type .eq. 2) then
         df(1,1)  = vel*(dlam(1)/lt+lam(1)*dlam(3))+
     &        krhog*(dlam(2)/lt+lam(1)*lam(2)*dlam(3))
         d2f(1,1) = vel*(d2lam(1)/lt+lam(1)*d2lam(3)+2*dlam(1)*dlam(3))+
     &        krhog*(d2lam(2)/lt+lam(1)*lam(2)*d2lam(3)+2*dlam(2)*dlam(3))
      else if (l_kr_type .eq. 3) then
         f1 = lam(1)/lt * (vel + krhog*lam(2))
         df(1,1) = (dlam(1)*lam(2)*(vel+krhog*lam(2))
     &        + lam(1)*(-vel + krhog*lam(1))*dlam(2))/(lt**2)     
         d2f(1,1) = 2.d0*krhog*dlam(1)*dlam(2)/lt
     &        - 2.d0*(vel+krhog*lam(2))*dlam(1)*(dlam(1)+dlam(2))/lt**2
     &        - 2.d0*krhog*lam(1)*dlam(2)*(dlam(1)+dlam(2))/lt**2
     &        + 2.d0*lam(1)*(vel+krhog*lam(2))*(dlam(1)+dlam(2))**2/lt**3
     &        + (vel+krhog*lam(2))*d2lam(1)/lt
     &        + krhog*lam(1)*d2lam(2)/lt
     &        - lam(1)*(vel+krhog*lam(2))*(d2lam(1)+d2lam(2))/lt**2
      end if

c     :: left hand side
      call get_st(sat,s_l,rhoval,l_sat_res,ncomps)
      if (l_kr_type .eq. 3) then
         sat(1) = min(one-1.d-6,sat(1))
         sat(1) = max(1.d-6,sat(1))
         sat(2) = one-sat(1)
      end if
      lam   = lbd_l
      dlam  = dlbd_l
      d2lam = d2lbd_l
      slr(1,2) = sat(1)
      lt       = lam(1)+lam(2)
      if (l_kr_type .eq. 2) then
         df(1,2)  = vel*(dlam(1)/lt+lam(1)*dlam(3))+
     &        krhog*(dlam(2)/lt+lam(1)*lam(2)*dlam(3))
         d2f(1,2) = vel*(d2lam(1)/lt+lam(1)*d2lam(3)+2*dlam(1)*dlam(3))+
     &        krhog*(d2lam(2)/lt+lam(1)*lam(2)*d2lam(3)+2*dlam(2)*dlam(3))
      elseif (l_kr_type .eq. 3) then
         f2 = lam(1)/lt * (vel + krhog*lam(2))
         df(1,2) = (dlam(1)*lam(2)*(vel+krhog*lam(2))
     &        + lam(1)*(-vel + krhog*lam(1))*dlam(2))/(lt**2)
         d2f(1,2) = 2.d0*krhog*dlam(1)*dlam(2)/lt
     &        - 2.d0*(vel+krhog*lam(2))*dlam(1)*(dlam(1)+dlam(2))/lt**2
     &        - 2.d0*krhog*lam(1)*dlam(2)*(dlam(1)+dlam(2))/lt**2
     &        + 2.d0*lam(1)*(vel+krhog*lam(2))*(dlam(1)+dlam(2))**2/lt**3
     &        + (vel+krhog*lam(2))*d2lam(1)/lt
     &        + krhog*lam(1)*d2lam(2)/lt
     &        - lam(1)*(vel+krhog*lam(2))*(d2lam(1)+d2lam(2))/lt**2
      end if

c     find arg min of flux if there is a minimum in the interval
      stmp = slr(1,:)
      dtmp = df(1,:)
      ncycle = 1
      if (dtmp(1)*dtmp(2) .lt. zero) then
         do nc = 1,100
            sat(1) = half*(stmp(1)+stmp(2))
            sat(2) = one-sat(1)
            call get_ld(lam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
            call get_dld(dlam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
            if (l_kr_type .eq. 2) then
               dftmp = vel*(d2lam(1)/lt+lam(1)*d2lam(3)+2*dlam(1)*dlam(3))+
     &              krhog*(d2lam(2)/lt+lam(1)*lam(2)*d2lam(3)+2*dlam(2)*dlam(3))
            else if (l_kr_type .eq. 3) then
               dftmp = (dlam(1)*lam(2)*(vel+krhog*lam(2))
     &              + lam(1)*(-vel + krhog*lam(1))*dlam(2))/(lt**2)
            end if
            if (dftmp*dtmp(1) .lt. zero) then
               dtmp(2)   = dftmp
               stmp(2) = sat(1)
            elseif (dftmp*dtmp(2) .lt. zero) then
               dtmp(1)   = dftmp
               stmp(1) = sat(1)
            end if
            if (dabs(stmp(1)-stmp(2)) .lt. 1.d-15) then
               goto 79
            end if
         end do
         
 79      ncycle = 2
         sat(1) = half*(stmp(1)+stmp(2))
         sat(2) = one -sat(1)
         call get_ld(lam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
         call get_dld(dlam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
         call get_d2ld(d2lam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
         slr(2,1) =slr(1,2)
         df(2,1)  = df(1,2)
         d2f(2,1) = d2f(1,2)
         slr(1,2) = sat(1)
         slr(2,2) = sat(1)
         lt     = lam(1)+lam(2)
         if (l_kr_type .eq. 2) then
            df(1,2)  = vel*(dlam(1)/lt+lam(1)*dlam(3))+
     &           krhog*(dlam(2)/lt+lam(1)*lam(2)*dlam(3))
            d2f(1,2) = vel*(d2lam(1)/lt+lam(1)*d2lam(3)+2*dlam(1)*dlam(3))+
     &           krhog*(d2lam(2)/lt+lam(1)*lam(2)*d2lam(3)+2*dlam(2)*dlam(3))
         elseif (l_kr_type .eq. 3) then
            f2 = lam(1)/lt * (vel + krhog*lam(2))
            df(1,2) = (dlam(1)*lam(2)*(vel+krhog*lam(2))
     &           + lam(1)*(-vel + krhog*lam(1))*dlam(2))/(lt**2)
c     
            d2f(1,2) = 2.d0*krhog*dlam(1)*dlam(2)/lt
     &           - 2.d0*(vel+krhog*lam(2))*dlam(1)*(dlam(1)+dlam(2))/lt**2
     &           - 2.d0*krhog*lam(1)*dlam(2)*(dlam(1)+dlam(2))/lt**2
     &           + 2.d0*lam(1)*(vel+krhog*lam(2))*(dlam(1)+dlam(2))**2/lt**3
     &           + (vel+krhog*lam(2))*d2lam(1)/lt
     &           + krhog*lam(1)*d2lam(2)/lt
     &           - lam(1)*(vel+krhog*lam(2))*(d2lam(1)+d2lam(2))/lt**2
         end if
         df(2,2)  = df(1,2)
         d2f(2,2) = d2f(2,2)
      end if
               
      do m = 1,ncycle
         eigtmp = max(dabs(df(m,1)),dabs(df(m,2)))
         if (sign(one,d2f(m,1))-sign(one,d2f(m,2)) .ne. zero .or. 
     &        ncycle .eq. 2) then
            do nc = 1,100
               sat(1) = half*(slr(m,1) + slr(m,2))
               sat(2) = one-sat(1)
               call get_ld(lam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
               call get_dld(dlam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
               call get_d2ld(d2lam,sat,rtmp,muval,ncomps,l_kr_type,l_kr_coef,zero)
               lt = lam(1)+lam(2)
               if (l_kr_type .eq. 2) then
                  dftmp = vel*(d2lam(1)/lt+lam(1)*d2lam(3)+2*dlam(1)*dlam(3))+
     &                 krhog*(d2lam(2)/lt+lam(1)*lam(2)*d2lam(3)+2*dlam(2)*dlam(3))
               elseif(l_kr_type .eq. 3) then
                  dftmp = 2.d0*krhog*dlam(1)*dlam(2)/lt
     &                 - 2.d0*(vel+krhog*lam(2))*dlam(1)*(dlam(1)+dlam(2))/lt**2
     &                 - 2.d0*krhog*lam(1)*dlam(2)*(dlam(1)+dlam(2))/lt**2
     &                 + 2.d0*lam(1)*(vel+krhog*lam(2))*(dlam(1)+dlam(2))**2/lt**3
     &                 + (vel+krhog*lam(2))*d2lam(1)/lt
     &                 + krhog*lam(1)*d2lam(2)/lt
     &                 - lam(1)*(vel+krhog*lam(2))*(d2lam(1)+d2lam(2))/lt**2
               end if
               if (sign(one,dftmp)-sign(one,d2f(m,1)) .eq. zero  ) then
                  d2f(m,1)  = dftmp
                  slr(m,1) = sat(1)
               else
                  d2f(m,2)  = dftmp
                  slr(m,2) = sat(1)
               end if                        
               if (dabs(slr(m,1)-slr(m,2)) .lt. 1.d-15) then
                  go to 89
               end if
            end do
            
 89         if (l_kr_type .eq. 2) then
               eigtmp = vel*(dlam(1)/lt+lam(1)*dlam(3))+
     &              krhog*(dlam(2)/lt+lam(1)*lam(2)*dlam(3))
            else if (l_kr_type .eq. 3) then
               eigtmp = ((vel + krhog*lam(2))*dlam(1)*lam(2) 
     &              + lam(1)*(-vel+krhog*lam(1))*dlam(2))/(lt**2)
            end if
         end if
         
c     Alternative for kr_type == 3
c     ### For calculation of the inflection point (which depends on vtot)
c     ### The value of f' at the inflection point will give the maximum speed
c     
c     ### use a newton iteration scheme to solve for the roots of 
c     ### g(x) := f(ipoint) - f(sat)/(ipoint - sat) - f'(ipoint) == 0
c$$$  if ((grav.eq.0.d0).and.(m.le.0.5d0)) then
c$$$  write(*,*) 'calling findipoint, but grav==0!',d2f(1),d2f(2),krhog,sat(1),lt
c$$$  endif
c$$$  call findipoint(ipoint,fgp,stmp,f1,df(1),muval,rhoval
c$$$  .                 ,ncomps,vel,krhog,kr_type,sat_res)
c$$$  fgtmp = fgp
c$$$  call findipoint(ipoint,fgp,stmp,f2,df(2),muval,rhoval
c$$$  .                 ,ncomps,vel,krhog,kr_type,sat_res)
c$$$  c     ### f'(ipoint) will be the fastest speed and is calculated in
c$$$  c     ### findipoint, so it should get returned properly without the
c$$$  c     ### extra function calls
c$$$  eigtmp = max(dabs(fgp),dabs(fgtmp))
         eigmax = max(eigmax,dabs(eigtmp/(one-l_sat_res)))
      end do 

      end 

      subroutine find_upwind_eig(s_l,s_r,ncomps,
     &     lbd_l,dlbd_l,lbd_r,dlbd_r,
     &     vel,krcoef,nkrcoef,krhog,eigmax)
c        
c ::: Estimate maximum eigenvalues
c        There are possibly two extrema for the eigenvalues.  
c        In those cases, we determine the minimum of the flux function,
c        and find the extrema in both directions.
c        The value is given by (dv1/ds)
c
      implicit none

c ::: Arguments
      integer ncomps
      integer nkrcoef
      REAL_T s_l(ncomps), s_r(ncomps)
      REAL_T lbd_l(2), lbd_r(2)
      REAL_T dlbd_l(3),  dlbd_r(3)
      REAL_T krcoef(nkrcoef)
      REAL_T vel
      REAL_T rhoval(ncomps)
      REAL_T krhog
      REAL_T eigmax

c ::: Internal variables
      integer l_kr_type
      REAL_T eigtmp, l_sat_res
      REAL_T lt
      REAL_T lam(2),dlam(3)
      REAL_T vlr(2),vtmp

      eigmax = zero
      l_kr_type = nint(krcoef(1))
      l_sat_res = krcoef(3)
      vlr(1) = lbd_r(1)/(lbd_r(1)+lbd_r(2))*(vel+lbd_r(2)*krhog)
      vlr(2) = lbd_l(1)/(lbd_l(1)+lbd_l(2))*(vel+lbd_l(2)*krhog)
      vtmp   = half*(vlr(1)+vlr(2))
      vtmp   = cvmgt(zero,vtmp,dabs(vtmp) .lt. 1.d-15)
      lam(1) = cvmgp(lbd_l(1),lbd_r(1),vtmp)
      dlam(1) = cvmgp(dlbd_l(1),dlbd_r(1),vtmp)

      vlr(1) = lbd_r(2)/(lbd_r(1)+lbd_r(2))*(vel-lbd_r(1)*krhog)
      vlr(2) = lbd_l(2)/(lbd_l(1)+lbd_l(2))*(vel-lbd_l(1)*krhog)
      vtmp   = half*(vlr(1)+vlr(2))
      vtmp   = cvmgt(zero,vtmp,dabs(vtmp) .lt. 1.d-15)
      lam(2) = cvmgp(lbd_l(2),lbd_r(2),vtmp)
      dlam(2) = cvmgp(dlbd_l(2),dlbd_r(2),vtmp)
      dlam(3) = cvmgp(dlbd_l(3),dlbd_r(3),vtmp)
      lt = lam(1) + lam(2)
      if (l_kr_type .eq. 2) then
         eigtmp = vel*(dlam(1)/lt+lam(1)*dlam(3))+
     &        krhog*(dlam(2)/lt+lam(1)*lam(2)*dlam(3))
      else if (l_kr_type .eq. 3) then
         eigtmp = ((vel + krhog*lam(2))*dlam(1)*lam(2) 
     &        + lam(1)*(-vel+krhog*lam(1))*dlam(2))/(lt*lt)
      end if
      eigmax = dabs(eigtmp/(one-l_sat_res))
      
      end 

      subroutine findipoint(ipoint,fgp,sat,fs1,dfs1,muval,rhoval,
     .     ncomps,vtot,gstar,kr_type,kr_coef,sat_res)
      implicit none
      integer ncomps,it,kr_type
      REAL_T ipoint,tol
      REAL_T sat(ncomps),lam(ncomps),dlam(ncomps),d2lam(ncomps)
      REAL_T s(ncomps),rhoval(ncomps),muval(ncomps),rtmp(ncomps)
      REAL_T fs1,dfs1,h,hprime,fg,fgp,fgpp
      REAL_T vtot,gstar,dels,lt,sat_res,kr_coef
      parameter (tol = 1.D-14)
c     ### Minimize the function f(x)-f(s1)/(x-s1) - f'(x) ==0
c     ### to find the saturation where there is an inflection point
c
c     xnew = xold - f'(x)/f''(x)
c
      it = 0
c     ### Start the iterator halfway between the two saturations
      s(1) = 0.5d0*(sat(1)+sat(2))
      rtmp(:) = one
      dels = 1
 100  if (abs(dels) .gt. tol) then
         if (it .lt. 20) then
c     ### Helper functions that evaluate the derivatives of lambda
c     ### that are necessary for calculating the derivative of the flux
c     ### functions
            call get_ld(lam,s,rtmp,muval,ncomps,kr_type,kr_coef,0.d0)
            call get_dld(dlam,s,rtmp,muval,ncomps,kr_type,kr_coef,0.d0)
            call get_d2ld(d2lam,s,rtmp,muval,ncomps,kr_type,kr_coef,0.d0)
            lt = lam(1)+lam(2)
c     ### Flux function evaluated at the new guess, s(1)
            fg = (lam(1)/lt)*(vtot+gstar*lam(2))
c     ### First derivative of the flux function evaluated at the new guess
c     
            fgp = ((vtot + gstar*lam(2))*dlam(1)*lam(2) 
     .           + lam(1)*(-vtot+gstar*lam(1))*dlam(2))/(lt**2)
c     ### Second derivative of the flux function evaluated at the new guess
c     ###
            fgpp= 2.d0*gstar*dlam(1)*dlam(2)/lt
     .           - 2.d0*(vtot+gstar*lam(2))*dlam(1)*(dlam(1)+dlam(2))/lt**2
     .           - 2.d0*gstar*lam(1)*dlam(2)*(dlam(1)+dlam(2))/lt**2
     .           + 2.d0*lam(1)*(vtot+gstar*lam(2))*(dlam(1)+dlam(2))**2/lt**3
     .           + (vtot+gstar*lam(2))*d2lam(1)/lt
     .           + gstar*lam(1)*d2lam(2)/lt
     .           - lam(1)*(vtot+gstar*lam(2))*(d2lam(1)+d2lam(2))/lt**2
            write(*,*) s(1), fg, fgp, fgpp, dels, it
            if (s(1).eq.sat(1)) then
               write(*,*) sat(1),s(1),fg,fgp,fgpp
               call bl_abort('findipoint: divide by zero!')
            elseif (s(1).lt.0.d0) then
               write(*,*) s(1), sat(1), fg, fgp, fgpp
               call bl_abort('negative saturations!')
            endif
         else
c     ### Max number of iterations was exceeded, throw error
            if (dels.lt.10e-8) then
               goto 110
            else 
               write(*,*) 'Newton', dels, s(1), sat(1), it, fg, fgp, fgpp
               call bl_abort('Newton iteration failed in find ipoint!')
            endif
         endif
         h = (fg - fs1)/(s(1)-sat(1)) - fgp
         hprime = fgp/(s(1)-sat(1)) - (fg - fs1)*(s(1)-sat(1))**(-2.d0)
     .        - fgpp
         if (hprime .eq. 0.d0) then
            goto 110
         else
            dels = -h/hprime
            s(1) = s(1) + dels
         endif

         if (s(1).gt.1.d0) then
            s(1) = 1.d0
            goto 110
         endif
         it = it+1
         goto 100
      endif
      
 110  ipoint = s(1)
c     ### One more time to get the correct derivative 
c     ### (make this a helper function too or rewrite)
      call get_ld(lam,s,rtmp,muval,ncomps,kr_type,kr_coef,0.d0)
      call get_dld(dlam,s,rtmp,muval,ncomps,kr_type,kr_coef,0.d0)
      call get_d2ld(d2lam,s,rtmp,muval,ncomps,kr_type,kr_coef,0.d0)
      lt = lam(1)+lam(2)
      fgp = ((vtot + gstar*lam(2))*dlam(1)*lam(2) 
     .     + lam(1)*(-vtot+gstar*lam(1))*dlam(2))/(lt**2) 
      end



c     =========================================================
c     Polymer Model
c     =========================================================

      subroutine get_aux_i(s,u,compcount,compidx,rhoval,ncomps,nphases)
c
c     ----------------------------------------------------------
c     Convert from primary to auxiliary variables for the 
c     polymer problem.
c     ----------------------------------------------------------
c 
      integer ncomps, nphases
      integer compcount(nphases)
      integer compidx(nphases,ncomps)
      REAL_T  s(ncomps)
      REAL_T  u(ncomps)
      REAL_T  sat
      REAL_T  rhoval(ncomps)

      integer k
      REAL_T  c
      
      do k = 1,nphases
         if (compcount(k) .eq. 2) then
            sat = s(compidx(k,1))/rhoval(compidx(k,1))
     &            + s(compidx(k,2))/rhoval(compidx(k,2))
            sat = dmin1(one,sat)
            sat = dmax1(zero,sat)
                  
            c = zero
            if (sat .gt. 1.d-10) then
               c = s(compidx(k,2))/(s(compidx(k,1))+s(compidx(k,2)))
            end if
         else
            sat = s(compidx(k,1))/rhoval(compidx(k,1))
         end if
      end do
      u(1) = sat
      u(2) = c*sat

      end

      subroutine eigen_polymer(vel,gstar,amuw,amul,athick,
     &     u,eigl,eigr,eigval,cstruc,ldefic)
c
c     ----------------------------------------------------------
c     eigenvalue solver
c     assumes components in order s,c   
c     ----------------------------------------------------------
c

c     arguments
      REAL_T vel
      REAL_T gstar
      REAL_T amuw, amul, athick
      REAL_T u(2)
      REAL_T eigl(2,2),eigr(2,2),eigval(2), cstruc(2,2)
      logical ldefic

c     local variables
      REAL_T s, c
      REAL_T amua, akrl, akra
      REAL_T amobl, amoba, amobt
      REAL_T vt, v
      REAL_T dmuadc, dvdc, dvds
      REAL_T dkrlds, dkrads, d2krlds2, d2krads2
      REAL_T dmoblds, dmobads, dmobtds
      REAL_T d2moblds, d2mobads, d2mobtds, d2vds2, d2vdsdc
      REAL_T diff
      REAL_T alen1, alen2
      REAL_T, parameter :: eps = 1.d-12
      
      vt = vel
      s  = u(1)
      c  = cvmgt(u(2)/s, zero, s .gt. zero)

      amua = amuw+athick*c
      akrl = (one-s)**2
      akra = s**2
      amobl=akrl/amul
      amoba=akra/amua
      amobt=amobl+amoba
      v=amoba*(vt+gstar*amobl)/amobt
      dmuadc=athick
      dvdc=-dmuadc*v*amobl/(amua*amobt)
      dkrlds=-2*(one-s)
      dkrads=2*s
      dmoblds=dkrlds/amul
      dmobads=dkrads/amua
      dmobtds=dmoblds+dmobads
      dvds=(dmobads*(vt+gstar*amobl)+amoba*gstar*dmoblds-v*dmobtds)/amobt
      d2krlds2=two
      d2krads2=two
      d2moblds=d2krlds2/amul
      d2mobads=d2krads2/amua
      d2mobtds=d2mobads+d2moblds
      d2vds2=(d2mobads*(vt+gstar*amobl)+2*dmobads*gstar*dmoblds
     &     +amoba*gstar*d2moblds-2*dvds*dmobtds-v*d2mobtds)/amobt
      d2vdsdc=-dmuadc*
     &     (dvds*amobl+v*(amoba*dmoblds-amobl*dmobads)/amobt)/(amua*amobt)

      eigval(1) = dvds
      eigval(2) = cvmgt(v/s, dvds, s .gt. 1.d-12)

      diff=eigval(2)-eigval(1)
      ldefic = (dabs(diff)-dabs(dvdc)) .le. eps
      if (.not.ldefic) then
         alen1 = sqrt(one+c*c)
         alen2 = sqrt(c*c +(c*dvdc+s*diff)**2) 
         eigr(1,1)=one/alen1
         eigr(2,1)=c/alen1
         eigr(1,2)=dvdc/alen2
         eigr(2,2)=(c*dvdc+s*diff)/alen2
         eigl(1,1)=(one+c*dvdc/(s*diff))*alen1
         eigl(2,1)=-c/(s*diff)*alen2
         eigl(1,2)=-dvdc/(s*diff)*alen1
         eigl(2,2)=one/(s*diff)*alen2
         cstruc(1,1)=(d2vds2*s - c*d2vdsdc)/s
         cstruc(1,2)=d2vdsdc/s
         cstruc(2,1)=-(s*diff+c*dvdc)/s**2
         cstruc(2,2)=dvdc/s**2
        
      elseif (s.gt.1.0d-12) then
         alen1 = sqrt(one+c*c)
         eigr(1,2)=one/alen1
         eigr(2,2)=c/alen1
         eigr(1,1)=zero
         eigr(2,1)=one
         eigl(2,1)=s
         eigl(1,1)=-c/s*alen1
         eigl(2,2)=zero
         eigl(1,2)=one/s*alen1        
         cstruc(1,2)=d2vdsdc/s
         cstruc(1,1)=(d2vds2*s - c*d2vdsdc)/s
         cstruc(2,2)=dvdc/s**2
         cstruc(2,1)=-(s*diff+c*dvdc)/s**2

      else
         eigr(1,1)=one
         eigr(2,1)=zero
         eigr(1,2)=zero
         eigr(2,2)=one
         eigl(1,1)=one
         eigl(2,1)=zero
         eigl(1,2)=zero
         eigl(2,2)=one
         cstruc(1,1)=d2vds2
         cstruc(1,2)=zero
         cstruc(2,1)=zero
         cstruc(2,2)=zero

      endif
      
      ldefic =ldefic .and. s.gt.1.d-10

      end


      subroutine cmpflx_polymer(um,up,ulft,urgt,eigvl,eigvr,flux,vel,gravity,
     &     amuw,amul,athick,rhoa,rhol,kappa,dx,dt)
c
c     ----------------------------------------------------------
c     This compute the flux w.r.t. to the auxiliary variables.
c     arguments
c     ----------------------------------------------------------
c
      REAL_T um(2), up(2)
      REAL_T ulft(2), urgt(2)
      REAL_T eigvl(2), eigvr(2)
      REAL_T flux(2)
      REAL_T vel
      REAL_T gravity
      REAL_T amuw, amul, athick
      REAL_T rhoa, rhol
      REAL_T kappa
      REAL_T dx,dt
      
c     internal variables
      REAL_T gstar,vtot

      vtot  = vel
      gstar = kappa*gravity*(rhoa-rhol)

      if (abs(gstar).gt. 1.d-15*abs(vtot))then
         call eoflux(ulft,urgt,flux,vtot,gstar,amuw,amul,athick)
      else
         if (vtot.gt. 1.d-15) then
            call evflux(ulft,flux,vtot,gstar,amuw,amul,athick)
         else
            call evflux(urgt,flux,vtot,gstar,amuw,amul,athick)
         endif
      endif

c ::: artificial viscosity.  Not used.
c$$$      sp1=eigvl(1)-eigvr(1)
c$$$      sp2=eigvl(2)-eigvr(2)
c$$$      satl = um(1) 
c$$$      satr = up(1) 
c$$$      cl = zero
c$$$      if (dabs(um(1)).gt.1.0d-10) then
c$$$        cl=um(2)/um(1)
c$$$      endif
c$$$      cr = zero
c$$$      if (dabs(up(1)).gt.1.0d-10) then
c$$$        cr=up(2)/up(1)
c$$$      endif
c$$$      dsmax=dmax1(dabs(cl-cr),dabs(satl-satr))
c$$$      dsmax=cvmgt(.2d0*dsmax,zero,dsmax.gt..2)
c$$$      dsmax=0.
c$$$      avisc=dmax1(zero,sp1,sp2,dx*dsmax/dt)
c$$$
c$$$      coefv = 0
c$$$      do ic=1,2
c$$$         flux(ic)=flux(ic)-coefv*avisc*(up(ic)- um(ic))
c$$$      end do

      end

      subroutine eoflux(ulft,urgt,fluxsv,vtot,grav,amuw,amul,athick)

      REAL_T ulft(2),urgt(2),fluxsv(2)
      REAL_T vtot, grav
      REAL_T amuw,amul,athick

      REAL_T fluxg(2),lleft(2),lright(2),kleft(2),kright(2),da(2)
      REAL_T fluxl(2),fluxr(2)
      REAL_T eigll(2,2),eigrr(2,2),eigrl(2,2),eiglr(2,2)
      REAL_T eigvalr(2),eigvall(2),alpha(2)
      REAL_T delu(2)
      REAL_T eiglex(2,2),eigrex(2,2),evalex(2)
      REAL_T spdotl(2,2),spdotr(2,2)
      REAL_T struc(2,2)
      logical ldef
      REAL_T uexp(2)
      integer k, l, ic
      integer isdef,itotdef
      REAL_T grlaml(2,2),grlamr(2,2)
      REAL_T defl, defr, def
      REAL_T dflxc1, dflxc2, dusq, spbar,du1,du2
      REAL_T ampl1, ampl2, viscpc
      REAL_T flag1, flag2,flag3, flag4, flgmin,flgmax
      REAL_T small,cdef,difmag

      itotdef = 0
      small   = 1.d-10
      cdef    = 0.1d0
      difmag  = half

      call evflux(ulft,fluxl,vtot,grav,amuw,amul,athick)

      call eigen_polymer(vtot,grav,amuw,amul,athick,ulft,eigll,eiglr,
     &   eigvall,grlaml,ldef)

      call evflux(urgt,fluxr,vtot,grav,amuw,amul,athick)

      call eigen_polymer(vtot,grav,amuw,amul,athick,urgt,eigrl,eigrr,
     &   eigvalr,grlamr,ldef)

      uexp(1)=0.5*(urgt(1)+ulft(1))
      uexp(2)=0.5*(urgt(2)+ulft(2))

      call eigen_polymer(vtot,grav,amuw,amul,athick,uexp,eiglex,
     &   eigrex,evalex,struc,ldef)

      do k=1,2
         do l=1,2
            spdotl(k,l) = grlaml(k,1)*eigrex(1,l) +
     &           grlaml(k,2)*eigrex(2,l)
            spdotr(k,l) = grlamr(k,1)*eigrex(1,l) +
     &           grlamr(k,2)*eigrex(2,l)
         end do
      end do

      do ic=1,2
         delu(ic)=urgt(ic)-ulft(ic)
         alpha(ic) = 0.
      end do

      do k=1,2
         do ic=1,2
            alpha(ic)=alpha(ic)+eiglex(ic,k)*delu(k)
         end do
      end do

      defl= dabs(alpha(1))*cdef*dabs(spdotl(1,1)-spdotl(2,1))
     &     +dabs(alpha(2))*cdef*dabs(spdotl(1,2)-spdotl(2,2))
     &     -dabs(eigvall(1)-eigvall(2))

      defr= dabs(alpha(1))*cdef*dabs(spdotr(1,1)-spdotr(2,1))
     &     +dabs(alpha(2))*cdef*dabs(spdotr(1,2)-spdotr(2,2))
     &     -dabs(eigvalr(1)-eigvalr(2))

      def = dmax1(defl,defr)

      if (def .gt. 0 .or. ldef)then
         isdef=0
      else
         isdef=1
      endif
      dflxc1 = fluxr(1) - fluxl(1)
      dflxc2 = fluxr(2) - fluxl(2)
      dusq = delu(1)**2 + delu(2)**2
      
      spbar = (delu(1)*dflxc1 + delu(2)*dflxc2)/dmax1(dusq,small)
      
      lleft(1)  = eigvall(1)
      lright(1) = eigvalr(1)
      lleft(2)  = eigvall(2)
      lright(2) = eigvalr(2)
      
      kleft(1)  = spdotl(1,1)
      kright(1) = spdotr(1,1)
      kleft(2)  = spdotl(2,2)
      kright(2) = spdotr(2,2)
      viscpc = zero
      du1=delu(1)
      du2=delu(2)

      ampl1 = alpha(1)
      ampl2 = alpha(2)

      if (isdef.eq.0)then
         flag1= dmin1(lleft(1),lright(1))
         flag2= dmax1(lleft(1),lright(1))
         flag3= dmin1(lleft(2),lright(2))
         flag4= dmax1(lleft(2),lright(2))
         flgmin=dmin1(flag1,flag3)
         flgmax=dmax1(flag2,flag4)
         
         if (flgmin*flgmax .lt. zero ) then
            viscpc = dmax1(dabs(flgmin),dabs(flgmax))
            viscpc = viscpc*difmag
            ampl1 = zero
            ampl2 = zero
         else
            ampl1 = cvmgp(one,zero,dusq - small)
            ampl2 = zero
         endif
      endif

      if (isdef.eq. 0)then
         
         eigrex(1,1) = cvmgp(du1,one,dusq-small)
         eigrex(2,1) = cvmgp(du2,zero,dusq-small)
         eigrex(1,2) = zero
         eigrex(2,2) = zero
         lleft(1)  = dmax1(eigvall(2),eigvall(1))
         lright(1) = dmin1(eigvalr(2),eigvalr(1))
         lleft(2)  = zero
         lright(2) = zero
         kleft(1)  = zero
         kleft(2)  = zero
         kright(1) = zero
         kright(2) = zero

      endif

      call flxinc(ampl1,lleft(1),lright(1),kleft(1),kright(1),
     &     spbar,da(1))

      call flxinc(ampl2,lleft(2),lright(2),kleft(2),kright(2),
     &     spbar,da(2))

      fluxg(1) = cvmgp(fluxl(1),fluxr(1),spbar)
      fluxg(2) = cvmgp(fluxl(2),fluxr(2),spbar)

      fluxg(1) = fluxg(1) + da(1)*eigrex(1,1) + da(2)*eigrex(1,2)
      fluxg(2) = fluxg(2) + da(1)*eigrex(2,1) + da(2)*eigrex(2,2)

      fluxg(1) = fluxg(1) - du1*viscpc
      fluxg(2) = fluxg(2) - du2*viscpc

      do  k = 1,2
        fluxsv(k) = fluxg(k)
      end do
      
      end
      
      subroutine cubicg(al,ar,kl,kr,amp,smalsp,dq,dqtot,numbr)
c
c     ----------------------------------------------------------
c     This subroutine fits a cubic spline for the local values of the
c     wave speed and its directional derivatives on the left and right
c     sides of the wave for the purpose of finding approximate values
c     for the sonic points (if any) on the wave curve. We assume that there
c     are no more than three sonic points possible in any wave family,
c     i. e., that the lack of convexity is no worse than that of the
c     Buckley-Leverett equations with gravity.
c
c     Inputs:
c     al,ar - wave speeds on the left and right sides.
c     kl,kr - da/du on the left and right side.
c     amp - amplitude of the wave connecting the left and right states,
c     amp .ge. 0 .
c     smalsp - a floor value for the wave speed, used to protect divides.
c
c     Outputs:
c
c     tl, tr = array containing the locations along the wave curve of
c     the left and right endpoints on intervals contributing to the
c     Engquist - Osher flux. 0 < tl(k) < tr(k) < amp,  a(t) < 0
c     between tl(k) and tr(k).
c     numbr = number of intervals of the above form.

c     NEW VERSION (PC - 12/11/86):
c     Outputs:
c     dqtot = difference of the indefinite integral of the cubic along
c     the wave curve
c     dq = array containing the differences of the indefinite integral
c     between sonic points along branches where the wave speed is negative.

c     The cubic spline formula for a is given by
c     a(t) = ar*t + al*(1-t) + t*(1-t)*(dr*t + dl*(1-t)) , 0 < t < 1
c     we find the locations of the extrema of this function (if any), by
c     differentiating wrt t and solving the resulting quadratic.For the
c     purposes of locating the sonic points, we assume that the wave speed
c     is piecewise linear with breaks at the endpoints and at any extrema
c     for which 0 < t < 1 . The sonic points are then located, the signs
c     of the local divergence of the characteristics at those points
c     determined, and we exit the routine.

c     First, the coefficients of the quadratic, and the determination of its
c     zeros.
c     ----------------------------------------------------------
c
      REAL_T al,ar,kl,kr,amp,smalsp
      REAL_T dq(4),dqtot
      integer numbr

       REAL_T tl(4),tr(4),aext(4),text(4),asnc(5),tsnc(5)
      REAL_T dela,dl,dr,deld,aa,bb,cc,disc
      REAL_T root1,root2, temp
      REAL_T delt, frac, npar

      integer kk, k, ktot

      al = cvmgp(al,-smalsp,abs(al) - smalsp)
      dela = (ar - al)
      dl = - dela + kl*amp
      dr = dela - kr*amp
      deld = dr - dl
      aa = -3*deld
      bb = 2*(deld - dl)
      cc = dela + dl

      disc = bb**2 - 4*aa*cc
      kk = 1
      if ((disc.lt.zero).or.(dabs(aa).lt.smalsp)) then
         text(1) = zero
         text(2) = one
         aext(1) = al
         aext(2) = ar
      else
         text(1) = zero
         aext(1) = al
         kk = 1

         disc = dsqrt(disc)
         root1 = (-bb - disc)*half/aa
         root2 = (-bb + disc)*half/aa

         if (root1.gt.root2) then
            temp = root1
            root1 = root2
            root2 = temp
         endif

         if ((root1.gt.zero).and.(root1.lt.one)) then
            kk = kk + 1
            text(kk) = root1
            aext(kk) = afcn(root1,al,ar,dl,dr)
         endif

         if ((root2.gt.zero).and.(root2.lt.one)) then
            kk = kk + 1
            text(kk) = root2
            aext(kk) = afcn(root2,al,ar,dl,dr)
            
         endif

         kk = kk + 1
         text(kk) = one
         aext(kk) = ar

      endif
      ktot = kk - 1
      kk = 1
      asnc(1) = al
      tsnc(1) = zero

c     Now, the calculation of the sonic points using the piecewise linear
c     fit to a.

      do  k=1,ktot
         if (aext(k)*aext(k+1).lt.zero) then
            kk = kk + 1
            delt = text(k+1) - text(k)
            frac = dabs(aext(k))/dmax1(dabs(aext(k)) + dabs(aext(k+1)),smalsp)
            tsnc(kk) = frac*delt + text(k)
            asnc(kk) = aext(k) + (aext(k+1) - aext(k))*frac

         endif
      end do

      kk = kk + 1
      asnc(kk) = ar
      tsnc(kk) = one

c     finally, the calculation of the negative branches for the E-O flux
c     integration.

      dqtot = qfcn(one,al,ar,dl,dr) - qfcn(zero,al,ar,dl,dr)

      if (kk.eq.2) then
         if (al.gt.zero) then
            tl(1) = zero
            tr(1) = zero
            numbr = 0
            dq(1) = zero
         else
            tl(1) = zero
            tr(1) = amp
            numbr = 1
            dq(1) = dqtot
         endif
      else
         npar = -1
         if (al.gt.zero) npar = 1
         numbr = 0

         do k = 2,kk
            npar = -npar
            if (npar.eq.1) then
               numbr = numbr + 1
               tl(numbr) = tsnc(k-1)*amp
               tr(numbr) = tsnc(k)*amp
               dq(numbr) = (qfcn(tsnc(k),al,ar,dl,dr) -
     &              qfcn(tsnc(k-1),al,ar,dl,dr))

            endif
         end do
      end if

      contains

        function afcn(tt,aal,aar,ddl,ddr) result(r)
        REAL_T tt,aal,aar,ddl,ddr,r
        r = aal + tt*((aar-aal) + (one-tt)*(ddr*tt + ddl*(one-tt)))
        end function afcn

        function qfcn(t,aal,aar,ddl,ddr) result(r)
        REAL_T t,aal,aar,ddl,ddr,r
        r = t*(aal*(one-half*t) + t*half*(aar+ddl)
     &       - t*t*(third*(2*ddl-ddr) + t*.25d0*(ddr-ddl)))
     &     
        end function qfcn

      end

      subroutine evflux(u,flux,vt,grav,amuw,amul,athick)
c
c     ----------------------------------------------------------
c     evaluate flux 
c     ----------------------------------------------------------
c
      REAL_T u(2)
      REAL_T flux(2)
      REAL_T vt, grav, amuw, amul, athick
      
      REAL_T s, c
      REAL_T amua, akrl, akra, amobl, amoba, amobt, v

      s=u(1)
      c=zero
      if (s.gt.zero) c=u(2)/s
      amua=amuw+athick*c
      akrl=(one-s)**2
      akra=s**2
      amobl=akrl/amul
      amoba=akra/amua
      amobt=amobl+amoba
      v=amoba*(vt+grav*amobl)/amobt

      flux(1)=v
      flux(2)=c*v

      end

      subroutine flxinc(amp,al,ar,kl,kr,spav,da)
      
      REAL_T amp,al,ar,kl,kr,spav,da

      REAL_T dq(4),absamp,smalsp,dqtot
      integer k,numbr
      absamp = dabs(amp)
      smalsp = 1.d-10

      call cubicg(al,ar,kl,kr,absamp,smalsp,dq,dqtot,numbr)
      da = zero
      if (numbr.gt.0) then
         do k = 1,numbr
            da = da + dq(k)
         end do
      end if

      if (spav.gt.zero) then
         da = amp*da
      else
         da = amp*(da - dqtot)
      endif

      end
