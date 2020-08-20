#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "Prob_PM_F.H" 
#include "ArrayLim.H"

#include "RockManager_F.H"


c ::: -----------------------------------------------------------
c ::: Initiate data in probdata.H 
c ::: -----------------------------------------------------------
      subroutine FORT_INITPARAMS(ncomps_in,nphases_in,model_in,rhoval_in,
     &     mu_in,phaseidx_in,gravity_in,gravity_dir_in)

      implicit none

#include "probdata.H"

      integer ncomps_in,nphases_in,model_in
      integer phaseidx_in(ncomps_in)
      REAL_T  rhoval_in(ncomps_in)
      REAL_T  mu_in(ncomps_in)
      REAL_T  gravity_in
      integer gravity_dir_in

      integer n
  
      ncomps   = ncomps_in
      nphases  = nphases_in
      model    = model_in
      do n = 1,ncomps
           rhoval(n) = rhoval_in(n)
            muval(n) = mu_in(n) 
         phaseidx(n) = phaseidx_in(n)
      end do
      gravity      = gravity_in
      gravity_dir  = gravity_dir_in
      end

c ::: -----------------------------------------------------------
c ::: Initiate tracers data in probdata.H 
c ::: -----------------------------------------------------------
      subroutine FORT_TCRPARAMS(ntracers_in)

      implicit none

#include "probdata.H"      

      integer ntracers_in
  
      ntracers = ntracers_in

      end

c ::: -----------------------------------------------------------
c ::: Initiate tracers data in probdata.H 
c ::: -----------------------------------------------------------
      subroutine FORT_AUXPARAMS(naux_in)

      implicit none

#include "probdata.H"      

      integer naux_in
  
      naux = naux_in

      end

c ::: -----------------------------------------------------------
c ::: Initiate data in probdata.H for get_kr function
c ::: -----------------------------------------------------------
      subroutine FORT_KR_INIT(sThresh_in,use_shifted_kr_in)

      implicit none

#include "probdata.H"      

      REAL_T sThresh_in
      integer use_shifted_kr_in

      sThresh = sThresh_in
      use_shifted_kr = use_shifted_kr_in
      end


c ::: -----------------------------------------------------------
c     Some analytical functions for kappa 
c       0: constant
c       1: layered
c       2: random
c       3: smooth
c ::: -----------------------------------------------------------
      subroutine getkappa(kappa_val,x,y,z,mediumtype,mediumval,nlayer)

      implicit none

      integer mediumtype,nlayer
      REAL_T x,y,z
      REAL_T mediumval(nlayer)
      REAL_T rdnum
      REAL_T ydiv(nlayer+1)
      REAL_T kappa_val
      integer i,j

      REAL_T, parameter :: alpha(4,3) = reshape(
     & (/ 0.02d0, 0.052d0, 0.063d0, 0.575d0,
     &    0.11d0, 0.33d0, 0.145d0, 0.042d0,
     &    0.07d0, 0.33d0, 0.0246d0, 0.01d0 /), (/4,3/) )

      REAL_T, parameter :: phi(4,3) = reshape(
     & (/ 0.356d0, 0.434d0, 0.744d0, 0.875d0,
     &    0.242d0, 0.222d0, 0.133d0, 0.248d0,
     &    0.125d0, 0.534d0, 0.988d0, 0.623d0 /), (/4,3/) )

#if 0
      alpha(1,1) = 0.02d0
      alpha(1,2) = 0.11d0
      alpha(1,3) = 0.07d0
      alpha(2,1) = 0.052d0
      alpha(2,2) = 0.33d0
      alpha(2,3) = 0.0956d0
      alpha(3,1) = 0.063d0
      alpha(3,2) = 0.145d0
      alpha(3,3) = 0.0246d0
      alpha(4,1) = 0.575d0
      alpha(4,2) = 0.042d0
      alpha(4,3) = 0.01d0

      phi(1,1) = 0.356d0
      phi(1,2) = 0.242d0
      phi(1,3) = 0.125d0
      phi(2,1) = 0.434d0
      phi(2,2) = 0.222d0
      phi(2,3) = 0.534d0
      phi(3,1) = 0.744d0
      phi(3,2) = 0.133d0
      phi(3,3) = 0.988d0
      phi(4,1) = 0.875d0
      phi(4,2) = 0.248d0
      phi(4,3) = 0.623d0
#endif

      if (mediumtype .eq. 0) then
c
c         constant kappa
c
          kappa_val = mediumval(1)


      elseif (mediumtype .eq. 1) then
c
c         layered kappa
c 

c         layered in the y-direction
         ydiv(nlayer+1) = one
         do i = 1,nlayer
            ydiv(i) = (dble(i-1))/dble(nlayer)
         end do 

         do i = 1,nlayer
            if (y .ge. ydiv(i) .and. y .lt. ydiv(i+1)) then
               kappa_val = mediumval(i)
            end if
         end do

         if (y .lt. ydiv(1)) then
            kappa_val = mediumval(1)
         else if (y .gt. ydiv(nlayer+1)) then
            kappa_val = mediumval(nlayer)
         end if

c        layered in the x-direction
c$$$         ydiv(nlayer+1) = 1.0
c$$$         do i = 1,nlayer
c$$$            ydiv(i) = (dble(i-1))/dble(nlayer)
c$$$         end do 
c$$$
c$$$         do i = 1,nlayer
c$$$            if (x .ge. ydiv(i) .and. x .lt. ydiv(i+1)) then
c$$$               kappa_val = mediumval(i)
c$$$            end if
c$$$         end do
c$$$
c$$$         if (x .lt. ydiv(1)) then
c$$$            kappa_val = mediumval(1)
c$$$         else if (x .gt. ydiv(nlayer+1)) then
c$$$            kappa_val = mediumval(nlayer)
c$$$         end if

      elseif  (mediumtype .eq. 2) then
c
c        random kappa
c
         kappa_val = mediumval(1)
         do i = 1,4
            do j= 1,3
               kappa_val = kappa_val + alpha(i,j)*500.d-8*
     &              SIN(Pi*dble(j+5)*(x + phi(i,j)))*SIN(Pi*dble(j+5)*(y + phi(i,j)))
            end do
	 end do
      elseif  (mediumtype .eq. 3) then      
c
c        smooth-varying kappa
c
         kappa_val = mediumval(1) * (1+ SIN(Pi*y)**4 + SIN(Pi*min(half,2*x)))

      elseif (mediumtype .eq. 4) then
         call blutilrand(rdnum)
         kappa_val = mediumval(1) * (1 + 1.d-8*rdnum);

      else
         call bl_abort('medium type not defined')
      endif

      end


c ::: ----------------------------------------------------------
c ::: calculate viscosity mu for single-phase 2-component case
c ::: ----------------------------------------------------------
      subroutine get_mu_1phase_2comp (s1,s2,mu)

      implicit none

      REAL_T  s1,s2,mu

      REAL_T  s(2)
      REAL_T  M, Mquart, c
      
      M      = one
      Mquart = M**0.25d0    
      s(1) = max(s1,zero)
      s(2) = max(s2,zero)
      c  = s(1)/(s(1) + s(2))
      c  = min(max(c,zero),one)
      mu = one/(one-c+Mquart*c)**4

      end


c ::: -----------------------------------------------------------
c ::: Determine viscosity for components for the different cases.
c ::: -----------------------------------------------------------
      subroutine get_mu(mu,s,muval,compcount,compidx,nphases,ncomps)

      implicit none

      integer nphases,ncomps
      integer compcount(nphases), compidx(nphases,ncomps)
      integer k,nc
      REAL_T  mu(nphases), muval(ncomps), s(ncomps)
      REAL_T  dummy, c

      if (nphases .eq. 1) then
         if (compcount(1) .eq. 1) then
            mu(1) = muval(compidx(1,1))
         else 
            dummy = zero
            do nc = 2,compcount(1) 
               dummy = dummy + s(compidx(1,nc))
            end do 
            call get_mu_1phase_2comp(s(compidx(1,1)),dummy,mu(1))
         end if

      elseif (nphases .eq. 2) then
         do k = 1,nphases
            mu(k) =muval(compidx(k,1))
            if (compcount(k) .eq. 2) then
               c = s(compidx(k,1))+s(compidx(k,2))
               if (c .gt. 1.0d-10) then
                  c = s(compidx(k,2))/c
                  mu(k) = mu(k) + c*muval(compidx(k,2))
               end if
            end if
         end do

      elseif (nphases .gt. 2) then
         do k = 1,nphases
            mu(k) = muval(compidx(k,1))
         end do
      end if                 
      
      end

c ::: -----------------------------------------------------------
c ::: Compute the capillary pressure.
c     CPL type: 
c       0 => None
c       1 => Leverett
c       2 => Brooks Corey
c       3 => vanGenuchten
c ::: -----------------------------------------------------------
      subroutine get_cpl(pc,s,rho,phi,kappa,sigma,ncomps,
     &     sres,cpl_type,cpl_coef)
      implicit none

      integer ncomps
      integer cpl_type
      REAL_T  pc
      REAL_T  s(ncomps)
      REAL_T  rho(ncomps)
      REAL_T  phi, kappa
      REAL_T  sigma
      REAL_T  sres,cpl_coef

      REAL_T  js,m
      REAL_T  sat(ncomps)

      call bl_abort('Should not be in get_cpl')

      call get_st(sat,s,rho,sres,ncomps)

      if (cpl_type .eq. 0) then 
         pc = zero

      elseif (cpl_type .eq. 1) then
         js = exp(5*(1-s(2)/rho(2)))
         pc = sigma*js*sqrt(phi/kappa)

      elseif (cpl_type .eq. 2) then
c        lambda = 2
         pc = one/(sigma*sqrt(sat(1)))

      elseif (cpl_type .eq. 3) then
c        sigma is user-specified alpha

c        p_c near zero sat approximated by linear function
         
c         if (sat(1) .lt. 0.05d0) then
c            pc = 4d1 - 4d2*sat(1)
c         else
c            pc = sqrt(one/(sat(1)*sat(1)) - one)/sigma
c         end if
         sat(1) = max(1.d-4,sat(1))
         m = cpl_coef            
         pc = (((sat(1)**(-one/m))-one)**(one-m))/sigma
      end if

      end

c ::: -----------------------------------------------------------
c ::: Compute the inverse of the capillary pressure.
c     CPL type: 
c       0 => None
c       1 => Leverett
c       2 => Corey
c       3 => vanGenuchten
c ::: -----------------------------------------------------------
      subroutine get_invcpl(s,pc,rho,phi,kappa,sigma,ncomps,
     &     sres,cpl_type,cpl_coef)

      implicit none

      integer ncomps
      integer cpl_type
      REAL_T  pc
      REAL_T  s(ncomps)
      REAL_T  rho(ncomps)
      REAL_T  phi, kappa
      REAL_T  sigma
      REAL_T  sres, cpl_coef

      REAL_T  m
      REAL_T  sat(ncomps)

      call bl_abort('Should not be in get_invcpl')

      if (cpl_type .eq. 0) then 
         sat = zero

      elseif (cpl_type .eq. 1) then
         if (ncomps .ne. 2) then
            print *,'cpl_type=1 requires ncomps=2'
            call bl_pd_abort()
         endif
         if (pc > zero) then
            sat(2) = one - 0.2d0*log(pc*sqrt(kappa/phi)/sigma)
         else
            sat(2) = zero
         end if
         if (ncomps > 1) then
            sat(1) = one - sat(2)
         end if

      elseif (cpl_type .eq. 2) then
c        lambda = 2
         if (pc > zero) then
            sat(1) = one/(sigma*pc)**2
         else
            sat(1) = one
         end if
         if (ncomps > 1) then
            sat(2) = one-sat(1)
         end if

      elseif (cpl_type .eq. 3) then
         
c        - sigma is user-specified alpha
c        - p_c near zero sat approximated by linear function
         
c         if (pc .gt. 20) then
c            sat(1) = zero
c            sat(2) = one
c         else
c            sat(1) = pc*sigma
c            sat(1) = sqrt(one/(one+sat(1)*sat(1)))
c         end if

         m = cpl_coef
         if (pc > zero) then
            sat(1) = (one + (pc*sigma)**(one/(one-m)))**(-m)
         else
            sat(1) = one
         end if
         if (ncomps > 1) then
            sat(2) = one-sat(1)
         end if

      end if

      if (cpl_type .ge. 0) then
         call get_invst(s,sat,rho,sres,ncomps)
      endif

      end

c ::: -----------------------------------------------------------
c ::: Compute the derivative of capillary pressure w.r.t. s
c     CPL type: 
c       0 => None
c       1 => Leverett
c       2 => Corey
c       3 => vanGenuchten
c ::: -----------------------------------------------------------
      subroutine get_dcplds(dpc,s,rho,phi,kappa,sigma,ncomps,
     &     sres,cpl_type,cpl_coef)

      implicit none

      integer ncomps
      integer cpl_type
      REAL_T  dpc
      REAL_T  s(ncomps)
      REAL_T  rho(ncomps)
      REAL_T  phi, kappa
      REAL_T  sigma
      REAL_T  sres,cpl_coef

      REAL_T  js,m,n
      REAL_T  sat(ncomps)

      call bl_abort('Should not bin in get_dcplds')

      call get_st(sat,s,rho,sres,ncomps)

      if (cpl_type .eq. 0) then 
         dpc = zero

      elseif (cpl_type .eq. 1) then
         js  = exp(5*(1-sat(2)))
         dpc = 0

      elseif (cpl_type .eq. 2) then

c        Brooks Corey model with lambda = 2

         dpc = -half/(sigma*sat(1)**(-1.5d0));

      elseif (cpl_type .eq. 3) then

c        van Genuchten model 
c        - sigma is user-specified alpha
         
         m = cpl_coef
         n = -one/m

c        rewrote to speed it up a little 
c        dpc = (one-m)*(sat(1)**n-one)**(-m)*n*(sat(1)**(n-one))/sigma
         sat(1) = max(1.d-3,sat(1))
         dpc = sat(1)**n
         dpc = n*(one-m)*(dpc-one)**(-m)*dpc/(sat(1)*sigma)
            
      end if

      end

c ::: -----------------------------------------------------------
c ::: Compute the relative permeability
c     kr type: 
c       0 => None
c       1 => linear
c       2 => quadratic
c       3 => vanGenuchten
c     Applicable only for two-phase at the moment
c ::: -----------------------------------------------------------
      subroutine get_kr(kr,s,rho,nc,kr_type,kr_coef,sres)

      implicit none

#include "probdata.H" 

      integer nc
      integer kr_type
      REAL_T  kr(nc)
      REAL_T  s(nc)
      REAL_T  rho(nc)
      REAL_T  kr_coef
      REAL_T  sres

      integer i
      REAL_T  alpha, alpha_scale
      REAL_T  m,n,twom
      REAL_T  sat(nc),satl

      REAL_T k_st, k_st1, k_sat, msat, mt, ds, dels, dels1, vgk, sa

      vgk(m,sa) = sqrt(sa)*(one-(one-sa**(one/m))**m)**2.d0

      call bl_abort('Should not be in get_kr')

      call get_st(sat,s,rho,sres,nc)

      if (kr_type .eq. 0) then 
         kr = half
      elseif (kr_type .eq. 1) then
         do i = 1,nc
            kr(i) = sat(i)
         end do
      elseif (kr_type .eq. 2) then
         do i = 1,nc
            kr(i) = sat(i)*sat(i)
         end do
      elseif (kr_type .eq. 3) then

         if (use_shifted_kr .ne. 1) then

            if (sThresh .lt. 0) then
c     The standard van Genuchten parameterization, no fits/shifts/etc
               m = kr_coef
               n = 1.d0/m
               kr(1) = vgk(m,sat(1))
            else
c
c     Try kr-smoothing: If se<sThresh, fit
c     cubic Hermite polynomial matching kr(s=1), kr(s=sThresh),
c     dkr/ds(sThresh) and setting dkr/ds(s=1) = constant
c     
               m = kr_coef
               if (sat(1) .le. sThresh) then
                  kr(1) = vgk(m,sat(1))
               else
                  k_st = vgk(m,sThresh)
                  ds = (one - sThresh)*1.d-3
                  k_st1 = vgk(m,sThresh+ds)
                  msat = (k_st - one)/(one - sThresh)
                  mt = -(k_st1-k_st)/ds
                  dels = one - sat(1)
                  dels1 = one - sThresh
                  kr(1) = one + dels*dels * msat/dels1 + 
     &                 dels*dels * (dels-dels1) * (mt-2*msat) 
     &                 / (dels1*dels1)
               endif
               
            endif

c     Shifted eval: Georges original implementation, activated when use_shifted_kr=1
         else
c
c     Van Genuchten-Gualem function
c     - typical well-behaved value is m = 0.5, while for ASCEM:   m = .28 
c     - the gradient at s = 1 is inf.  We modify the original function by 
c       introducing a factor alpha to obtain finite gradient at s = 1.  
c       It is then scaled to obtain kr = 1 at s = 1.
c
            m      = kr_coef
            n      = 1.d0/m
            twom   = 2.d0*m
            alpha  = 0.97d0
            alpha_scale = dsqrt(alpha)*(one-(one-alpha**n)**m)**2
            kr(1) = sqrt(alpha*sat(1))*(one-(one-(alpha*sat(1))**n)**m)**2
            kr(1) = kr(1)/alpha_scale
         endif

c     # no scaling for the gas permeability function, add checks in the derivatives
c     # to avoid blow-up
c     # In order to write the gas permeability function in terms of s_eff, we write
c     # the gas saturation in terms of s_eff.  This way we know that when get_kr
c     # is called in other functions, we are using consistent saturations. 
         
         if (nc > 1) then
c            satl = sat(1)*(one-sres) + sres
c            if (satl .eq. one) then
c               kr(2) = zero
c            else
c               kr(2) = ((one-satl)**third)*(one-satl**n)**twom
c            end if
            kr(2) = one - kr(1)
         end if
      end if

      end

c ::: ---------------------------------------------------------------
c ::: Compute the saturation, taking residual saturation into account
c     Applicable only for two-phase at the moment
c ::: ---------------------------------------------------------------
      subroutine get_st(st,s,rho,sres,ncomps)

      implicit none

      integer ncomps
      REAL_T  st(ncomps)
      REAL_T  s(ncomps)
      REAL_T  rho(ncomps)
      REAL_T  sres

      call bl_abort('Should not be in get_st')
      s(1) = min(rho(1),s(1))
      s(1) = max(zero,s(1))

      st(1) = (s(1)/rho(1)-sres)/(one-sres)
      st(1) = min(one,st(1))
      st(1) = max(zero,st(1))

      if (ncomps > 1) then
         st(2) = one - st(1)
      end if

      end

c ::: ---------------------------------------------------------------
c ::: Compute the inverse of saturation function.
c     Applicable only for two-phase at the moment
c ::: ---------------------------------------------------------------
      subroutine get_invst(s,st,rho,sres,ncomps)

      implicit none

      integer ncomps
      REAL_T  st(ncomps)
      REAL_T  s(ncomps)
      REAL_T  rho(ncomps)
      REAL_T  sres

      integer n

      call bl_abort('Should not be in get_invst')

      s(1) = st(1)*(one-sres)+sres
      s(1) = max(sres,s(1))
      s(1) = min(one,s(1))
      if (ncomps > 1) then
         s(2) = one-s(1)
      end if

      do n = 1,ncomps
         s(n) = s(n)*rho(n)
      end do

      end

c ::: -----------------------------------------------------------
c ::: Compute lambda = kr/mu
c     kr type: 
c       0 => None
c       1 => linear
c       2 => quadratic
c       3 => vanGenuchten
c     Applicable only for two-phase at the moment
c ::: -----------------------------------------------------------
      subroutine get_ld(kr,s,rho,muval,ncomps,kr_type,kr_coef,sres)

      implicit none

      integer ncomps
      integer kr_type
      REAL_T  kr(ncomps)
      REAL_T  s(ncomps)
      REAL_T  rho(ncomps)
      REAL_T  muval(ncomps)
      REAL_T  kr_coef
      REAL_T  sres

      integer i

      call bl_abort('Should not be in get_ld')

      if (kr_type .ge. 0) then
         call get_kr(kr,s,rho,ncomps,kr_type,kr_coef,sres)
         
         do i = 1,ncomps
            kr(i) = kr(i)/muval(i)
         end do
      endif
      end

c ::: -----------------------------------------------------------
c ::: Compute the first derivatives of ld
c     kr type: 
c       0 => None
c       1 => linear
c       2 => quadratic : d(ld1), d(ld1ld2), d(1/ldT)
c       3 => vanGenuchten : d(ld1), d(ld2)
c     Applicable only for two-phase at the moment
c ::: -----------------------------------------------------------
      subroutine get_dld(dld,s,rho,muval,ncomps,kr_type,kr_coef,sres)

      implicit none

      integer ncomps
      integer kr_type
      REAL_T  dld(3)
      REAL_T  s(ncomps)
      REAL_T  rho(ncomps)
      REAL_T  muval(ncomps)
      REAL_T  kr_coef
      REAL_T  sres,tmp,tmpa
      REAL_T  alpha, alphs, alpha_scale

      REAL_T  m
      REAL_T  sat(ncomps),satl

      call get_st(sat,s,rho,sres,ncomps)

      dld(:) = zero
      if (kr_type .eq. 1 .and. ncomps .gt. 1) then
         dld(1) = 1/muval(1)
         dld(2) = (1-2*sat(1))/(muval(1)*muval(2))
         dld(3) = sat(1)/muval(1)+sat(2)/muval(2)
         dld(3) = -(1/muval(1)-1/muval(2))/dld(3)**2

      elseif (kr_type .eq. 2 .and. ncomps .gt. 1) then
         dld(1) = 2*sat(1)/muval(1)
         dld(2) = (2*sat(1)-6*sat(1)*sat(1)+4*sat(1)**3)/(muval(1)*muval(2))
         dld(3) = sat(1)*sat(1)/muval(1)+sat(2)*sat(2)/muval(2)
         dld(3) = -2*(sat(1)/muval(1)-sat(2)/muval(2))/dld(3)**2

      elseif (kr_type .eq. 3) then
c
c     van Genuchten function
c     - please refer to description in get_kr
c     - dld(1) and dld(2) holds the derivative of lambda1 and lambda2 wrt s
c     - The derivative of the permeability function for water is ill-behaved
c       when the water saturation is 0.  This should not happen but   
c       just in case we explicitly test for this case and 
c       set the derivative values to zero.
c
          m     = kr_coef
          alpha = 0.97d0
          alphs = alpha*sat(1)
          tmpa  = one - alphs**(one/m)
          satl  = sat(1)*(one-sres) + sres
          tmp   = one - satl**(one/m)

          alpha_scale = dsqrt(alpha)*(one-(one-alpha**(one/m))**m)**2

          if (sat(1).eq.zero) then
             dld(1) = zero
             dld(2) = zero
          elseif (satl.eq.one) then
c     # If the saturation is 1, the derivative of lambda2 wrt s should be 0
c     # the functional form blows up, so we enforce that condition here
             dld(2) = zero
             dld(1) = two*(alphs**(-half+(one/m)))*(tmpa**(-one+m))*(one-tmpa**m) +
     &            (half*(one-tmpa**m)**two)/(sqrt(alphs))
             dld(1) = dld(1)/alpha_scale
          else 
             dld(1) = two*(alphs**(-half+(one/m)))*(tmpa**(-one+m))*(one-tmpa**m) +
     &            (half*(one-tmpa**m)**two)/(sqrt(alphs))
             dld(1) = dld(1)/alpha_scale
             dld(2) = -tmp**(two*m)/(3.0*(one-satl)**(two3rd)) - 
     &            two*((one-satl)**(third))
     &            *(satl**(-one+one/m))*(tmp**(-one+two*m))
             dld(2) = dld(2)*(one-sres)
          endif
          dld(2) = -alpha*dld(1)/muval(2)
          dld(1) = alpha*dld(1)/muval(1)
c          dld(2) = dld(2)/muval(2)
      end if

      end

      subroutine get_d2ld(d2ld,s,rho,muval,ncomps,kr_type,kr_coef,sres)
      implicit none
c
c     -----------------------------------------------------------
c     Compute the second derivatives of ld: 
c     kr type: 
c       0 => None
c       1 => linear
c       2 => quadratic : d2(ld1),d2(ld1ld2), d2(1/ldT)
c       3 => vanGenuchten : d2(ld1),d2(ld2)
c     Applicable only for two-phase at the moment
c     -----------------------------------------------------------
c
      integer ncomps
      integer kr_type
      REAL_T  d2ld(3)
      REAL_T  s(ncomps)
      REAL_T  rho(ncomps)
      REAL_T  muval(ncomps)
      REAL_T  sres,kr_coef

      REAL_T  m,tmp
      REAL_T  sat(ncomps),satl,lt,s1
      REAL_T  alpha,alphs,alpha_scale
      REAL_T  tmpa

      call get_st(sat,s,rho,sres,ncomps)
      satl = sat(1)*(one-sres) + sres

      if (kr_type == 0) then 
         d2ld = zero

      elseif (kr_type .eq. 1) then
         d2ld(1) = zero
         d2ld(2) = -2/(muval(1)*muval(2))
         d2ld(3) = sat(1)/muval(1)+sat(2)/muval(2)
         d2ld(3) = 2*(1/muval(1)-1/muval(2))**2/d2ld(3)**3

      elseif (kr_type .eq. 2) then
         d2ld(1) = 2/muval(1)
         d2ld(2) = (2-12*sat(1)+12*sat(1)*sat(1))/(muval(1)*muval(2))
         lt      = sat(1)*sat(1)/muval(1)+sat(2)*sat(2)/muval(2)
         d2ld(3) = -2*(sat(1)/muval(1)-sat(2)/muval(2))/lt**2
         d2ld(3) = -2*(1/muval(1)+1/muval(2))/lt**2
     .        + lt*d2ld(3)**2

      elseif (kr_type .eq. 3) then
c
c     van Genuchten function
c     - please refer to description in get_kr
c     - d2ld(1) and d2ld(2) holds 2nd derivatives of lambda1 and lambda2 wrt s
c     - The 2nd derivatives of the permeability function for water is ill-behaved
c       when the water saturation is 0.  This should not happen but   

c       just in case we explicitly test for this case and 
c       set the derivative values to zero.
c
c     ### Eventually would like to use the 'realistic' m=.28
         m     = kr_coef
         alpha = .97d0
         alphs = alpha*sat(1)
         tmpa  = (1-alphs**(1/m))
         tmp   = 1-satl**(1/m)
         s1    = satl
         
         alpha_scale = dsqrt(alpha)*(one-(one-alpha**(one/m))**m)**2
         
         if (sat(1).eq.zero) then 
            d2ld(1) = zero
            d2ld(2) = zero
         elseif(satl.eq.one) then
            d2ld(2) = zero
            d2ld(1) = 2*alphs**(-1.5d0+2/m)*(tmpa**(-2+2*m))
     .           - (2/m)*(-1+m)*alphs**(-1.5d0+(2/m))*tmpa**(-2+m)*(1-(tmpa**m))
     .           + alphs**(-1.5d0+(1/m))*tmpa**(-1+m)*(1-tmpa**m)
     .           + (-1+2/m)*alphs**(-1.5d0+1/m)*tmpa**(-1+m)*(1-tmpa**m)
     .           - .25d0*(1-tmpa**m)**2/alphs**1.5d0
            d2ld(1) = d2ld(1)/alpha_scale
         else
            d2ld(1) = 2*alphs**(-1.5d0+2/m)*(tmpa**(-2+2*m))
     .           - (2/m)*(-1+m)*alphs**(-1.5d0+(2/m))*tmpa**(-2+m)*(1-(tmpa**m))
     .           + alphs**(-1.5d0+(1/m))*tmpa**(-1+m)*(1-tmpa**m)
     .           + (-1+2/m)*alphs**(-1.5d0+1/m)*tmpa**(-1+m)*(1-tmpa**m)
     .           - .25d0*(1-tmpa**m)**2/alphs**1.5d0
            d2ld(1) = d2ld(1)/alpha_scale
            d2ld(2) = -2*tmp**(2*m)/(9*(1-s1)**(5.d0/3.d0)) 
     .           + (2/m)*((-1+2*m)*(1-s1)**(third)
     .           *(s1**(-2+(2/m)))*(tmp**(-2+2*m)))
     .           - 2*(-1+(1/m))*((1-s1)**(third))
     .           *s1**(-2+(1/m))*(tmp**(-1+2*m))
     .           + 4*s1**(-1+(1/m))*tmp**(-1+2*m)/(3*(1-s1)**(two3rd))
            d2ld(2) = d2ld(2)*(one-sres)**2
c         	write(*,*) 'second derivative check', d2ld(1),d2ld(2),sat(1)
         endif
         d2ld(2) = -d2ld(1)*alpha**2/muval(2)
         d2ld(1) = d2ld(1)*alpha**2/muval(1)
c         d2ld(2) = d2ld(2)/muval(2)
c         write(*,*) 'second derivative check', d2ld(1),d2ld(2),sat(1)
      end if

      end

      subroutine FORT_FIND_INV_FLUX(sol,vel,nc,vtot,gstar,muval,
     &     ncomps,rmID,matID,time)
c     given a velocity of a phase, find the corresponding saturation 

      implicit none

      integer nc,ncomps,rmID,matID,n
      REAL_T time
      REAL_T sol,muval(ncomps)
      REAL_T vel,vtot,gstar,Sr,kr

      REAL_T ld(2),stmp,rtmp(2),mtmp(2)
      REAL_T s1,s2,f1,f2,ftmp
      REAL_T, parameter :: eps = 1.0d-14

      if (ncomps < 2) then
         mtmp(1) = muval(1)
         mtmp(2) = 0.018
         call bl_abort('FORT_FIND_INV_FLUX only works for 2-phase flow!')
      end if

      call ROCK_MANAGER_RESIDSAT(rmID,matID,time,Sr,1)

      s1 = Sr
      s2 = Sr
      f2 = 1e20
      rtmp = one

c     first find the minimum for vT = 0
      do n = 1,100
         s1 = s1 + 0.01
         call ROCK_MANAGER_RELPERM(rmID,s1,matID,time,kr,1)
         ld(1) = kr / muval(1)
         ld(2) = (1.d0 - kr) / muval(2)
         if (nc .eq. 1) then
            f1 = ld(1)/(ld(1)+ld(2))*(vtot + ld(2)*gstar) 
         else
            f1 = ld(2)/(ld(1)+ld(2))*(vtot - ld(1)*gstar) 
         end if     
         if (f1 < f2) then
            s2 = s1
            f2 = f1
         end if
      end do

      s1 = Sr + eps
      if (abs(vel) .lt. 1e-15) then
         sol  = zero
      else
      call ROCK_MANAGER_RELPERM(rmID,s1,matID,time,kr,1)
      ld(1) = kr / muval(1)
      ld(2) = (1.d0 - kr) / muval(2)
      if (nc .eq. 1) then
         f1 = vel - ld(1)/(ld(1)+ld(2))*(vtot + ld(2)*gstar) 
      else
         f1 = vel - ld(2)/(ld(1)+ld(2))*(vtot - ld(1)*gstar) 
      end if     
 60   s1 = s2
      call ROCK_MANAGER_RELPERM(rmID,s1,matID,time,kr,1)
      ld(1) = kr / muval(1)
      ld(2) = (1.d0 - kr) / muval(2)
      if (nc .eq. 1) then
         f2 = vel - ld(1)/(ld(1)+ld(2))*(vtot + ld(2)*gstar) 
      else
         f2 = vel - ld(2)/(ld(1)+ld(2))*(vtot - ld(1)*gstar) 
      end if
      if (sign(one,f1)-sign(one,f2) .ne. zero) then
         do n = 1,100
            stmp = half*(s1+s2)
            call ROCK_MANAGER_RELPERM(rmID,stmp,matID,time,kr,1)
            ld(1) = kr / muval(1)
            ld(2) = (1.d0 - kr) / muval(2)
            if (nc .eq. 1) then
               ftmp = vel - ld(1)/(ld(1)+ld(2))*(vtot + ld(2)*gstar) 
            else
               ftmp = vel - ld(2)/(ld(1)+ld(2))*(vtot - ld(1)*gstar) 
            end if

            if (ftmp*f1 .lt. zero) then
               f2 = ftmp
               s2 = stmp
            elseif (ftmp*f2 .lt. zero) then
               f1 = ftmp
               s1 = stmp
            endif

            if (dabs(s1-s2) .lt. eps) then
               goto 70
            endif
         enddo
      endif

 70   sol = half*(s1+s2)
      end if
      end

      subroutine FORT_FIND_INV_RFLUX(sol,vel,gstar,muval,
     &     ncomps,rmID,matID,time)
c     given the inflow velocity, find the saturation based on richard's equation 

      implicit none

      integer ncomps,rmID,matID,n
      REAL_T  sol,muval(ncomps)
      REAL_T  vel,gstar,time,Sr,kr

      REAL_T s1,s2,stmp,f1,f2,ftmp
      REAL_T, parameter :: eps = 1.0d-14

      call ROCK_MANAGER_RESIDSAT(rmID,matID,time,Sr,1)

      s1 = Sr
      s2 = Sr
      f2 = 1e20

      s1 = Sr
      call ROCK_MANAGER_RELPERM(rmID,s1,matID,time,kr,1)
      f1 = vel - kr*gstar/muval(1)

      s2 = one
      call ROCK_MANAGER_RELPERM(rmID,s2,matID,time,kr,1)
      f2 = vel - kr*gstar/muval(1)

      if (sign(one,f1)-sign(one,f2) .ne. zero) then
         do n = 1,100
            stmp = half*(s1+s2)
            call ROCK_MANAGER_RELPERM(rmID,stmp,matID,time,kr,1)
            ftmp = vel - kr*gstar/muval(1)

            if (ftmp*f1 .lt. zero) then
               f2 = ftmp
               s2 = stmp
            elseif (ftmp*f2 .lt. zero) then
               f1 = ftmp
               s1 = stmp
            endif

            if (dabs(s1-s2) .lt. eps) then
               goto 80
            endif
         enddo
      endif

 80   sol = half*(s1+s2)

      end

c ::: ----------------------------------------------------------
c ::: Get component index and component count
c ::: ----------------------------------------------------------
      subroutine get_compidx(compcount,compidx,phaseidx,nphases,ncomps)
      
      integer nphases,ncomps
      integer compcount(nphases), compidx(nphases,ncomps)
      integer phaseidx(ncomps)
      integer k,l

      do k = 1,nphases
         compcount(k) = 0
         do l = 1,ncomps
            compidx(k,l) = 0
         end do
      end do

      do k = 1,ncomps
         compcount(phaseidx(k)+1) = compcount(phaseidx(k)+1) + 1
         compidx(phaseidx(k)+1,compcount(phaseidx(k)+1)) = k
      end do

      end

