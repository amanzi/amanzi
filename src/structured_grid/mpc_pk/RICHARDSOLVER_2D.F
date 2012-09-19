#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "RICHARDSOLVER_F.H"
#include "ArrayLim.H"

#define SDIM 2

      subroutine FORT_RS_GXPA(e,DIMS(e),
     &     c,DIMS(c),lo,hi,dx,a,dir,nc)
      implicit none      
      integer dir, nc
      integer DIMDEC(e)
      integer DIMDEC(c)
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  e(DIMV(e),nc)
      REAL_T  c(DIMV(c),nc)
      REAL_T  a(nc)
      
      integer i, j, n
      REAL_T tmult

      if (dir.eq.0) then
         tmult = 1.d0/dx(dir+1)
         do n=1,nc
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  e(i,j,n) = -( tmult*(c(i,j,n) - c(i-1,j,n)) + a(n))
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         tmult = 1.d0/dx(dir+1)
         do n=1,nc
            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  e(i,j,n) = -( tmult*(c(i,j,n) - c(i,j-1,n)) + a(n))
               enddo
            enddo
         enddo
      else
         print *, 'FORT_RS_GXPA::bad dir: ',dir
         call bl_pd_abort()
      endif
      end     


      subroutine FORT_RS_XMULTYZ(X,DIMS(X),Y,DIMS(Y),Z,DIMS(Z),lo,hi,nc)
      implicit none      
      integer DIMDEC(X)
      integer DIMDEC(Y)
      integer DIMDEC(Z)
      integer lo(SDIM),hi(SDIM)
      REAL_T  X(DIMV(X),nc)
      REAL_T  Y(DIMV(Y),nc)
      REAL_T  Z(DIMV(Z),nc)
      integer nc
      
      integer i, j, n

      do n=1,nc
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               X(i,j,n) = X(i,j,n) * Y(i,j,n) * Z(i,j,n)
            enddo
         enddo
      enddo
      end     

      subroutine FORT_RS_CTE_UPW(coef,DIMS(coef),
     &     lam,DIMS(lam),vel,DIMS(vel),lo,hi,dir,nc,dirbclo,dirbchi,do_upwind)
      implicit none      
      integer DIMDEC(coef)
      integer DIMDEC(lam)
      integer DIMDEC(vel)
      integer lo(SDIM),hi(SDIM)
      REAL_T  coef(DIMV(coef),0:nc-1)
      REAL_T  lam(DIMV(lam),0:nc-1)
      REAL_T  vel(DIMV(vel))
      REAL_T  mult
      integer dir, nc, dirbclo, dirbchi, do_upwind
      integer i, j, n

      if (dir.eq.0) then
         if (do_upwind.eq.1) then
            do n=0,nc-1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     coef(i,j,n) = merge(lam(i-1,j,n),lam(i,j,n),vel(i,j)>=zero)
                  enddo
               enddo
            enddo
         else
            do n=0,nc-1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     coef(i,j,n) = 0.5d0 * (lam(i-1,j,n)+lam(i,j,n))
                  enddo
               enddo
            enddo
         endif
         if (dirbclo.eq.1) then 
            do n=0,nc-1
               i = lo(1)
               do j = lo(2),hi(2)
                  coef(i,j,n) = lam(i-1,j,n)
               enddo
            enddo
         endif
         if (dirbchi.eq.1) then 
            do n=0,nc-1
               i = hi(1)+1
               do j = lo(2),hi(2)
                  coef(i,j,n) = lam(i,j,n)
               enddo
            enddo
         endif
      else if (dir.eq.1) then
         if (do_upwind.eq.1) then
            do n=0,nc-1
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     coef(i,j,n) = merge(lam(i,j-1,n),lam(i,j,n),vel(i,j)>=zero)
                  enddo
               enddo
            enddo
         else
            do n=0,nc-1
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     coef(i,j,n) = 0.5d0 * (lam(i,j-1,n)+lam(i,j,n))
                  enddo
               enddo
            enddo
         endif
         if (dirbclo.eq.1) then 
            do n=0,nc-1
               j = lo(2)
               do i = lo(1),hi(1)
                  coef(i,j,n) = lam(i,j-1,n)
               enddo
            enddo
         endif
         if (dirbchi.eq.1) then 
            do n=0,nc-1
               j = hi(2)+1
               do i = lo(1),hi(1)
                  coef(i,j,n) = lam(i,j,n)
               enddo
            enddo
         endif
      else
         print *, 'FORT_CALC_RS_H::bad dir: ',dir
         call bl_pd_abort()
      endif
      end     

      subroutine FORT_RS_PDOTRES(Res,DIMS(Res),RSn,DIMS(RSn),
     &     RSnp1,DIMS(RSnp1),Porositynp1,DIMS(Porositynp1),
     &     Porosityn,DIMS(Porosityn),dt,lo,hi,nc)
      implicit none      
      integer DIMDEC(Res)
      integer DIMDEC(RSn)
      integer DIMDEC(RSnp1)
      integer DIMDEC(Porositynp1)
      integer DIMDEC(Porosityn)
      integer lo(SDIM),hi(SDIM)
      integer nc
      REAL_T  Res(DIMV(Res),nc)
      REAL_T  RSn(DIMV(RSn),nc)
      REAL_T  RSnp1(DIMV(RSnp1),nc)
      REAL_T  Porositynp1(DIMV(Porositynp1))
      REAL_T  Porosityn(DIMV(Porosityn))
      REAL_T  rho(nc)
      REAL_T  dt
      
      integer i, j, n

c     On input, Res contains Div(rho.DarcyVel)
      do n=1,nc
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Res(i,j,n) = (Porositynp1(i,j)*RSnp1(i,j,n) - Porosityn(i,j)*RSn(i,j,n))/dt + Res(i,j,n)
            enddo
         enddo
      enddo
      end


      subroutine FORT_RS_ALTUP(rhosat,DIMS(rhosat),p,DIMS(p),
     &     dp,DIMS(dp),kappa,DIMS(kappa),coef,DIMS(coef),n_coef,
     &     ls_factor,lo,hi,sTOL)

      implicit none

#include "probdata.H"

      integer DIMDEC(rhosat)
      integer DIMDEC(p)
      integer DIMDEC(dp)
      integer DIMDEC(kappa)
      integer DIMDEC(coef)
      integer n_coef
      REAL_T  ls_factor
      integer lo(SDIM),hi(SDIM)
      REAL_T  rhosat(DIMV(rhosat),ncomps)
      REAL_T  p(DIMV(p),ncomps)
      REAL_T  dp(DIMV(dp),ncomps)
      REAL_T  kappa(DIMV(kappa))
      REAL_T  coef(DIMV(coef),n_coef)
      REAL_T  sTOL

      integer i,j
      integer l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res, dpcdse, tphi, pcap_new
      REAL_T  rs_old, s_old, s_new, ds, p_new, pcap_old, rs_new
      logical is_saturated
      REAL_T  saturated_s_thresh
      parameter (saturated_s_thresh = 1.d0 - 1.d-8)

      tphi = one
      if (ncomps.ne.1) then
         call bl_pd_abort('Richard solver only works for single comp')
      endif

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            l_cpl_type = nint(coef(i,j,1))
            l_cpl_coef = coef(i,j,2)
            l_sigma    = coef(i,j,3)
            l_sat_res  = coef(i,j,4)


c     Get s^{k+1,*} = s^{k}  -  (ds/dp)^k . fac.dp,  ds/dp = - 1 / [(dpc/dse).(dse/ds)] 
c     If s^{k+1,*} < sTOL, then s^{k+1} = s^{k+1,*}, p^{k+1} from s^{k+1}
c       else p^{k+1} = p^{k} - fac.dp, s^{k+1} from p^{k+1}
c
c     Note: This routine modifies only dp, and only if s^{k}<sTOL

            rs_old = rhosat(i,j,1) 
            s_old = rs_old / rhoval(1)
            is_saturated = (s_old .ge. saturated_s_thresh)

            if (.not. is_saturated) then
               call get_dcplds(dpcdse,rs_old,rhoval,tphi,kappa(i,j),
     &              l_sigma,ncomps,zero,l_cpl_type,l_cpl_coef)
               ds = -ls_factor*(one-l_sat_res)*dp(i,j,1)/dpcdse
               s_new = s_old - ds
               s_new = MIN(one,MAX(zero,s_new))
               ds = s_new - s_old
               is_saturated = (s_new .ge. sTOL)
            endif
            
            if (.not. is_saturated) then
               call get_cpl(pcap_old,rs_old,rhoval(1),tphi,
     &              kappa(i,j),l_sigma,ncomps,l_sat_res,l_cpl_type)

               rs_new = rhoval(1) * s_new

               call get_cpl(pcap_new,rs_new,rhoval(1),tphi,
     &              kappa(i,j),l_sigma,ncomps,l_sat_res,l_cpl_type)

               if (abs(dp(i,j,1)-pcap_new+pcap_old).gt.1.d-2) then
c                  print *,i,j,dp(i,j,1)-pcap_new+pcap_old
               endif

               dp(i,j,1) = pcap_new - pcap_old

            endif
         end do 
      end do
      end
