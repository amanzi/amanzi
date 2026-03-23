
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "Prob_PM_F.H"
#include "PorousMedia_F.H"
#include "ArrayLim.H"

#define SDIM 2


c ::: ----------------------------------------------------------
c ::: Volume-weight average the fine grid data onto the coarse
c ::: grid.  Overlap is given in coarse grid coordinates.      
c ::: crse      =  coarse grid data
c ::: nvar	= number of components in arrays
c ::: fine      = fine grid data
c ::: cv        = coarse grid volume array
c ::: fv        = fine grid volume array
c ::: lo,hi     = index limits of overlap (crse grid)
c ::: ratios    = IntVect refinement ratio
c ::: ----------------------------------------------------------

      subroutine FORT_AVGDOWN (
     &     crse,DIMS(crse),nvar,
     &     fine,DIMS(fine),
     &     cv,DIMS(cv),
     &     fv,DIMS(fv),
     &     lo,hi,ratios)

      integer  DIMDEC(crse)
      integer  DIMDEC(cv)
      integer  DIMDEC(fine)
      integer  DIMDEC(fv)
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T     cv(DIMV(cv))
      REAL_T   fine(DIMV(fine),nvar)
      REAL_T     fv(DIMV(fv))

      integer  i, j, n, ic, jc, ioff, joff
      integer  lratx,lraty

      lratx = ratios(1)
      lraty = ratios(2)

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = zero
            end do
         end do
c
c     sum fine data
c
         do joff = 0, lraty-1
            do jc = lo(2), hi(2)
               j = jc*lraty + joff
               do ioff = 0, lratx-1
                  do ic = lo(1), hi(1)
                     i = ic*lratx + ioff
                     crse(ic,jc,n) = crse(ic,jc,n) +
     &                    fv(i,j)*fine(i,j,n)
                  end do
               end do
            end do
         end do
c
c     divide out by volume weight
c
         do ic = lo(1), hi(1)
            do jc = lo(2), hi(2)
               crse(ic,jc,n) = crse(ic,jc,n)/cv(ic,jc)
            end do
         end do
      end do

      end

c ::: ---------------------------------------------------------
c ::: SUMMASS
c :::            MASS = sum{ vol(i,j)*rho(i,j) }
c :::
c ::: IPUTS / OUTPUTS:
c :::  ho        => density field
c :::  IMS(rho)  => index limits of rho aray
c :::  lo,hi     => index limits of grid interior
c :::  dx	 => cell size
c :::  mass     <=  total mass
c ::: ---------------------------------------------------------
      subroutine FORT_SUMMASS(rho,DIMS(rho),DIMS(grid),dx,mass)

      integer DIMDEC(rho)
      integer DIMDEC(grid)
      REAL_T  mass, dx(2)
      REAL_T  rho(DIMV(rho))

      integer i, j
      REAL_T  vol

      vol = dx(1)*dx(2)

      mass = zero
      do i = ARG_L1(grid), ARG_H1(grid)
         do j = ARG_L2(grid), ARG_H2(grid)
            mass = mass + vol*rho(i,j)
         end do
      end do

      end

c ::: ----------------------------------------------------------
c ::: This routine fills an edge-centered fab from a cell-centered
c ::: fab using simple linear interpolation.
c :::
c ::: INPUTS / OUTPUTS:
c :::  lo,hi      => index limits of the of the cell-centered fab
c :::  DIMS(cfab) => index limits of the cell-centered fab
c :::  cfab       => cell-centered data
c :::  DIMS(efab) => index limits of the edge-centered fab
c :::  efab       => edge-centered fab to fill
c :::  nc         => Number of components in the fab to fill
c :::  dir        => direction data needs to be shifted to get to edges
c ::: ----------------------------------------------------------

      subroutine FORT_CEN2EDG(lo, hi, 
     &     DIMS(cfab), cfab,
     &     DIMS(efab), efab, nc, dir,
     &     isharm)
      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1)+1, hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i-1,j,n))
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2)+1, hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i,j-1,n))
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1)+1, hi(1)
                     if((cfab(i,j,n) * cfab(i-1,j,n)).gt.0.d0)then
                        efab(i,j,n)
     &                       = 2.d0*(cfab(i,j,n) * cfab(i-1,j,n))/
     &                       (cfab(i,j,n) + cfab(i-1,j,n))
                     else
                        efab(i,j,n)=0.d0
                     endif
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2)+1, hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i,j-1,n)).gt.0.d0)then
                        efab(i,j,n)
     &                       = 2.d0*(cfab(i,j,n) * cfab(i,j-1,n))/
     &                       (cfab(i,j,n) + cfab(i,j-1,n))
                     else
                        efab(i,j,n)=0.d0
                     endif
                  end do
               end do
            end do
         end if
      end if
      end

c ::: -------------------------------------------------------------------
c ::: Do linear in dir, pc transverse to dir, leave alone the fine values
c ::: lining up with coarse edges--assume these have been set to hold the 
c ::: values you want to interpolate to the rest.
c ::: -------------------------------------------------------------------
      subroutine EDGE_INTERP(flo, fhi, nc, ratio, dir,
     &     fine, fine_l0, fine_l1, fine_h0, fine_h1)
      implicit none
      integer flo(0:2-1), fhi(0:2-1), nc, ratio(0:2-1), dir
      integer fine_l0, fine_l1, fine_h0, fine_h1
      DOUBLE PRECISION fine(fine_l0:fine_h0,fine_l1:fine_h1,nc)
      integer i,j,n,P,M
      DOUBLE PRECISION val, df

      if (dir.eq.0) then
         do n=1,nc
            do j=flo(1),fhi(1),ratio(1)
               do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                  df = fine(i+ratio(dir),j,n)-fine(i,j,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                        fine(i+M,P,n) = val
                     enddo
                  enddo                     
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=flo(1),fhi(1)-ratio(dir),ratio(1)
               do i=flo(0),fhi(0)
                  df = fine(i,j+ratio(dir),n)-fine(i,j,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                        fine(P,j+M,n) = val
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      end

c ::: ---------------------------------------------------------------
c ::: For edge-based data, fill fine values with piecewise-constant 
c ::: interp of coarse data.
c ::: Operate only on faces that overlap--ie, only fill the fine faces  
c ::: that make up each coarse face, leave the in-between faces alone.
c ::: ----------------------------------------------------------------

      subroutine PC_EDGE_INTERP(lo, hi, nc, ratio, dir,
     &     crse, crse_l0, crse_l1, crse_h0, crse_h1,
     &     fine, fine_l0, fine_l1, fine_h0, fine_h1)
      implicit none
      integer lo(2),hi(2), nc, ratio(0:2-1), dir
      integer crse_l0, crse_l1, crse_h0, crse_h1
      integer fine_l0, fine_l1, fine_h0, fine_h1
      DOUBLE PRECISION crse(crse_l0:crse_h0,crse_l1:crse_h1,nc)
      DOUBLE PRECISION fine(fine_l0:fine_h0,fine_l1:fine_h1,nc)
      integer i,j,ii,jj,n,L
      if (dir.eq.0) then
         do n=1,nc
            do j=lo(2),hi(2)
               jj = ratio(1)*j
               do i=lo(1),hi(1)
                  ii = ratio(0)*i
                  do L=0,ratio(1)-1
                     fine(ii,jj+L,n) = crse(i,j,n)
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=lo(2),hi(2)
               jj = ratio(1)*j
               do i=lo(1),hi(1)
                  ii = ratio(0)*i
                  do L=0,ratio(0)-1
                     fine(ii+L,jj,n) = crse(i,j,n)
                  enddo
               enddo
            enddo
         enddo
      endif

      end



      subroutine FORT_MAX_TRACDIFF_DT(lo, hi, psv, DIMS(psv),
     &     psadx, DIMS(psadx), psady, DIMS(psady),
     &     delta, max_dt)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(psv)
      integer DIMDEC(psadx)
      integer DIMDEC(psady)
      REAL_T  delta(SDIM), max_dt
      REAL_T    psv(DIMV(psv))
      REAL_T  psadx(DIMV(psadx))
      REAL_T  psady(DIMV(psady))
      REAL_T  denom, vol, dx, dy
      integer i, j

      dx = delta(1)
      dy = delta(2)
      vol = dx * dy

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            denom = (psadx(i,j) + psadx(i+1,j)) * dy
     &           +  (psady(i,j) + psady(i,j+1)) * dx
            max_dt = MIN(max_dt, vol*psv(i,j)/denom)
         enddo
      enddo
      end

c ::: -----------------------------------------------------------
c ::: Compute the inverse of the capillary pressure. 
c ::: Assume capillary pressure is positive.  Otherwise, we have 
c ::: saturated condition.
c :::    It calls get_invcpl in EOS_F.F. 
c ::: -----------------------------------------------------------

      subroutine FORT_MK_INV_CPL(pc,DIMS(pc),s,DIMS(s),
     &     phi,DIMS(phi),kappa,DIMS(kappa),
     &     coef,DIMS(coef),n_coef)

      implicit none

#include "probdata.H"

      integer DIMDEC(pc)
      integer DIMDEC(s)
      integer DIMDEC(phi)
      integer DIMDEC(kappa)
      integer DIMDEC(coef)
      integer n_coef
      REAL_T  pc(DIMV(pc))
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  phi(DIMV(phi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  coef(DIMV(coef),n_coef)

      integer i,j
      integer l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res


      call bl_abort('Should not be here in FORT_MK_INV_CPL')

      do j =  ARG_L2(s),ARG_H2(s)
         do i = ARG_L1(s),ARG_H1(s)
            l_cpl_type = nint(coef(i,j,1))
            l_cpl_coef = coef(i,j,2)
            l_sigma    = coef(i,j,3)
            l_sat_res  = coef(i,j,4)

            if (pc(i,j) .gt. zero) then
               call get_invcpl(s(i,j,:),pc(i,j),rhoval,phi(i,j),
     &              kappa(i,j),l_sigma,ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
            else
               s(i,j,1) = rhoval(1)
               if (ncomps .gt. 1) then
                  s(i,j,2) = zero
               end if
            end if
         end do 
      end do

      end

c ::: -----------------------------------------------------------
c ::: Compute the capillary pressure, p_c.
c :::    It calls get_cpl in EOS_F.F. 
c ::: -----------------------------------------------------------

      subroutine FORT_MK_CPL(pc,DIMS(pc),s,DIMS(s),
     &     phi,DIMS(phi),kappa,DIMS(kappa),
     &     coef,DIMS(coef),n_coef,lo,hi,bc)

      implicit none

#include "probdata.H"

      integer DIMDEC(pc)
      integer DIMDEC(s)
      integer DIMDEC(phi)
      integer DIMDEC(kappa)
      integer DIMDEC(coef)
      integer n_coef
      integer lo(SDIM),hi(SDIM)
      integer bc(SDIM,2)
      REAL_T  pc(DIMV(pc))
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  phi(DIMV(phi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  coef(DIMV(coef),n_coef)

      integer i,j
      integer l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res, dpc

      call bl_abort('Should not be here in FORT_MK_CPL')

      do j =  ARG_L2(s),ARG_H2(s)
         do i = ARG_L1(s),ARG_H1(s)
            l_cpl_type = nint(coef(i,j,1))
            l_cpl_coef = coef(i,j,2)
            l_sigma    = coef(i,j,3)
            l_sat_res  = coef(i,j,4)
            if (i .eq. lo(1)-1 .and. bc(1,1) .eq. HOEXTRAP) then
               dpc = pc(i+2,j)-pc(i+1,j)
               pc(i,j) = pc(i+1,j) - dpc/2.0
            elseif (i .eq. hi(1)+1 .and. bc(1,2) .eq. HOEXTRAP) then
               dpc = pc(i-1,j)-pc(i-2,j)
               pc(i,j) = pc(i-1,j) + dpc/2.0
            elseif (j .eq. lo(2)-1 .and. bc(2,1) .eq. HOEXTRAP) then
               dpc = pc(i,j+2)-pc(i,j+1)
               pc(i,j) = pc(i,j+1) - dpc/2.0
            elseif (j .eq. hi(2)+1 .and. bc(2,2) .eq. HOEXTRAP) then
               dpc = pc(i,j-1)-pc(i,j-2)
               pc(i,j) = pc(i,j-1) + dpc/2.0
            else
               call get_cpl(pc(i,j),s(i,j,:),rhoval(1),phi(i,j),
     &              kappa(i,j),l_sigma,ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
            endif
         end do 
      end do
      end


c ::: -----------------------------------------------------------
c ::: Compute the lambda = kr/mu for all components
c       -  Calls get_ld in EOS_F.F. 
c ::: -----------------------------------------------------------

      subroutine FORT_MK_LAMBDA(lbd,DIMS(lbd),s,DIMS(s),
     &     kr,DIMS(kr),n_kr)

#include "probdata.H"

      integer DIMDEC(lbd)
      integer DIMDEC(s)
      integer DIMDEC(kr)
      integer n_kr
      REAL_T  lbd(DIMV(lbd),ncomps)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  kr(DIMV(kr),n_kr)

      integer i,j
      integer l_kr_type
      REAL_T  l_kr_coef, l_sat_res

      call bl_abort('Should not be here in FORT_MK_LAMBDA')

      do j = ARG_L2(lbd),ARG_H2(lbd)
         do i = ARG_L1(lbd),ARG_H1(lbd)
            l_kr_type = nint(kr(i,j,1))
            l_kr_coef = kr(i,j,2)
            l_sat_res = kr(i,j,3)
            
            call get_ld(lbd(i,j,1:ncomps),s(i,j,1:ncomps),
     &           rhoval,muval,ncomps,
     &           l_kr_type,l_kr_coef,l_sat_res)

         end do 
      end do

      end
c ::: -----------------------------------------------------------
c ::: Compute the derivative of lambda = kr/mu for all components
c       - Calls get_dld in EOS_F.F. 
c ::: -----------------------------------------------------------

      subroutine FORT_MK_DLAMBDA(lbd,DIMS(lbd),s,DIMS(s),
     &     kr,DIMS(kr),n_kr)

#include "probdata.H"

      integer DIMDEC(lbd)
      integer DIMDEC(s)
      integer DIMDEC(kr)
      integer n_kr
      REAL_T  lbd(DIMV(lbd),3)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  kr(DIMV(kr),n_kr)

      integer i,j
      integer l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T  ld_tmp(3), s_tmp(2)

      do j = ARG_L2(lbd),ARG_H2(lbd)
         do i = ARG_L1(lbd),ARG_H1(lbd)
            l_kr_type = nint(kr(i,j,1))
            l_kr_coef = kr(i,j,2)
            l_sat_res = kr(i,j,3)

            s_tmp  = zero
            s_tmp(1:ncomps) = s(i,j,1:ncomps)

            call get_dld(ld_tmp,s_tmp,rhoval,muval,ncomps,
     &           l_kr_type,l_kr_coef,l_sat_res)

            lbd(i,j,:) = ld_tmp

         end do 
      end do

      end

c :: ------------------------------------------------------------
c :: RICHARD_AJAC
c      Compute the coefficient for the analytical jacobian in the  
c      equilibrium richard equation.  DON'T USE AT THE MOMENT. 
c
c     dfx(i,j,1) = lbdx(i,j,1)*dpcdn(i-1,j)
c     dfx(i,j,2) = lbdx(i,j,1)*dpcdn(i+1,j)
c     dfx(i,j,3) = (lbdx(i,j,1)+lbdx(i+1,j,1))*dpcdn(i,j)
c     dfy(i,j,1) = lbdy(i,j,1)*dpcdn(i-1,j)
c     dfy(i,j,2) = lbdy(i,j,1)*dpcdn(i+1,j)
c     dfy(i,j,3) = (lbdy(i,j,1)+lbdy(i+1,j,1))*dpcdn(i,j)
c
c :: ------------------------------------------------------------         
      subroutine FORT_RICHARD_AJAC(s, DIMS(s), 
     &     bx,DIMS(bx),by,DIMS(by),
     &     umac, DIMS(umac), vmac, DIMS(vmac),
     &     kpx, DIMS(kpx), kpy, DIMS(kpy),
     &     lbd,DIMS(lbd),pc,DIMS(pc),
     &     phi,DIMS(phi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     cplcoef,DIMS(cplcoef),ncplcoef,
     &     lo,hi,domlo,domhi,dx,bc, 
     &     inflow_bc_lo,inflow_bc_hi,do_upwind)

      implicit none

#include "probdata.H"

      integer DIMDEC(s)
      integer DIMDEC(bx)
      integer DIMDEC(by)
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(kpx)
      integer DIMDEC(kpy)
      integer DIMDEC(pc)
      integer DIMDEC(lbd)
      integer DIMDEC(phi)
      integer DIMDEC(kappa)
      integer DIMDEC(krcoef)
      integer DIMDEC(cplcoef)
      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      integer bc(SDIM,2)
      integer inflow_bc_lo(SDIM),inflow_bc_hi(SDIM)
      integer nkrcoef,ncplcoef
      REAL_T  dx(SDIM)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  bx(DIMV(bx),3)
      REAL_T  by(DIMV(by),3)
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))  
      REAL_T  kpx(DIMV(kpx))
      REAL_T  kpy(DIMV(kpy))
      REAL_T  pc(DIMV(pc))
      REAL_T  lbd(DIMV(lbd),ncomps)
      REAL_T  cplcoef(DIMV(cplcoef),ncplcoef)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  phi(DIMV(phi))
      REAL_T  kappa(DIMV(kappa))
      integer do_upwind

      integer i, j,l_kr_type,l_cpl_type
      REAL_T l_kr_coef, l_kr_sres, l_cpl_coef,l_sigma,l_cpl_sres
      REAL_T  gravx, gravy
      REAL_T  multx(2), multy(2), ldx(2), ldy(2), dldx(4), dldy(4)

      REAL_T  sat(ncomps), stmp(ncomps), rtmp(ncomps), dltmp(3)
      REAL_T, allocatable :: dpcdn(:,:), dlddn(:,:,:)

      allocate(dpcdn(DIMV(s)))
      allocate(dlddn(DIMV(s),ncomps))

      rtmp(:) = one
      do j = ARG_L2(s),ARG_H2(s)
         do i = ARG_L1(s),ARG_H1(s)

            l_kr_type  = nint(krcoef(i,j,1))
            l_kr_coef  = krcoef(i,j,2)
            l_kr_sres  = krcoef(i,j,3)

            l_cpl_type  = nint(cplcoef(i,j,1))
            l_cpl_coef  = cplcoef(i,j,2)
            l_sigma     = cplcoef(i,j,3)
            l_cpl_sres  = cplcoef(i,j,4)

            stmp = zero
            stmp(1:ncomps) = s(i,j,1:ncomps) 

            call get_st(sat,stmp,rhoval,l_cpl_sres,ncomps)

            if (one-sat(1) .lt. 1.d-6) then
               sat(1) = one - 1.d-6
               if (ncomps > 1) then
                  sat(2) = 1.d-6
               end if
            end if

            call get_dcplds(dpcdn(i,j),sat,rtmp,phi(i,j),
     &           kappa(i,j),l_sigma,ncomps,zero,l_cpl_type,l_cpl_coef)
            call get_dld(dltmp,sat,rtmp,muval,ncomps,
     &           l_kr_type,l_kr_coef,l_kr_sres)
         
            dpcdn(i,j) = dpcdn(i,j)/(rhoval(1)*(one-l_cpl_sres))
            dlddn(i,j,1) = dltmp(1)/(rhoval(1)*(one-l_kr_sres))

         end do
      end do

      if (gravity_dir.eq.0) then
         gravx = -rhoval(1)*gravity
         gravy = 0.d0
      else if (gravity_dir.eq.1) then
         gravx = 0.d0
         gravy = -rhoval(1)*gravity
      endif
      bx = zero
      by = zero
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            if (do_upwind.eq.1) then
               ldx(1)  = merge(lbd(i-1,j,1),lbd(i,j,1),umac(i,j)>=zero)
               ldx(2)  = merge(lbd(i,j,1),lbd(i+1,j,1),umac(i+1,j)>=zero)
               ldy(1)  = merge(lbd(i,j-1,1),lbd(i,j,1),vmac(i,j)>=zero)
               ldy(2)  = merge(lbd(i,j,1),lbd(i,j+1,1),vmac(i,j+1)>=zero)
               dldx(1) = merge(dlddn(i-1,j,1),zero,umac(i,j)>=zero)
               dldx(2) = merge(zero,dlddn(i+1,j,1),umac(i+1,j)>=zero) 
               dldx(3) = merge(zero,dlddn(i,j,1),umac(i,j)>=zero)
               dldx(4) = merge(dlddn(i,j,1),zero,umac(i+1,j)>=zero) 
               dldy(1) = merge(dlddn(i,j-1,1),zero,vmac(i,j)>=zero)
               dldy(2) = merge(zero,dlddn(i,j+1,1),vmac(i,j+1)>=zero) 
               dldy(3) = merge(zero,dlddn(i,j,1),vmac(i,j)>=zero)
               dldy(4) = merge(dlddn(i,j,1),zero,vmac(i,j+1)>=zero)

               ldx(1)  = kpx(i,j)*ldx(1)
               ldx(2)  = kpx(i+1,j)*ldx(2)
               ldy(1)  = kpy(i,j)*ldy(1)
               ldy(2)  = kpy(i,j+1)*ldy(2)
               dldx(1) = kpx(i,j)*dldx(1)
               dldx(2) = kpx(i+1,j)*dldx(2)
               dldx(3) = kpx(i,j)*dldx(3)
               dldx(4) = kpx(i+1,j)*dldx(4)
               dldy(1) = kpy(i,j)*dldy(1)
               dldy(2) = kpy(i,j+1)*dldy(2)
               dldy(3) = kpy(i,j)*dldy(3)
               dldy(4) = kpy(i,j+1)*dldy(4)
            else
               ldx(1)  = kpx(i,j)*half*(lbd(i,j,1)+lbd(i-1,j,1))
               ldx(2)  = kpx(i+1,j)*half*(lbd(i,j,1)+lbd(i+1,j,1))
               ldy(1)  = kpy(i,j)*half*(lbd(i,j,1)+lbd(i,j-1,1))
               ldy(2)  = kpy(i,j+1)*half*(lbd(i,j,1)+lbd(i,j+1,1))
               dldx(1) = kpx(i,j)*half*dlddn(i-1,j,1)
               dldx(2) = kpx(i+1,j)*half*dlddn(i+1,j,1)
               dldx(3) = kpx(i,j)*half*dlddn(i,j,1)
               dldx(4) = kpx(i+1,j)*half*dlddn(i,j,1)
               dldy(1) = kpy(i,j)*half*dlddn(i,j-1,1)
               dldy(2) = kpy(i,j+1)*half*dlddn(i,j+1,1)
               dldy(3) = kpy(i,j)*half*dlddn(i,j,1)
               dldy(4) = kpy(i,j+1)*half*dlddn(i,j,1)
            end if

            multx = one
            multy = one
            if (i .eq. domlo(1)) then
               if (bc(1,1) .eq. EXT_DIR) then
                  multx(1) = two 
                  ldx(1) = kpx(i,j)*lbd(i-1,j,1)
                  dldx(1) = kpx(i,j)*dlddn(i-1,j,1)
                  dldx(3) = zero
               elseif (bc(1,1) .eq. FOEXTRAP) then
                  multx(1) = zero
                  if (inflow_bc_lo(1) .eq. 1) then
                     ldx(1)  = zero
                     dldx(1) = zero
                     dldx(3) = zero
                  end if
               end if
            end if
            if (i .eq. domhi(1)) then
               if (bc(1,2) .eq. EXT_DIR) then
                  multx(2) = two
                  ldx(2) = kpx(i+1,j)*lbd(i+1,j,1)
                  dldx(2) = kpx(i+1,j)*dlddn(i+1,j,1)
                  dldx(4) = zero
               elseif (bc(1,2) .eq. FOEXTRAP) then
                  multx(2) = zero
                  if (inflow_bc_hi(1) .eq. 1) then
                     ldx(2) = zero
                     dldx(2) = zero
                     dldx(4) = zero
                  end if
               end if
            end if
            if (j .eq. domlo(2)) then
               if (bc(2,1) .eq. EXT_DIR) then
                  multy(1) = two
                  ldy(1) = kpy(i,j)*lbd(i,j-1,1)
                  dldy(1) = kpy(i,j)*dlddn(i,j-1,1)
                  dldy(3) = zero
               elseif (bc(2,1) .eq. FOEXTRAP) then
                  multy(1) = zero
                  if (inflow_bc_lo(2) .eq. 1) then
                     ldy(1) = zero
                     dldy(1) = zero
                     dldy(3) = zero
                  end if
               end if
            end if
            if (j .eq. domhi(2)) then
               if (bc(2,2) .eq. EXT_DIR) then
                  multy(2) = two 
                  ldy(2) =  kpy(i,j+1)*lbd(i,j+1,1)
                  dldy(2) =  kpy(i,j+1)*dlddn(i,j+1,1)
                  dldy(4) = zero
               elseif (bc(2,2) .eq. FOEXTRAP) then
                  multy(2) = zero
                  if (inflow_bc_hi(2) .eq. 1) then
                     ldy(2) = zero
                     dldy(2) = zero
                     dldy(4) = zero
                  end if
               end if
            end if

            bx(i,j,2)   = (-dldx(1)*(-multx(1)*(pc(i,j)-pc(i-1,j))/dx(1) - gravx) -
     &           ldx(1)*(multx(1)*dpcdn(i-1,j)/dx(1)))/dx(1)
            bx(i,j,1) = (dldx(2)*(-multx(2)*(pc(i+1,j)-pc(i,j))/dx(1) - gravx) -
     &           ldx(2)*(multx(2)*dpcdn(i+1,j)/dx(1)))/dx(1)
            bx(i,j,3)   = (dldx(4)*(-multx(2)*(pc(i+1,j)-pc(i,j))/dx(1) - gravx) - 
     &           dldx(3)*(-multx(1)*(pc(i,j)-pc(i-1,j))/dx(1) - gravx) +
     &           (ldx(1)+ldx(2))*dpcdn(i,j)/dx(1))/dx(1)
            
            by(i,j,2)   = (-dldy(1)*(-multy(1)*(pc(i,j)-pc(i,j-1))/dx(2) - gravy) -
     &           ldy(1)*(multy(1)*dpcdn(i,j-1)/dx(2)))/dx(2)
            by(i,j,1) = (dldy(2)*(-multy(2)*(pc(i,j+1)-pc(i,j))/dx(2) - gravy) -
     &           ldy(2)*(multy(2)*dpcdn(i,j+1)/dx(2)))/dx(2)
            by(i,j,3)   = (dldy(4)*(-multy(2)*(pc(i,j+1)-pc(i,j))/dx(2) - gravy) - 
     &           dldy(3)*(-multy(1)*(pc(i,j)-pc(i,j-1))/dx(2) - gravy) +
     &           (ldy(1)+ldy(2))*dpcdn(i,j)/dx(2))/dx(2)
         end do
      end do

      deallocate(dpcdn)
      deallocate(dlddn)

      end

c :: ------------------------------------------------------------
c :: RICHARD_NJAC
c      Compute the coefficient for the numerical jacobian in the  
c      equilibrium richard's equation (w.r.t to n).
c
c     dfx(i,j,1) = dDF(i,j)/dn(i-1,j)
c     dfx(i,j,2) = dDF(i,j)/dn(i+1,j)
c     dfx(i,j,3) = dDF(i,j)/dn(i,j) 
c     dfy(i,j,1) = dDF(i,j)/dn(i,j-1)
c     dfy(i,j,2) = dDF(i,j)/dn(i,j+1)
c     dfy(i,j,3) = dDF(i,j)/dn(i,j)
c
c :: ------------------------------------------------------------    
      subroutine FORT_RICHARD_NJAC(s, DIMS(s),
     &     bx,DIMS(bx),by,DIMS(by),
     &     umac, DIMS(umac), vmac, DIMS(vmac),
     &     kpx, DIMS(kpx), kpy, DIMS(kpy),
     &     lbd,DIMS(lbd),pc,DIMS(pc),
     &     phi,DIMS(phi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     cplcoef,DIMS(cplcoef),ncplcoef,
     &     lo,hi,domlo,domhi,dx, bc, 
     &     inflow_bc_lo,inflow_bc_hi,deps,
     &     do_upwind)

      implicit none

#include "probdata.H"

      integer DIMDEC(s)
      integer DIMDEC(bx)
      integer DIMDEC(by)
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(kpx)
      integer DIMDEC(kpy)
      integer DIMDEC(pc)
      integer DIMDEC(lbd)
      integer DIMDEC(phi)
      integer DIMDEC(kappa)
      integer DIMDEC(krcoef)
      integer DIMDEC(cplcoef)
      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      integer bc(SDIM,2)
      integer inflow_bc_lo(SDIM),inflow_bc_hi(SDIM)
      integer nkrcoef,ncplcoef
      REAL_T  dx(SDIM)
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  bx(DIMV(bx),3)
      REAL_T  by(DIMV(by),3)
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  kpx(DIMV(kpx))
      REAL_T  kpy(DIMV(kpy))
      REAL_T  pc(DIMV(pc))
      REAL_T  lbd(DIMV(lbd),ncomps)
      REAL_T  cplcoef(DIMV(cplcoef),ncplcoef)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  phi(DIMV(phi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  deps
      integer do_upwind

      integer i, j,n,l_kr_type,l_cpl_type
      REAL_T l_kr_coef, l_kr_sres, l_cpl_coef,l_sigma,l_cpl_sres
      REAL_T  gravx, gravy, rgx, rgy
      REAL_T  multx(2), multy(2), ldx(2), ldy(2)
      REAL_T  comp(ncomps)
      REAL_T  ldtmp(ncomps),ldn(3),pcn(3),rgtmpx(3),rgtmpy(3), dcx(3), dcy(3)

      if (gravity_dir.eq.0) then
         gravx = -rhoval(1)*gravity
         gravy = 0.d0
      else if (gravity_dir.eq.1) then
         gravx = 0.d0
         gravy = -rhoval(1)*gravity
      endif

      bx = zero
      by = zero
      do j=lo(2), hi(2)
         do i=lo(1), hi(1)

            if (do_upwind.eq.1) then
               ldx(1) = merge(kpx(i,j)*lbd(i-1,j,1),kpx(i,j)*lbd(i,j,1),umac(i,j)>=zero)
               ldx(2) = merge(kpx(i+1,j)*lbd(i,j,1),kpx(i+1,j)*lbd(i+1,j,1),umac(i+1,j)>=zero)
               ldy(1) = merge(kpy(i,j)*lbd(i,j-1,1),kpy(i,j)*lbd(i,j,1),vmac(i,j)>=zero)
               ldy(2) = merge(kpy(i,j+1)*lbd(i,j,1),kpy(i,j+1)*lbd(i,j+1,1),vmac(i,j+1)>=zero)
            else
               ldx(1) = kpx(i,j)*half*(lbd(i,j,1)+lbd(i-1,j,1))
               ldx(2) = kpx(i+1,j)*half*(lbd(i,j,1)+lbd(i+1,j,1))
               ldy(1) = kpy(i,j)*half*(lbd(i,j,1)+lbd(i,j-1,1))
               ldy(2) = kpy(i,j+1)*half*(lbd(i,j,1)+lbd(i,j+1,1))
            end if
            multx = one
            multy = one
            if (i .eq. domlo(1)) then
               if (bc(1,1) .eq. EXT_DIR) then
                  multx(1) = two 
                  ldx(1) = kpx(i,j)*lbd(i-1,j,1)
               elseif (bc(1,1) .eq. FOEXTRAP) then
                  multx(1) = zero
                  if (inflow_bc_lo(1) .eq. 1) then
                     ldx(1) = zero
                  end if
               end if
            end if
            if (i .eq. domhi(1)) then
               if (bc(1,2) .eq. EXT_DIR) then
                  multx(2) = two
                  ldx(2) = kpx(i+1,j)*lbd(i+1,j,1)
               elseif (bc(1,2) .eq. FOEXTRAP) then
                  multx(2) = zero
                  if (inflow_bc_hi(1) .eq. 1) then
                     ldx(2) = zero
                  end if
               end if
            end if
            if (j .eq. domlo(2)) then
               if (bc(2,1) .eq. EXT_DIR) then
                  multy(1) = two
                  ldy(1) = kpy(i,j)*lbd(i,j-1,1)
               elseif (bc(2,1) .eq. FOEXTRAP) then
                  multy(1) = zero
                  if (inflow_bc_lo(2) .eq. 1) then
                     ldy(1) = zero
                  end if
               end if
            end if
            if (j .eq. domhi(2)) then
               if (bc(2,2) .eq. EXT_DIR) then
                  multy(2) = two 
                  ldy(2) = kpy(i,j+1)*lbd(i,j+1,1)
               elseif (bc(2,2) .eq. FOEXTRAP) then
                  multy(2) = zero
                  if (inflow_bc_hi(2) .eq. 1) then
                     ldy(2) = zero
                  end if
               end if
            end if

            rgx = (ldx(2)*(-multx(2)*(pc(i+1,j)-pc(i,j))/dx(1) - gravx) -
     &           ldx(1)*(-multx(1)*(pc(i,j)-pc(i-1,j))/dx(1) - gravx))/dx(1) 
            rgy = (ldy(2)*(-multy(2)*(pc(i,j+1)-pc(i,j))/dx(2) - gravy) - 
     &           ldy(1)*(-multy(1)*(pc(i,j)-pc(i,j-1))/dx(2) - gravy))/dx(2)           

c           x-direction
            do n = -1,1
               ldn = lbd(i-1:i+1,j,1)
               pcn = pc(i-1:i+1,j)
               comp(1) = (one+deps)*s(i+n,j,1)
               dcx(n+2) = deps*s(i+n,j,1)
               if (comp(1)  > rhoval(1)) then
                  comp(1) = (one-deps)*s(i+n,j,1)
                  dcx(n+2) = - dcx(n+2)
               end if
               if (ncomps > 1) then 
                  comp(2) = rhoval(2)*(one - comp(1)/rhoval(1))
               end if
               l_kr_type  = nint(krcoef(i+n,j,1))
               l_kr_coef  = krcoef(i+n,j,2)
               l_kr_sres  = krcoef(i+n,j,3)
               call get_ld(ldtmp,comp,rhoval,muval,ncomps,l_kr_type,l_kr_coef,l_kr_sres)
               ldn(n+2) = ldtmp(1)
               l_cpl_type  = nint(cplcoef(i+n,j,1))
               l_cpl_coef  = cplcoef(i+n,j,2)
               l_sigma     = cplcoef(i+n,j,3)
               l_cpl_sres  = cplcoef(i+n,j,4)
               call get_cpl(pcn(n+2),comp,rhoval,phi(i+n,j),kappa(i+n,j),l_sigma,ncomps,
     &              l_cpl_sres,l_cpl_type,l_cpl_coef)
               if (do_upwind.eq.1) then
                  ldx(1) = merge(ldx(1),
     &                 merge(kpx(i,j)*ldn(1),kpx(i,j)*ldn(2),umac(i,j)>=zero),ldx(1)==zero)
                  ldx(2) = merge(ldx(2),
     &                 merge(kpx(i+1,j)*ldn(2),kpx(i+1,j)*ldn(3),umac(i+1,j)>=zero),ldx(2)==zero)
               else
                  ldx(1) = merge(ldx(1),kpx(i,j)*half*(ldn(1)+ldn(2)),ldx(1)==zero)
                  ldx(2) = merge(ldx(2),kpx(i+1,j)*half*(ldn(2)+ldn(3)),ldx(2)==zero)
               end if
               if (i .eq. domlo(1) .and. bc(1,1) .eq. EXT_DIR) then
                  ldx(1) = merge(ldx(1),kpx(i,j)*ldn(1),ldx(1)==zero)
               elseif (i .eq. domhi(1) .and. bc(1,2) .eq. EXT_DIR) then
                  ldx(2) = merge(ldx(2),kpx(i+1,j)*ldn(3),ldx(2)==zero)
               end if

               rgtmpx(n+2) = (ldx(2)*(-multx(2)*(pcn(3)-pcn(2))/dx(1) - gravx) -
     &           ldx(1)*(-multx(1)*(pcn(2)-pcn(1))/dx(1) - gravx))/dx(1) 
            end do

c           y-direction
            do n = -1,1
               ldn = lbd(i,j-1:j+1,1)
               pcn = pc(i,j-1:j+1)
               comp(1) = (one+deps)*s(i,j+n,1)
               dcy(n+2) = deps*s(i,j+n,1)
               if (comp(1)  > rhoval(1)) then
                  comp(1) = (one-deps)*s(i,j+n,1)
                  dcy(n+2) = - dcy(n+2)
               end if  
               if (ncomps > 1) then 
                  comp(2) = rhoval(2)*(one - comp(1)/rhoval(1))
               end if
               l_kr_type  = nint(krcoef(i,j+n,1))
               l_kr_coef  = krcoef(i,j+n,2)
               l_kr_sres  = krcoef(i,j+n,3)
               call get_ld(ldtmp,comp,rhoval,muval,ncomps,l_kr_type,l_kr_coef,l_kr_sres)
               ldn(n+2) = ldtmp(1)
               l_cpl_type  = nint(cplcoef(i,j+n,1))
               l_cpl_coef  = cplcoef(i,j+n,2)
               l_sigma     = cplcoef(i,j+n,3)
               l_cpl_sres  = cplcoef(i,j+n,4)
               call get_cpl(pcn(n+2),comp,rhoval,phi(i,j+n),kappa(i,j+n),l_sigma,ncomps,
     &              l_cpl_sres,l_cpl_type,l_cpl_coef)
               if (do_upwind.eq.1) then
                  ldy(1) = merge(ldy(1),merge(kpy(i,j)*ldn(1),kpy(i,j)*ldn(2),vmac(i,j)>=zero),ldy(1)==zero)
                  ldy(2) = merge(ldy(2),merge(kpy(i,j+1)*ldn(2),kpy(i,j+1)*ldn(3),vmac(i,j+1)>=zero),ldy(2)==zero)
               else
                  ldy(1) = merge(ldy(1),kpy(i,j)*half*(ldn(1)+ldn(2)),ldy(1)==zero)
                  ldy(2) = merge(ldy(2),kpy(i,j+1)*half*(ldn(2)+ldn(3)),ldy(2)==zero)
               end if
               if (j .eq. domlo(2) .and. bc(2,1) .eq. EXT_DIR) then
                  ldy(1) = merge(ldy(1),kpy(i,j)*ldn(1),ldy(1)==zero)
               elseif (j .eq. domhi(2) .and. bc(2,2) .eq. EXT_DIR) then
                  ldy(2) = merge(ldy(2),kpy(i,j+1)*ldn(3),ldy(2)==zero)
               end if

               rgtmpy(n+2) =  (ldy(2)*(-multy(2)*(pcn(3)-pcn(2))/dx(2) - gravy) - 
     &           ldy(1)*(-multy(1)*(pcn(2)-pcn(1))/dx(2) - gravy))/dx(2)                
            end do

            bx(i,j,2) = (rgtmpx(1) - rgx)/dcx(1)
            bx(i,j,1) = (rgtmpx(3) - rgx)/dcx(3)
            bx(i,j,3) = (rgtmpx(2) - rgx)/dcx(2)
            by(i,j,2) = (rgtmpy(1) - rgy)/dcy(1)
            by(i,j,1) = (rgtmpy(3) - rgy)/dcy(3)
            by(i,j,3) = (rgtmpy(2) - rgy)/dcy(2)

         end do
      end do
      
      end

c :: ------------------------------------------------------------
c :: RICHARD_NJAC2
c      Compute the coefficient for the numerical jacobian in the  
c      equilibrium richard's equation (w.r.t to p).
c
c     dfx(i,j,1) = dDF(i,j)/dp(i-1,j)
c     dfx(i,j,2) = dDF(i,j)/dp(i+1,j)
c     dfx(i,j,3) = dDF(i,j)/dp(i,j) 
c     dfy(i,j,1) = dDF(i,j)/dp(i,j-1)
c     dfy(i,j,2) = dDF(i,j)/dp(i,j+1)
c     dfy(i,j,3) = dDF(i,j)/dp(i,j)
c
c :: ------------------------------------------------------------    
      subroutine FORT_RICHARD_NJAC2(
     &     bx,DIMS(bx),by,DIMS(by),
     &     umac, DIMS(umac), vmac, DIMS(vmac),
     &     kpx, DIMS(kpx), kpy, DIMS(kpy),
     &     lbd,DIMS(lbd),p,DIMS(p),
     &     phi,DIMS(phi),kappa,DIMS(kappa),
     &     krcoef,DIMS(krcoef),nkrcoef,
     &     cplcoef,DIMS(cplcoef),ncplcoef,
     &     lo,hi,domlo,domhi,dx, bc, 
     &     inflow_bc_lo,inflow_bc_hi,deps,do_upwind)

      implicit none

#include "probdata.H"

      integer DIMDEC(bx)
      integer DIMDEC(by)
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(kpx)
      integer DIMDEC(kpy)
      integer DIMDEC(p)
      integer DIMDEC(lbd)
      integer DIMDEC(phi)
      integer DIMDEC(kappa)
      integer DIMDEC(krcoef)
      integer DIMDEC(cplcoef)
      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      integer bc(SDIM,2)
      integer inflow_bc_lo(SDIM),inflow_bc_hi(SDIM)
      integer nkrcoef,ncplcoef
      REAL_T  dx(SDIM)
      REAL_T  bx(DIMV(bx),3)
      REAL_T  by(DIMV(by),3)
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  kpx(DIMV(kpx))
      REAL_T  kpy(DIMV(kpy))
      REAL_T  p(DIMV(p))
      REAL_T  lbd(DIMV(lbd),ncomps)
      REAL_T  cplcoef(DIMV(cplcoef),ncplcoef)
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  phi(DIMV(phi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  deps
      integer do_upwind

      integer i, j,n,l_kr_type,l_cpl_type
      REAL_T l_kr_coef, l_kr_sres, l_cpl_coef,l_sigma,l_cpl_sres
      REAL_T  gravx, gravy, rgx, rgy
      REAL_T  multx(2), multy(2), ldx(2), ldy(2)
      REAL_T  comp(ncomps),cap
      REAL_T  ldtmp(ncomps),ldn(3),pn(3),rgtmpx(3),rgtmpy(3), dcx(3), dcy(3)

      if (gravity_dir.eq.0) then
         gravx = -rhoval(1)*gravity
         gravy = 0.d0
      else if (gravity_dir.eq.1) then
         gravx = 0.d0
         gravy = -rhoval(1)*gravity
      endif

      bx = zero
      by = zero
      do j=lo(2), hi(2)
         do i=lo(1), hi(1)

            if (do_upwind.eq.1) then
               ldx(1) = merge(kpx(i,j)*lbd(i-1,j,1),kpx(i,j)*lbd(i,j,1),umac(i,j)>=zero)
               ldx(2) = merge(kpx(i+1,j)*lbd(i,j,1),kpx(i+1,j)*lbd(i+1,j,1),umac(i+1,j)>=zero)
               ldy(1) = merge(kpy(i,j)*lbd(i,j-1,1),kpy(i,j)*lbd(i,j,1),vmac(i,j)>=zero)
               ldy(2) = merge(kpy(i,j+1)*lbd(i,j,1),kpy(i,j+1)*lbd(i,j+1,1),vmac(i,j+1)>=zero)
            else
               ldx(1) = kpx(i,j)*half*(lbd(i,j,1)+lbd(i-1,j,1))
               ldx(2) = kpx(i+1,j)*half*(lbd(i,j,1)+lbd(i+1,j,1))
               ldy(1) = kpy(i,j)*half*(lbd(i,j,1)+lbd(i,j-1,1))
               ldy(2) = kpy(i,j+1)*half*(lbd(i,j,1)+lbd(i,j+1,1))
            end if
            multx = one
            multy = one
            if (i .eq. domlo(1)) then
               if (bc(1,1) .eq. EXT_DIR) then
                  multx(1) = two 
                  ldx(1) = kpx(i,j)*lbd(i-1,j,1)
               elseif (bc(1,1) .eq. FOEXTRAP) then
                  multx(1) = zero
                  if (inflow_bc_lo(1) .eq. 1) then
                     ldx(1) = zero
                  end if
               end if
            end if
            if (i .eq. domhi(1)) then
               if (bc(1,2) .eq. EXT_DIR) then
                  multx(2) = two
                  ldx(2) = kpx(i+1,j)*lbd(i+1,j,1)
               elseif (bc(1,2) .eq. FOEXTRAP) then
                  multx(2) = zero
                  if (inflow_bc_hi(1) .eq. 1) then
                     ldx(2) = zero
                  end if
               end if
            end if
            if (j .eq. domlo(2)) then
               if (bc(2,1) .eq. EXT_DIR) then
                  multy(1) = two
                  ldy(1) = kpy(i,j)*lbd(i,j-1,1)
               elseif (bc(2,1) .eq. FOEXTRAP) then
                  multy(1) = zero
                  if (inflow_bc_lo(2) .eq. 1) then
                     ldy(1) = zero
                  end if
               end if
            end if
            if (j .eq. domhi(2)) then
               if (bc(2,2) .eq. EXT_DIR) then
                  multy(2) = two 
                  ldy(2) = kpy(i,j+1)*lbd(i,j+1,1)
               elseif (bc(2,2) .eq. FOEXTRAP) then
                  multy(2) = zero
                  if (inflow_bc_hi(2) .eq. 1) then
                     ldy(2) = zero
                  end if
               end if
            end if

            rgx = (ldx(2)*(multx(2)*(p(i+1,j)-p(i,j))/dx(1) - gravx) -
     &           ldx(1)*(multx(1)*(p(i,j)-p(i-1,j))/dx(1) - gravx))/dx(1) 
            rgy = (ldy(2)*(multy(2)*(p(i,j+1)-p(i,j))/dx(2) - gravy) - 
     &           ldy(1)*(multy(1)*(p(i,j)-p(i,j-1))/dx(2) - gravy))/dx(2)           

c           x-direction
            do n = -1,1 
               l_kr_type  = nint(krcoef(i+n,j,1))
               l_kr_coef  = krcoef(i+n,j,2)
               l_kr_sres  = krcoef(i+n,j,3)
               l_cpl_type  = nint(cplcoef(i+n,j,1))
               l_cpl_coef  = cplcoef(i+n,j,2)
               l_sigma     = cplcoef(i+n,j,3)
               l_cpl_sres  = cplcoef(i+n,j,4)
               ldn = lbd(i-1:i+1,j,1)
               pn  = p(i-1:i+1,j)
               pn(n+2)  = pn(n+2) + deps/l_sigma
               dcx(n+2) = deps/l_sigma
               cap = -pn(n+2)
               call get_invcpl(comp,cap,rhoval,phi(i+n,j),kappa(i+n,j),l_sigma,
     &              ncomps,l_cpl_sres,l_cpl_type,l_cpl_coef)
               call get_ld(ldtmp,comp,rhoval,muval,ncomps,l_kr_type,l_kr_coef,l_kr_sres)
               ldn(n+2) = ldtmp(1)
               if (do_upwind.eq.1) then
                  ldx(1) = merge(ldx(1),merge(kpx(i,j)*ldn(1),kpx(i,j)*ldn(2),umac(i,j)>=zero),ldx(1)==zero)
                  ldx(2) = merge(ldx(2),merge(kpx(i+1,j)*ldn(2),kpx(i+1,j)*ldn(3),umac(i+1,j)>=zero),ldx(2)==zero)
               else
                  ldx(1) = merge(ldx(1),kpx(i,j)*half*(ldn(1)+ldn(2)),ldx(1)==zero)
                  ldx(2) = merge(ldx(2),kpx(i+1,j)*half*(ldn(2)+ldn(3)),ldx(2)==zero)
               end if
               if (i .eq. domlo(1) .and. bc(1,1) .eq. EXT_DIR) then
                  ldx(1) = merge(ldx(1),kpx(i,j)*ldn(1),ldx(1)==zero)
               elseif (i .eq. domhi(1) .and. bc(1,2) .eq. EXT_DIR) then
                  ldx(2) = merge(ldx(2),kpx(i+1,j)*ldn(3),ldx(2)==zero)
               end if

               rgtmpx(n+2) = (ldx(2)*(multx(2)*(pn(3)-pn(2))/dx(1) - gravx) -
     &           ldx(1)*(multx(1)*(pn(2)-pn(1))/dx(1) - gravx))/dx(1) 
            end do

c           y-direction
            do n = -1,1
               l_kr_type  = nint(krcoef(i,j+n,1))
               l_kr_coef  = krcoef(i,j+n,2)
               l_kr_sres  = krcoef(i,j+n,3)
               l_cpl_type = nint(cplcoef(i,j+n,1))
               l_cpl_coef = cplcoef(i,j+n,2)
               l_sigma    = cplcoef(i,j+n,3)
               l_cpl_sres = cplcoef(i,j+n,4)
               ldn = lbd(i,j-1:j+1,1)
               pn  = p(i,j-1:j+1)
               pn(n+2)  = pn(n+2) + deps/l_sigma
               dcy(n+2) = deps/l_sigma
               cap = -pn(n+2)
               call get_invcpl(comp,cap,rhoval,phi(i,j+n),kappa(i,j+n),l_sigma,
     &              ncomps,l_cpl_sres,l_cpl_type,l_cpl_coef)
               call get_ld(ldtmp,comp,rhoval,muval,ncomps,l_kr_type,l_kr_coef,l_kr_sres)
               ldn(n+2) = ldtmp(1)
               if (do_upwind.eq.1) then
                  ldy(1) = merge(ldy(1),merge(kpy(i,j)*ldn(1),kpy(i,j)*ldn(2),vmac(i,j)>=zero),ldy(1)==zero)
                  ldy(2) = merge(ldy(2),merge(kpy(i,j+1)*ldn(2),kpy(i,j+1)*ldn(3),vmac(i,j+1)>=zero),ldy(2)==zero)
               else
                  ldy(1) = merge(ldy(1),kpy(i,j)*half*(ldn(1)+ldn(2)),ldy(1)==zero)
                  ldy(2) = merge(ldy(2),kpy(i,j+1)*half*(ldn(2)+ldn(3)),ldy(2)==zero)
               end if
               if (j .eq. domlo(2) .and. bc(2,1) .eq. EXT_DIR) then
                  ldy(1) = merge(ldy(1),kpy(i,j)*ldn(1),ldy(1)==zero)
               elseif (j .eq. domhi(2) .and. bc(2,2) .eq. EXT_DIR) then
                  ldy(2) = merge(ldy(2),kpy(i,j+1)*ldn(3),ldy(2)==zero)
               end if
               rgtmpy(n+2) = (ldy(2)*(multy(2)*(pn(3)-pn(2))/dx(2) - gravy) - 
     &           ldy(1)*(multy(1)*(pn(2)-pn(1))/dx(2) - gravy))/dx(2)                

            end do

            bx(i,j,2) = (rgtmpx(1) - rgx)/dcx(1)
            bx(i,j,1) = (rgtmpx(3) - rgx)/dcx(3)
            bx(i,j,3) = (rgtmpx(2) - rgx)/dcx(2)
            by(i,j,2) = (rgtmpy(1) - rgy)/dcy(1)
            by(i,j,1) = (rgtmpy(3) - rgy)/dcy(3) 
            by(i,j,3) = (rgtmpy(2) - rgy)/dcy(2)
         end do
      end do
      end

c ::: -----------------------------------------------------------
c ::: Compute dn/dp = -dn/dpc
c ::: -----------------------------------------------------------

      subroutine FORT_RICHARD_ALPHA(alpha,DIMS(alpha),s,DIMS(s),
     &     phi,DIMS(phi),kappa,DIMS(kappa),
     &     coef,DIMS(coef),n_coef,lo,hi)

      implicit none

#include "probdata.H"

      integer DIMDEC(alpha)
      integer DIMDEC(s)
      integer DIMDEC(phi)
      integer DIMDEC(kappa)
      integer DIMDEC(coef)
      integer n_coef
      integer lo(SDIM),hi(SDIM)
      REAL_T  alpha(DIMV(alpha))
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  phi(DIMV(phi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  coef(DIMV(coef),n_coef)

      integer i,j
      integer l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res, dpcdn
      REAL_T  rtmp(ncomps), stmp(ncomps), sat(ncomps)

      rtmp = one
      do j = lo(2)-1,hi(2)+1
         do i = lo(1)-1,hi(1)+1
            l_cpl_type = nint(coef(i,j,1))
            l_cpl_coef = coef(i,j,2)
            l_sigma    = coef(i,j,3)
            l_sat_res  = coef(i,j,4)

            stmp = zero
            stmp(1:ncomps) = s(i,j,1:ncomps) 

            call get_st(sat,stmp,rhoval,l_sat_res,ncomps)

            if (dabs(one-sat(1)) .lt. 1.d-6) then
               sat(1) = one - 1.d-6
            end if
            if (ncomps > 1) then
               sat(2) = 1.d-6
            end if

            call get_dcplds(dpcdn,sat,rtmp,phi(i,j),
     &           kappa(i,j),l_sigma,ncomps,zero,l_cpl_type,l_cpl_coef)
         
            alpha(i,j) =  - phi(i,j)/dpcdn*rhoval(1)*(one-l_sat_res)
         end do 
      end do

      end

c
c    ----------------------------------------------------------
c    RICHARD_VEL_BC
c      Set the richard's velocity boundary condition
c      res = res + \nabla vel at boundary
c    ----------------------------------------------------------
c  
      subroutine FORT_RICHARD_VELBC(rg,DIMS(rg),
     &     ux,DIMS(ux),uy,DIMS(uy),
     &     lo,hi,domlo,domhi,dx,
     &     inflow_bc_lo,inflow_bc_hi,mult)

      implicit none
      
      integer DIMDEC(rg)
      integer DIMDEC(ux)
      integer DIMDEC(uy)
      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      integer inflow_bc_lo(SDIM),inflow_bc_hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  rg(DIMV(rg))
      REAL_T  ux(DIMV(ux))
      REAL_T  uy(DIMV(uy))
      REAL_T  mult
      
      integer i, j
      if (lo(1) .eq. domlo(1) .and. inflow_bc_lo(1) .eq. 1) then
           do j = lo(2),hi(2)
              rg(lo(1),j) = rg(lo(1),j) - mult*ux(lo(1),j)/dx(1)
           end do
      end if
      if (hi(1) .eq. domhi(1) .and. inflow_bc_hi(1) .eq. 1) then
           do j = lo(2),hi(2)
              rg(hi(1),j) = rg(hi(1),j) + mult*ux(hi(1)+1,j)/dx(1)
           end do
        end if
      if (lo(2) .eq. domlo(2) .and. inflow_bc_lo(2) .eq. 1) then
           do i = lo(1),hi(1)
              rg(i,lo(2)) = rg(i,lo(2)) - mult*uy(i,lo(2))/dx(2)
           end do
      end if
      if (hi(2) .eq. domhi(2) .and. inflow_bc_hi(2) .eq. 1) then
           do i = lo(1),hi(1)
              rg(i,hi(2)) = rg(i,hi(2)) + mult*uy(i,hi(2)+1)/dx(2)
           end do
      end if
   
      end     

c ::: ----------------------------------------------------------
c ::: Average the edge-centered velocities to a cell-centered quantity.
c ::: ----------------------------------------------------------
      subroutine FORT_AVG_UMAC(umac,DIMS(umac),vmac,DIMS(vmac),
     &                         vel,DIMS(vel),lo,hi,upwind_flag)

      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(vel)
      integer lo(SDIM),hi(SDIM),upwind_flag
      REAL_T umac(DIMV(umac))
      REAL_T vmac(DIMV(vmac))
      REAL_T  vel(DIMV(vel),SDIM)

      integer i,j

      if (upwind_flag .eq. 0) then
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vel(i,j,1) = half * (umac(i,j) + umac(i+1,j))
               vel(i,j,2) = half * (vmac(i,j) + vmac(i,j+1))
            end do
         end do
      else
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (umac(i,j) * umac(i+1,j) .le. 0.d0) then 
                  vel(i,j,1) = half * (umac(i,j) + umac(i+1,j))
               else
                  vel(i,j,1) = merge(umac(i,j),umac(i+1,j),umac(i,j).gt.0.d0)
               endif

               if (vmac(i,j) * vmac(i,j+1) .le. 0.d0) then 
                  vel(i,j,2) = half * (vmac(i,j) + vmac(i,j+1))
               else
                  vel(i,j,2) = merge(vmac(i,j),vmac(i,j+1),vmac(i,j).gt.0.d0)
               endif
            end do
         end do
      endif
      end

c ::: -------------------------------------------------------------
c ::: Copy the edge-centered velocities to a cell-centered quantity.
c ::: -------------------------------------------------------------
      subroutine FORT_CPY_UMAC(umac,DIMS(umac),vmac,DIMS(vmac),
     &                         vel,DIMS(vel),lo,hi,shift)

      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(vel)
      integer lo(SDIM),hi(SDIM)
      integer shift
      REAL_T umac(DIMV(umac))
      REAL_T vmac(DIMV(vmac))
      REAL_T  vel(DIMV(vel),SDIM)

      integer i,j


      do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          vel(i,j,1) = umac(i+shift,j) 
          vel(i,j,2) = vmac(i,j+shift)
        end do
      end do

      end


c ::: ---------------------------------------------------------
c ::: Define linear pressure profile
c ::: ---------------------------------------------------------
      
      subroutine FORT_LINEAR_PRESSURE(lo,hi,p,DIMS(p),ncomps,
     &     dx,problo,probhi,ref_val,ref_loc,gradp)

      implicit none

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(p)
      integer ncomps
      REAL_T  p(DIMV(p))
      REAL_T  rhoval(ncomps)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM)
      REAL_T  probhi(SDIM)
      REAL_T  ref_val
      REAL_T  ref_loc(SDIM)
      REAL_T  gradp(SDIM)
      
      integer i,j
      REAL_T  dp(SDIM), x, y

      do j = lo(2),hi(2)
         y = MAX(problo(2), MIN(probhi(2), (dble(j)+half)*dx(2) + problo(2)))
         dp(2) = (y-ref_loc(2))*gradp(2)
         do i = lo(1),hi(1)
            x = MAX(problo(1), MIN(probhi(1), (dble(i)+half)*dx(1) + problo(1)))
            dp(1) = (x-ref_loc(1))*gradp(1)
            p(i,j) = ref_val + dp(1) + dp(2)
         end do
      end do 
      end


c ::: ---------------------------------------------------------
c ::: Define hydraulic head
c ::: ---------------------------------------------------------
 
      subroutine FORT_HYD_HEAD(lo, hi, ph, DIMS(ph),
     &     rhog, gdir, dx, plo, nc)

      implicit none

      integer DIMDEC(ph)
      integer nc
      integer lo(SDIM), hi(SDIM)
      REAL_T  ph(DIMV(ph),nc)
      REAL_T  dx(SDIM)
      REAL_T  plo(SDIM)
      REAL_T  rhog(nc)
      integer gdir

      integer i,j,n
      REAL_T  x, y, fac

      if (gdir.eq.0) then 
         do n = 1,nc
            fac = 1.d0 / rhog(n)
            do j = ARG_L2(ph), ARG_H2(ph)
               do i = ARG_L1(ph), ARG_H1(ph)
                  x = (dble(i)+half)*dx(1) + plo(1)
                  ph(i,j,n) = x  +  ph(i,j,n) * fac
               enddo
            enddo
         enddo
      else if (gdir.eq.1) then 
         do n = 1,nc
            fac = 1.d0 / rhog(n)
            do j = ARG_L2(ph), ARG_H2(ph)
               y = (dble(j)+half)*dx(2) + plo(2)
               do i = ARG_L1(ph), ARG_H1(ph)
                  ph(i,j,n) = y  +  ph(i,j,n) * fac
               enddo
            enddo
         enddo
      else
         do n = 1,nc
            fac = 1.d0 / rhog(n)
            do j = ARG_L2(ph), ARG_H2(ph)
               do i = ARG_L1(ph), ARG_H1(ph)
                  ph(i,j,n) = ph(i,j,n) * fac
               enddo
            enddo
         enddo
      endif
      end


c ::: ------------------------------------------------------------------
c ::: Define approximate steady state condition based on inflow velocity
c ::: ------------------------------------------------------------------
      
      subroutine FORT_STEADYSTATE(rho,DIMS(rho),rhoval,muval,ncomps,
     &     kappa,DIMS(kappa),matID,DIMS(matID),rmID,time,inflow_vel,
     &     nc,grav,lo,hi)

      implicit none

      integer nc, rmID
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(rho)
      integer DIMDEC(kappa)
      integer DIMDEC(matID)
      integer ncomps
      REAL_T  rho(DIMV(rho),ncomps)
      REAL_T  kappa(DIMV(kappa))
      integer matID(DIMV(matID))
      REAL_T  rhoval(ncomps)
      REAL_T  muval(ncomps)
      REAL_T  inflow_vel
      REAL_T  grav, time
      
      integer i,j,l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T  vtot, gstar, sol

      vtot  = zero
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (ncomps == 1) then
               gstar = -kappa(i,j)*rhoval(1)*grav
               call FORT_FIND_INV_RFLUX(sol,inflow_vel,gstar,muval,
     &              ncomps,rmID,matID,time)
            else
               gstar = -kappa(i,j)*(rhoval(1)-rhoval(2))*grav
               call FORT_FIND_INV_FLUX(sol,inflow_vel,nc,vtot,gstar,muval,
     &              ncomps,rmID,matID,time)
            end if

            rho(i,j,1) = max(rho(i,j,1),rhoval(1)*sol)

            if (ncomps > 1) then
               rho(i,j,2) = min(rho(i,j,2),rhoval(2)*(1.0-sol))
            endif

         end do
      end do 

      end

#define NCROSS 1
#define DXY 1

      subroutine FORT_TENSORDIFF(lo, hi, a,DIMS(a), Dcc, DIMS(Dcc),
     &                      u, DIMS(u), bx, DIMS(bx), b1x, DIMS(b1x),
     &                      v, DIMS(v), by, DIMS(by), b1y, DIMS(b1y),
     &                      bc)

      implicit none

      integer nc, lo(SDIM), hi(SDIM), bc(SDIM,2)
      integer DIMDEC(a)
      integer DIMDEC(Dcc)
      integer DIMDEC(u)
      integer DIMDEC(bx)
      integer DIMDEC(b1x)
      integer DIMDEC(v)
      integer DIMDEC(by)
      integer DIMDEC(b1y)
      REAL_T  a(DIMV(a),2)
      REAL_T  Dcc(DIMV(Dcc))
      REAL_T  u(DIMV(u))
      REAL_T  bx(DIMV(bx))
      REAL_T  b1x(DIMV(b1x),NCROSS)
      REAL_T  v(DIMV(v))
      REAL_T  by(DIMV(by))
      REAL_T  b1y(DIMV(b1y),NCROSS)
      
      integer i,j
      REAL_T u00, u01, u10, u11
      REAL_T v00, v01, v10, v11
      REAL_T uface, vface, velInv, vel, aLface, aTface, Dface

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            uface = u(i,j)
            if (i.eq.lo(1) .and. bc(1,1).ne.INT_DIR) then
               v00 = v(i  ,j  )
               v01 = v(i  ,j+1)
            else
               v00 = v(i-1,j  )
               v01 = v(i-1,j+1)
            endif
            if (i.eq.hi(1)+1 .and. bc(1,2).ne.INT_DIR) then
               v10 = v(i-1,j  )
               v11 = v(i-1,j+1)
            else
               v10 = v(i  ,j  )
               v11 = v(i  ,j+1)
            endif
            vface = 0.25d0*(v00+v01+v10+v11)
            aLface = 0.5d0*(a(i,j,1)+a(i-1,j,1))
            aTface = 0.5d0*(a(i,j,2)+a(i-1,j,2))
            Dface = 0.5d0*(Dcc(i,j)+Dcc(i-1,j))
            vel = sqrt(uface*uface + vface*vface)
            if (vel .ne. 0.d0) then
               velInv = 1.d0 / vel
            else
               velInv = 0.d0
            endif
            bx(i,j) = Dface + aTface*vel + (aLface - aTface)*uface*uface*velInv
            b1x(i,j,DXY) = (aLface - aTface)*uface*vface*velInv
         enddo
      enddo

      do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)

            if (j.eq.lo(2) .and. bc(2,1).ne.INT_DIR) then
               u00 = u(i  ,j  )
               u10 = u(i+1,j  )
            else
               u00 = u(i  ,j-1)
               u10 = u(i+1,j-1)
            endif
            if (j.eq.hi(2)+1 .and. bc(2,2).ne.INT_DIR) then
               u01 = u(i  ,j-1)
               u11 = u(i+1,j-1)
            else
               u01 = u(i  ,j  )
               u11 = u(i+1,j  )
            endif
            uface = 0.25d0*(u00+u10+u01+u11)
            vface = v(i,j)
            aLface = 0.5d0*(a(i,j,1)+a(i,j-1,1))
            aTface = 0.5d0*(a(i,j,2)+a(i,j-1,2))
            Dface = 0.5d0*(Dcc(i,j)+Dcc(i,j-1))
            vel = sqrt(uface*uface + vface*vface)
            if (vel .ne. 0.d0) then
               velInv = 1.d0 / vel
            else
               velInv = 0.d0
            endif
            by(i,j) = Dface + aTface*vel + (aLface - aTface)*vface*vface*velInv
            b1y(i,j,DXY) = (aLface - aTface)*uface*vface*velInv
         enddo
      enddo
      end

c :: ----------------------------------------------------------
c :: VISCEXTRAP
c ::             extrapolate viscous terms to ghost zones
c ::
c :: INPUTS / OUTPUTS:
c ::  visc        <=> array containing visc terms
c ::  DIMS(visc)   => index limits for visc
c ::  lo,hi        => index limits for grid interior
c ::  ncomp        => number of components
c :: ----------------------------------------------------------
c ::
      subroutine FORT_VISCEXTRAP (visc,DIMS(visc),lo,hi,ncomp)

      integer DIMDEC(visc)
      integer lo(SDIM), hi(SDIM)
      integer ncomp
      REAL_T  visc(DIMV(visc),ncomp)

      integer i, j, n

      do n = 1, ncomp
c        ::::: left side
         i = lo(1)-1
         do j = lo(2), hi(2)
c           visc(i,j,n) = two*visc(i+1,j,n) - visc(i+2,j,n)
            visc(i,j,n) = visc(i+1,j,n)
         end do

c        ::::: right side
         i = hi(1)+1
         do j = lo(2), hi(2)
c           visc(i,j,n) = two*visc(i-1,j,n) - visc(i-2,j,n)
            visc(i,j,n) = visc(i-1,j,n)
         end do

c        ::::: bottom side
         j = lo(2)-1
         do i = lo(1), hi(1)
c           visc(i,j,n) = two*visc(i,j+1,n) - visc(i,j+2,n)
            visc(i,j,n) = visc(i,j+1,n)
         end do

c        ::::: top side
         j = hi(2)+1
         do i = lo(1), hi(1)
c           visc(i,j,n) = two*visc(i,j-1,n) - visc(i,j-2,n)
            visc(i,j,n) = visc(i,j-1,n)
         end do

c        ::::: corners
         i = lo(1)-1
         j = lo(2)-1
c        visc(i,j,n) = two*visc(i+1,j+2,n) - visc(i+2,j+2,n)
         visc(i,j,n) = visc(i+1,j+1,n)

         j = hi(2)+1
c        visc(i,j,n) = two*visc(i+1,j-1,n) - visc(i+2,j-2,n)
         visc(i,j,n) = visc(i+1,j-1,n)

         i = hi(1)+1
c        visc(i,j,n) = two*visc(i-1,j-1,n) - visc(i-2,j-2,n) 
         visc(i,j,n) = visc(i-1,j-1,n)

         j = lo(2)-1
c        visc(i,j,n) = two*visc(i-1,j+1,n) - visc(i-2,j+2,n)
         visc(i,j,n) = visc(i-1,j+1,n)
      end do
       
      end

      subroutine FORT_INITKEDGE (cdat, DIMS(cdat), edat, DIMS(edat), 
     &     lo, hi, dir)
      implicit none
      integer DIMDEC(cdat)
      integer DIMDEC(edat)
      integer lo(SDIM), hi(SDIM)
      integer dir
c     FIXME: add support for vector kappa
      REAL_T  cdat(DIMV(cdat),SDIM)
      REAL_T  edat(DIMV(edat))

      integer i,j,nc

      nc = dir+1
      if (dir .eq. 0) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               edat(i,j) = 2.d0 / 
     &              (1.d0/cdat(i,j,nc) + 1.d0/cdat(i-1,j,nc))
            enddo
         enddo
      else
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
               edat(i,j) = 2.d0 / 
     &              (1.d0/cdat(i,j,nc) + 1.d0/cdat(i,j-1,nc))
            enddo
         enddo
      endif
      end

