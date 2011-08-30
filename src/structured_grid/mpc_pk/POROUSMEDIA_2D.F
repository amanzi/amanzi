c
c $Id: POROUSMEDIA_2D.F,v 1.65 2011-08-15 17:53:23 gpau Exp $
c

#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "PROB_PM_F.H"
#include "POROUSMEDIA_F.H"
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
c ::: SUMMASSSQ
c :::             MASS = sum{ vol(i,j)*rho(i,j)*rho(i,j) }
c :::
c ::: INPUTS / OUTPUTS:
c :::  rho        => density field
c :::  DIMS(rho)  => index limits of rho aray
c :::  lo,hi      => index limits of grid interior
c :::  ds   	 => cell size
c :::  mass      <=  total mass
c ::: ----------------------------------------------------------
       subroutine FORT_SUMMASSSQ(rho,DIMS(rho),DIMS(grid),dx,mass)

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, dx(SDIM)
       REAL_T  rho(DIMV(rho))

       integer i, j
       REAL_T  vol

       vol = dx(1)*dx(2)

       mass = zero
       do i = ARG_L1(grid), ARG_H1(grid)
          do j = ARG_L2(grid), ARG_H2(grid)
             mass = mass + vol*rho(i,j)*rho(i,j)
          end do
       end do

       end

c ::: ----------------------------------------------------------
c ::: SUMPRESS
c :::             SUM = sum{ press(i,j) * surroundingVolume(i,j) }
c :::
c ::: INPUTS / OUTPUTS:
c :::  press      => pressure field on nodes
c :::  DIMS(press)=> index limits of press aray
c :::  vol        => volume field (zero-ed out under finer grids)
c :::  DIMS(vol)  => index limits of vol aray
c :::  lo,hi      => index limits of grid interior in nodal space
c :::  total_press<= sum of volume-weighted nodal pressures
c :::  total_vol  <= sum of volumes
c ::: ----------------------------------------------------------

       subroutine FORT_SUMPRESS(press,DIMS(press),vol,DIMS(vol),lo,hi,
     &                          total_press,total_vol)

       integer lo(2),hi(2)
       integer DIMDEC(press)
       integer DIMDEC(vol)
       REAL_T  press(DIMV(press))
       REAL_T  vol(DIMV(vol))
       REAL_T  total_press,total_vol

       integer i, j
       REAL_T  surroundingVolume

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             surroundingVolume = vol(i,j) + vol(i-1,j) + vol(i,j-1) + vol(i-1,j-1)
             total_press = total_press + press(i,j) * surroundingVolume
             total_vol   = total_vol   +              surroundingVolume

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

c ::: ---------------------------------------------------------------
c ::: For edge-based data, fill fine values with gradient-interpolated 
c ::: coarse data.
c ::: Operate only on faces that overlap--ie, only fill the fine faces  
c ::: that make up each coarse face, leave the in-between faces alone.
c ::: ----------------------------------------------------------------

      subroutine CV_EDGE_INTERP(lo, hi, domlo, domhi,
     &     nc, ratio, dir, dx,
     &     crse, crse_l0, crse_l1, crse_h0, crse_h1,
     &     fine, fine_l0, fine_l1, fine_h0, fine_h1)

      implicit none
      integer lo(2),hi(2), domlo(2),domhi(2)
      integer nc, ratio(0:2-1), dir
      integer crse_l0, crse_l1, crse_h0, crse_h1
      integer fine_l0, fine_l1, fine_h0, fine_h1
      DOUBLE PRECISION dx(0:2-1)
      DOUBLE PRECISION crse(crse_l0:crse_h0,crse_l1:crse_h1,nc)
      DOUBLE PRECISION fine(fine_l0:fine_h0,fine_l1:fine_h1,nc)

      integer i,j,ii,jj,n,L,idx
      DOUBLE PRECISION grad

      if (dir.eq.0) then
         do n=1,nc
            do j=lo(2)+1,hi(2)-1
               jj = ratio(1)*j
               do i=lo(1)+1,hi(1)-1
                  ii = ratio(0)*i
                  
c                 approximate d/dy
                  if (j .eq. domlo(2)) then
                     grad = (crse(i,j+1,n) - crse(i,j,n))
                  elseif (j .eq. domhi(2)) then
                     grad = (crse(i,j,n) - crse(i,j-1,n))
                  else
                     grad = (crse(i,j+1,n)-crse(i,j-1,n))/(2.0)
                  endif
                  do L=0,ratio(1)/2-1
                     idx = -ratio(1)/2 + L 
                     fine(ii,jj+L,n)=(crse(i,j,n)+0.25*dble(idx)*grad)
                     idx = ratio(1)/2 + L
                     fine(ii,jj+idx,n)=(crse(i,j,n)+0.25*dble(L+1)*grad)
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=lo(2)+1,hi(2)-1
               jj = ratio(1)*j
               do i=lo(1)+1,hi(1)-1
                  ii = ratio(0)*i

c                 approximate d/dx
                  if (i .eq. domlo(1)) then
                     grad = (crse(i+1,j,n) - crse(i,j,n))
                  elseif (i .eq. domhi(1)) then
                     grad = (crse(i,j,n) - crse(i-1,j,n))
                  else
                     grad = (crse(i+1,j,n)-crse(i-1,j,n))/(2.0)
                  endif
                  do L=0,ratio(0)/2-1
                     idx = -ratio(0)/2 + L 
                     fine(ii+L,jj,n) = (crse(i,j,n)+0.25*dble(idx)*grad)
                     idx = ratio(0)/2 + L
                     fine(ii+idx,jj,n)= (crse(i,j,n)+0.25*dble(L+1)*grad)
                  enddo
               enddo
            enddo
         enddo
      endif

      end


c ::: -------------------------------------------------------------------
c ::: Do linear in dir, pc transverse to dir, leave alone the fine values
c ::: lining up with coarse edges--assume these have been set to hold the 
c ::: values you want to interpolate to the rest.
c ::: -------------------------------------------------------------------

      subroutine CVF_EDGE_INTERP(clo, chi, 
     &     nc, ratio, cdx, fdx,
     &     crsex, crsex_l0, crsex_l1, crsex_h0, crsex_h1,
     &     crsey, crsey_l0, crsey_l1, crsey_h0, crsey_h1,
     &     finex, finex_l0, finex_l1, finex_h0, finex_h1,
     &     finey, finey_l0, finey_l1, finey_h0, finey_h1)

      implicit none
      integer clo(0:2-1), chi(0:2-1), nc, ratio(0:2-1)
      integer crsex_l0, crsex_l1, crsex_h0, crsex_h1
      integer crsey_l0, crsey_l1, crsey_h0, crsey_h1
      integer finex_l0, finex_l1, finex_h0, finex_h1
      integer finey_l0, finey_l1, finey_h0, finey_h1
      DOUBLE PRECISION cdx(0:2-1), fdx(0:2-1)
      DOUBLE PRECISION crsex(crsex_l0:crsex_h0,crsex_l1:crsex_h1,nc)
      DOUBLE PRECISION crsey(crsey_l0:crsey_h0,crsey_l1:crsey_h1,nc)
      DOUBLE PRECISION finex(finex_l0:finex_h0,finex_l1:finex_h1,nc)
      DOUBLE PRECISION finey(finey_l0:finey_h0,finey_l1:finey_h1,nc)
      integer i,j,n,ii,jj
      DOUBLE PRECISION df, r(3), s(3)

      do n=1,nc
         do j=clo(1),chi(1)
            jj = ratio(1)*j
            do i=clo(0),chi(0)
               ii = ratio(0)*i

               df = (crsex(i+1,j,n)-crsex(i,j,n))/cdx(0)+
     &              (crsey(i,j+1,n)-crsey(i,j,n))/cdx(1)
               
               r(1) = (finex(ii+1,jj,n)-finex(ii,jj,n))/fdx(0) + 
     &                (finey(ii,jj+1,n)-finey(ii,jj,n))/fdx(1) - 0.25*df
               r(2) = (finex(ii+2,jj,n)-finex(ii+1,jj,n))/fdx(0)  + 
     &                (finey(ii+1,jj+1,n)-finey(ii+1,jj,n))/fdx(1)- 0.25*df
               r(3) = (finex(ii+1,jj+1,n)-finex(ii,jj+1,n))/fdx(0)+ 
     &                (finey(ii,jj+2,n)-finey(ii,jj+1,n))/fdx(1)  - 0.25*df

               r = r*fdx(0)*fdx(0)

               s(1) = (-0.5*r(2)-r(1)-0.5*r(3))
               s(2) = (-0.75*r(2)-0.5*r(1)-0.25*r(3))
               s(3) = (-0.25*r(2)-0.5*r(1)-0.75*r(3))
               
               finex(ii+1,jj,n) = finex(ii+1,jj,n) - (s(2)-s(1))/fdx(0)
               finex(ii+1,jj+1,n) = finex(ii+1,jj+1,n) - (0.-s(3))/fdx(0)
               finey(ii,jj+1,n) = finey(ii,jj+1,n) - (s(3)-s(1))/fdx(1)
               finey(ii+1,jj+1,n) = finey(ii+1,jj+1,n) - (0.-s(2))/fdx(1)
            enddo
         enddo
      enddo

      end


c ::: ---------------------------------------------------------------
c ::: For edge-based data, fill coarse values with average of fine data 
c ::: one the same edge
c ::: ----------------------------------------------------------------

      subroutine EDGE_AVGDOWN(lo, hi, nc, ratio, dir,
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
                  crse(i,j,n) = zero
                  do L=0,ratio(1)-1
                     crse(i,j,n) = crse(i,j,n) + fine(ii,jj+L,n) 
                  enddo
                  crse(i,j,n) = crse(i,j,n)/ratio(1)
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do i=lo(1),hi(1)
               ii = ratio(0)*i            
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  crse(i,j,n) = zero
                  do L=0,ratio(0)-1
                     crse(i,j,n) = crse(i,j,n) + fine(ii+L,jj,n) 
                  enddo
                   crse(i,j,n) = crse(i,j,n)/ratio(0)
               enddo
            enddo
         enddo
      endif

      end

c ::: ----------------------------------------------------------
c ::: SPECTEMPVISC
c :::    compute variable diffusivity constant
c ::: ----------------------------------------------------------

      subroutine FORT_SPECTEMPVISC(lo,hi,
     &     S, DIMS(S),
     &     B, DIMS(B),
     &     phi, DIMS(phi),
     &     visc_coef,
     &     ncomp, diffuse)

      implicit none

      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(S)
      integer DIMDEC(B)
      integer DIMDEC(phi)
      integer ncomp, diffuse
      REAL_T  S(DIMV(S),ncomp)
      REAL_T  B(DIMV(B),ncomp)
      REAL_T  phi(DIMV(phi))
      REAL_T  visc_coef(ncomp)

      integer i, j, n
      REAL_T  rho

      if (diffuse .ne. -1) then
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)
               do n=1,ncomp
                  B(i,j,n) = visc_coef(n)
               end do
            end do
         end do         

      else

         do j=lo(2), hi(2)
            do i=lo(1), hi(1)
               rho = zero
               do n=1,ncomp
                  rho = rho + S(i,j,n) 
               end do
               do n=1,ncomp
                  B(i,j,n) = rho*phi(i,j)*visc_coef(n) 
               end do
            end do
         end do

      end if

      end

c ::: ----------------------------------------------------------
c ::: CALCEOS
c :::   compute \phi * (1-U(t)) / dt
c ::: ----------------------------------------------------------

      subroutine FORT_CALCEOS(coef,DIMS(coef),
     &     s, DIMS(s), rho, DIMS(rho), phi, DIMS(phi),
     &     ncomp,lo,hi,dt)

      implicit none

      integer DIMDEC(coef)
      integer DIMDEC(s)
      integer DIMDEC(rho)
      integer DIMDEC(phi)
      integer lo(SDIM),hi(SDIM)
      integer ncomp
      REAL_T  coef(DIMV(coef))
      REAL_T  s(DIMV(s),ncomp)
      REAL_T  rho(DIMV(rho),ncomp)
      REAL_T  phi(DIMV(phi))
      REAL_T  dt

      integer i, j, n
      REAL_T  s_sum

      do j=lo(2), hi(2)
         do i=lo(1), hi(1)
            s_sum = zero
            do n = 1,ncomp
               s_sum = s_sum +s(i,j,n)/rho(i,j,n)
            end do
            
            coef(i,j) = -phi(i,j)*(one-s_sum)/dt
         end do
      end do

      end

      subroutine  FORT_RECOMP_UPDATE(lo, hi,
     &                               update, DIMS(update),
     &                               phi,    DIMS(phi),
     &                               xflux,  DIMS(xflux),
     &                               yflux,  DIMS(yflux),
     &                               vol,    DIMS(vol),
     &                               nc)
      implicit none
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(update)
      integer DIMDEC(phi)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(vol)
      REAL_T update(DIMV(update),nc)
      REAL_T phi(DIMV(phi))
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T vol(DIMV(vol))

      integer i, j, n
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            do n=1,nc
               update(i,j,n)=-((xflux(i+1,j,n)-xflux(i,j,n))
     &              + (yflux(i,j+1,n)-yflux(i,j,n)))/vol(i,j)/phi(i,j)
            end do
         end do
      end do
      end

c :: ----------------------------------------------------------
c :: CHEM_DUMMY
c ::    a dummy chemical reaction.
c :: ----------------------------------------------------------
      subroutine FORT_CHEM_DUMMY(s,DIMS(s),dt,nc)
      implicit none
      integer nc
      integer DIMDEC(s)
      REAL_T  s(DIMV(s),nc), dt

      integer i,j,r
      REAL_T  deltan,rate

      rate = 1.d-8
      do i = ARG_L1(s),ARG_H1(s)
         do j =  ARG_L2(s),ARG_H2(s)
            deltan = rate*dt*s(i,j,1)*s(i,j,2)
            r = s(i,j,1)/(s(i,j,1)+s(i,j,2))
            s(i,j,1) = s(i,j,1) - deltan
            s(i,j,2) = s(i,j,2) - deltan
            s(i,j,3) = 1000 - s(i,j,1) - s(i,j,2)
         end do 
      end do

      end

c ::: -----------------------------------------------------------
c ::: Compute the capillary pressure, p_c.
c :::    It calls get_cpl in EOS_F.F. 
c ::: -----------------------------------------------------------

      subroutine FORT_MK_INV_CPL(pc,DIMS(pc),s,DIMS(s),
     &     phi,DIMS(phi),kappa,DIMS(kappa),
     &     coef,DIMS(coef),n_coef,lo,hi)

      implicit none

#include "probdata.H"

      integer DIMDEC(pc)
      integer DIMDEC(s)
      integer DIMDEC(phi)
      integer DIMDEC(kappa)
      integer DIMDEC(coef)
      integer n_coef
      integer lo(SDIM),hi(SDIM)
      REAL_T  pc(DIMV(pc))
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  phi(DIMV(phi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  coef(DIMV(coef),n_coef)

      integer i,j
      integer l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res

      do j = lo(2)-1,hi(2)+1
         do i = lo(1)-1,hi(1)+1
            l_cpl_type = nint(coef(i,j,1))
            l_cpl_coef = coef(i,j,2)
            l_sigma    = coef(i,j,3)
            l_sat_res  = coef(i,j,4)

            call get_invcpl(s(i,j,:),pc(i,j),rhoval,phi(i,j),
     &           kappa(i,j),l_sigma,ncomps,l_sat_res,l_cpl_type,l_cpl_coef)

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

      do j = lo(2)-1,hi(2)+1
         do i = lo(1)-1,hi(1)+1
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
               call get_cpl(pc(i,j),s(i,j,:),rhoval,phi(i,j),
     &              kappa(i,j),l_sigma,ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
            endif
         end do 
      end do

      end

c ::: -----------------------------------------------------------
c ::: Smooth the capillary pressure, p_c.
c ::: -----------------------------------------------------------

      subroutine FORT_SMOOTH_CPL(pc,DIMS(pc),
     &     coef,DIMS(coef),n_coef,lo,hi)

      implicit none

#include "probdata.H"

      integer DIMDEC(pc)
      integer DIMDEC(coef)
      integer n_coef
      integer lo(SDIM),hi(SDIM)
      REAL_T  pc(DIMV(pc))
      REAL_T  coef(DIMV(coef),n_coef)

      integer i,j
      integer l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res
      REAL_T, allocatable :: pctmp(:,:)
      
      allocate(pctmp(DIMV(pc)))
      pctmp = pc

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            l_cpl_type = nint(coef(i,j,1))
            l_cpl_coef = coef(i,j,2)
            l_sigma    = coef(i,j,3)
            l_sat_res  = coef(i,j,4)

            pc(i,j) = pctmp(i,j)+pctmp(i-1,j)+pctmp(i+1,j)+pctmp(i,j-1)+pctmp(i,j+1)
            pc(i,j) = pc(i,j)/5.d0
            
         end do 
      end do

      deallocate(pctmp)

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

c ::: ----------------------------------------------------------
c ::: GETDIFFUSE_CPL
c       Compute the diffusive coefficient for the capillary solve.
c       - Only solve for water ( the first component)
c ::: ----------------------------------------------------------         
      subroutine FORT_GETDIFFUSE_CPL(lbd, DIMS(lbd), 
     &     dfx, DIMS(dfx), dfy, DIMS(dfy),
     &     lo, hi, domlo, domhi, bc, ncomp)

      implicit none

#include "probdata.H"

      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM),domhi(SDIM)
      integer DIMDEC(lbd)
      integer DIMDEC(dfx)
      integer DIMDEC(dfy)
      integer ncomp
      integer bc(SDIM,2)
      REAL_T  lbd(DIMV(lbd),ncomps)
      REAL_T  dfx(DIMV(dfx),ncomp)
      REAL_T  dfy(DIMV(dfy),ncomp)

      integer i, j
      REAL_T  lambda(ncomps),H_coef,tmp(2)

c ::: diffusion coefficient in x-dir
      do i=lo(1), hi(1)+1
         do j=lo(2), hi(2)
            lambda = lbd(i-1,j,:)
            tmp(1) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
            lambda = lbd(i,j,:)
            tmp(2) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
            H_coef = half*(tmp(1)+tmp(2))
c            H_coef = one/tmp(1) + one/tmp(2)
c            H_coef = two/H_coef
            if (i .eq. domlo(1)) then
               if (bc(1,1) .eq. EXT_DIR) then
                  lambda = lbd(i,j,:)
               elseif (bc(1,1) .eq. SEEPAGE) then
                  lambda = lbd(i,j,:)
               end if
               H_coef  = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
            elseif (i .eq. domhi(1)+1) then
               if (bc(1,2) .eq. EXT_DIR) then
                  lambda = lbd(i,j,:)
               elseif (bc(1,2) .eq. SEEPAGE) then
                  lambda = lbd(i,j,:)
               end if       
               H_coef  = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
            endif

            dfx(i,j,1) = H_coef*rhoval(1)
         end do
      end do


c ::: diffusion coefficient in y-dir
      do j=lo(2), hi(2)+1
         do i=lo(1), hi(1)
            lambda = lbd(i,j-1,:)
            tmp(1) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
            lambda = lbd(i,j,:)
            tmp(2) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
            H_coef = half*(tmp(1)+tmp(2))
c            H_coef = one/tmp(1) + one/tmp(2)
c            H_coef = two/H_coef
            if (j .eq. domlo(2)) then
               if (bc(2,1) .eq. EXT_DIR) then
                  lambda = lbd(i,j,:)
               elseif (bc(2,1) .eq. SEEPAGE) then
                  lambda = lbd(i,j,:)
               end if
               H_coef = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
            elseif (j .eq. domhi(2)+1) then
               if (bc(2,2) .eq. EXT_DIR) then
                  lambda = lbd(i,j,:)
               elseif (bc(2,2) .eq. SEEPAGE) then
                  lambda = lbd(i,j,:)
               end if       
               H_coef = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
            end if

            dfy(i,j,1) = H_coef*rhoval(1)

         end do
      end do

      end

c :: ----------------------------------------------------------
c :: GETDIFFUSE_CPL_dp
c      Compute the diffusive coefficient for the jacobian in 
c      capillary solve.  
c      - Only valid for water (the first component)
c :: ----------------------------------------------------------         
      subroutine FORT_GETDIFFUSE_CPL_dp(s, h_cc, DIMS(s), 
     &     lbd, DIMS(lbd), dfx, DIMS(dfx), dfy, DIMS(dfy), 
     &     phi, DIMS(phi), kappa, DIMS(kappa),
     &     coef,DIMS(coef),n_coef,
     &     lo, hi, domlo, domhi, bc, ncomp)

      implicit none

#include "probdata.H"

      integer DIMDEC(s)
      integer DIMDEC(lbd)
      integer DIMDEC(dfx)
      integer DIMDEC(dfy)
      integer DIMDEC(phi)
      integer DIMDEC(kappa)
      integer DIMDEC(coef)
      integer n_coef
      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM),domhi(SDIM)
      integer bc(SDIM,2)
      integer ncomp
      REAL_T  s(DIMV(s),ncomps)
      REAL_T  h_cc(DIMV(s))
      REAL_T  lbd(DIMV(lbd),ncomps)
      REAL_T  dfx(DIMV(dfx),ncomp)
      REAL_T  dfy(DIMV(dfy),ncomp)
      REAL_T  phi(DIMV(phi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  coef(DIMV(coef),n_coef)

      integer i, j
      REAL_T  lambda(ncomps),H_coef,dpcds
      integer l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res

      do i = ARG_L1(s),ARG_H1(s)
         do j = ARG_L2(s),ARG_H2(s)
            lambda = lbd(i,j,1:2)
            h_cc(i,j) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
            if (h_cc(i,j) .gt. 1.d-12) then
               l_cpl_type = nint(coef(i,j,1))
               l_cpl_coef = coef(i,j,2)
               l_sigma    = coef(i,j,3)
               l_sat_res  = coef(i,j,4)
               call get_dcplds(dpcds,s(i,j,:),rhoval,phi(i,j),
     &              kappa(i,j),l_sigma,ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
               h_cc(i,j) = h_cc(i,j)*dpcds*rhoval(1)
            else
               h_cc(i,j) = zero
            end if
         end do
      end do

c ::: diffusion coefficient in x-dir
      do i=lo(1), hi(1)+1
         do j=lo(2), hi(2)

            H_coef = half*(h_cc(i-1,j)+h_cc(i,j))
c            H_coef = one/h_cc(i-1,j) + one/h_cc(i,j)
c            H_coef = two/H_coef
            if (i .eq. domlo(1)) then
               if (bc(1,1) .eq. EXT_DIR) then
                  H_coef = h_cc(i,j)
               elseif (bc(1,1) .eq. SEEPAGE) then
                  H_coef = h_cc(i,j)
               end if
            elseif (i .eq. domhi(1)+1) then
               if (bc(1,2) .eq. EXT_DIR) then
                  H_coef = h_cc(i,j)
               elseif (bc(1,2) .eq. SEEPAGE) then
                  H_coef = h_cc(i,j)
               end if  
            end if
            dfx(i,j,1) = H_coef
         end do
      end do


c ::: diffusion coefficient in y-dir
      do j=lo(2), hi(2)+1
         do i=lo(1), hi(1)

            H_coef = half*(h_cc(i,j-1)+h_cc(i,j))
c            H_coef = one/h_cc(i,j-1) + one/h_cc(i,j)
c            H_coef = two/H_coef
            if (j .eq. domlo(2)) then
               if (bc(2,1) .eq. EXT_DIR) then
                  H_coef = h_cc(i,j)
               elseif (bc(2,1) .eq. SEEPAGE) then
                  H_coef = h_cc(i,j)
               end if
               
            elseif (j .eq. domhi(2)+1) then
               if (bc(2,2) .eq. EXT_DIR) then
                  H_coef = h_cc(i,j)
               elseif (bc(2,2) .eq. SEEPAGE) then
                  H_coef = h_cc(i,j)
               end if 
            end if

            dfy(i,j,1) = H_coef

         end do
      end do

      end

c ::: ---------------------------------------------------------------
c ::: RICHARD_COEF
c       Compute the diffusive coefficient for the richard equation.
c       - Only solve for water ( the first component)
c ::: ---------------------------------------------------------------    
      subroutine FORT_RICHARD_COEF(lbd, DIMS(lbd), 
     &     dfx, DIMS(dfx), dfy, DIMS(dfy),
     &     umac, DIMS(umac), vmac, DIMS(vmac),
     &     lo, hi, domlo, domhi, bc, 
     &     inflow_bc_lo, inflow_bc_hi, nc, do_upwind)

      implicit none

#include "probdata.H"

      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM),domhi(SDIM)
      integer DIMDEC(lbd)
      integer DIMDEC(dfx)
      integer DIMDEC(dfy)
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer nc
      integer bc(SDIM,2)
      integer inflow_bc_lo(SDIM), inflow_bc_hi(SDIM)
      REAL_T  lbd(DIMV(lbd),ncomps)
      REAL_T  dfx(DIMV(dfx))
      REAL_T  dfy(DIMV(dfy))
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      logical do_upwind

      integer i, j
      REAL_T  H_coef

c ::: diffusion coefficient in x-dir
      do i=lo(1), hi(1)+1
         do j=lo(2), hi(2)
            if (do_upwind) then
               H_coef = cvmgp(lbd(i-1,j,nc),lbd(i,j,nc),umac(i,j))
            else
               H_coef = half*(lbd(i-1,j,nc)+lbd(i,j,nc))
            end if
            if (i .eq. domlo(1)) then
               if (bc(1,1) .eq. EXT_DIR) then
                  H_coef = lbd(i-1,j,nc)
               elseif (bc(1,1) .eq. FOEXTRAP .and. inflow_bc_lo(1) .eq. 1) then
                  H_coef = zero
               end if
            elseif (i .eq. domhi(1)+1) then
               if (bc(1,2) .eq. EXT_DIR) then
                  H_coef = lbd(i,j,nc)
               elseif (bc(1,2) .eq. FOEXTRAP .and. inflow_bc_hi(1) .eq. 1) then
                  H_coef = zero
               end if    
            endif
            dfx(i,j) = H_coef 
         end do
      end do


c ::: diffusion coefficient in y-dir
      do j=lo(2), hi(2)+1
         do i=lo(1), hi(1)
            if (do_upwind) then
               H_coef = cvmgp(lbd(i,j-1,nc),lbd(i,j,nc),vmac(i,j))
            else
               H_coef = half*(lbd(i,j-1,nc)+lbd(i,j,nc))
            end if
            if (j .eq. domlo(2)) then
               if (bc(2,1) .eq. EXT_DIR) then
                  H_coef = lbd(i,j-1,nc) 
               elseif (bc(2,1) .eq. FOEXTRAP .and. inflow_bc_lo(2) .eq. 1) then
                  H_coef = zero
               end if    
            elseif (j .eq. domhi(2)+1) then
               if (bc(2,2) .eq. EXT_DIR) then
                  H_coef = lbd(i,j,nc) 
               elseif (bc(2,2) .eq. FOEXTRAP .and. inflow_bc_hi(2) .eq. 1) then
                  H_coef = zero
               end if    
            endif
            dfy(i,j) = H_coef
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
      logical do_upwind

      integer i, j,l_kr_type,l_cpl_type
      REAL_T l_kr_coef, l_kr_sres, l_cpl_coef,l_sigma,l_cpl_sres
      REAL_T  grav
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
            call get_dld(dltmp,sat,rtmp,muval, 
     &           l_kr_type,l_kr_coef,l_kr_sres)
         
            dpcdn(i,j) = dpcdn(i,j)/(rhoval(1)*(one-l_cpl_sres))
            dlddn(i,j,1) = dltmp(1)/(rhoval(1)*(one-l_kr_sres))

         end do
      end do

      grav = -rhoval(1)*gravity
      bx = zero
      by = zero
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            if (do_upwind) then
               ldx(1)  = cvmgp(lbd(i-1,j,1),lbd(i,j,1),umac(i,j))
               ldx(2)  = cvmgp(lbd(i,j,1),lbd(i+1,j,1),umac(i+1,j))
               ldy(1)  = cvmgp(lbd(i,j-1,1),lbd(i,j,1),vmac(i,j))
               ldy(2)  = cvmgp(lbd(i,j,1),lbd(i,j+1,1),vmac(i,j+1))
               dldx(1) = cvmgp(dlddn(i-1,j,1),zero,umac(i,j))
               dldx(2) = cvmgp(zero,dlddn(i+1,j,1),umac(i+1,j)) 
               dldx(3) = cvmgp(zero,dlddn(i,j,1),umac(i,j))
               dldx(4) = cvmgp(dlddn(i,j,1),zero,umac(i+1,j)) 
               dldy(1) = cvmgp(dlddn(i,j-1,1),zero,vmac(i,j))
               dldy(2) = cvmgp(zero,dlddn(i,j+1,1),vmac(i,j+1)) 
               dldy(3) = cvmgp(zero,dlddn(i,j,1),vmac(i,j))
               dldy(4) = cvmgp(dlddn(i,j,1),zero,vmac(i,j+1))

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

            bx(i,j,2)   = (-dldx(1)*(-multx(1)*(pc(i,j)-pc(i-1,j))/dx(1)) -
     &           ldx(1)*(multx(1)*dpcdn(i-1,j)/dx(1)))/dx(1)
            bx(i+1,j,1) = (dldx(2)*(-multx(2)*(pc(i+1,j)-pc(i,j))/dx(1)) -
     &           ldx(2)*(multx(2)*dpcdn(i+1,j)/dx(1)))/dx(1)
            bx(i,j,3)   = (dldx(4)*(-multx(2)*(pc(i+1,j)-pc(i,j))/dx(1)) - 
     &           dldx(3)*(-multx(1)*(pc(i,j)-pc(i-1,j))/dx(1)) +
     &           (ldx(1)+ldx(2))*dpcdn(i,j)/dx(1))/dx(1)
            
            by(i,j,2)   = (-dldy(1)*(-multy(1)*(pc(i,j)-pc(i,j-1))/dx(2) - grav) -
     &           ldy(1)*(multy(1)*dpcdn(i,j-1)/dx(2)))/dx(2)
            by(i,j+1,1) = (dldy(2)*(-multy(2)*(pc(i,j+1)-pc(i,j))/dx(2) - grav) -
     &           ldy(2)*(multy(2)*dpcdn(i,j+1)/dx(2)))/dx(2)
            by(i,j,3)   = (dldy(4)*(-multy(2)*(pc(i,j+1)-pc(i,j))/dx(2) - grav) - 
     &           dldy(3)*(-multy(1)*(pc(i,j)-pc(i,j-1))/dx(2) - grav) +
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
      logical do_upwind

      integer i, j,n,l_kr_type,l_cpl_type
      REAL_T l_kr_coef, l_kr_sres, l_cpl_coef,l_sigma,l_cpl_sres
      REAL_T  grav, rgx, rgy
      REAL_T  multx(2), multy(2), ldx(2), ldy(2)
      REAL_T  comp(ncomps)
      REAL_T  ldtmp(ncomps),ldn(3),pcn(3),rgtmpx(3),rgtmpy(3), dcx(3), dcy(3)
c
c     gravity is nonzero only in the y direction
c
      grav = -rhoval(1)*gravity
      bx = zero
      by = zero
      do j=lo(2), hi(2)
         do i=lo(1), hi(1)

            if (do_upwind) then
               ldx(1) = cvmgp(kpx(i,j)*lbd(i-1,j,1),kpx(i,j)*lbd(i,j,1),umac(i,j))
               ldx(2) = cvmgp(kpx(i+1,j)*lbd(i,j,1),kpx(i+1,j)*lbd(i+1,j,1),umac(i+1,j))
               ldy(1) = cvmgp(kpy(i,j)*lbd(i,j-1,1),kpy(i,j)*lbd(i,j,1),vmac(i,j))
               ldy(2) = cvmgp(kpy(i,j+1)*lbd(i,j,1),kpy(i,j+1)*lbd(i,j+1,1),vmac(i,j+1))
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

            rgx = (ldx(2)*(-multx(2)*(pc(i+1,j)-pc(i,j))/dx(1)) -
     &           ldx(1)*(-multx(1)*(pc(i,j)-pc(i-1,j))/dx(1)))/dx(1) 
            rgy =  (ldy(2)*(-multy(2)*(pc(i,j+1)-pc(i,j))/dx(2) - grav) - 
     &           ldy(1)*(-multy(1)*(pc(i,j)-pc(i,j-1))/dx(2) - grav))/dx(2)           

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
               if (do_upwind) then
                  ldx(1) = cvmgz(ldx(1),cvmgp(kpx(i,j)*ldn(1),kpx(i,j)*ldn(2),umac(i,j)),ldx(1))
                  ldx(2) = cvmgz(ldx(2),cvmgp(kpx(i+1,j)*ldn(2),kpx(i+1,j)*ldn(3),umac(i+1,j)),ldx(2))
               else
                  ldx(1) = cvmgz(ldx(1),kpx(i,j)*half*(ldn(1)+ldn(2)),ldx(1))
                  ldx(2) = cvmgz(ldx(2),kpx(i+1,j)*half*(ldn(2)+ldn(3)),ldx(2))
               end if
               if (i .eq. domlo(1) .and. bc(1,1) .eq. EXT_DIR) then
                  ldx(1) = cvmgz(ldx(1),kpx(i,j)*ldn(1),ldx(1))
               elseif (i .eq. domhi(1) .and. bc(1,2) .eq. EXT_DIR) then
                  ldx(2) = cvmgz(ldx(2),kpx(i+1,j)*ldn(3),ldx(2))
               end if

               rgtmpx(n+2) = (ldx(2)*(-multx(2)*(pcn(3)-pcn(2))/dx(1)) -
     &           ldx(1)*(-multx(1)*(pcn(2)-pcn(1))/dx(1)))/dx(1) 
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
               if (do_upwind) then
                  ldy(1) = cvmgz(ldy(1),cvmgp(kpy(i,j)*ldn(1),kpy(i,j)*ldn(2),vmac(i,j)),ldy(1))
                  ldy(2) = cvmgz(ldy(2),cvmgp(kpy(i,j+1)*ldn(2),kpy(i,j+1)*ldn(3),vmac(i,j+1)),ldy(2))
               else
                  ldy(1) = cvmgz(ldy(1),kpy(i,j)*half*(ldn(1)+ldn(2)),ldy(1))
                  ldy(2) = cvmgz(ldy(2),kpy(i,j+1)*half*(ldn(2)+ldn(3)),ldy(2))
               end if
               if (j .eq. domlo(2) .and. bc(2,1) .eq. EXT_DIR) then
                  ldy(1) = cvmgz(ldy(1),kpy(i,j)*ldn(1),ldy(1))
               elseif (j .eq. domhi(2) .and. bc(2,2) .eq. EXT_DIR) then
                  ldy(2) = cvmgz(ldy(2),kpy(i,j+1)*ldn(3),ldy(2))
               end if

               rgtmpy(n+2) =  (ldy(2)*(-multy(2)*(pcn(3)-pcn(2))/dx(2) - grav) - 
     &           ldy(1)*(-multy(1)*(pcn(2)-pcn(1))/dx(2) - grav))/dx(2)                
            end do

            bx(i,j,2)   = (rgtmpx(1) - rgx)/dcx(1)
            bx(i+1,j,1) = (rgtmpx(3) - rgx)/dcx(3)
            bx(i,j,3)   = (rgtmpx(2) - rgx)/dcx(2)
            by(i,j,2)   = (rgtmpy(1) - rgy)/dcy(1)
            by(i,j+1,1) = (rgtmpy(3) - rgy)/dcy(3)
            by(i,j,3)   = (rgtmpy(2) - rgy)/dcy(2)

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
     &     lbd,DIMS(lbd),pc,DIMS(pc),
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
      logical do_upwind

      integer i, j,n,l_kr_type,l_cpl_type
      REAL_T l_kr_coef, l_kr_sres, l_cpl_coef,l_sigma,l_cpl_sres
      REAL_T  grav, rgx, rgy
      REAL_T  multx(2), multy(2), ldx(2), ldy(2)
      REAL_T  comp(ncomps),cap
      REAL_T  ldtmp(ncomps),ldn(3),pcn(3),rgtmpx(3),rgtmpy(3), dcx(3), dcy(3)
c
c     gravity is nonzero only in the y direction
c
      grav = -rhoval(1)*gravity
      bx = zero
      by = zero
      do j=lo(2), hi(2)
         do i=lo(1), hi(1)

            if (do_upwind) then
               ldx(1) = cvmgp(kpx(i,j)*lbd(i-1,j,1),kpx(i,j)*lbd(i,j,1),umac(i,j))
               ldx(2) = cvmgp(kpx(i+1,j)*lbd(i,j,1),kpx(i+1,j)*lbd(i+1,j,1),umac(i+1,j))
               ldy(1) = cvmgp(kpy(i,j)*lbd(i,j-1,1),kpy(i,j)*lbd(i,j,1),vmac(i,j))
               ldy(2) = cvmgp(kpy(i,j+1)*lbd(i,j,1),kpy(i,j+1)*lbd(i,j+1,1),vmac(i,j+1))
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

            rgx = (ldx(2)*(-multx(2)*(pc(i+1,j)-pc(i,j))/dx(1)) -
     &           ldx(1)*(-multx(1)*(pc(i,j)-pc(i-1,j))/dx(1)))/dx(1) 
            rgy = (ldy(2)*(-multy(2)*(pc(i,j+1)-pc(i,j))/dx(2) - grav) - 
     &           ldy(1)*(-multy(1)*(pc(i,j)-pc(i,j-1))/dx(2) - grav))/dx(2)           

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
               pcn = pc(i-1:i+1,j)
               cap = pc(i+n,j) + deps/l_sigma
               dcx(n+2) = deps/l_sigma
               call get_invcpl(comp,cap,rhoval,phi(i+n,j),kappa(i+n,j),l_sigma,
     &              ncomps,l_cpl_sres,l_cpl_type,l_cpl_coef)
               call get_ld(ldtmp,comp,rhoval,muval,ncomps,l_kr_type,l_kr_coef,l_kr_sres)
               ldn(n+2) = ldtmp(1)
               pcn(n+2) = cap
               if (do_upwind) then
                  ldx(1) = cvmgz(ldx(1),cvmgp(ldn(1),ldn(2),umac(i,j)),ldx(1))
                  ldx(2) = cvmgz(ldx(2),cvmgp(ldn(2),ldn(3),umac(i+1,j)),ldx(2))
               else
                  ldx(1) = cvmgz(ldx(1),half*(ldn(1)+ldn(2)),ldx(1))
                  ldx(2) = cvmgz(ldx(2),half*(ldn(2)+ldn(3)),ldx(2))
               end if
               if (i .eq. domlo(1) .and. bc(1,1) .eq. EXT_DIR) then
                  ldx(1) = cvmgz(ldx(1),ldn(1),ldx(1))
               elseif (i .eq. domhi(1) .and. bc(1,2) .eq. EXT_DIR) then
                  ldx(2) = cvmgz(ldx(2),ldn(3),ldx(2))
               end if

               rgtmpx(n+2) = (ldx(2)*(-multx(2)*(pcn(3)-pcn(2))/dx(1)) -
     &           ldx(1)*(-multx(1)*(pcn(2)-pcn(1))/dx(1)))/dx(1) 
            end do

c           y-direction
            do n = -1,1
               l_kr_type  = nint(krcoef(i,j+n,1))
               l_kr_coef  = krcoef(i,j+n,2)
               l_kr_sres  = krcoef(i,j+n,3)
               l_cpl_type  = nint(cplcoef(i,j+n,1))
               l_cpl_coef  = cplcoef(i,j+n,2)
               l_sigma     = cplcoef(i,j+n,3)
               l_cpl_sres  = cplcoef(i,j+n,4)
               ldn = lbd(i,j-1:j+1,1)
               pcn = pc(i,j-1:j+1)
               cap = pc(i,j+n) + deps/l_sigma
               dcy(n+2) = deps/l_sigma
               call get_invcpl(comp,cap,rhoval,phi(i,j+n),kappa(i,j+n),l_sigma,
     &              ncomps,l_cpl_sres,l_cpl_type,l_cpl_coef)
               call get_ld(ldtmp,comp,rhoval,muval,ncomps,l_kr_type,l_kr_coef,l_kr_sres)
               ldn(n+2) = ldtmp(1)
               pcn(n+2) = cap
               if (do_upwind) then
                  ldy(1) = cvmgz(ldy(1),cvmgp(ldn(1),ldn(2),vmac(i,j)),ldy(1))
                  ldy(2) = cvmgz(ldy(2),cvmgp(ldn(2),ldn(3),vmac(i,j+1)),ldy(2))
               else
                  ldy(1) = cvmgz(ldy(1),half*(ldn(1)+ldn(2)),ldy(1))
                  ldy(2) = cvmgz(ldy(2),half*(ldn(2)+ldn(3)),ldy(2))
               end if
               if (j .eq. domlo(2) .and. bc(2,1) .eq. EXT_DIR) then
                  ldy(1) = cvmgz(ldy(1),ldn(1),ldy(1))
               elseif (j .eq. domhi(2) .and. bc(2,2) .eq. EXT_DIR) then
                  ldy(2) = cvmgz(ldy(2),ldn(3),ldy(2))
               end if
               rgtmpy(n+2) = (ldy(2)*(-multy(2)*(pcn(3)-pcn(2))/dx(2) - grav) - 
     &           ldy(1)*(-multy(1)*(pcn(2)-pcn(1))/dx(2) - grav))/dx(2)                

            end do

            bx(i,j,2)   = (rgtmpx(1) - rgx)/dcx(1)
            bx(i+1,j,1) = (rgtmpx(3) - rgx)/dcx(3)
            bx(i,j,3)   = (rgtmpx(2) - rgx)/dcx(2)
            by(i,j,2)   = (rgtmpy(1) - rgy)/dcy(1)
            by(i,j+1,1) = (rgtmpy(3) - rgy)/dcy(3)
            by(i,j,3)   = (rgtmpy(2) - rgy)/dcy(2)
         end do
      end do

      end

c ::: -----------------------------------------------------------
c ::: Compute the capillary pressure, p_c.
c :::    It calls get_cpl in EOS_F.F. 
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
         
            alpha(i,j) = phi(i,j)/dpcdn*rhoval(1)*(one-l_sat_res)
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
     &     lo,hi,domlo,domhi,dx,
     &     inflow_bc_lo,inflow_bc_hi,
     &     inflow_vel_lo,inflow_vel_hi, mult)

      implicit none
      
      integer DIMDEC(rg)
      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      integer inflow_bc_lo(SDIM),inflow_bc_hi(SDIM)
      REAL_T  inflow_vel_lo(SDIM),inflow_vel_hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  rg(DIMV(rg))
      REAL_T  mult
      
      integer i, j

      if (lo(1) .eq. domlo(1) .and. inflow_bc_lo(1) .eq. 1) then
           do j = lo(2),hi(2)
              rg(lo(1),j) = rg(lo(1),j) - mult*inflow_vel_lo(1)/dx(1)
           end do
      end if
      if (hi(1) .eq. domhi(1) .and. inflow_bc_hi(1) .eq. 1) then
           do j = lo(2),hi(2)
              rg(hi(1),j) = rg(hi(1),j) + mult*inflow_vel_hi(1)/dx(1)
           end do
      end if
      if (lo(2) .eq. domlo(2) .and. inflow_bc_lo(2) .eq. 1) then
           do i = lo(1),hi(1)
              rg(i,lo(2)) = rg(i,lo(2)) - mult*inflow_vel_lo(2)/dx(2)
           end do
      end if
      if (hi(2) .eq. domhi(2) .and. inflow_bc_hi(2) .eq. 1) then
           do i = lo(1),hi(1)
              rg(i,hi(2)) = rg(i,hi(2)) + mult*inflow_vel_hi(2)/dx(2)
           end do
      end if

      end     

c ::: ----------------------------------------------------------
c ::: compute divergence of velocity
c ::: ----------------------------------------------------------
      subroutine FORT_DIV_UMAC(s,DIMS(s),
     &     ux,DIMS(ux),uy,DIMS(uy),lo,hi,dx)
      
      integer lo(2), hi(2)
      integer DIMDEC(ux)
      integer DIMDEC(uy)
      integer DIMDEC(s)
      REAL_T  ux(DIMV(ux))
      REAL_T  uy(DIMV(uy))
      REAL_T  s(DIMV(s))
      REAL_T  dx(2)

      integer i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            s(i,j) = s(i,j) + (ux(i+1,j)-ux(i,j))/dx(1) +
     &              (uy(i,j+1)-uy(i,j))/dx(2)
         end do
      end do

      end


c ::: -----------------------------------------------------------
c ::: This routine compute the phase velocity of component 1
c ::: -----------------------------------------------------------
      subroutine FORT_UPHASE(uphase,DIMS(uphase),vphase,DIMS(vphase),
     &     umac,DIMS(umac),vmac,DIMS(vmac),
     &     pc,DIMS(pc),lbd,DIMS(lbd),
     &     kpx,DIMS(kpx),kpy,DIMS(kpy),
     &     sbc,dom_lo,dom_hi,dx,lo,hi)

      implicit none

#include "probdata.H"

      integer DIMDEC(uphase), DIMDEC(vphase)
      integer DIMDEC(umac), DIMDEC(vmac)
      integer DIMDEC(kpx), DIMDEC(kpy)
      integer DIMDEC(pc)
      integer DIMDEC(lbd)
      integer lo(SDIM),hi(SDIM)
      integer dom_lo(SDIM),dom_hi(SDIM)
      integer sbc(SDIM,2)
      REAL_T  dx(SDIM)
      REAL_T  uphase(DIMV(uphase))
      REAL_T  vphase(DIMV(vphase))
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  kpx(DIMV(kpx)),kpy(DIMV(kpy))  
      REAL_T  pc(DIMV(pc))
      REAL_T  lbd(DIMV(lbd),ncomps)

      integer i,j
      REAL_T  dpc
      REAL_T  ld(ncomps),le(ncomps)

c     phase velocity in the x-direction
      do i = lo(1),hi(1)+1
         do j = lo(2),hi(2)
            if (model .eq. 2) then
               dpc = (pc(i,j) - pc(i-1,j))/dx(1)
               ld = lbd(i,j,:)
               le = lbd(i-1,j,:)
               ld(1) = half*(ld(1) + le(1))
               ld(2) = half*(ld(2) + le(2))

               if (i .eq. dom_lo(1)) then
                  if (sbc(1,1) .eq. EXT_DIR .or. sbc(1,1) .eq. SEEPAGE) then
                     dpc = dpc/2d0
                     ld = lbd(i-1,j,:)
                  endif
               else if (i .eq. dom_hi(1)+1) then
                  if (sbc(1,2) .eq. EXT_DIR .or. sbc(1,2) .eq. SEEPAGE) then
                     dpc = dpc/2d0
                     ld = lbd(i,j,:)
                  endif
               end if
               uphase(i,j) = ld(1)/(ld(1)+ld(2))*(umac(i,j) + kpx(i,j)*ld(2)*dpc)
            end if
         end do
      end do

      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)
               
            if (model .eq. 2) then  
               dpc = (pc(i,j) - pc(i,j-1))/dx(2)
               ld = lbd(i,j,:)
               le = lbd(i,j-1,:)
               ld(1) = half*(ld(1)+le(1))
               ld(2) = half*(ld(2)+le(2))

               if (j .eq. dom_lo(2)) then
                  if (sbc(2,1) .eq. EXT_DIR .or. sbc(2,1) .eq. SEEPAGE) then 
                     dpc = dpc/2.d0
                     ld = lbd(i,j-1,:)
                  endif
               elseif (j .eq. dom_hi(2)+1) then
                  if (sbc(2,2) .eq. EXT_DIR .or. sbc(2,2) .eq. SEEPAGE) then
                     dpc = dpc/2.d0
                     ld = lbd(i,j,:)
                  endif
               end if
               vphase(i,j) = ld(1)/(ld(1)+ld(2))*(vmac(i,j) 
     &              + kpy(i,j)*ld(2)*dpc)
            end if
         end do
      end do

      end

c :: ------------------------------------------------------------
c :: UPHASE_P
c      Compute the phase velocity based on the phase pressure  
c
c :: ------------------------------------------------------------    
      subroutine FORT_UPHASE_P(ux,DIMS(ux),uy,DIMS(uy),
     &     lbd,DIMS(lbd),p,DIMS(p),
     &     lo,hi,domlo,domhi,dx,bc, 
     &     inflow_bc_lo,inflow_bc_hi,
     &     inflow_vel_lo,inflow_vel_hi)

      implicit none

#include "probdata.H"

      integer DIMDEC(ux)
      integer DIMDEC(uy)
      integer DIMDEC(p)
      integer DIMDEC(lbd)
      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      integer bc(SDIM,2)
      integer inflow_bc_lo(SDIM),inflow_bc_hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  ux(DIMV(ux))
      REAL_T  uy(DIMV(uy))
      REAL_T  p(DIMV(p))
      REAL_T  lbd(DIMV(lbd),ncomps)
      REAL_T  inflow_vel_lo(SDIM),inflow_vel_hi(SDIM)
      logical do_upwind
      integer i, j
      REAL_T  grav, mult, ld, dp
c
c     gravity is nonzero only in the y direction
c
      grav = -rhoval(1)*gravity
      ux = zero
      uy = zero
      do_upwind = .true.
      do j=lo(2), hi(2)
         do i=lo(1), hi(1)+1
            dp = - (p(i,j)-p(i-1,j))/dx(1)
            mult = one
            if (do_upwind) then
               ld   =  cvmgp(lbd(i-1,j,1),lbd(i,j,1),dp)
            else
               ld   = half*(lbd(i,j,1)+lbd(i-1,j,1))
            end if
            if (i .eq. domlo(1)) then
               if (bc(1,1) .eq. EXT_DIR) then
                  mult = two 
                  ld   = lbd(i-1,j,1)
               elseif (bc(1,1) .eq. FOEXTRAP) then
                  mult = zero
                  if (inflow_bc_lo(1) .eq. 1) then
                     ld = zero
                     ux(i,j) = inflow_vel_lo(1)
                  end if
               end if
            end if
            if (i .eq. domhi(1)+1) then
               if (bc(1,2) .eq. EXT_DIR) then
                  mult = two
                  ld   = lbd(i,j,1)
               elseif (bc(1,2) .eq. FOEXTRAP) then
                  mult = zero
                  if (inflow_bc_hi(1) .eq. 1) then
                     ld = zero
                     ux(i,j) = inflow_vel_hi(1)
                  end if
               end if
            end if

            ux(i,j) = ux(i,j) + ld*mult*dp

         end do
      end do
            
      do j=lo(2), hi(2)+1
         do i=lo(1), hi(1)
            dp = - ((p(i,j)-p(i,j-1))/dx(2) - grav)
            mult = one
            if (do_upwind) then
               ld   =  cvmgp(lbd(i,j-1,1),lbd(i,j,1),dp)
            else
               ld   = half*(lbd(i,j,1)+lbd(i,j-1,1))
            end if
            if (j .eq. domlo(2)) then
               if (bc(2,1) .eq. EXT_DIR) then
                  mult = two
                  ld = lbd(i,j-1,1)
               elseif (bc(2,1) .eq. FOEXTRAP) then
                  mult = zero
                  if (inflow_bc_lo(2) .eq. 1) then
                     ld = zero
                     uy(i,j) = inflow_vel_lo(2)
                  end if
               end if
            end if
            if (j .eq. domhi(2)+1) then
               if (bc(2,2) .eq. EXT_DIR) then
                  mult = two 
                  ld = lbd(i,j,1)
               elseif (bc(2,2) .eq. FOEXTRAP) then
                  mult = zero
                  if (inflow_bc_hi(2) .eq. 1) then
                     ld = zero
                     uy(i,j) = inflow_vel_hi(2)
                  end if
               end if
            end if

            uy(i,j) = uy(i,j) - ld*(mult*(p(i,j)-p(i,j-1))/dx(2) - grav)  

         end do
      end do

      end   



c ::: ----------------------------------------------------------
c ::: Average the edge-centered velocities to a cell-centered quantity.
c ::: ----------------------------------------------------------
      subroutine FORT_AVG_UMAC(umac,DIMS(umac),vmac,DIMS(vmac),
     &                         vel,DIMS(vel),lo,hi)

      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(vel)
      integer lo(SDIM),hi(SDIM)
      REAL_T umac(DIMV(umac))
      REAL_T vmac(DIMV(vmac))
      REAL_T  vel(DIMV(vel),SDIM)

      integer i,j

      do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          vel(i,j,1) = half * (umac(i,j) + umac(i+1,j))
          vel(i,j,2) = half * (vmac(i,j) + vmac(i,j+1))
        end do
      end do

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
c ::: Define hydrostatic boundary condition
c ::: ---------------------------------------------------------
      
      subroutine FORT_HYDRO(rho,DIMS(rho),rhoval,ncomps,
     &     cplcoef,DIMS(cplcoef),ncplcoef,dx,wt_loc,grav)

      implicit none

      integer DIMDEC(rho)
      integer DIMDEC(cplcoef)
      integer ncomps, ncplcoef
      REAL_T  rho(DIMV(rho),ncomps)
      REAL_T  cplcoef(DIMV(cplcoef),ncplcoef)
      REAL_T  rhoval(ncomps)
      REAL_T  dx(SDIM)
      REAL_T  wt_loc
      REAL_T  grav
      
      integer i,j,k,l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res
      REAL_T  y
      REAL_T  pc,phi,kappa,satlim(2),stmp(2)

c     Gauss Legendre integration
      REAL_T, parameter :: qd(5) = (/ 0.0d0,
     &                                0.906179845938664d0*half,
     &                                0.538469310105683d0*half,
     &                               -0.906179845938664d0*half,
     &                               -0.538469310105683d0*half /)


      REAL_T, parameter :: wt(5) = (/ 0.568888888888889d0*half,
     &                                0.236926885056189d0*half,
     &                                0.478628670499366d0*half,
     &                                0.236926885056189d0*half,
     &                                0.478628670499366d0*half /)

c     dummy
      phi   = zero
      kappa = zero

      do j = ARG_L2(rho), ARG_H2(rho)
         do i = ARG_L1(rho), ARG_H1(rho)
            rho(i,j,:) = zero
            do k = 1,5
               y = (dble(j)+half+qd(k))*dx(2)
               if (y .le. wt_loc) then
                  rho(i,j,1) = rho(i,j,1) + rhoval(1)*wt(k)
                  rho(i,j,2) = rho(i,j,2) + rhoval(2)*wt(k)
               else
                  l_cpl_type = nint(cplcoef(i,j,1))
                  l_cpl_coef = cplcoef(i,j,2)
                  l_sigma    = cplcoef(i,j,3)
                  l_sat_res  = cplcoef(i,j,4)
                  satlim(1) = 1.1*l_sat_res*rhoval(1)
                  satlim(2) = (one - 1.1*l_sat_res)*rhoval(2)
                  pc = (rhoval(1)-rhoval(2))*grav*(y-wt_loc)
                  call get_invcpl(stmp,pc,rhoval,phi,kappa,l_sigma,
     &              ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
                  if (stmp(1) .lt. satlim(1)) then
                     stmp = satlim
                  end if
                  rho(i,j,1) = rho(i,j,1) + stmp(1)*wt(k)
                  rho(i,j,2) = rho(i,j,2) + stmp(2)*wt(k)
               end if
            end do 
         end do
      end do 

      end

c ::: ------------------------------------------------------------------
c ::: Define approximate steady state condition based on inflow velocity
c ::: ------------------------------------------------------------------
      
      subroutine FORT_STEADYSTATE(rho,DIMS(rho),rhoval,muval,ncomps,
     &     kappa,DIMS(kappa),krcoef,DIMS(krcoef),nkrcoef,
     &     dx,inflow_vel,nc,grav)

      implicit none

      integer nc
      integer DIMDEC(rho)
      integer DIMDEC(kappa)
      integer DIMDEC(krcoef)
      integer ncomps, nkrcoef
      REAL_T  rho(DIMV(rho),ncomps)
      REAL_T  kappa(DIMV(kappa))
      REAL_T  krcoef(DIMV(krcoef),nkrcoef)
      REAL_T  rhoval(ncomps)
      REAL_T  muval(ncomps)
      REAL_T  dx(SDIM)
      REAL_T  inflow_vel
      REAL_T  grav
      
      integer i,j,l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T  vtot, gstar, sol
      REAL_T  rtmp(2)

      vtot  = zero
      do j = ARG_L2(rho), ARG_H2(rho)
         do i = ARG_L1(rho), ARG_H1(rho)
            l_kr_type = nint(krcoef(i,j,1))
            l_kr_coef = krcoef(i,j,2)
            l_sat_res = krcoef(i,j,3)

            rtmp = 1.2
            rtmp(1:ncomps) = rhoval(1:ncomps)
            gstar = -kappa(i,j)*(rtmp(1)-rtmp(2))*grav

            call FORT_FIND_INV_FLUX(sol,inflow_vel,nc,vtot,gstar,muval,
     &           ncomps,l_kr_type,l_kr_coef)

            sol = sol*(1.0-l_sat_res)+l_sat_res
            
            rho(i,j,1) = rhoval(1)*sol

            if (ncomps > 1) then
               rho(i,j,2) = rhoval(2)*(1.0-sol)
            endif

         end do
      end do 

      end


c ::: ---------------------------------------------------------
c ::: Define Rockhold boundary condition
c ::: ---------------------------------------------------------
      
      subroutine FORT_ROCKHOLD(rho,DIMS(rho),rhoval,ncomps,
     &     depth,pressure,nz,cplcoef,DIMS(cplcoef),ncplcoef,dx,grav,probhi)


      implicit none

      integer DIMDEC(rho)
      integer DIMDEC(cplcoef)
      integer ncomps, ncplcoef, nz
      REAL_T  rho(DIMV(rho),ncomps)
      REAL_T  cplcoef(DIMV(cplcoef),ncplcoef)
      REAL_T  rhoval(ncomps)
      REAL_T  dx(SDIM)
      REAL_T  probhi(SDIM)
      REAL_T  depth(nz),pressure(nz)
      REAL_T  grav
      
      integer i,j,k,l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res,phi,kappa
      REAL_T  y,dz
      integer zlo,zhi
      REAL_T  pc,stmp(2)

c     Gauss Legendre integration
      REAL_T, parameter :: qd(5) = (/ 0.0d0,
     &                                0.906179845938664d0*half,
     &                                0.538469310105683d0*half,
     &                               -0.906179845938664d0*half,
     &                               -0.538469310105683d0*half /)


      REAL_T, parameter :: wt(5) = (/ 0.568888888888889d0*half,
     &                                0.236926885056189d0*half,
     &                                0.478628670499366d0*half,
     &                                0.236926885056189d0*half,
     &                                0.478628670499366d0*half /)
      phi = zero
      kappa = zero
      dz = depth(nz)/dble(nz-1)
      do j = ARG_L2(rho), ARG_H2(rho)
         do i = ARG_L1(rho), ARG_H1(rho)
            rho(i,j,:) = zero
            do k = 1,5
               y = probhi(2) - (dble(j)+half+qd(k))*dx(2)

               if (y .gt. depth(nz)) then
                  rho(i,j,1) = rho(i,j,1) + rhoval(1)*wt(k)
               else if (y .lt. depth(1)) then
                  pc = -pressure(1)*grav*rhoval(1)
                  l_cpl_type = nint(cplcoef(i,j,1))
                  l_cpl_coef = cplcoef(i,j,2)
                  l_sigma    = cplcoef(i,j,3)
                  l_sat_res  = cplcoef(i,j,4)
                  call get_invcpl(stmp,pc,rhoval,phi,kappa,l_sigma,
     &              ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
                  rho(i,j,1) = rho(i,j,1) + stmp(1)*wt(k)
                  rho(i,j,2) = rho(i,j,2) + stmp(2)*wt(k)
               else
                  zlo = nint(y/dz+0.5)
                  zhi = nint(y/dz+0.5)+1
                  if (zlo .lt. 1 .or. zhi .lt. 1) then
                     print *, zlo, zhi,y
                  end if
                  pc = pressure(zlo) + (y-depth(zlo))/(depth(zhi)-depth(zlo))*(pressure(zhi)-pressure(zlo))
                  pc = -pc*grav*rhoval(1)
                  l_cpl_type = nint(cplcoef(i,j,1))
                  l_cpl_coef = cplcoef(i,j,2)
                  l_sigma    = cplcoef(i,j,3)
                  l_sat_res  = cplcoef(i,j,4)
                  call get_invcpl(stmp,pc,rhoval,phi,kappa,l_sigma,
     &              ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
                  rho(i,j,1) = rho(i,j,1) + stmp(1)*wt(k)
                  rho(i,j,2) = rho(i,j,2) + stmp(2)*wt(k)
               end if
            end do 
         end do
      end do 
      end

c ::: ---------------------------------------------------------
c ::: Define Rockhold boundary condition
c ::: ---------------------------------------------------------
      
c$$$      subroutine FORT_ROCKHOLD(rho,DIMS(rho),rhoval,ncomps,
c$$$     &     cplcoef,DIMS(cplcoef),ncplcoef,dx,wt_loc,grav,
c$$$     &     vgk,vga,vgn,ts,tr,nmat,q,h0,nn,nz,ns)
c$$$
c$$$      implicit none
c$$$
c$$$      integer DIMDEC(rho)
c$$$      integer DIMDEC(cplcoef)
c$$$      integer ncomps, ncplcoef
c$$$      REAL_T  rho(DIMV(rho),ncomps)
c$$$      REAL_T  cplcoef(DIMV(cplcoef),ncplcoef)
c$$$      REAL_T  rhoval(ncomps)
c$$$      REAL_T  dx(SDIM)
c$$$      REAL_T  wt_loc
c$$$      REAL_T  grav
c$$$      
c$$$      integer i,j,k,l_cpl_type
c$$$      REAL_T  l_cpl_coef, l_sigma, l_sat_res
c$$$      REAL_T  y
c$$$      REAL_T  pc,phi,kappa,satlim(2),stmp(2)
c$$$      integer nn,nmat,nz,ns
c$$$      integer mtype(nz)
c$$$      integer i,iter
c$$$      REAL_T  q,h0,hlo,hhi,tmpk,cond
c$$$      REAL_T  z(nz),vgk(nz),vga(nz),vgn(nz),ts(nz),tr(nz),
c$$$     &     h(nz),thb(nz),tht(nz),tt(nz),hk(nz),gam(ns),sl(ns),
c$$$     &     htmp(ns+1),uk(ns+1)
c$$$
c$$$      REAL_T, allocatable :: z(:),vgk(:),
c$$$      
c$$$
c$$$      hlo = 0
c$$$      hhi = 25
c$$$    
c$$$      do i = 1,nmat
c$$$        iter = 0
c$$$        tmpk = cond(hhi,vga(i),vgn(i),vgk(i))
c$$$        if ( tmpk .ge. q ) then
c$$$  100     hhi = hhi*2.d0
c$$$          iter = iter + 1
c$$$          tmpk = cond(hhi,vga(i),vgn(i),vgk(i))
c$$$          if( tmpk .ge. q ) then
c$$$              goto 100
c$$$           end if
c$$$        end if
c$$$      end do
c$$$      
c$$$      
c$$$      do i = 1,nn
c$$$
c$$$        read(1,*) z(i),mtype(i)
c$$$        if( i .gt. 1 ) then
c$$$          if( z(i) .ge. z(i-1) ) then
c$$$            write(*,*)' '
c$$$            write(*,*)' program terminating...'
c$$$            write(*,*)' z must be in decreasing order.'
c$$$            stop
c$$$          endif
c$$$        endif
c$$$      end do
c$$$
c$$$      dh = (hhi - hlo)/float(ns)
c$$$      h(1) = h0
c$$$      jj = mtype(1)
c$$$
c$$$      call slope(ns,hlo,dh,vga(jj),vgn(jj),vgk(jj),htmp,sl,uk)
c$$$      call gamma(ns,q,uk,sl,htmp,gam)
c$$$
c$$$      do j = 2,nn
c$$$        if( mtype(j) .ne. mtype(j-1) ) then
c$$$          jj = mtype(j)
c$$$          call slope(ns,hlo,dh,vga(jj),vgn(jj),vgk(jj),htmp,sl,uk)
c$$$          call gamma(ns,q,uk,sl,htmp,gam)
c$$$        endif
c$$$        i = max(1,int((h(j-1) - hlo)/dh + 1.))   
c$$$        tst = gam(i)*exp(sl(i)*h(j-1)) - 1.d0
c$$$        if( tst .eq. 0.d0 ) then
c$$$          h(j) = h(j-1)
c$$$        else
c$$$          z0 = z(j-1)
c$$$          dz = z(j) - z(j-1)
c$$$          if( h(j-1) .lt. htmp(1) ) then
c$$$            hnew = g0(q,vgk(jj),dz,h(j-1))
c$$$            if( hnew .le. htmp(1) ) then
c$$$              h(j) = hnew
c$$$            else
c$$$              z0 = z0 + f0(htmp(1),q,vgk(jj)) - f0(h(j-1),q,vgk(jj))
c$$$              dz = z(j) - z0
c$$$              hnew = g(sl(1),gam(1),htmp(1),dz)
c$$$              if( hnew .le. htmp(2) ) then
c$$$                h(j) = hnew
c$$$              else
c$$$                do kk = 2,ns
c$$$                  z0 = z0 + f(sl(kk),gam(kk),htmp(kk))
c$$$     &                    - f(sl(kk-1),gam(kk-1),htmp(kk-1))
c$$$                  dz = z(j) - z0
c$$$                  hnew = g(sl(kk),gam(kk),htmp(kk),dz)
c$$$                  if( hnew .le. htmp(kk) ) then
c$$$                    h(j) = hnew
c$$$                    goto 600
c$$$                  endif
c$$$                end do
c$$$              endif
c$$$            endif
c$$$          else
c$$$            hnew = g(sl(i),gam(i),h(j-1),dz)
c$$$            if( hnew .ge. htmp(i) .and. hnew .lt. htmp(i+1)) then
c$$$              h(j) = hnew
c$$$            elseif( hnew .lt. htmp(i) ) then
c$$$              if( i .eq. 1 ) then
c$$$                dz = dz + f(sl(1),gam(1),htmp(1)) 
c$$$     &                  - f(sl(1),gam(1),h(j-1))
c$$$                h(j) = g0(q,vgk(jj),dz,htmp(1))
c$$$              else
c$$$                k = i
c$$$                do kk = k,2,-1 
c$$$                  if( kk .eq. k ) then
c$$$                    z0 = z0 + f(sl(kk),gam(kk),htmp(kk))
c$$$     &                      - f(sl(kk),gam(kk),h(j-1))
c$$$                  else
c$$$                    z0 = z0 + f(sl(kk),gam(kk),htmp(kk))
c$$$     &                      - f(sl(kk),gam(kk),htmp(kk+1))
c$$$                  endif
c$$$                  dz = z(j) - z0
c$$$                  hnew = g(sl(kk-1),gam(kk-1),htmp(kk),dz)
c$$$                  if( hnew .ge. htmp(kk-1) ) then
c$$$                    h(j) = hnew
c$$$                    goto 600
c$$$                  endif
c$$$                end do
c$$$                z0 = z0 + f(sl(1),gam(1),htmp(1)) 
c$$$     &                  - f(sl(1),gam(1),h(j-1))
c$$$                dz = z(j) - z0
c$$$                h(j) = g0(q,vgk(jj),dz,htmp(1)) 
c$$$              endif
c$$$            elseif( hnew .ge. htmp(i+1) ) then
c$$$              k = i
c$$$              do kk = k,ns-2
c$$$                if( kk .eq. k ) then
c$$$                  z0 = z0 + f(sl(kk),gam(kk),htmp(kk+1))
c$$$     &                    - f(sl(kk),gam(kk),h(j-1))
c$$$                else
c$$$                  z0 = z0 + f(sl(kk),gam(kk),htmp(kk+1))
c$$$     &                    - f(sl(kk),gam(kk),htmp(kk))
c$$$                endif
c$$$                dz = z(j) - z0
c$$$                hnew = g(sl(kk+1),gam(kk+1),htmp(kk+1),dz)
c$$$                if( hnew .lt. htmp(kk+2) ) then
c$$$                  h(j) = hnew
c$$$                  goto 600
c$$$                endif
c$$$              end do
c$$$            endif
c$$$          endif
c$$$  600   endif
c$$$      end do
c$$$
c$$$      contains 
c$$$        subroutine gamma(ns,q,uk,sl,htmp,gam)
c$$$        implicit none
c$$$        integer ns
c$$$        REAL_T q, uk(ns+1),sl(ns),htmp(ns+1),gam(ns)
c$$$        do i = 1,ns
c$$$           gam(i) = q/uk(i)*dexp(-sl(i)*htmp(i))
c$$$        end do
c$$$        end
c$$$      
c$$$        subroutine get_theta(nz,nn,mtype,ts,tr,vga,vgn,h,thb,tht)
c$$$        implicit none
c$$$        integer i,jb,jt,nz,nn
c$$$        integer mtype(nz)
c$$$        REAL_T  vga(nz),vgn(nz),ts(nz),tr(nz),
c$$$     &       h(nz),thb(nz),tht(nz),theta
c$$$        do i = 1,nn
c$$$           jb = mtype(i)
c$$$           thb(i) = theta(ts(jb),tr(jb),vga(jb),vgn(jb),h(i))
c$$$           tht(i) = thb(i)
c$$$           if( i .lt. nn ) then
c$$$              if( mtype(i+1) .ne. mtype(i) ) then
c$$$                 jt = mtype(i+1)
c$$$                 tht(i) = theta(ts(jt),tr(jt),vga(jt),vgn(jt),h(i))
c$$$              endif
c$$$           endif
c$$$        end do
c$$$        end
c$$$       
c$$$        subroutine slope(ns,h1,dh,vga,vgn,vgk,htmp,sl,uk)
c$$$        implicit none
c$$$        integer ns,i
c$$$        REAL_T h1,dh,vga,vgn,vgk,htmp(ns+1),sl(ns),uk(ns+1)
c$$$        do i = 1,ns+1
c$$$           htmp(i) = h1 + float(i-1)*dh
c$$$           uk(i) = cond(htmp(i),vga,vgn,vgk)
c$$$        end do
c$$$        do i = 1,ns
c$$$           sl(i) = dlog(uk(i)/uk(i+1))/(htmp(i+1) - htmp(i))
c$$$           sl(i) = max(sl(i),1.d-15)
c$$$        end do
c$$$        end
c$$$
c$$$        REAL_T function cond(h,vga,vgn,vgk)
c$$$        implicit none
c$$$        REAL_T h, vga, vgn, vgk, t1,vgm
c$$$        if( h .le. 0.d0 ) then
c$$$           cond = vgk
c$$$        else 
c$$$           vgm = 1.d0 - 1.d0/vgn
c$$$           t1 = 1.d0 + (vga*h)**vgn
c$$$           cond = (vgk*(1.d0 - (vga*h)**(vgn - 1.d0)
c$$$     &          *t1**(-vgm))**2.d0)/t1**(vgm*0.5d0)
c$$$        endif
c$$$        end function cond
c$$$
c$$$        REAL_T function theta(ts,tr,vga,vgn,h)
c$$$        implicit none
c$$$        REAL_T ts,tr,vga,vgn,h,vgm
c$$$        vgm = 1.d0 - 1.d0/vgn
c$$$        if( h .ge. 0.d0 ) then
c$$$           theta = ((1.d0 + (vga*h)**vgn)**(-vgm))*(ts - tr) + tr
c$$$        else
c$$$           theta = ts
c$$$        endif
c$$$        end function theta
c$$$      
c$$$        REAL_T function g(a,gam,ho,z)
c$$$        implicit none
c$$$        REAL_T a,gam,ho,z
c$$$        g = -1.d0/a*dlog(gam - (gam*dexp(a*ho) - 1.d0)*dexp(a*(z-ho)))
c$$$        return
c$$$        end function g
c$$$       
c$$$        REAL_T function g0(q,vgk,dz,h0)
c$$$        implicit none
c$$$        REAL_T q,vgk,dz,h0
c$$$        g0 = (q/vgk - 1.d0)*dz + h0
c$$$        end function g0
c$$$
c$$$        REAL_T function f(a,g,h)
c$$$        implicit none
c$$$        REAL_T a,g,h
c$$$        f = 1.d0/a*dlog(abs(g*dexp(a*h) - 1.d0)) - h
c$$$        end function f
c$$$
c$$$        REAL_T function f0(h,q,vgk)
c$$$        implicit none
c$$$        REAL_T h,q,vgk
c$$$        f0 = h/(q/vgk - 1.d0)
c$$$        end function f0
c$$$
c$$$      end

c :::
c ::: Fill ghost cells if nghost > 1 
c :::
      subroutine FORT_PATCH_GHOST(p,DIMS(p),ncomp,face,domlo,domhi)

      implicit none

      integer DIMDEC(p)
      integer face
      integer ncomp
      integer domlo(SDIM)
      integer domhi(SDIM)
      REAL_T  p(DIMV(p),ncomp)

c     Local variables
      integer i,j,n

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

c     y-faces
      if (face .eq. YLO) then
         do j = ARG_L2(p),domlo(2)-2
            do n = 1,ncomp
               p(:,j,n) = p(:,domlo(2)-1,n)
            end do
         end do

      elseif (face .eq. YHI) then
         do j = domhi(2)+2,ARG_H2(p)
            do n = 1,ncomp
               p(:,j,n) = p(:,domhi(2)+1,n)
            end do
         end do
      end if

      
c     x-faces
      if (face .eq. XLO) then
         do i = ARG_L1(p),domlo(1)-2
            do n = 1,ncomp
               p(i,:,n) = p(domlo(1)-1,:,n)
            end do
         end do

      elseif (face .eq. XHI) then
         do i = domhi(1)+2,ARG_H1(p)
            do n = 1,ncomp
               p(i,:,n) = p(domhi(1)+1,:,n)
            end do
         end do
      end if

#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end
