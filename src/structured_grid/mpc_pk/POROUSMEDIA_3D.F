c
c $Id: POROUSMEDIA_3D.F,v 1.35 2011-06-07 21:08:16 gpau Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "POROUSMEDIA_F.H"
#include "ArrayLim.H"

#define SDIM 3

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

      integer  i, j, k, n, ic, jc, kc, ioff, joff, koff
      integer  lratx, lraty, lratz

      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = zero
               end do
            end do
         end do
c         
c     sum fine data
c
!$omp parallel private(i,j,k,ic,jc,kc,ioff,joff,koff)
         do koff = 0, lratz-1
!$omp do
            do kc = lo(3),hi(3)
               k = kc*lratz + koff
               do joff = 0, lraty-1
                  do jc = lo(2), hi(2)
                     j = jc*lraty + joff
                     do ioff = 0, lratx-1
                        do ic = lo(1), hi(1)
                           i = ic*lratx + ioff
                           crse(ic,jc,kc,n) = crse(ic,jc,kc,n) + fv(i,j,k)*fine(i,j,k,n)
                        end do
                     end do
                  end do
               end do
            end do
!$omp end do nowait
         end do
!$omp end parallel
c
c     divide out by volume weight
c
!$omp parallel do private(ic,jc,kc)
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = crse(ic,jc,kc,n)/cv(ic,jc,kc)
               end do
            end do
         end do
!$omp end parallel do

      end do

      end

c ::: ----------------------------------------------------------
c ::: SUMMASS
c :::             MASS = sum{ vol(i,j,k)*rho(i,j,k) }
c :::
c ::: INPUTS / OUTPUTS:
c :::  rho        => density field
c :::  DIMS(rho)  => index limits of rho aray
c :::  lo,hi      => index limits of grid interior
c :::  delta	 => cell size
c :::  mass      <=  total mass
c ::: ----------------------------------------------------------

       subroutine FORT_SUMMASS(rho,DIMS(rho),DIMS(grid),delta,mass)

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM)
       REAL_T  rho(DIMV(rho))

       integer i, j, k
       REAL_T  vol

       vol = delta(1)*delta(2)*delta(3)

       mass = zero

       do k = ARG_L3(grid), ARG_H3(grid)
          do j = ARG_L2(grid), ARG_H2(grid)
             do i = ARG_L1(grid), ARG_H1(grid)
                mass = mass + rho(i,j,k)
             end do
          end do
       end do

       mass = mass*vol

       end


c ::: ----------------------------------------------------------
c ::: SUMMASSSQ
c :::             MASS = sum{ vol(i,j)*rho(i,j)*rho(i,j) }
c :::
c ::: INPUTS / OUTPUTS:
c :::  rho        => density field
c :::  DIMS(rho)  => index limits of rho aray
c :::  lo,hi      => index limits of grid interior
c :::  delta	 => cell size
c :::  mass      <=  total mass
c :::  tmp        => temp column array
c ::: ----------------------------------------------------------

       subroutine FORT_SUMMASSSQ(rho,DIMS(rho),DIMS(grid),delta,mass)

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM)
       REAL_T  rho(DIMV(rho))

       integer i, j, k
       REAL_T  vol

       vol = delta(1)*delta(2)*delta(3)

       mass = zero

       do k = ARG_L3(grid), ARG_H3(grid)
          do j = ARG_L2(grid), ARG_H2(grid)
             do i = ARG_L1(grid), ARG_H1(grid)
                mass = mass + rho(i,j,k)*rho(i,j,k)
             end do
          end do
       end do

       mass = mass*vol

       end

c ::: ----------------------------------------------------------
c ::: SUMPRESS
c :::             SUM = sum{ press(i,j,k) * surroundingVolume(i,j,k) }
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

       integer lo(3),hi(3)
       integer DIMDEC(press)
       integer DIMDEC(vol)
       REAL_T  press(DIMV(press))
       REAL_T  vol(DIMV(vol))
       REAL_T  total_press,total_vol

       integer i, j, k
       REAL_T  surroundingVolume

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                surroundingVolume = vol(i,j  ,k  ) + vol(i-1,j  ,k  ) + 
     &               vol(i,j-1,k  ) + vol(i-1,j-1,k  ) +
     &               vol(i,j  ,k-1) + vol(i-1,j  ,k-1) + 
     &               vol(i,j-1,k-1) + vol(i-1,j-1,k-1)
                total_press = total_press + press(i,j,k) * surroundingVolume
                total_vol   = total_vol   +                surroundingVolume

             end do
          end do
       end do

       end

c ::: -----------------------------------------------------------------------
c ::: This routine fills an edge-centered fab from a cell-centered
c ::: fab using simple linear interpolation.
c ::: INPUTS / OUTPUTS:
c ::: lo,hi      => index limits of the region of the edge-centered fab
c                   to be filled
c ::: DIMS(cfab) => index limits of the cell-centered fab
c ::: cfab       => cell-centered data
c ::: DIMS(efab) => index limits of the edge-centered fab
c ::: efab       => edge-centered fab to fill
c ::: nc         => Number of components in the fab to fill
c ::: dir        => direction data needs to be shifted to get to edges
c ::: -----------------------------------------------------------------------

      subroutine FORT_CEN2EDG(lo, hi, 
     &     DIMS(cfab), cfab,
     &     DIMS(efab), efab, nc, dir, isharm
     &     )
      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)
      integer i,j,k,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1)+1, hi(1)
                        efab(i,j,k,n) =
     &                       half*(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2)+1, hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) =
     &                       half*(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3)+1, hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) =
     &                       half*(cfab(i,j,k,n) + cfab(i,j,k-1,n))
                     end do
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1)+1, hi(1)
                        if((cfab(i,j,k,n) * cfab(i-1,j,k,n)) .gt.0.d0)
     &                       then
                           efab(i,j,k,n) =
     &                          2.d0*(cfab(i,j,k,n) * cfab(i-1,j,k,n))/
     &                          (cfab(i,j,k,n) + cfab(i-1,j,k,n))
                        else
                           efab(i,j,k,n)=0.d0
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2)+1, hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j-1,k,n)).gt.0.d0)
     &                       then
                           efab(i,j,k,n) =
     &                          2.d0*(cfab(i,j,k,n) * cfab(i,j-1,k,n))/
     &                          (cfab(i,j,k,n) + cfab(i,j-1,k,n))
                        else
                           efab(i,j,k,n)=0.d0
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3)+1, hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j,k-1,n)).gt.0.d0)
     &                       then
                           efab(i,j,k,n) =
     &                          2.d0*(cfab(i,j,k,n) * cfab(i,j,k-1,n))/
     &                          (cfab(i,j,k,n) + cfab(i,j,k-1,n))
                        else
                           efab(i,j,k,n)=0.d0
                        endif
                     end do
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

      subroutine FORT_EDGE_INTERP(flo, fhi, nc, ratio, dir,
     &     fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)
      implicit none
      integer flo(0:3-1), fhi(0:3-1), nc, ratio(0:3-1), dir
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION
     &     fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,n,P,M,L
      DOUBLE PRECISION val, df, dtmp

      dtmp = 1.0d0 / dble(ratio(dir))

      if (dir.eq.0) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                     df = fine(i+ratio(dir),j,k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) + df*dble(M)*dtmp
                        do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(i+M,P,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1)-ratio(dir),ratio(1)
                  do i=flo(0),fhi(0)
                     df = fine(i,j+ratio(dir),k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) + df*dble(M)*dtmp
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(P,j+M,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=flo(2),fhi(2)-ratio(dir),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0),ratio(0)
                     df = fine(i,j,k+ratio(dir),n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) + df*dble(M)*dtmp
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                              fine(P,L,k+M,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      endif
      end

c ::: ------------------------------------------------------------------------------------
c ::: For edge-based data, fill fine values with piecewise-constant interp of coarse data.
c ::: Operate only on faces that overlap--ie, only fill the fine faces that make up each
c ::: coarse face, leave the in-between faces alone.
c ::: ------------------------------------------------------------------------------------

      subroutine FORT_PC_CF_EDGE_INTERP(lo, hi, nc, ratio, dir,
     &     crse, crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2,
     &     fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)
      implicit none
      integer lo(3),hi(3), nc, ratio(0:3-1), dir
      integer crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION
     &     crse(crse_l0:crse_h0,crse_l1:crse_h1,crse_l2:crse_h2,nc)
      DOUBLE PRECISION
     &     fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,L, P

      if (dir.eq.0) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(1)-1
                           fine(ii,jj+L,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(1)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj+P,kk,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
      end

c ::: ---------------------------------------------------------------
c ::: For edge-based data, fill coarse values with average of fine data 
c ::: one the same edge
c ::: ----------------------------------------------------------------

      subroutine FORT_EDGE_AVGDOWN(lo, hi, nc, ratio, dir,
     &     crse, crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2,
     &     fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)

      implicit none
      integer lo(SDIM),hi(SDIM), nc, ratio(0:SDIM-1), dir
      integer crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION
     &     crse(crse_l0:crse_h0,crse_l1:crse_h1,crse_l2:crse_h2,nc)
      DOUBLE PRECISION
     &     fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,L, P
      DOUBLE PRECISION dtmp

      if (dir.eq.0) then
         dtmp = 1.0d0 / (ratio(1)*ratio(2))
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     crse(i,j,k,n) = zero
                     do P=0,ratio(2)-1
                        do L=0,ratio(1)-1
                           crse(i,j,k,n) = crse(i,j,k,n) + fine(ii,jj+L,kk+P,n) 
                        enddo
                     enddo
                     crse(i,j,k,n) = crse(i,j,k,n)*dtmp
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         dtmp = 1.0d0 / (ratio(0)*ratio(2))
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     crse(i,j,k,n) = zero
                     do P=0,ratio(2)-1
                        do L=0,ratio(0)-1
                           crse(i,j,k,n) = crse(i,j,k,n) + fine(ii+L,jj,kk+P,n) 
                        enddo
                     enddo
                     crse(i,j,k,n) = crse(i,j,k,n)*dtmp
                  enddo
               enddo
            enddo
         enddo
      else
         dtmp = 1.0d0 / (ratio(0)*ratio(1))
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     crse(i,j,k,n) = zero
                     do P=0,ratio(1)-1
                        do L=0,ratio(0)-1
                           crse(i,j,k,n) = crse(i,j,k,n) + fine(ii+L,jj+P,kk,n) 
                        enddo
                     enddo
                     crse(i,j,k,n) = crse(i,j,k,n)*dtmp
                  enddo
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
      REAL_T S(DIMV(S),ncomp)
      REAL_T B(DIMV(B),ncomp)
      REAL_T phi(DIMV(phi))
      REAL_T visc_coef(ncomp)
      
      integer i, j, k, n
      REAL_T rho

      if (diffuse .ne. -1) then
         do k = lo(3),hi(3)
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  do n=1,ncomp
                     B(i,j,k,n) = visc_coef(n)
                  end do
               end do
            end do
         end do         

      else
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  rho = zero
                  do n=1,ncomp
                     rho = rho + S(i,j,k,n) 
                  end do
                  do n=1,ncomp
                     B(i,j,k,n) = rho*phi(i,j,k)*visc_coef(n)
                  end do
               end do
            end do
         end do
      end if

      end
      
      subroutine  FORT_RECOMP_UPDATE(lo, hi,
     &                               update, DIMS(update),
     &                               phi,    DIMS(phi),
     &                               xflux,  DIMS(xflux),
     &                               yflux,  DIMS(yflux),
     &                               zflux,  DIMS(zflux),
     &                               vol,    DIMS(vol),
     &                               nc)
      implicit none
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(update)
      integer DIMDEC(phi)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(vol)
      REAL_T update(DIMV(update),nc)
      REAL_T phi(DIMV(phi))
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T zflux(DIMV(zflux),nc)
      REAL_T vol(DIMV(vol))

      integer i, j, k, n

      do n = 1,nc      
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  update(i,j,k,n)=-((xflux(i+1,j,k,n)-xflux(i,j,k,n))
     &                 + (yflux(i,j+1,k,n)-yflux(i,j,k,n))
     &                 + (zflux(i,j,k+1,n)-zflux(i,j,k,n)))
     &                 /vol(i,j,k)/phi(i,j,k)
               end do
            end do
         end do
      end do
      end


c ::: ----------------------------------------------------------
c ::: CHEM_DUMMY
c :::    a dummy chemical reaction.
c ::: ----------------------------------------------------------
      subroutine FORT_CHEM_DUMMY(s,DIMS(s),lo,hi,dt,nc)
      implicit none
      integer lo(SDIM),hi(SDIM), nc
      integer DIMDEC(s)
      REAL_T  s(DIMV(s),nc), dt

      integer i,j,k,r
      REAL_T  deltan,rate

      rate = 1.d-8

      do k = ARG_L3(s),ARG_H3(s)
         do j =  ARG_L2(s),ARG_H2(s)
            do i = ARG_L1(s),ARG_H1(s)
               deltan = rate*dt*s(i,j,k,1)*s(i,j,k,2)
               r = s(i,j,k,1)/(s(i,j,k,1)+s(i,j,k,2))
               s(i,j,k,1) = s(i,j,k,1) - deltan
               s(i,j,k,2) = s(i,j,k,2) - deltan
               s(i,j,k,3) = 1000.0d0 - s(i,j,k,1) - s(i,j,k,2)
            end do
         end do 
      end do

      end


c ::: -----------------------------------------------------------
c ::: Compute the capillary pressure, p_c
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

      integer i,j,k
      integer l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res, dpc

      if (ncomps .ne. nphases .and. ncomps .ne. 2) then
         call bl_abort('FORT_MK_CPL: only for ncomps = nphases = 2.')
      end if

!$omp parallel do private(i,j,k)
      do k = ARG_L3(pc),ARG_H3(pc)
         do j = ARG_L2(pc),ARG_H2(pc)
            do i = ARG_L1(pc),ARG_H1(pc)
               l_cpl_type = nint(coef(i,j,k,1))
               l_cpl_coef = coef(i,j,k,2)
               l_sigma    = coef(i,j,k,3)
               l_sat_res  = coef(i,j,k,4)
               if (i .eq. lo(1)-1 .and. bc(1,1) .eq. HOEXTRAP) then
                  dpc = pc(i+2,j,k)-pc(i+1,j,k)
                  pc(i,j,k) = pc(i+1,j,k) - dpc/2.0
               elseif (i .eq. hi(1)+1 .and. bc(1,2) .eq. HOEXTRAP) then
                  dpc = pc(i-1,j,k)-pc(i-2,j,k)
                  pc(i,j,k) = pc(i-1,j,k) + dpc/2.0
               elseif (j .eq. lo(2)-1 .and. bc(2,1) .eq. HOEXTRAP) then
                  dpc = pc(i,j+2,k)-pc(i,j+1,k)
                  pc(i,j,k) = pc(i,j+1,k) - dpc/2.0
               elseif (j .eq. hi(2)+1 .and. bc(2,2) .eq. HOEXTRAP) then
                  dpc = pc(i,j-1,k)-pc(i,j-2,k)
                  pc(i,j,k) = pc(i,j-1,k) + dpc/2.0
               elseif (k .eq. lo(3)-1 .and. bc(3,1) .eq. HOEXTRAP) then
                  dpc = pc(i,j,k+2)-pc(i,j,k+1)
                  pc(i,j,k) = pc(i,j,k+1) - dpc/2.0
               elseif (j .eq. hi(3)+1 .and. bc(3,2) .eq. HOEXTRAP) then
                  dpc = pc(i,j,k-1)-pc(i,j,k-2)
                  pc(i,j,k) = pc(i,j,k-1) + dpc/2.0
               else
                  call get_cpl(pc(i,j,k),s(i,j,k,1:2),rhoval,phi(i,j,k),
     &                 kappa(i,j,k),l_sigma,ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
               end if
            end do
         end do 
      end do
!$omp end parallel do

      end

c ::: -----------------------------------------------------------
c ::: Compute the lambda = kr/mu for all components
c :::    It calls get_ld in EOS_F.F. 
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

      integer i,j,k
      integer l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T  tlbd(ncomps)

!$omp parallel do private(i,j,k,tlbd)
      do k = ARG_L3(lbd),ARG_H3(lbd)
         do j = ARG_L2(lbd),ARG_H2(lbd)
            do i = ARG_L1(lbd),ARG_H1(lbd)
               l_kr_type = nint(kr(i,j,k,1))
               l_kr_coef = kr(i,j,k,2)
               l_sat_res = kr(i,j,k,3)
               call get_ld(tlbd,s(i,j,k,:),rhoval,muval,ncomps,
     &              l_kr_type,l_kr_coef,l_sat_res)
               lbd(i,j,k,:) = tlbd
            end do
         end do 
      end do
!$omp end parallel do

      end

c ::: -----------------------------------------------------------
c ::: Compute the derivative of lambda = kr/mu for all components
c :::    It calls get_ld in EOS_F.F. 
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

      integer i,j,k
      integer l_kr_type
      REAL_T  l_kr_coef, l_sat_res
      REAL_T  tlbd(3)

!$omp parallel do private(i,j,k,tlbd)
      do k = ARG_L3(lbd),ARG_H3(lbd)
         do j = ARG_L2(lbd),ARG_H2(lbd)
            do i = ARG_L1(lbd),ARG_H1(lbd)
               l_kr_type = nint(kr(i,j,k,1))
               l_kr_coef = kr(i,j,k,2)
               l_sat_res = kr(i,j,k,3)
               call get_dld(tlbd,s(i,j,k,:),rhoval,muval,ncomps,
     &              l_kr_type,l_kr_coef,l_sat_res)
               lbd(i,j,k,:) = tlbd
            end do
         end do 
      end do
!$omp end parallel do

      end

c ::: ----------------------------------------------------------
c ::: GETDIFFUSE_CPL
c :::    Get the diffusive coefficient for the capillary solve.
c :::    We only solve for water
c ::: ----------------------------------------------------------         
      subroutine FORT_GETDIFFUSE_CPL(lbd, DIMS(lbd), 
     &     dfx, DIMS(dfx), dfy, DIMS(dfy), dfz, DIMS(dfz),
     &     lo, hi, domlo, domhi, bc, ncomp)

      implicit none

#include "probdata.H"

      integer lo(SDIM),hi(SDIM)
      integer domlo(SDIM),domhi(SDIM)
      integer DIMDEC(lbd)
      integer DIMDEC(dfx)
      integer DIMDEC(dfy)
      integer DIMDEC(dfz)
      integer ncomp
      integer bc(SDIM,2)
      REAL_T  lbd(DIMV(lbd),ncomps)
      REAL_T  dfx(DIMV(dfx),ncomp)
      REAL_T  dfy(DIMV(dfy),ncomp)
      REAL_T  dfz(DIMV(dfz),ncomp)

      integer i, j, k
      REAL_T  lambda(ncomps),H_coef,tmp(2)

      if (model .ne. 2) then
         call bl_abort('FORT_GETDIFFUSE_CPL: valid only for model = 2.')
      end if

!$omp parallel private(i,j,k,lambda,H_coef,tmp)
c ::: diffusion coefficient in x-dir
!$omp do
      do k = lo(3), hi(3)
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)+1
               lambda = lbd(i-1,j,k,:)
               tmp(1) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               lambda = lbd(i,j,k,:)
               tmp(2) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
c               H_coef = half*(tmp(1)+tmp(2))
               H_coef = one/tmp(1) + one/tmp(2)
               H_coef = two/H_coef

               if (i .eq. domlo(1)) then
                  if (bc(1,1) .eq. EXT_DIR) then
                     lambda = lbd(i,j,k,:)
                  elseif (bc(1,1) .eq. SEEPAGE) then
                     lambda = lbd(i,j,k,:)
                  end if
                  H_coef  = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               elseif (i .eq. domhi(1)) then
                  if (bc(1,2) .eq. EXT_DIR) then
                     lambda = lbd(i,j,k,:)
                  elseif (bc(1,2) .eq. SEEPAGE) then
                     lambda = lbd(i,j,k,:)
                  end if       
                  H_coef  = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               endif

               dfx(i,j,k,1) = H_coef*rhoval(1)

            end do
         end do
      end do
!$omp end do nowait

c ::: diffusion coefficient in y-dir
!$omp do
      do k = lo(3), hi(3)
         do j=lo(2), hi(2)+1
            do i=lo(1), hi(1)
               lambda = lbd(i,j-1,k,:)
               tmp(1) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               lambda = lbd(i,j,k,:)
               tmp(2) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
c               H_coef = half*(tmp(1)+tmp(2))
               H_coef = one/tmp(1) + one/tmp(2)
               H_coef = two/H_coef
               if (j .eq. domlo(2)) then
                  if (bc(2,1) .eq. EXT_DIR) then
                     lambda = lbd(i,j,k,:)
                  elseif (bc(1,1) .eq. SEEPAGE) then
                     lambda = lbd(i,j,k,:)
                  end if
                  H_coef  = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               elseif (j .eq. domhi(2)) then
                  if (bc(2,2) .eq. EXT_DIR) then
                     lambda = lbd(i,j,k,:)
                  elseif (bc(2,2) .eq. SEEPAGE) then
                     lambda = lbd(i,j,k,:)
                  end if       
                  H_coef  = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               endif

               dfy(i,j,k,1) = H_coef*rhoval(1)

            end do
         end do
      end do
!$omp end do nowait

c ::: diffusion coefficient in z-dir
!$omp do
      do k = lo(3), hi(3)+1
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)
               lambda = lbd(i,j,k-1,:)
               tmp(1) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               lambda = lbd(i,j,k,:)
               tmp(2) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
c               H_coef = half*(tmp(1)+tmp(2))
               H_coef = one/tmp(1) + one/tmp(2)
               H_coef = two/H_coef
               if (k .eq. domlo(3)) then
                  if (bc(3,1) .eq. EXT_DIR) then
                     lambda = lbd(i,j,k,:)
                  elseif (bc(3,1) .eq. SEEPAGE) then
                     lambda = lbd(i,j,k,:)
                  end if
                  H_coef  = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               elseif (k .eq. domhi(3)) then
                  if (bc(3,2) .eq. EXT_DIR) then
                     lambda = lbd(i,j,k,:)
                  elseif (bc(3,2) .eq. SEEPAGE) then
                     lambda = lbd(i,j,k,:)
                  end if       
                  H_coef  = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               endif

               dfz(i,j,k,1) = H_coef*rhoval(1)

            end do
         end do
      end do
!$omp end do
!$omp end parallel

      end

c ::: ----------------------------------------------------------
c ::: GETDIFFUSE_CPL_dp
c :::    get the diffusive coefficient for the capillary term with dp term
c ::: ----------------------------------------------------------         
      subroutine FORT_GETDIFFUSE_CPL_dp(s, h_cc, DIMS(s), lbd, DIMS(lbd),
     &     dfx, DIMS(dfx), dfy, DIMS(dfy), dfz, DIMS(dfz),
     &     phi, DIMS(phi), kappa, DIMS(kappa),    
     &     coef,DIMS(coef),n_coef,
     &     lo, hi, domlo, domhi, bc, ncomp)

      implicit none

#include "probdata.H"

      integer DIMDEC(s)
      integer DIMDEC(lbd)
      integer DIMDEC(dfx)
      integer DIMDEC(dfy)
      integer DIMDEC(dfz)
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
      REAL_T  dfz(DIMV(dfz),ncomp)
      REAL_T  phi(DIMV(phi))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  coef(DIMV(coef),n_coef)

      integer i, j, k
      REAL_T  lambda(ncomps),H_coef,dpcds
      integer l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res
      if (model .ne. 2) then
         call bl_abort('FORT_GETDIFFUSE_CPL: valid only for model = 2.')
      end if

!$omp parallel private(i,j,k,lambda,H_coef,tmp,dpcds)
!$omp do
      do k = ARG_L3(s),ARG_H3(s)
         do j = ARG_L2(s),ARG_H2(s)
            do i = ARG_L1(s),ARG_H1(s)
               lambda = lbd(i,j,k,1:2)
               h_cc(i,j,k) = lambda(1)*lambda(2)/(lambda(1)+lambda(2))
               if (h_cc(i,j,k) .gt. 1.d-12) then
                  l_cpl_type = nint(coef(i,j,k,1))
                  l_cpl_coef = coef(i,j,k,2)
                  l_sigma    = coef(i,j,k,3)
                  l_sat_res  = coef(i,j,k,4)
                  call get_dcplds(dpcds,s(i,j,k,:),rhoval,phi(i,j,k),
     &              kappa(i,j,k),l_sigma,ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
                  h_cc(i,j,k) = h_cc(i,j,k)*dpcds*rhoval(1)
               end if
            end do
         end do
      end do
!$omp end do nowait

c     diffusion coefficient in x-dir
!$omp do
      do k = lo(3), hi(3)
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)+1

c               H_coef = half*(h_cc(i-1,j,k)+h_cc(i,j,k))
               H_coef = one/h_cc(i-1,j,k) + one/h_cc(i,j,k)
               H_coef = two/H_coef
               if (i .eq. domlo(1)) then
                  if (bc(1,1) .eq. EXT_DIR) then
                     H_coef = h_cc(i,j,k)
                  elseif (bc(1,1) .eq. SEEPAGE) then
                     H_coef = h_cc(i,j,k)
                  end if
               elseif (i .eq. domhi(1)) then
                  if (bc(1,2) .eq. EXT_DIR) then
                     H_coef = h_cc(i,j,k)
                  elseif (bc(1,2) .eq. SEEPAGE) then
                     H_coef = h_cc(i,j,k)
                  end if  
               end if
               dfx(i,j,k,1) = H_coef 
            end do
         end do
      end do
!$omp end do nowait

c     diffusion coefficient in y-dir
!$omp do
      do k = lo(3), hi(3)
         do j=lo(2), hi(2)+1
            do i=lo(1), hi(1)
c               H_coef = half*(h_cc(i,j-1,k)+h_cc(i,j,k))
               H_coef = one/h_cc(i,j-1,k) + one/h_cc(i,j,k)
               H_coef = two/H_coef
               if (j .eq. domlo(2)) then
                  if (bc(2,1) .eq. EXT_DIR) then
                     H_coef = h_cc(i,j,k)
                  elseif (bc(2,1) .eq. SEEPAGE) then
                     H_coef = h_cc(i,j,k)
                  end if
               elseif (j .eq. domhi(2)) then
                  if (bc(2,2) .eq. EXT_DIR) then
                     H_coef = h_cc(i,j,k)
                  elseif (bc(2,2) .eq. SEEPAGE) then
                     H_coef = h_cc(i,j,k)
                  end if  
               end if
               dfy(i,j,k,1) = H_coef
            end do
         end do
      end do
!$omp end do nowait

c     diffusion coefficient in z-dir
!$omp do
      do k = lo(3), hi(3)+1
         do j=lo(2), hi(2)
            do i=lo(1), hi(1)
               H_coef = half*(h_cc(i,j,k-1)+h_cc(i,j,k))
               H_coef = one/h_cc(i,j,k-1) + one/h_cc(i,j,k)
               H_coef = two/H_coef
               if (k .eq. domlo(3)) then
                  if (bc(3,1) .eq. EXT_DIR) then
                     H_coef = h_cc(i,j,k)
                  elseif (bc(3,1) .eq. SEEPAGE) then
                     H_coef = h_cc(i,j,k)
                  end if
               elseif (k .eq. domhi(3)) then
                  if (bc(3,2) .eq. EXT_DIR) then
                     H_coef = h_cc(i,j,k)
                  elseif (bc(3,2) .eq. SEEPAGE) then
                     H_coef = h_cc(i,j,k)
                  end if  
               end if
               dfz(i,j,k,1) = H_coef
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      end


c ::: -----------------------------------------------------------
c ::: This routine compute the phase velocity of component 1
c ::: -----------------------------------------------------------
      subroutine FORT_UPHASE(uphase,DIMS(uphase),
     &     vphase,DIMS(vphase),wphase,DIMS(wphase),
     &     umac,DIMS(umac),vmac,DIMS(vmac),wmac,DIMS(wmac),
     &     pc,DIMS(pc),lbd,DIMS(lbd),
     &     kpx,DIMS(kpx),kpy,DIMS(kpy),kpz,DIMS(kpz),
     &     sbc,dom_lo,dom_hi,dx,lo,hi)

      implicit none

#include "probdata.H"

      integer DIMDEC(uphase), DIMDEC(vphase), DIMDEC(wphase)
      integer DIMDEC(umac), DIMDEC(vmac), DIMDEC(wmac)
      integer DIMDEC(kpx), DIMDEC(kpy), DIMDEC(kpz)
      integer DIMDEC(pc)
      integer DIMDEC(lbd)
      integer lo(SDIM),hi(SDIM)
      integer dom_lo(SDIM),dom_hi(SDIM)
      integer sbc(SDIM,2)
      REAL_T  dx(SDIM)
      REAL_T  uphase(DIMV(uphase))
      REAL_T  vphase(DIMV(vphase))
      REAL_T  wphase(DIMV(wphase))
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  wmac(DIMV(wmac))
      REAL_T  kpx(DIMV(kpx)),kpy(DIMV(kpy)),kpz(DIMV(kpz))
      REAL_T  pc(DIMV(pc))
      REAL_T  lbd(DIMV(lbd),ncomps)

      integer i,j,k
      REAL_T  dpc
      REAL_T  ld(ncomps),le(ncomps)

!$omp parallel private(i,j,k,dpc,ld,le)
c     phase velocity in the x-direction
!$omp do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1

               if (model .eq. 2) then
                  dpc = (pc(i,j,k) - pc(i-1,j,k))/dx(1)
                  ld = lbd(i,j,k,:)
                  le = lbd(i-1,j,k,:)
                  ld(1) = half*(ld(1) + le(1))
                  ld(2) = half*(ld(2) + le(2))
                  
                  if (i .eq. dom_lo(1)) then
                     if (sbc(1,1) .eq. EXT_DIR .or.
     &                   sbc(1,1) .eq. SEEPAGE) then
                        dpc = half*dpc
                        ld = lbd(i-1,j,k,:)
                     endif
                  else if (i .eq. dom_hi(1)+1) then
                     if (sbc(1,2) .eq. EXT_DIR .or.
     &                   sbc(1,2) .eq. SEEPAGE) then
                        dpc = half*dpc
                        ld = lbd(i,j,k,:)
                     endif
                  end if
                  uphase(i,j,k) = ld(1)/(ld(1)+ld(2))*(umac(i,j,k) + kpx(i,j,k)*ld(2)*dpc)
               end if
            end do
         end do
      end do
!$omp end do nowait

c     phase velocity in the y-direction
!$omp do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)

               if (model .eq. 2) then
                  dpc = (pc(i,j,k) - pc(i,j-1,k))/dx(2)
                  ld = lbd(i,j,k,:)
                  le = lbd(i,j-1,k,:)
                  ld(1) = half*(ld(1) + le(1))
                  ld(2) = half*(ld(2) + le(2))
                  
                  if (j .eq. dom_lo(2)) then
                     if (sbc(2,1) .eq. EXT_DIR .or.
     &                   sbc(2,1) .eq. SEEPAGE) then
                        dpc = half*dpc
                        ld = lbd(i,j-1,k,:)
                     endif
                  else if (j .eq. dom_hi(2)+1) then
                     if (sbc(2,2) .eq. EXT_DIR .or.
     &                   sbc(2,2) .eq. SEEPAGE) then
                        dpc = half*dpc
                        ld = lbd(i,j,k,:)
                     endif
                  end if
                  vphase(i,j,k) = ld(1)/(ld(1)+ld(2))*(vmac(i,j,k) + kpy(i,j,k)*ld(2)*dpc)
               end if
            end do
         end do
      end do
!$omp end do

c     phase velocity in the z-direction
!$omp do
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               if (model .eq. 2) then
                  dpc = (pc(i,j,k) - pc(i,j,k-1))/dx(3)
                  ld = lbd(i,j,k,:)
                  le = lbd(i,j,k-1,:)
                  ld(1) = half*(ld(1) + le(1))
                  ld(2) = half*(ld(2) + le(2))
                  
                  if (k .eq. dom_lo(3)) then
                     if (sbc(3,1) .eq. EXT_DIR .or.
     &                   sbc(3,1) .eq. SEEPAGE) then
                        dpc = half*dpc
                        ld = lbd(i,j,k-1,:)
                     endif
                  else if (k .eq. dom_hi(3)+1) then
                     if (sbc(3,2) .eq. EXT_DIR .or.
     &                   sbc(3,2) .eq. SEEPAGE) then
                        dpc = half*dpc
                        ld = lbd(i,j,k,:)
                     endif
                  end if
                  wphase(i,j,k) = ld(1)/(ld(1)+ld(2))*(wmac(i,j,k) + kpz(i,j,k)*ld(2)*dpc)
               end if
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      end

c ::: ----------------------------------------------------------
c ::: compute divergence of velocity
c ::: ----------------------------------------------------------
      subroutine FORT_DIV_UMAC(s,DIMS(s),
     &     ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz),lo,hi,dx)
      
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(ux)
      integer DIMDEC(uy)
      integer DIMDEC(uz)
      integer DIMDEC(s)
      REAL_T  ux(DIMV(ux))
      REAL_T  uy(DIMV(uy))
      REAL_T  uz(DIMV(uz))
      REAL_T  s(DIMV(s))
      REAL_T  dx(SDIM)

      integer i,j,k

      do k = lo(3),hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k) = s(i,j,k) + (ux(i+1,j,k)-ux(i,j,k))/dx(1) +
     &              (uy(i,j+1,k)-uy(i,j,k))/dx(2) + 
     &              (uz(i,j,k+1)-uz(i,j,k))/dx(3) 
            end do
         end do
      end do

      end

c ::: ----------------------------------------------------------
c ::: Average the edge-centered velocities to a cell-centered quantity.
c ::: ----------------------------------------------------------
      subroutine FORT_AVG_UMAC(umac,DIMS(umac),vmac,DIMS(vmac),
     &                         wmac,DIMS(wmac),vel,DIMS(vel),lo,hi)

      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(wmac)
      integer DIMDEC(vel)
      integer lo(SDIM),hi(SDIM)
      REAL_T umac(DIMV(umac))
      REAL_T vmac(DIMV(vmac))
      REAL_T wmac(DIMV(wmac))
      REAL_T  vel(DIMV(vel),SDIM)

      integer i,j,k

!$omp parallel private(i,j,k)
!$omp do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vel(i,j,k,1) = half * (umac(i,j,k) + umac(i+1,j,k))
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vel(i,j,k,2) = half * (vmac(i,j,k) + vmac(i,j+1,k))
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vel(i,j,k,3) = half * (wmac(i,j,k) + wmac(i,j,k+1))
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      end

c ::: -------------------------------------------------------------
c ::: Copy the edge-centered velocities to a cell-centered quantity.
c ::: -------------------------------------------------------------
      subroutine FORT_CPY_UMAC(umac,DIMS(umac),vmac,DIMS(vmac),
     &                         wmac,DIMS(wmac),vel,DIMS(vel),lo,hi,shift)

      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(wmac)
      integer DIMDEC(vel)
      integer lo(SDIM),hi(SDIM)
      integer shift
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  wmac(DIMV(wmac))
      REAL_T  vel(DIMV(vel),SDIM)

      integer i,j,k

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vel(i,j,k,1) = umac(i+shift,j,k) 
            end do
        end do
      end do

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vel(i,j,k,2) = vmac(i,j+shift,k)
            end do
        end do
      end do

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vel(i,j,k,3) = wmac(i,j,k+shift)
            end do
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
      
      integer i,j,k,kk,l_cpl_type
      REAL_T  l_cpl_coef, l_sigma, l_sat_res
      REAL_T  z
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

      do k = ARG_L3(rho), ARG_H3(rho)
         do j = ARG_L2(rho), ARG_H2(rho)
            do i = ARG_L1(rho), ARG_H1(rho)
               rho(i,j,k,1) = zero
               do kk = 1,5
                  z = (dble(k)+half+qd(kk))*dx(3)
                  if (z .le. wt_loc) then
                     rho(i,j,k,1) = rho(i,j,k,1) + rhoval(1)*wt(kk)

                  else
                     l_cpl_type = nint(cplcoef(i,j,k,1))
                     l_cpl_coef = cplcoef(i,j,k,2)
                     l_sigma    = cplcoef(i,j,k,3)
                     l_sat_res  = cplcoef(i,j,k,4)

                     satlim(1) = l_sat_res*rhoval(1)
                     satlim(2) = (one - l_sat_res)*rhoval(2)

                     pc = (rhoval(1)-rhoval(2))*grav*(z-wt_loc)
                     call get_invcpl(stmp,pc,rhoval,phi,kappa,l_sigma,
     &                    ncomps,l_sat_res,l_cpl_type,l_cpl_coef)
                     
                     if (stmp(1) .lt. satlim(1)) then
                        stmp = satlim
                     end if

                     rho(i,j,k,1) = rho(i,j,k,1) + stmp(1)*wt(kk)
                     rho(i,j,k,2) = rho(i,j,k,2) + stmp(2)*wt(kk)

                  end if 
               end do 
            end do  
         end do
      end do 

      end

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
      integer i,j,k,n

#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5

c     z-faces
      if (face .eq. ZLO) then
         do k = ARG_L3(p),domlo(3)-2
            do n = 1,ncomp
               p(:,:,k,n) = p(:,:,domlo(3)-1,n)
            end do
         end do

      elseif (face .eq. ZHI) then
         do k = domhi(3)+2,ARG_H3(p)
            do n = 1,ncomp
               p(:,:,k,n) = p(:,:,domhi(3)+1,n)
            end do
         end do
      end if

c     y-faces
      if (face .eq. YLO) then
         do j = ARG_L2(p),domlo(2)-2
            do n = 1,ncomp
               p(:,j,:,n) = p(:,domlo(2)-1,:,n)
            end do
         end do

      elseif (face .eq. YHI) then
         do j = domhi(2)+2,ARG_H2(p)
            do n = 1,ncomp
               p(:,j,:,n) = p(:,domhi(2)+1,:,n)
            end do
         end do
      end if

      
c     x-faces
      if (face .eq. XLO) then
         do i = ARG_L1(p),domlo(1)-2
            do n = 1,ncomp
               p(i,:,:,n) = p(domlo(1)-1,:,:,n)
            end do
         end do

      elseif (face .eq. XHI) then
         do i = domhi(1)+2,ARG_H1(p)
            do n = 1,ncomp
               p(i,:,:,n) = p(domhi(1)+1,:,:,n)
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
