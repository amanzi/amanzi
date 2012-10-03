#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "POROUS_F.H"
#include "ArrayLim.H"

#define SDIM 3

c ::: ----------------------------------------------------------
c ::: Define the permeability, kappa.
c :::    Kappa is defined with 3 ghost cells on each side to 
c :::    acommodate the reaction module.  However, kpedge are
c :::    only defined for edges of the valid region.
c ::: ---------------------------------------------------------
      subroutine FORT_INITKAPPA (kappa,DIMS(kappa),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),kappaz,DIMS(kpz),
     &     lo,hi,dx,prob_hi,currlev,maxlev,
     &     mediumtype,mediumval_x,mediumval_y,mediumval_z,nlayer,fratio)
 
      implicit none

      integer DIMDEC(kappa),DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer lo(SDIM), hi(SDIM)
      integer currlev,maxlev
      integer mediumtype,nlayer
      integer fratio
      REAL_T  kappa(DIMV(kappa))
      REAL_T  kappax(DIMV(kpx)),kappay(DIMV(kpy)),kappaz(DIMV(kpz))
      REAL_T  dx(SDIM),prob_hi(SDIM)
      REAL_T  mediumval_x(nlayer), mediumval_y(nlayer)
      REAL_T  mediumval_z(nlayer)
      integer i, j, k, m, n, ratio, ratio2
      integer imin, imax, jmin, jmax, kmin, kmax
      REAL_T  x, y, z, vp, wp, kt
      REAL_T  dxfine(SDIM), xprobhi1, xprobhi2, xprobhi3, xratio2
      REAL_T  k1, k2, k3

      ratio    = fratio**(maxlev-currlev)
      ratio2   = ratio*ratio
      xratio2  = one / dble(ratio2)
      xprobhi1 = one / prob_hi(1)
      xprobhi2 = one / prob_hi(2)
      xprobhi3 = one / prob_hi(3)

      do i = 1,SDIM
         dxfine(i) = dx(i)/(ratio*prob_hi(i))
      end do

      imin = lo(1) - 1
      jmin = lo(2) - 1
      kmin = lo(3) - 1
      imax = hi(1) + 1
      jmax = hi(2) + 1
      kmax = hi(3) + 1

      if (mediumtype .eq. 5) then
         call bl_abort('FORT_INITKAPPA: mediumtype cannot be 5')
      end if

!$omp parallel private(i,j,k,m,n,x,y,z,k1,k2,k3,vp,wp,kt)
!$omp do
      do k = ARG_L3(kappa), ARG_H3(kappa)
         z = (dble(k)+half)*dx(3)*xprobhi3
         do j = ARG_L2(kappa), ARG_H2(kappa)
            y = (dble(j)+half)*dx(2)*xprobhi2
            do i = ARG_L1(kappa), ARG_H1(kappa)
               x = (dble(i)+half)*dx(1)*xprobhi1
               call getkappa(k1,x,z,y,mediumtype,mediumval_x,nlayer)
               call getkappa(k2,x,z,y,mediumtype,mediumval_y,nlayer)
               call getkappa(k3,x,z,y,mediumtype,mediumval_z,nlayer)
               kappa(i,j,k) = 3*(one/k1+one/k2+one/k3)
            end do
         end do
      end do
!$omp end do nowait

c ::: kpedge in x-dir
!$omp do
      do k = kmin,kmax
         z = (dble(k))*dx(3)*xprobhi3
         do j = jmin,jmax        
            y = (dble(j))*dx(2)*xprobhi2
            do i = imin,imax+1
               x = (dble(i))*dx(1)*xprobhi1
               kappax(i,j,k) = zero
               do m = 1,ratio
                  vp = y + (2*m-1)*dxfine(2)*half
                  do n = 1,ratio
                     wp = z + (2*n-1)*dxfine(3)*half
                     call getkappa(k1,x-dxfine(1)*half,wp,vp,mediumtype,
     &                    mediumval_x,nlayer)
                     call getkappa(k2,x+dxfine(1)*half,wp,vp,mediumtype,
     &                    mediumval_x,nlayer)
                     kt = one/k1 + one/k2
                     kappax(i,j,k) = kappax(i,j,k) + two/kt
                  end do
               end do
               kappax(i,j,k) = kappax(i,j,k)*xratio2
            end do
         end do
      end do
!$omp end do nowait

c ::: kpedge in y-dir
!$omp do
      do k = kmin,kmax
         z = (dble(k))*dx(3)*xprobhi3
         do j = jmin,jmax+1 
            y = (dble(j))*dx(2)*xprobhi2
            do i = imin,imax         
               x = (dble(i))*dx(1)*xprobhi1
               kappay(i,j,k) = zero
               do m = 1,ratio
                  vp = x + (2*m-1)*dxfine(1)*half
                  do n = 1,ratio
                     wp = z + (2*n-1)*dxfine(3)*half
                     call getkappa(k1,vp,wp,y-dxfine(2)*half,mediumtype,
     &                    mediumval_y,nlayer)
                     call getkappa(k2,vp,wp,y+dxfine(2)*half,mediumtype,
     &                    mediumval_y,nlayer)
                     kt = one/k1 + one/k2
                     kappay(i,j,k) = kappay(i,j,k) + two/kt
                  end do
               end do
               kappay(i,j,k) = kappay(i,j,k)*xratio2
            end do
         end do
      end do
!$omp end do nowait

c ::: kpedge in z-dir
!$omp do
      do k = kmin,kmax+1
         z = (dble(k))*dx(3)*xprobhi3
         do j = jmin,jmax 
            y = (dble(j))*dx(2)*xprobhi2
            do i = imin,imax       
               x = (dble(i))*dx(1)*xprobhi1
               kappaz(i,j,k) = zero
               do m = 1,ratio
                  vp = x + (2*m-1)*dxfine(1)*half
                  do n = 1,ratio
                     wp = y + (2*n-1)*dxfine(2)*half 
                     call getkappa(k1,vp,z-dxfine(3)*half,wp,mediumtype,
     &                    mediumval_z,nlayer)
                     call getkappa(k2,vp,z+dxfine(3)*half,wp,mediumtype,
     &                    mediumval_z,nlayer)
                     kt = one/k1 + one/k2
                     kappaz(i,j,k) = kappaz(i,j,k) + two/kt
                  end do
               end do
               kappaz(i,j,k) = kappaz(i,j,k)*xratio2
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      end

      subroutine FORT_INITKAPPA2 (kappadata, DIMS(kappadata), kappa,DIMS(kappa),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),kappaz,DIMS(kpz),
     &     lo,hi,currlev,maxlev,fratio)

      integer DIMDEC(kappadata)
      integer DIMDEC(kappa),DIMDEC(kpx),DIMDEC(kpy),DIMDEC(kpz)
      integer lo(SDIM), hi(SDIM)
      integer currlev,maxlev
      integer fratio
      REAL_T  kappadata(DIMV(kappadata),SDIM)
      REAL_T  kappa(DIMV(kappa))
      REAL_T  kappax(DIMV(kpx)),kappay(DIMV(kpy)),kappaz(DIMV(kpz))
      integer i, j, k, m, n, l, ratio, ratio2,dir
      integer imin, imax, jmin, jmax, kmin, kmax
      integer fi,fj,fk
      REAL_T  k1, k2, kl1,kl2

      ratio = fratio**(maxlev-currlev)
      ratio2 = ratio*ratio

      imin = lo(1) - 1
      jmin = lo(2) - 1
      kmin = lo(3) - 1
      imax = hi(1) + 1
      jmax = hi(2) + 1
      kmax = hi(3) + 1

      do k = lo(3)-3,hi(3)+3
         fk = ratio*k
         do j = lo(2)-3,hi(2)+3
            fj = ratio*j
            do i = lo(1)-3,hi(1)+3        
               fi = ratio*i
               kappa(i,j,k) = zero
               do l = 0,ratio-1
                  do m = 0,ratio-1
                     do n = 0,ratio-1
                        do dir = 1,SDIM
                           kappa(i,j,k) = kappa(i,j,k) + 
     &                          one/kappadata(fi+l,fj+m,fk+n,dir)
                        end do
                     end do
                  end do
               end do
               kappa(i,j,k) = dble(SDIM*ratio*ratio*ratio)/kappa(i,j,k)
            end do
         end do
      end do

c     ::: kpedge in x-dir
      do i = imin,imax+1
         fi = ratio*i
         do k = kmin,kmax
            fk = ratio*k
            do j = jmin,jmax         
               kappax(i,j,k) = 0
               fj = ratio*j
               k1 = 0
               k2 = 0 
               do n = 0,ratio-1
                  do m = 0,ratio-1
                     kl1 = 0
                     kl2 = 0
                     do l = 0,ratio-1                    
                        kl1 = kl1 + one/kappadata(fi+l,fj+m,fk+n,1)
                        kl2 = kl2 + one/kappadata(fi-1-l,fj+m,fk+n,1)
                     end do
                     k1 = k1 + dble(ratio)/kl1
                     k2 = k2 + dble(ratio)/kl2
                  end do
               end do
               k1 = k1/dble(ratio2)
               k2 = k2/dble(ratio2)          
               kappax(i,j,k) = two/(one/k1 + one/k2)
            end do
         end do
      end do

c     ::: kpedge in y-dir
      do j = jmin,jmax+1
         fj = ratio*j
         do k = kmin,kmax
            fk = ratio*k
            do i = imin,imax         
               kappay(i,j,k) = 0
               fi = ratio*i
               k1 = 0
               k2 = 0
               do l = 0,ratio-1
                  do m = 0,ratio-1
                     kl1 = 0
                     kl2 = 0
                     do n = 0,ratio-1    
                        kl1 = kl1 + one/kappadata(fi+m,fj+l,fk+n,2)
                        kl2 = kl2 + one/kappadata(fi+m,fj-1-l,fk+n,2)   
                     end do
                     k1 = k1 + dble(ratio)/kl1
                     k2 = k2 + dble(ratio)/kl2
                  end do
               end do
               k1 = k1/dble(ratio2)
               k2 = k2/dble(ratio2)          
               kappay(i,j,k) = two/(one/k1 + one/k2)    
            end do
         end do
      end do

c     ::: kpedge in z-dir
      do k = kmin,kmax+1
         fk = ratio*k
         do j = jmin,jmax
            fj = ratio*j
            do i = imin,imax         
               kappaz(i,j,k) = 0
               fi = ratio*i
               k1 = 0
               k2 = 0
               do l = 0,ratio-1
                  do m = 0,ratio-1
                     kl1 = 0
                     kl2 = 0
                     do n = 0,ratio-1                    
                        kl1 = kl1 + one/kappadata(fi+m,fj+n,fk+l,3)
                        kl2 = kl2 + one/kappadata(fi+m,fj+n,fk-1-l,3)   
                     end do
                     k1 = k1 + dble(ratio)/kl1
                     k2 = k2 + dble(ratio)/kl2
                  end do
               end do
               k1 = k1/dble(ratio2)
               k2 = k2/dble(ratio2)          
               kappaz(i,j,k) = two/(one/k1 + one/k2) 
            end do
         end do
      end do

      end 

c ::: -----------------------------------------------------------
c ::: Define variable porosity, phi
c ::: ----------------------------------------------------------- 
      subroutine FORT_INITPHI (phi,DIMS(phi),
     &     domhi,dx,mediumtype,mediumval,nlayer)

      integer DIMDEC(phi)
      integer domhi(SDIM)
      integer mediumtype,nlayer
      REAL_T  phi(DIMV(phi))
      REAL_T  dx(SDIM)
      REAL_T  mediumval(nlayer)

      integer i,j,k
      REAL_T  x(3),xmax(3)

      if (mediumtype .eq. 5) then 
         call bl_abort('FORT_INITKAPPA: mediumtype cannot be 5')
      end if

c     dummy for 2D problem
      x    = zero
      xmax = zero
      do i = 1,SDIM
         xmax(i) = dble(domhi(i))*dx(i)
      end do

      do k = ARG_L3(phi),ARG_H3(phi)
         x(3) = (dble(k)+half)*dx(3) 
         do j = ARG_L2(phi), ARG_H2(phi)
            x(2) = (dble(j)+half)*dx(2) 
            do i = ARG_L1(phi), ARG_H1(phi)         
               x(1) = (dble(i)+half)*dx(1)
               !
               ! TODO - this doesn't have right number of arguments!!!
               !
               call getkappa(phi(i,j,k),x,xmax,mediumtype,mediumval,nlayer)
            end do
         end do
      end do

      end

c ::: -----------------------------------------------------------
c ::: Define variable porosity, phi
c ::: ----------------------------------------------------------- 
      subroutine FORT_INITPHI2 (phidata,DIMS(phidata),
     &     phi,DIMS(phi),lo,hi,ratio)

      integer DIMDEC(phi)
      integer DIMDEC(phidata)
      integer lo(SDIM), hi(SDIM)
      integer ratio
      REAL_T  phi(DIMV(phi))
      REAL_T  phidata(DIMV(phidata))

      integer i, j, k,l,m,n,fi,fj,fk

      do k = lo(3),hi(3)
         fk = ratio*k
         do j = lo(2),hi(2)
            fj = ratio*j
            do i = lo(1),hi(1)        
               fi = ratio*i
               phi(i,j,k) = zero
               do l = 0,ratio-1
                  do m = 0,ratio-1
                     do n = 0,ratio-1
                        phi(i,j,k) = phi(i,j,k) + phidata(fi+l,fj+m,fk+n)
                     end do
                  end do
               end do
               phi(i,j,k) = phi(i,j,k)/dble(ratio*ratio*ratio)
            end do
         end do
      end do

      end


c ::: -----------------------------------------------------------
c ::: Average from finest level to current level
c ::: ----------------------------------------------------------- 
      subroutine FORT_INITKR (krdata,DIMS(krdata),
     &     kr,DIMS(kr),nkr,currlev,maxlev,fratio)

      implicit none 

      integer DIMDEC(kr)
      integer DIMDEC(krdata)
      integer nkr,currlev,maxlev
      integer fratio
      REAL_T  kr(DIMV(kr),nkr)
      REAL_T  krdata(DIMV(krdata),nkr)

      integer i, j, k,l,m,n,fi,fj,fk,ratio,nc

      if (fratio .le. 0) then
         ratio = one
      else
         ratio = fratio**(maxlev-currlev)
      end if

      do nc = 1,nkr
         do k = ARG_L3(kr), ARG_H3(kr)
            fk = ratio*k
            do j = ARG_L2(kr), ARG_H2(kr)
               fj = ratio*j
               do i = ARG_L1(kr), ARG_H1(kr)        
                  fi = ratio*i
                  kr(i,j,k,nc) = zero
                  do l = 0,ratio-1
                     do m = 0,ratio-1
                        do n = 0,ratio-1
                           kr(i,j,k,nc) = kr(i,j,k,nc) + 
     &                          krdata(fi+l,fj+m,fk+n,nc)
                        end do
                     end do
                  end do
                  kr(i,j,k,nc) = kr(i,j,k,nc)/dble(ratio*ratio*ratio)
               end do
            end do
         end do
      end do

      end

      subroutine FORT_MK_MACCOEF (coefx,DIMS(coefx),
     &     coefy,DIMS(coefy),coefz,DIMS(coefz),
     &     kpx,DIMS(kpx),kpy,DIMS(kpy),kpz,DIMS(kpz),
     &     lbd,DIMS(lbd),lo,hi,domlo,domhi,bc)
c 
c     ----------------------------------------------------------
c	Construct the cell-centered coefficient of grad p to be 
c        averaged onto edges and used in the MAC solve.
c     ----------------------------------------------------------
c 
      implicit none

#include "probdata.H" 

c ::: arguments

      integer DIMDEC(coefx)
      integer DIMDEC(coefy)
      integer DIMDEC(coefz)
      integer DIMDEC(kpx)
      integer DIMDEC(kpy)
      integer DIMDEC(kpz)
      integer DIMDEC(lbd)
      integer lo(SDIM), hi(SDIM)
      integer domlo(SDIM),domhi(SDIM)
      integer bc(SDIM,2)
      REAL_T  coefx (DIMV(coefx))
      REAL_T  coefy (DIMV(coefy))
      REAL_T  coefz (DIMV(coefz))
      REAL_T  kpx (DIMV(kpx))
      REAL_T  kpy (DIMV(kpy))
      REAL_T  kpz (DIMV(kpz))
      REAL_T  lbd (DIMV(lbd),ncomps)
      
c ::: internal variable
      integer i, j, k
      REAL_T  ld(ncomps)

!$omp parallel private(i,j,k,ld)
c ::: kpedge in x-dir
!$omp do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)  
            do i = lo(1),hi(1)+1

               ld = lbd(i,j,k,:)
               coefx(i,j,k) = ld(1) + ld(2)
               ld = lbd(i-1,j,k,:)
               coefx(i,j,k) = half*(coefx(i,j,k) + ld(1) + ld(2))

               if (i .eq. domlo(1)) then
                  if (bc(1,1) .eq. EXT_DIR) then
                     ld = lbd(i-1,j,k,:)
                  elseif (bc(1,1) .eq. SEEPAGE) then
                     ld = lbd(i,j,k,:)
                  end if
                  coefx(i,j,k) = ld(1) + ld(2)
               elseif (i .eq. domhi(1)) then
                  if (bc(1,2) .eq. EXT_DIR) then
                     ld = lbd(i,j,k,:)
                  elseif (bc(1,2) .eq. SEEPAGE) then
                     ld = lbd(i-1,j,k,:)
                  endif
                  coefx(i,j,k) = ld(1) + ld(2)  
               end if
               coefx(i,j,k) = kpx(i,j,k)*coefx(i,j,k)

            end do
         end do
      end do
!$omp end do nowait

c ::: kpedge in y-dir
!$omp do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)  

               ld = lbd(i,j,k,:)
               coefy(i,j,k) = ld(1) + ld(2) 
               ld = lbd(i,j-1,k,:)
               coefy(i,j,k) = half*(coefy(i,j,k) + ld(1) + ld(2))

               if (j .eq. domlo(2)) then
                  if (bc(2,1) .eq. EXT_DIR) then
                     ld = lbd(i,j-1,k,:)
                  elseif (bc(2,1) .eq. SEEPAGE) then
                     ld = lbd(i,j,k,:)
                  endif
                  coefy(i,j,k) = ld(1) + ld(2)                     
               elseif (j .eq. domhi(2)) then
                  if (bc(2,2) .eq. EXT_DIR) then
                     ld = lbd(i,j,k,:)
                  elseif (bc(2,2) .eq. SEEPAGE) then
                     ld = lbd(i,j-1,k,:)
                  endif
                  coefy(i,j,k) = ld(1) + ld(2)
               end if
               coefy(i,j,k) = kpy(i,j,k)*coefy(i,j,k)
            end do
         end do
      end do
!$omp end do nowait

c ::: kpedge in z-dir
!$omp do
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1) 

               ld = lbd(i,j,k,:)
               coefz(i,j,k) = ld(1) + ld(2)
               ld = lbd(i,j,k-1,:)
               coefz(i,j,k) = half*(coefz(i,j,k) + ld(1) + ld(2))

               if (k .eq. domlo(3)) then
                  if (bc(3,1) .eq. EXT_DIR) then
                     ld = lbd(i,j,k,:)
                  elseif (bc(3,1) .eq. SEEPAGE) then
                     ld = lbd(i,j,k-1,:)
                  endif
                  coefz(i,j,k) = ld(1) + ld(2)                    
               elseif (k .eq. domhi(3)) then
                  if (bc(3,2) .eq. EXT_DIR) then
                     ld = lbd(i,j,k,:)
                  elseif (bc(3,2) .eq. SEEPAGE) then
                     ld = lbd(i,j,k-1,:)
                  endif
                  coefz(i,j,k) = ld(1) + ld(2)  
               end if
               coefz(i,j,k) = kpz(i,j,k)*coefz(i,j,k)
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      end

c ::: ----------------------------------------------------------------------
c ::: This routine scale the permeability according to some layer structures
c ::: ----------------------------------------------------------------------
      subroutine FORT_RESCALE(kappa,DIMS(kappa),kpval,nkpval,domlo,domhi)

      integer DIMDEC(kappa)
      integer nkpval
      integer domlo(SDIM),domhi(SDIM)
      REAL_T  kappa(DIMV(kappa))
      REAL_T  kpval(nkpval)

      integer i,j,k,n
      REAL_T  y,ydiv(nkpval+1)

      ydiv(nkpval+1) = dble(domhi(3))
      do n = 1,nkpval
         ydiv(n) = dble(domlo(3)) + 
     &        dble(domhi(3)-domlo(3))*dble(n-1)/dble(nkpval)
      end do 
      
      do i=ARG_L1(kappa),ARG_H1(kappa)
         do j=ARG_L2(kappa),ARG_H2(kappa)
            do k = ARG_L3(kappa),ARG_H3(kappa)
               y = dble(k) + half
               do n = 1,nkpval
                  if (y .ge. ydiv(n) .and. y .lt. ydiv(n+1)) then
                     kappa(i,j,k) = kappa(i,j,k)*kpval(n)/kpval(1)
                  end if
               end do
               
               if (y .lt. ydiv(1)) then
                  kappa(i,j,k) = kappa(i,j,k)
               else if (y .ge. ydiv(nkpval+1)) then
                  kappa(i,j,k) = kappa(i,j,k)*kpval(nkpval)/kpval(1)
               end if
            end do
         end do
      end do

      end

c ::: --------------------------------------------------------------------------
c ::: This routine create random permeability according to some layer structures
c ::: --------------------------------------------------------------------------
      subroutine FORT_PHIRAND(kappa,DIMS(kappa),kpval,nkpval,dkp,
     &     domlo,domhi,iuc)

      integer DIMDEC(kappa)
      integer nkpval,iuc
      integer domlo(SDIM),domhi(SDIM)
      REAL_T  dkp
      REAL_T  kappa(DIMV(kappa))
      REAL_T  kpval(nkpval)

      integer i,j,k,n
      REAL_T  phirnd
      REAL_T  y,ydiv(nkpval+1)

      call blutilinitrand(iuc)

      ydiv(nkpval+1) = dble(domhi(3))
      do n = 1,nkpval
         ydiv(n) = dble(domlo(3)) + 
     &        dble(domhi(3)-domlo(3))*dble(n-1)/dble(nkpval)
      end do 
      
      do i=ARG_L1(kappa),ARG_H1(kappa)
         do j=ARG_L2(kappa),ARG_H2(kappa)
            do k = ARG_L3(kappa),ARG_H3(kappa)
               y = dble(k) + half
               call blutilrand(phirnd)
               do n = 1,nkpval
                  if (y .ge. ydiv(n) .and. y .lt. ydiv(n+1)) then
                     kappa(i,j,k) = kpval(n)*(1+2*dkp*(0.5-phirnd))
                  end if
               end do

               if (y .lt. ydiv(1)) then
                  kappa(i,j,k) = kpval(1)*(1+2*dkp*(0.5-phirnd))
               else if (y .ge. ydiv(nkpval+1)) then
                  kappa(i,j,k) = kpval(nkpval)*(1+2*dkp*(0.5-phirnd))
               end if
            end do
         end do
      end do

      end
