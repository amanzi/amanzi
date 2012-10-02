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

#define SDIM 2
c ::: ----------------------------------------------------------
c ::: Define the permeability, kappa.
c :::    Kappa is defined with 3 ghost cells on each side to 
c :::    acommodate the reaction module.  However, kpedge are
c :::    only defined for edges of the valid region.
c ::: ----------------------------------------------------------
      subroutine FORT_INITKAPPA (kappa,DIMS(kappa),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),
     &     lo,hi,dx,prob_hi,currlev,maxlev,
     &     mediumtype,mediumval_x,mediumval_y,nlayer,fratio)

      implicit none

      integer DIMDEC(kappa),DIMDEC(kpx),DIMDEC(kpy)
      integer lo(SDIM), hi(SDIM)
      integer currlev,maxlev
      integer mediumtype,nlayer
      integer fratio
      REAL_T  kappa(DIMV(kappa))
      REAL_T  kappax(DIMV(kpx)),kappay(DIMV(kpy))
      REAL_T  dx(SDIM),prob_hi(SDIM)
      REAL_T  mediumval_x(nlayer),mediumval_y(nlayer)

      integer i,j,k,ratio
      integer imin,imax,jmin,jmax
      REAL_T  x,y,z,vp,kt
      REAL_T  k1, k2
      REAL_T  dxfine(SDIM)

      ratio = fratio**(maxlev-currlev)
      dxfine(1) = dx(1)/dble(ratio)/prob_hi(1)
      dxfine(2) = dx(2)/dble(ratio)/prob_hi(2)

      imin = lo(1) -1
      jmin = lo(2) -1
      imax = hi(1) +1
      jmax = hi(2) +1

      if (mediumtype .eq. 5) then 
         call bl_abort('FORT_INITKAPPA: mediumtype cannot be 5')
      end if
      
c     dummy for 2D problem
      z = zero

      do j = ARG_L2(kappa), ARG_H2(kappa)
         y = (dble(j)+half)*dx(2) / prob_hi(2)
         do i = ARG_L1(kappa), ARG_H1(kappa)         
            x = (dble(i)+half)*dx(1) / prob_hi(1)
            call getkappa(k1,x,y,z,mediumtype,mediumval_x,nlayer)
            call getkappa(k2,x,y,z,mediumtype,mediumval_y,nlayer)
            kappa(i,j) = two*(one/k1+one/k2)
         end do
      end do

c ::: kpedge in x-dir
      do i = imin, imax+1
         x = (dble(i))*dx(1)/prob_hi(1)
         do j = jmin,jmax         
            kappax(i,j) = zero
            y = (dble(j))*dx(2)/prob_hi(2)
            do k = 1,ratio
               vp = y + (2*k-1)*dxfine(2)*half
               call getkappa(k1,x-dxfine(1)*half,vp,z,mediumtype,mediumval_x,nlayer)
               call getkappa(k2,x+dxfine(1)*half,vp,z,mediumtype,mediumval_x,nlayer)
               kt = one/k1 + one/k2
               kappax(i,j) = kappax(i,j) + two/kt
            end do
            kappax(i,j) = kappax(i,j)/dble(ratio)
         end do
      end do

c ::: kpedge in y-dir
      do j = jmin,jmax+1
         y = (dble(j))*dx(2)/prob_hi(2)
         do i = imin,imax        
            kappay(i,j) = zero
            x = (dble(i))*dx(1)/prob_hi(1)
            do k = 1,ratio
               vp = x + (2*k-1)*dxfine(1)*half
               call getkappa(k1,vp,y-dxfine(2)*half,z,mediumtype,mediumval_y,nlayer)
               call getkappa(k2,vp,y+dxfine(2)*half,z,mediumtype,mediumval_y,nlayer)
               kt = one/k1 + one/k2
               kappay(i,j) = kappay(i,j) + two/kt
            end do
            kappay(i,j) = kappay(i,j)/dble(ratio)
         end do
      end do
      
      end

c ::: ----------------------------------------------------------
c ::: Define the permeability, kappa.
c :::    Kappa is defined with 3 ghost cells on each side to 
c :::    acommodate the reaction module.  However, kpedge are
c :::    only defined for edges of the valid region.
c ::: ----------------------------------------------------------
      subroutine FORT_INITKAPPA2 (kappadata, DIMS(kappadata), 
     &     kappa,DIMS(kappa),
     &     kappax,DIMS(kpx),kappay,DIMS(kpy),
     &     lo,hi,currlev,maxlev,fratio)


      integer DIMDEC(kappadata)
      integer DIMDEC(kappa),DIMDEC(kpx),DIMDEC(kpy)
      integer lo(SDIM), hi(SDIM)
      integer currlev,maxlev
      integer fratio
      REAL_T  kappadata(DIMV(kappadata),SDIM)
      REAL_T  kappa(DIMV(kappa))
      REAL_T  kappax(DIMV(kpx)),kappay(DIMV(kpy))

      integer i,j,k,l,ratio,n
      integer fi,fj
      integer imin,imax,jmin,jmax
      REAL_T  k1, k2, kl1, kl2

      ratio = fratio**(maxlev-currlev)

      imin = lo(1) -1
      jmin = lo(2) -1
      imax = hi(1) +1
      jmax = hi(2) +1

      do j = lo(2)-3, hi(2)+3
         fj = ratio*j
         do i = lo(1)-3, hi(1)+3         
            fi = ratio*i
            kappa(i,j) = zero
            do k = 0,ratio-1
               do l = 0,ratio-1
                  do n = 1,SDIM
                     kappa(i,j) = kappa(i,j) + one/kappadata(fi+k,fj+l,n)
                  end do
               end do
            end do
            kappa(i,j) = dble(SDIM*ratio*ratio)/kappa(i,j)
         end do
      end do

#if 1
c    ::: kpedge in x-dir
      do i = imin, imax+1
         fi = ratio*i
         do j = jmin, jmax  
            fj = ratio*j
            kappax(i,j) = zero            
            do k = 0,ratio-1
               kl1 = 0
               kl2 = 0
               do l = 0,ratio-1
                  kl1 = kl1 + one/kappadata(fi+l,fj+k,1)
                  kl2 = kl2 + one/kappadata(fi-1-l,fj+k,1)
               end do
               kappax(i,j) = kappax(i,j) + two*ratio/(kl1 + kl2)
            enddo
            kappax(i,j) = kappax(i,j)/DBLE(ratio)
         end do
      end do
         
c     ::: kpedge in y-dir
      do j = jmin, jmax+1
         fj = ratio*j
         do i = imin, imax         
            fi = ratio*i
            kappay(i,j) = zero
            do k = 0,ratio-1   
               kl1 = 0
               kl2 = 0
               do l = 0,ratio-1
                  kl1 = kl1 + one/kappadata(fi+k,fj+l,2)
                  kl2 = kl2 + one/kappadata(fi+k,fj-1-l,2)
               end do
               kappay(i,j) = kappay(i,j) + two*ratio/(kl1 + kl2)
            end do
            kappay(i,j) = kappay(i,j)/DBLE(ratio)
         end do
      end do
#else
c    ::: kpedge in x-dir
      do i = imin, imax+1
         fi = ratio*i
         do j = jmin, jmax  
            fj = ratio*j
            kappax(i,j) = zero            
            k1 = 0
            k2 = 0
            do k = 0,ratio-1
               kl1 = 0
               kl2 = 0
               do l = 0,ratio-1
                  kl1 = kl1 + one/kappadata(fi+l,fj+k,1)
                  kl2 = kl2 + one/kappadata(fi-1-l,fj+k,1)
               end do
               k1 = k1 + dble(ratio)/kl1
               k2 = k2 + dble(ratio)/kl2
            end do
            k1 = k1/dble(ratio)
            k2 = k2/dble(ratio)          
            kappax(i,j) = two/(one/k1 + one/k2)
         end do
      end do
         
c     ::: kpedge in y-dir
      do j = jmin, jmax+1
         fj = ratio*j
         do i = imin, imax         
            fi = ratio*i
            kappay(i,j) = zero
            k1 = 0
            k2 = 0
            do k = 0,ratio-1   
               kl1 = 0
               kl2 = 0
               do l = 0,ratio-1
                  kl1 = kl1 + one/kappadata(fi+k,fj+l,2)
                  kl2 = kl2 + one/kappadata(fi+k,fj-1-l,2)
               end do
               k1 = k1 + dble(ratio)/kl1
               k2 = k2 + dble(ratio)/kl2
            end do
            k1 = k1/dble(ratio)
            k2 = k2/dble(ratio)          
            kappay(i,j) = two/(one/k1 + one/k2)
         end do
      end do
#endif
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

      integer i,j
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

      do j = ARG_L2(phi), ARG_H2(phi)
         x(2) = (dble(j)+half)*dx(2) 
         do i = ARG_L1(phi), ARG_H1(phi)         
            x(1) = (dble(i)+half)*dx(1) 
            call getkappa(phi(i,j),x,xmax,mediumtype,mediumval,nlayer)
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

      integer   i,j,k,l,fi,fj

      do j = lo(2), hi(2)
         fj = ratio*j
         do i = lo(1), hi(1)         
            fi = ratio*i
            phi(i,j) = zero
            do k = 0,ratio-1
               do l = 0,ratio-1
                  phi(i,j) = phi(i,j) + phidata(fi+k,fj+l)
               end do
            end do
            phi(i,j) = phi(i,j)/dble(ratio*ratio)

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

      integer   i,j,k,l,ratio,fi,fj,n

      if (fratio .le. 0) then
         ratio = one
      else
         ratio = fratio**(maxlev-currlev)
      end if
      do n = 1,nkr
         do j = ARG_L2(kr), ARG_H2(kr)
            fj = ratio*j
            do i = ARG_L1(kr), ARG_H1(kr)
               fi = ratio*i
               kr(i,j,n) = zero
               do k = 0,ratio-1
                  do l = 0,ratio-1
                     kr(i,j,n) = kr(i,j,n) + krdata(fi+k,fj+l,n)
                  end do
               end do
               kr(i,j,n) = kr(i,j,n)/dble(ratio*ratio)
            end do
         end do
      end do
      end

c ::: -----------------------------------------------------------
c ::: Construct the cell-centered coefficient of grad p to be 
c ::: averaged onto edges and used in the MAC solve.
c ::: ----------------------------------------------------------- 
      subroutine FORT_MK_MACCOEF (coefx,DIMS(coefx),
     &     coefy,DIMS(coefy),
     &     kpx,DIMS(kpx),kpy,DIMS(kpy),
     &     lbd,DIMS(lbd),
     &     lo,hi,domlo,domhi,bc)

      implicit none

#include "probdata.H" 

      integer DIMDEC(coefx)
      integer DIMDEC(coefy)
      integer DIMDEC(kpx)
      integer DIMDEC(kpy)
      integer DIMDEC(lbd)
      integer lo(SDIM), hi(SDIM)
      integer domlo(SDIM),domhi(SDIM)
      integer bc(SDIM,2)
      REAL_T  coefx(DIMV(coefx))
      REAL_T  coefy(DIMV(coefy))
      REAL_T  kpx(DIMV(kpx))
      REAL_T  kpy(DIMV(kpy))
      REAL_T  lbd(DIMV(lbd),ncomps)

      integer i, j
      REAL_T  ld(ncomps)

c ::: kpedge in x-dir
      do i = lo(1),hi(1)+1
         do j = lo(2),hi(2)  

            ld = lbd(i,j,:)
            coefx(i,j) = ld(1) + ld(2)
            ld = lbd(i-1,j,:)
            coefx(i,j) = half*(coefx(i,j) + ld(1) + ld(2))
            
            if (i .eq. domlo(1)) then
               if (bc(1,1) .eq. EXT_DIR) then
                  ld = lbd(i-1,j,:) 
               elseif (bc(1,1) .eq. SEEPAGE) then
                  ld = lbd(i,j,:)
               endif
               coefx(i,j) = ld(1) + ld(2)                    
            elseif (i .eq. domhi(1)+1) then
               if (bc(1,2) .eq. EXT_DIR) then
                  ld = lbd(i,j,:)
               elseif (bc(1,2) .eq. SEEPAGE) then
                  ld = lbd(i-1,j,:)
               endif
               coefx(i,j) = ld(1) + ld(2)
            end if
            coefx(i,j) = kpx(i,j)*coefx(i,j)
         end do
      end do

c ::: kpedge in y-dir
      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1) 
            ld = lbd(i,j,:)
            coefy(i,j) = ld(1) + ld(2)
            ld = lbd(i,j-1,:)
            coefy(i,j) = half*(coefy(i,j) + ld(1) + ld(2))

            if (j .eq. domlo(2)) then
               if (bc(2,1) .eq. EXT_DIR) then
                  ld = lbd(i,j-1,:)
               elseif (bc(2,1) .eq. SEEPAGE) then
                  ld = lbd(i,j,:)
               endif
               coefy(i,j) = ld(1) + ld(2)                     
            elseif (j .eq. domhi(2)+1) then
               if (bc(2,2) .eq. EXT_DIR) then
                  ld = lbd(i,j,:)
               elseif (bc(2,2) .eq. SEEPAGE) then
                  ld = lbd(i,j-1,:)
               endif
               coefy(i,j) = ld(1) + ld(2)     
            end if
            coefy(i,j) = kpy(i,j)*coefy(i,j)
         end do
      end do
      
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

      integer i,j,n
      REAL_T  y,ydiv(nkpval+1)

      ydiv(nkpval+1) = dble(domhi(2))
      do n = 1,nkpval
         ydiv(n) = dble(domlo(2)) + 
     &        dble(domhi(2)-domlo(2))*dble(n-1)/dble(nkpval)
      end do 
      
      do i=ARG_L1(kappa),ARG_H1(kappa)
         do j=ARG_L2(kappa),ARG_H2(kappa)
            y = dble(j) + half
            do n = 1,nkpval
               if (y .ge. ydiv(n) .and. y .lt. ydiv(n+1)) then
                  kappa(i,j) = kappa(i,j)*kpval(n)/kpval(1)
               end if
            end do

            if (y .lt. ydiv(1)) then
               kappa(i,j) = kappa(i,j)
            else if (y .ge. ydiv(nkpval+1)) then
               kappa(i,j) = kappa(i,j)*kpval(nkpval)/kpval(1)
            end if
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

      integer i,j,n
      REAL_T  phirnd
      REAL_T  y,ydiv(nkpval+1)

      call blutilinitrand(iuc)

      ydiv(nkpval+1) = dble(domhi(2))
      do n = 1,nkpval
         ydiv(n) = dble(domlo(2)) + 
     &        dble(domhi(2)-domlo(2))*dble(n-1)/dble(nkpval)
      end do 
      
      do i=ARG_L1(kappa),ARG_H1(kappa)
         do j=ARG_L2(kappa),ARG_H2(kappa)
            y = dble(j) + half
            call blutilrand(phirnd)
            do n = 1,nkpval
               if (y .ge. ydiv(n) .and. y .lt. ydiv(n+1)) then
                  kappa(i,j) = kpval(n)*(1+2*dkp*(0.5-phirnd))
               end if
            end do

            if (y .lt. ydiv(1)) then
               kappa(i,j) = kpval(1)*(1+2*dkp*(0.5-phirnd))
            else if (y .ge. ydiv(nkpval+1)) then
               kappa(i,j) = kpval(nkpval)*(1+2*dkp*(0.5-phirnd))
            end if
         end do
      end do

      end
