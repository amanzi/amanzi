#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "CHEMDRIVER_F.H" 
#include "ArrayLim.H"

#define SDIM BL_SPACEDIM

c :: ----------------------------------------------------------
c :: Initialize Chem Modules
c :: ----------------------------------------------------------
      subroutine FORT_INIT_CHEM_TOUGH

      call INITIALIZE_CHEM_MODULES

      end

c :: ----------------------------------------------------------
c :: Chemistry based on TOUGH codes.
c :: ----------------------------------------------------------
      subroutine FORT_CHEM_TOUGH(s,DIMS(s),press,DIMS(press),
     &     cr,DIMS(cr),fc,DIMS(fc),
     &     nc,lo,hi,dt,nd,idx_cr,mweight,
     &     kappa,DIMS(kappa),phi,DIMS(phi))

#include "probdata.H"

      integer lo(SDIM),hi(SDIM), nc, nd
      integer DIMDEC(s)
      integer DIMDEC(press)
      integer DIMDEC(kappa)
      integer DIMDEC(phi)
      integer DIMDEC(cr)
      integer DIMDEC(fc)
      integer idx_cr(4)
      REAL_T  s(DIMV(s),nc)
      REAL_T  press(DIMV(press))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  phi(DIMV(phi))
      REAL_T  cr(DIMV(cr),nd)
      REAL_T  fc(DIMV(fc))
      REAL_T  dt
      REAL_T  mweight(nc)
      REAL_T  fctotal

      integer i,j,n
      integer(8) ielem
      integer(8) dim_pth, dim_ut, dim_c, dim_pre, dim_amin, dim_rad, dim_pfug
      REAL_T  pth(20)
      REAL_T  ut(nc-1), c(idx_cr(1)), pre(idx_cr(2))
      REAL_T  amin(idx_cr(3)), rad(idx_cr(3))
      REAL_T  pfug(idx_cr(4))
      REAL_T  deltex,timetot

      dim_pth  = 20
      dim_ut   = nc-1
      dim_c    = idx_cr(1)
      dim_pre  = idx_cr(2)
      dim_amin = idx_cr(3)
      dim_rad  = idx_cr(3)
      dim_pfug = idx_cr(4)
      deltex   = dt
      timetot  = 0
      pth      = 0

      do i = lo(1),hi(1)
         do j = lo(2),hi(2)
       
            ielem  = (i-1)*j
c           pressure in Pa
            pth(1) = 1.02d5 + press(i,j)*1.0d5
c           temperature in C, set to 25.	
            pth(2) = 25.0
c           saturation of liquid, set to 1.
            pth(3) = 1.0
c	    rock bulk density in kg/m3, set to fix value 2530
            pth(4) = 2530.0
c           permeability in m^2
            pth(5) = kappa(i,j)*9.869233d-8
            pth(6) = kappa(i,j)*9.869233d-8
            pth(7) = kappa(i,j)*9.869233d-8
            pth(8) = phi(i,j)
c           density of water in kg/L
            pth(9) = 0.99716

c           conc in mol/L
            do n = 1,nc-1
               ut(n) = max(zero,s(i,j,n)/mweight(n))
            end do

c           mineral abundance in mol/dm3 = mol/L
            do n = 1,dim_pre
c              pre(n)  = cr(i,j,n+dim_c)
               pre(n) = s(i,j,nc+n-1)
            end do
            
            do n = 1,dim_c
               c(n) = cr(i,j,n)
            end do

c           mineral abundance in mol/dm3 = mol/L
            do n = 1,dim_pre
c               pre(n)  = cr(i,j,n+dim_c)
               pre(n) = s(i,j,nc+n-1)
            end do
            do n = 1,dim_amin
               amin(n) = cr(i,j,n+dim_c+dim_pre)
            end do
            do n = 1,dim_rad
		rad(n) = cr(i,j,n+dim_c+dim_pre+dim_amin)	
            end do
            do n = 1,dim_pfug
               pfug(n) = cr(i,j,n+dim_c+dim_pre+dim_amin+dim_rad)
            end do

            call FORT_CHEM_2ND(ielem,
     &           pth, dim_pth, ut, dim_ut, c, dim_c,
     &           pre, dim_pre, amin, dim_amin, rad, dim_rad, 
     &           pfug, dim_pfug, nc, deltex, fc(i,j))

            s(i,j,1) = one
            do n = 2,nc-1
               s(i,j,n) = (ut(n))*mweight(n)
               s(i,j,1) = s(i,j,1) - s(i,j,n)/rhoval(n)
            end do
            s(i,j,1) = s(i,j,1)*rhoval(1)
            do n = 1,dim_c
c               cr(i,j,n)=max(zero,c(n))
               cr(i,j,n) = c(n)
            end do
            
            do n = 1,dim_pre 
               cr(i,j,n+dim_c) = pre(n)
c              s(i,j,nc+n-1)   = pre(n)*mweight(nc)
c               s(i,j,nc+n-1)   = s(i,j,nc+n-1) + pre(n)
               s(i,j,nc+n-1)   = pre(n)
            end do
            
            do n = 1,dim_amin
c               cr(i,j,n+dim_c+dim_pre)= amin(n)
            end do
            do n = 1,dim_rad
c               cr(i,j,n+dim_c+dim_pre+dim_amin)=rad(n)	
            end do
            do n = 1,dim_pfug
c               cr(i,j,n+dim_c+dim_pre+dim_amin+dim_rad)=pfug(n)
            end do
c            press(i,j) = (pth(1) - 1.02d5)/1.0d5
c           permeability in m^2
c            kappa(i,j) = pth(5)/9.869233d-8
         end do
      end do

      end

c :: ----------------------------------------------------------
c :: Determine if aqueous components are saturated
c :: ----------------------------------------------------------
      subroutine FORT_SAT_TOUGH(s,DIMS(s),satv,DIMS(satv),
     &     press,DIMS(press),cr,DIMS(cr),
     &     nc,lo,hi,nd,idx_cr,mweight,
     &     kappa,DIMS(kappa),phi,DIMS(phi))

#include "probdata.H"

      integer lo(SDIM),hi(SDIM),nc,nd,idx_sat
      integer DIMDEC(s)
      integer DIMDEC(satv)
      integer DIMDEC(press)
      integer DIMDEC(kappa)
      integer DIMDEC(phi)
      integer DIMDEC(cr)
      integer idx_cr(4)
      REAL_T  s(DIMV(s),nc)
      REAL_T  satv(DIMV(satv))
      REAL_T  press(DIMV(press))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  phi(DIMV(phi))
      REAL_T  cr(DIMV(cr),nd)
      REAL_T  dt
      REAL_T  mweight(nc)

      integer i,j,n,NEQ
      integer(8) dim_pth, dim_ut, dim_c, dim_pre, dim_amin, dim_rad, dim_pfug
      REAL_T  pth(20)
      REAL_T  ut(nc-1), c(idx_cr(1)), pre(idx_cr(2))
      REAL_T  amin(idx_cr(3)), rad(idx_cr(3))
      REAL_T  pfug(idx_cr(4))
      REAL_T  st
      REAL_T  Y(50)

      dim_ut   = nc-1
      dim_c    = idx_cr(1)
      dim_pre  = idx_cr(2)
      dim_amin = idx_cr(3)
      dim_rad  = idx_cr(3)
      dim_pfug = idx_cr(4)
      pth = 0

      NEQ = ncomps
      
      do i = lo(1),hi(1)
         do j = lo(2),hi(2)

            pth(1) = 1.02d5 + press(i,j)*1.0d5
            pth(2) = 25.0
            pth(3) = 1.0
            pth(4) = 2530.0
            pth(5) = kappa(i,j)*9.869233d-8
            pth(6) = kappa(i,j)*9.869233d-8
            pth(7) = kappa(i,j)*9.869233d-8
            pth(8) = phi(i,j)
            pth(9) = 0.99716

c           conc in mol/L
            do n = 1,ncomps-1
               ut(n) = max(zero,s(i,j,n)/mweight(n))
               Y(n)  = ut(n)
            end do

            do n = 1,dim_c
               c(n) = cr(i,j,n)
            end do

            do n = 1,dim_pre
               pre(n) = s(i,j,ncomps+n-1)
               Y(ncomps-1+n) = pre(n)
            end do

            do n = 1,dim_amin
               amin(n) = cr(i,j,n+dim_c+dim_pre)
            end do
            do n = 1,dim_rad
		rad(n) = cr(i,j,n+dim_c+dim_pre+dim_amin)	
            end do
            do n = 1,dim_pfug
               pfug(n) = cr(i,j,n+dim_c+dim_pre+dim_amin+dim_rad)
            end do
            
            call Pre_Chem(pth,ut,c,pre,rad,amin,pfug)
            
            call MINERAL_SAT(st,Y,NEQ)

            satv(i,j) = st

         end do
      end do

      end

c :: ----------------------------------------------------------
c :: Initializes the secondary variables.
c     For efficiency: if s_current == s_previous, then
c     cr_current = cr_previous.
c :: ----------------------------------------------------------
      subroutine FORT_INIT_SECONDARY_TOUGH(s,DIMS(s),press,DIMS(press),
     &     cr,DIMS(cr),nc,lo,hi,nd,idx_cr,mweight,
     &     kappa,DIMS(kappa),phi,DIMS(phi))

c      implicit none

#include "probdata.H"

      integer lo(SDIM),hi(SDIM), nc, nd
      integer DIMDEC(s)
      integer DIMDEC(press)
      integer DIMDEC(kappa)
      integer DIMDEC(phi)
      integer DIMDEC(cr)
      integer idx_cr(4)
      REAL_T  s(DIMV(s),nc)
      REAL_T  press(DIMV(press))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  phi(DIMV(phi))
      REAL_T  cr(DIMV(cr),nd)
      REAL_T  dt
      REAL_T  mweight(nc)

      logical dif_from_prev
      integer i,j,n
      integer(8) dim_pth, dim_ut, dim_c, dim_pre, dim_amin, dim_rad, dim_pfug
      REAL_T  pth(20)
      REAL_T  ut(nc-1), c(idx_cr(1)), pre(idx_cr(2))
      REAL_T  amin(idx_cr(3)), rad(idx_cr(3))
      REAL_T  pfug(idx_cr(4))
      REAL_T  s_prev(nc)

      dim_pth  = 20
      dim_ut   = nc-1
      dim_c    = idx_cr(1)
      dim_pre  = idx_cr(2)
      dim_amin = idx_cr(3)
      dim_rad  = idx_cr(3)
      dim_pfug = idx_cr(4)
      pth      = 0

      do i = lo(1),hi(1)
         do j = lo(2),hi(2)
     
            dif_from_prev = .true.
            if (i .eq. lo(1) .and. j .eq. lo(2)) then
               do n = 1,nc
                  s_prev(nc) = s(i,j,n)
               end do
            else
               do n = 1,nc
                  if (dabs(s(i,j,n) - s_prev(n)) .le. 1.d-15) then
                     dif_from_prev = .false.
                  end if
               end do
            end if

            if (dif_from_prev) then          
c              pressure in Pa
               pth(1) = 1.02d5 + press(i,j)*1.0d5
c              temperature in C, set to 25.	
               pth(2) = 25.0
c              saturation of liquid, set to 1.
               pth(3) = 1.0
c	       rock bulk density in kg/m3, set to fix value 2530
               pth(4) = 2530.0
c              permeability in m^2
               pth(5) = kappa(i,j)*9.869233d-8
               pth(6) = kappa(i,j)*9.869233d-8
               pth(7) = kappa(i,j)*9.869233d-8
               pth(8) = phi(i,j)
c              density of water in kg/L
               pth(9) = 0.99716

c              conc in mol/L
               do n = 1,nc-1
                  ut(n) = max(zero,s(i,j,n)/mweight(n))
               end do

c              mineral abundance in mol/dm3 = mol/L
               do n = 1,dim_pre
                  pre(n) = s(i,j,nc+n-1)
               end do

               do n = 1,dim_ut
                  c(n) = ut(n)/pth(9)
               end do
               do n = dim_ut+1,dim_c
                  c(n) = 0.1
               end do
            
               do n = 1,dim_amin
                  amin(n) = rhosec(n+dim_c+dim_pre)
               end do

               do n = 1,dim_rad
                  rad(n) = rhosec(n+dim_c+dim_pre+dim_amin)	
               end do

               do n = 1,dim_pfug
                  pfug(n) = rhosec(n+dim_c+dim_pre+dim_amin+dim_rad)
               end do

               call INIT_SECONDARY(pth, dim_pth, ut, dim_ut, c, dim_c,
     &              pre, dim_pre, amin, dim_amin, rad, dim_rad, 
     &              pfug, dim_pfug)
            end if

            do n = 1,dim_c
               cr(i,j,n) = c(n)
            end do
            
            do n = 1,dim_pre 
               cr(i,j,n+dim_c) = pre(n)
            end do

            do n = 1,dim_amin
               cr(i,j,n+dim_c+dim_pre)= amin(n)
            end do
            do n = 1,dim_rad
               cr(i,j,n+dim_c+dim_pre+dim_amin)=rad(n)	
            end do
            do n = 1,dim_pfug
               cr(i,j,n+dim_c+dim_pre+dim_amin+dim_rad)=pfug(n)
            end do

         end do
      end do

      end
