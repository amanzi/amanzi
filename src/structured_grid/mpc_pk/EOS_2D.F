#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "PROB_PM_F.H" 
#include "ArrayLim.H"

#define SDIM 2

c ::: ----------------------------------------------------------
c ::: CHEM
c :::    A + B <=> C
c :::    Simple 3-component chemistry module.
c ::: ----------------------------------------------------------
      subroutine FORT_CHEM(s,DIMS(s),nc,lo,hi,dt,press,DIMS(press),
     &     kappa,DIMS(kappa),phi,DIMS(phi),cr,DIMS(cr),nd,idx_cr,mweight)
      implicit none

#include "probdata.H"

      integer lo(SDIM),hi(SDIM), nc, nd
      integer DIMDEC(s)
      REAL_T  s(DIMV(s),nc)
      REAL_T  dt

c     unused 
      integer DIMDEC(press)
      integer DIMDEC(kappa)
      integer DIMDEC(phi)
      integer DIMDEC(cr)
      integer idx_cr(4)
      REAL_T  press(DIMV(press))
      REAL_T  kappa(DIMV(kappa))
      REAL_T  phi(DIMV(phi))
      REAL_T  cr, mweight

      integer i,j,n,r,eq
      REAL_T  deltan(3), rate_f, rate_b, Keq, a,b,c, x(2)

      external FEX, JEX
      integer  NEQ, ITOL,ITASK,ISTATE,IOPT,IWORK(33),LRW,LIW,IPAR, MF
      REAL_T   Y(3), T, TOUT, RTOL, ATOL(3), RPAR(2), RWORK(67)

      
      rate_f = 1.d1
      rate_b = 1.d1

      if (eq .ne. 1) then

         RPAR(1) = rate_f
         RPAR(2) = rate_b

         NEQ = 3
         do i = ARG_L1(s),ARG_H1(s)
            do j = ARG_L2(s),ARG_H2(s)
               Y(1)   = max(zero,s(i,j,2))/18.d0
               Y(2)   = max(zero,s(i,j,3))/18.d0
	       if (model .eq. 1) then
c		  Y(3) = 0
		  Y(3) = max(zero,s(i,j,4))/32.d0
               else	
                  Y(3) = max(zero,s(i,j,4))/32.d0
  	       end if
               T      = zero
               TOUT   = dt
               ITOL   = 1
               RTOL   = 1.d-12
               ATOL   = 1.d-12
               ITASK  = 1
               ISTATE = 1
               IOPT   = 0
               LRW    = 67
               LIW    = 33
               MF     = 21

               
 99            CALL DVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     &              IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)   

               if (ISTATE .ne. 2) then
                  print *, "ISTATE = ", ISTATE
                  ISTATE = 1
                  goto 99
               end if

               s(i,j,2) = Y(1)*18.d0
               s(i,j,3) = Y(2)*18.d0
	       if (model .eq. 1) then
c                  s(i,j,4) = s(i,j,4) + Y(3)*32.0
		  s(i,j,4) = Y(3)*32.d0
 	       else
                  s(i,j,4) = Y(3)*32.d0
	       end if
               s(i,j,1) = rhoval(1) - s(i,j,2)-s(i,j,3)
               if (model .ne. 1) then
                  s(i,j,1) = s(i,j,1) - s(i,j,4)
               end if

            end do
         end do

      else

         Keq = rate_f/rate_b
         do i = ARG_L1(s),ARG_H1(s)
            do j = ARG_L2(s),ARG_H2(s)
               Y(1)   = max(zero,s(i,j,2))/18.d0
               Y(2)   = max(zero,s(i,j,3))/18.d0
               Y(3)   = max(zero,s(i,j,4))/32.d0
               a =  Keq
               b = -Keq*(Y(1)+Y(2))-1
               c =  Keq*(Y(1)*Y(2))-Y(3)
               x(1) = (-b + sqrt(b**2-4*a*c))/(2*a)
               x(2) = (-b - sqrt(b**2-4*a*c))/(2*a)

               if ((Y(1)-x(1)) .ge. 0 .and. (Y(2)-x(1)) .ge. 0 .and. 
     &             (Y(3)+x(1)) .ge. 0) then
                  s(i,j,2) = (Y(1)-x(1))*18.d0
                  s(i,j,3) = (Y(2)-x(1))*18.d0        
                  s(i,j,4) = (Y(3)+x(1))*32.d0
               else
                  s(i,j,2) = (Y(1)-x(2))*18.d0
                  s(i,j,3) = (Y(2)-x(2))*18.d0        
                  s(i,j,4) = (Y(3)+x(2))*32.d0
               end if

               s(i,j,1) = rhoval(1) - s(i,j,2)-s(i,j,3)-s(i,j,4)
            end do
         end do
      end if
      end

      subroutine FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
      integer NEQ, IPAR
      REAL_T RPAR(2), T, Y(NEQ), YDOT(NEQ)

      REAL_T rate_f, rate_b
      rate_f = RPAR(1)
      rate_b = RPAR(2)
       
      YDOT(1) = -rate_f*Y(1)*Y(2) + rate_b*Y(3)
      YDOT(2) = -rate_f*Y(1)*Y(2) + rate_b*Y(3)
      YDOT(3) = rate_f*Y(1)*Y(2) - rate_b*Y(3) 
  
      end

      subroutine JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      integer NEQ, ML,MU, NRPD, IPAR
      REAL_T PD(NRPD,NEQ), RPAR(2), T, Y(NEQ)
      REAL_T rate_f, rate_b
      rate_f = RPAR(1)
      rate_b = RPAR(2)

      PD(1,1) = -rate_f*Y(2)
      PD(1,2) = -rate_f*Y(1)
      PD(1,3) = rate_b
      PD(2,1) = -rate_f*Y(2)
      PD(2,2) = -rate_f*Y(1)
      PD(2,3) = rate_b
      PD(3,1) = rate_f*Y(2)
      PD(3,2) = rate_f*Y(1)
      PD(3,3) = -rate_b    
 
      end


      subroutine calc_rho(rho,s,p,rhoval,ncomp,nphases)
c ::: -----------------------------------------------------------
c ::: Compute the density.  Assume the first component is water
c ::: the second component is compressible gas. 
c ::: -----------------------------------------------------------

      integer ncomp
      integer nphases
      REAL_T  rho(ncomp)
      REAL_T  s(ncomp)
      REAL_T  rhoval(ncomp)
      REAL_T  p(nphases)

      integer i
      REAL_T  grad(ncomp)

      if (ncomp .ne. nphases) then
         call bl_abort('limited at present to cases where ncomp=nphases.')
      end if

      grad(1) = zero
      grad(2) = 0.808d-6*1.0d8
      grad(2) = zero

      do i = 1,ncomp
         rho(i) = rhoval(i) + grad(i)*p(i)
      end do

      end
      
      subroutine calc_drhodn(drhodn,rho,ncomp)
c ::: -----------------------------------------------------------
c ::: Compute the derivative of density.  First component is water
c ::: second component is compressible gas. 
c ::: -----------------------------------------------------------

      integer ncomp
      REAL_T  drhodn(ncomp)
      REAL_T  rho(ncomp)

      REAL_T  grad(ncomp)

      grad(1) = zero
      grad(2) = 0.808d-6*1.0d8
      grad(2) = zero

      drhodn(1) = zero
      drhodn(2) = -grad(2)/rho(2)

      end

      subroutine calc_1dpressure(p,y,p0,rho,grav,ncomp)
c ::: -----------------------------------------------------------
c ::: Compute the pressure,given the distribution of s 
c ::: -----------------------------------------------------------
      integer ncomp
      REAL_T  p,y,p0
      REAL_T  grav
      REAL_T  rho(ncomp)
      REAL_T  grad(ncomp)

      grad(1) = zero
      grad(2) = 0.808d-6*1.0d8
      grad(2) = zero

      p = (1/grad(2))*((rho(2)+grad(2)*p0)*
     &     exp(-grad(2)*grav*y) - rho(2))
      p = p0 - rho(2)*grav*y
      end
