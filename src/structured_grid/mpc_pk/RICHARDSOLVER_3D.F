#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "RICHARDSOLVER_F.H"
#include "ArrayLim.H"

#define SDIM 3

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
      
      integer i, j, k, n
      REAL_T tmult

      if (dir.eq.0) then
         tmult = 1.d0/dx(dir+1)
         do n=1,nc
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     e(i,j,k,n) = -( tmult*(c(i,j,k,n) - c(i-1,j,k,n)) + a(n))
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         tmult = 1.d0/dx(dir+1)
         do n=1,nc
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     e(i,j,k,n) = -( tmult*(c(i,j,k,n) - c(i,j-1,k,n)) + a(n))
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.2) then
         tmult = 1.d0/dx(dir+1)
         do n=1,nc
            do k = lo(3),hi(3)+1
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     e(i,j,k,n) = -( tmult*(c(i,j,k,n) - c(i,j,k-1,n)) + a(n))
                  enddo
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
      
      integer i, j, k, n

      do n=1,nc
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  X(i,j,k,n) = X(i,j,k,n) * Y(i,j,k,n) * Z(i,j,k,n)
               enddo
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
      integer i, j, k, n

      if (dir.eq.0) then
         if (do_upwind.eq.1) then
            do n=0,nc-1
               do k = lo(3),hi()
                  do j = lo(2),hi(2)
                     do i = lo(1),hi(1)+1
                        coef(i,j,k,n) = merge(lam(i-1,j,k,n),lam(i,j,k,n),vel(i,j,k)>=zero)
                     enddo
                  enddo
               enddo
            enddo
         else
            do n=0,nc-1
               do k = lo(3),hi(3)
                  do j = lo(2),hi(2)
                     do i = lo(1),hi(1)+1
                        coef(i,j,k,n) = 0.5d0 * (lam(i-1,j,k,n)+lam(i,j,k,n))
                     enddo
                  enddo
               enddo
            enddo
         endif
         if (dirbclo.eq.1) then 
            do n=0,nc-1
               i = lo(1)
               do k = lo(3),hi(3)
                  do j = lo(2),hi(2)
                     coef(i,j,k,n) = lam(i-1,j,k,n)
                  enddo
               enddo
            enddo
         endif
         if (dirbchi.eq.1) then 
            do n=0,nc-1
               i = hi(1)+1
               do k = lo(3),hi(3)
                  do j = lo(2),hi(2)
                     coef(i,j,k,n) = lam(i,j,k,n)
                  enddo
               enddo
            enddo
         endif
      else if (dir.eq.1) then
         if (do_upwind.eq.1) then
            do n=0,nc-1
               do k = lo(3),hi(3)
                  do j = lo(2),hi(2)+1
                     do i = lo(1),hi(1)
                        coef(i,j,k,n) = merge(lam(i,j-1,k,n),lam(i,j,k,n),vel(i,j,k)>=zero)
                     enddo
                  enddo
               enddo
            enddo
         else
            do n=0,nc-1
               do k = lo(3),hi(3)
                  do j = lo(2),hi(2)+1
                     do i = lo(1),hi(1)
                        coef(i,j,k,n) = 0.5d0 * (lam(i,j-1,k,n)+lam(i,j,k,n))
                     enddo
                  enddo
               enddo
            enddo
         endif
         if (dirbclo.eq.1) then 
            do n=0,nc-1
               j = lo(2)
               do k = lo(3),hi(3)
                  do i = lo(1),hi(1)
                     coef(i,j,k,n) = lam(i,j-1,k,n)
                  enddo
               enddo
            enddo
         endif
         if (dirbchi.eq.1) then 
            do n=0,nc-1
               j = hi(2)+1
               do k = lo(3),hi(3)
                  do i = lo(1),hi(1)
                     coef(i,j,k,n) = lam(i,j,k,n)
                  enddo
               enddo
            enddo
         endif
      else if (dir.eq.2) then
         if (do_upwind.eq.1) then
            do n=0,nc-1
               do k = lo(3),hi(3)+1
                  do j = lo(2),hi(2)
                     do i = lo(1),hi(1)
                        coef(i,j,k,n) = merge(lam(i,j,k-1,n),lam(i,j,k,n),vel(i,j,k)>=zero)
                     enddo
                  enddo
               enddo
            enddo
         else
            do n=0,nc-1
               do k = lo(3),hi(3)+1
                  do j = lo(2),hi(2)
                     do i = lo(1),hi(1)
                        coef(i,j,k,n) = 0.5d0 * (lam(i,j,k-1,n)+lam(i,j,k,n))
                     enddo
                  enddo
               enddo
            enddo
         endif
         if (dirbclo.eq.1) then 
            do n=0,nc-1
               k = lo(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     coef(i,j,k,n) = lam(i,j,k-1,n)
                  enddo
               enddo
            enddo
         endif
         if (dirbchi.eq.1) then 
            do n=0,nc-1
               k = hi(3)+1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     coef(i,j,k,n) = lam(i,j,k,n)
                  enddo
               enddo
            enddo
         endif
      else
         print *, 'FORT_CALC_RS_H::bad dir: ',dir
         call bl_pd_abort()
      endif
      end     

      subroutine FORT_RS_PDOTRES(Res,DIMS(Res),RSn,DIMS(RSn),
     &     RSnp1,DIMS(RSnp1),Porosity,DIMS(Porosity),
     &     rho,dt,lo,hi,nc)
      implicit none      
      integer DIMDEC(Res)
      integer DIMDEC(RSn)
      integer DIMDEC(RSnp1)
      integer DIMDEC(Porosity)
      integer lo(SDIM),hi(SDIM)
      integer nc
      REAL_T  Res(DIMV(Res),nc)
      REAL_T  RSn(DIMV(RSn),nc)
      REAL_T  RSnp1(DIMV(RSnp1),nc)
      REAL_T  Porosity(DIMV(Porosity))
      REAL_T  rho(nc)
      REAL_T  dt
      
      integer i, j, k, n

c     On input, Res contains Div(DarcyVel)
      do n=1,nc
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  Res(i,j,k,n) = Porosity(i,j,k)*(RSnp1(i,j,k,n) - RSn(i,j,k,n))/dt
     &                 + rho(n)*Res(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
      end
