#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "RICHARDSOLVER_F.H"
#include "ArrayLim.H"
#include "RockManager_F.H"

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

      subroutine FORT_RS_SATURATEDRES(Res,DIMS(Res),Pn,DIMS(Pn),
     &     Pnp1,DIMS(Pnp1),SpecificStorage,DIMS(SpecificStorage),
     &     Porosity,DIMS(Porosity),Src,DIMS(Src),dt,lo,hi,nc)
      implicit none      
      integer DIMDEC(Res)
      integer DIMDEC(Pn)
      integer DIMDEC(Pnp1)
      integer DIMDEC(Porosity)
      integer DIMDEC(SpecificStorage)
      integer DIMDEC(Src)
      integer lo(SDIM),hi(SDIM)
      integer nc
      REAL_T  Res(DIMV(Res),nc)
      REAL_T  Pn(DIMV(Pn),nc)
      REAL_T  Pnp1(DIMV(Pnp1),nc)
      REAL_T  Porosity(DIMV(Porosity))
      REAL_T  SpecificStorage(DIMV(SpecificStorage))
      REAL_T  Src(DIMV(Src),nc)
      REAL_T  rho(nc)
      REAL_T  dt
      
      integer i, j, n

c     On input, Res contains Div(rho.DarcyVel)
      do n=1,nc
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Res(i,j,n) = SpecificStorage(i,j)*(Pnp1(i,j,n) - Pn(i,j,n))/dt + Res(i,j,n) - Src(i,j,n)
            enddo
         enddo
      enddo
      end

      subroutine FORT_RS_RICHARDRES(Res,DIMS(Res),RSn,DIMS(RSn),
     &     RSnp1,DIMS(RSnp1),Porosity,DIMS(Porosity),Src,DIMS(Src),
     &     dt,lo,hi,nc)
      implicit none      
      integer DIMDEC(Res)
      integer DIMDEC(RSn)
      integer DIMDEC(RSnp1)
      integer DIMDEC(Porosity)
      integer DIMDEC(Src)
      integer lo(SDIM),hi(SDIM)
      integer nc
      REAL_T  Res(DIMV(Res),nc)
      REAL_T  RSn(DIMV(RSn),nc)
      REAL_T  RSnp1(DIMV(RSnp1),nc)
      REAL_T  Porosity(DIMV(Porosity))
      REAL_T  Src(DIMV(Src),nc)
      REAL_T  dt
      
      integer i, j, n

c     On input, Res contains Div(rho.DarcyVel)
      do n=1,nc
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Res(i,j,n) = Porosity(i,j)*(RSnp1(i,j,n) - RSn(i,j,n))/dt + Res(i,j,n) - Src(i,j,n)
            enddo
         enddo
      enddo
      end

      subroutine FORT_RS_ALTUP(rhosat,DIMS(rhosat),
     &     del_p,DIMS(del_p),alpha,DIMS(alpha),phi,DIMS(phi),
     &     matID,DIMS(matID),ls_factor,sThreshVS,rmID,time,lo,hi)

      implicit none

#include "probdata.H"

      integer rmID,lo(SDIM),hi(SDIM)
      REAL_T ls_factor,sThreshVS,time
      integer DIMDEC(rhosat)
      integer DIMDEC(p)
      integer DIMDEC(del_p)
      integer DIMDEC(alpha)
      integer DIMDEC(phi)
      integer DIMDEC(matID)
      REAL_T  rhosat(DIMV(rhosat))
      REAL_T   del_p(DIMV(del_p))
      REAL_T   alpha(DIMV(alpha))
      REAL_T     phi(DIMV(phi))
      integer  matID(DIMV(matID))

      integer i,j
      REAL_T dsdpc, del_s, s_k, s_kp1, rInv, pcap_k, pcap_kp1, lsInv

      rInv = 1/rhoval(1)
      lsInv = 1 / ls_factor
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            s_k = rhosat(i,j)*rInv
            call ROCK_MANAGER_DSDPCAP(rmID,s_k,matID(i,j),time,dsdpc,1)
            del_s = - ls_factor * del_p(i,j) * dsdpc * rInv
            s_kp1 = s_k - del_s
            if (s_kp1 .lt. sThreshVS) then
               call ROCK_MANAGER_PCAP(rmID,s_k  ,matID(i,j),time,pcap_k  ,1)
               call ROCK_MANAGER_PCAP(rmID,s_kp1,matID(i,j),time,pcap_kp1,1)
               del_p(i,j) = (pcap_kp1 - pcap_k) * lsInv
            endif
         end do 
      end do

      end
