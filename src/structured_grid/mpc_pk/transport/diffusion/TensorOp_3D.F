#include <REAL.H>

#include "TensorOp_F.H"
#include <LO_BCTYPES.H>
#include "ArrayLim.H"

      subroutine FORT_TOAPPLY(
     $     u, DIMS(u),
     $     alpha, beta,
     $     a, DIMS(a),
     $     muX, DIMS(muX),
     $     mu1X, DIMS(mu1X),
     $     muY, DIMS(muY),
     $     mu1Y, DIMS(mu1Y),
     $     muZ, DIMS(muZ),
     $     mu1Z, DIMS(mu1Z),
     $     out, DIMS(out),
     $     maskn,DIMS(maskn),
     $     maske,DIMS(maske),
     $     maskw,DIMS(maskw),
     $     masks,DIMS(masks),
     $     maskt,DIMS(maskt),
     $     maskb,DIMS(maskb),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     trandert,DIMS(trandert),
     $     tranderb,DIMS(tranderb),
     $     lo,hi,h
     $     )
      implicit none
      REAL_T alpha, beta
      integer DIMDEC(u)
      integer DIMDEC(a)
      integer DIMDEC(muX)
      integer DIMDEC(mu1X)
      integer DIMDEC(muY)
      integer DIMDEC(mu1Y)
      integer DIMDEC(muZ)
      integer DIMDEC(mu1Z)
      integer DIMDEC(out)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(maskt)
      integer DIMDEC(maskb)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer DIMDEC(trandert)
      integer DIMDEC(tranderb)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      REAL_T u(DIMV(u))
      REAL_T a(DIMV(a))
      REAL_T muX(DIMV(muX))
      REAL_T mu1X(DIMV(mu1X))
      REAL_T muY(DIMV(muY))
      REAL_T mu1Y(DIMV(mu1Y))
      REAL_T muZ(DIMV(muZ))
      REAL_T mu1Z(DIMV(mu1Z))
      REAL_T out(DIMV(out))
      integer  maskn(DIMV(maskn))
      integer  maske(DIMV(maske))
      integer  maskw(DIMV(maskw))
      integer  masks(DIMV(masks))
      integer  maskt(DIMV(maskt))
      integer  maskb(DIMV(maskb))

      REAL_T trandern(DIMV(trandern))
      REAL_T trandere(DIMV(trandere))
      REAL_T tranderw(DIMV(tranderw))
      REAL_T tranders(DIMV(tranders))
      REAL_T trandert(DIMV(trandert))
      REAL_T tranderb(DIMV(tranderb))

      integer i,j
      REAL_T dudye,dudyw,dudyn,dudys
      REAL_T dudxe,dudxw,dudxn,dudxs
      REAL_T xfluxw,xfluxe,yfluxs,yfluxn
      REAL_T hx, hxx, hx4
      REAL_T hy, hyy, hy4

      print *,'TensorOp_3D not ready yet'
      end


      subroutine FORT_TOFLUX(
     $     u, DIMS(u),
     $     muX, DIMS(muX),
     $     mu1X, DIMS(mu1X),
     $     muY, DIMS(muY),
     $     mu1Y, DIMS(mu1Y),
     $     muZ, DIMS(muZ),
     $     mu1Z, DIMS(mu1Z),
     $     xflux,DIMS(xflux),
     $     yflux,DIMS(yflux),
     $     zflux,DIMS(zflux),
     $     maskn,DIMS(maskn),
     $     maske,DIMS(maske),
     $     maskw,DIMS(maskw),
     $     masks,DIMS(masks),
     $     maskt,DIMS(maskt),
     $     maskb,DIMS(maskb),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     trandert,DIMS(trandert),
     $     tranderb,DIMS(tranderb),
     $     lo,hi,h
     $     )
      implicit none
      integer DIMDEC(u)
      integer DIMDEC(muX)
      integer DIMDEC(mu1X)
      integer DIMDEC(muY)
      integer DIMDEC(mu1Y)
      integer DIMDEC(muZ)
      integer DIMDEC(mu1Z)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(maskt)
      integer DIMDEC(maskb)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer DIMDEC(trandert)
      integer DIMDEC(tranderb)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      REAL_T u(DIMV(u))
      REAL_T muX(DIMV(muX))
      REAL_T mu1X(DIMV(mu1X))
      REAL_T muY(DIMV(muY))
      REAL_T mu1Y(DIMV(mu1Y))
      REAL_T muZ(DIMV(muZ))
      REAL_T mu1Z(DIMV(mu1Z))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T zflux(DIMV(zflux))
      integer  maskn(DIMV(maskn))
      integer  maske(DIMV(maske))
      integer  maskw(DIMV(maskw))
      integer  masks(DIMV(masks))
      integer  maskt(DIMV(maskt))
      integer  maskb(DIMV(maskb))

      REAL_T trandern(DIMV(trandern))
      REAL_T trandere(DIMV(trandere))
      REAL_T tranderw(DIMV(tranderw))
      REAL_T tranders(DIMV(tranders))
      REAL_T trandert(DIMV(trandert))
      REAL_T tranderb(DIMV(tranderb))

      integer i,j
      REAL_T hx,hy,hx4,hy4
      REAL_T dudxe,dudxw,dudxn,dudxs
      REAL_T dudye,dudyw,dudyn,dudys

      print *,'TensorOp_3D not ready yet'
      end

      subroutine FORT_TOGSRB (
     $     u, DIMS(u),
     $     rhs, DIMS(rhs),
     $     alpha, beta,
     $     a, DIMS(a),
     $     muX, DIMS(muX),
     $     mu1X, DIMS(mu1X),
     $     muY, DIMS(muY),
     $     mu1Y, DIMS(mu1Y),
     $     muZ, DIMS(muZ),
     $     mu1Z, DIMS(mu1Z),
     $     maskn,DIMS(maskn),
     $     fn, DIMS(fn),
     $     ftn, DIMS(ftn),
     $     maske,DIMS(maske),
     $     fe, DIMS(fe),
     $     fte, DIMS(fte),
     $     maskw,DIMS(maskw),
     $     fw, DIMS(fw),
     $     ftw, DIMS(ftw),
     $     masks,DIMS(masks),
     $     fs, DIMS(fs),
     $     fts, DIMS(fts),
     $     maskt,DIMS(maskt),
     $     ft, DIMS(ft),
     $     ftt, DIMS(ftt),
     $     maskb,DIMS(maskb),
     $     fb, DIMS(fb),
     $     ftb, DIMS(ftb),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     trandert,DIMS(trandert),
     $     tranderb,DIMS(tranderb),
     $     lo,hi,h,phaseflag
     $     )

      implicit none
      REAL_T alpha, beta
      integer DIMDEC(u)
      integer DIMDEC(rhs)
      integer DIMDEC(a)
      integer DIMDEC(muX)
      integer DIMDEC(mu1X)
      integer DIMDEC(muY)
      integer DIMDEC(mu1Y)
      integer DIMDEC(muZ)
      integer DIMDEC(mu1Z)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(maskt)
      integer DIMDEC(maskb)
      integer DIMDEC(fn)
      integer DIMDEC(fe)
      integer DIMDEC(fw)
      integer DIMDEC(fs)
      integer DIMDEC(ft)
      integer DIMDEC(fb)
      integer DIMDEC(ftn)
      integer DIMDEC(fte)
      integer DIMDEC(ftw)
      integer DIMDEC(fts)
      integer DIMDEC(ftt)
      integer DIMDEC(ftb)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer DIMDEC(trandert)
      integer DIMDEC(tranderb)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      integer phaseflag

      REAL_T u(DIMV(u))
      REAL_T rhs(DIMV(rhs))
      REAL_T a(DIMV(a))
      REAL_T muX(DIMV(muX))
      REAL_T mu1X(DIMV(mu1X))
      REAL_T muY(DIMV(muY))
      REAL_T mu1Y(DIMV(mu1Y))
      REAL_T muZ(DIMV(muZ))
      REAL_T mu1Z(DIMV(mu1Z))
      integer maskn(DIMV(maskn))
      integer maske(DIMV(maske))
      integer maskw(DIMV(maskw))
      integer masks(DIMV(masks))
      integer maskt(DIMV(maskt))
      integer maskb(DIMV(maskb))
      REAL_T fn(DIMV(fn))
      REAL_T fe(DIMV(fe))
      REAL_T fw(DIMV(fw))
      REAL_T fs(DIMV(fs))
      REAL_T ft(DIMV(ft))
      REAL_T fb(DIMV(fb))
      REAL_T ftn(DIMV(ftn))
      REAL_T fte(DIMV(fte))
      REAL_T ftw(DIMV(ftw))
      REAL_T fts(DIMV(fts))
      REAL_T ftt(DIMV(ftt))
      REAL_T ftb(DIMV(ftb))
      REAL_T trandern(DIMV(trandern))
      REAL_T trandere(DIMV(trandere))
      REAL_T tranderw(DIMV(tranderw))
      REAL_T tranders(DIMV(tranders))
      REAL_T trandert(DIMV(trandert))
      REAL_T tranderb(DIMV(tranderb))

      integer i,j,modx,mody
      integer istart,jstart
      REAL_T dudye,dudyw,dudyn,dudys
      REAL_T dudxe,dudxw,dudxn,dudxs
      REAL_T xfluxw,xfluxe,yfluxs,yfluxn
      REAL_T diag, oper
      REAL_T hx, hxx, hx4
      REAL_T hy, hyy, hy4
      REAL_T dw, de, ds, dn
      REAL_T dtw, dte, dts, dtn

      print *,'TensorOp_3D not ready yet'
      return

      end

      subroutine FORT_TOAPPLYBC (
     $     flagden, flagbc, maxorder,
     $     u,   DIMS(u),
     $     cdir, bct, bcl,
     $     bcval, DIMS(bcval),
     $     mask,  DIMS(mask),
     $     maskp, DIMS(maskp),
     $     maskm, DIMS(maskm),
     $     den,   DIMS(den),
     $     dent,  DIMS(dent),
     $     exttd, DIMS(exttd),
     $     tander,DIMS(tander),
     $     lo, hi, nc,
     $     h
     $     )
c
c     If the boundary is of Neumann type, set the ghost cell value to
c     that of the outermost point in the valid data (2nd order accurate)
c     and then fill the "den" array with the value "1"
c     
c     
c     If flagbc==1:
c     
c     If the boundary is of Dirichlet type, construct a polynomial
c     interpolation through the boundary location and internal points
c     (at locations x(-1:len-2) that generates the ghost cell value (at
c     location xInt).  Then fill the ghost cell with the interpolated value.
c     If flagden==1, load the "den" array with the interpolation
c     coefficient corresponding to outermost point in the valid region
c     ( the coef(0) corresponding to the location x(0) )
c      
c     Note: 
c     The bc type = LO_REFLECT_ODD is a special type of boundary condition.

c ::: the old FORT_APPLYBC put values out in ghost cells.  The new
c ::: FORT_APPLYBC makes a distinction between normal and tangential derivs.
c ::: the normal derivatives are still evaluated with ghost cells ( and
c ::: the den array for adjusting the diagonal element).  Tangential
c ::: derivatives are NOT computed from the ghost cells (note: we are
c ::: discussing here the tangential derivates which are centered on the
c ::: external cells walls of the rectangular domain.  Internal tangential
c ::: derivatives whose stencil extends outside the rectangular domain have
c ::: still ANOTHER calculational trick, one sided derivatives, which is 
c ::: implemented in the apply operator).  For these tangential derivatives,
c ::: if some of the cells in the stencil are uncovered, the derivative is
c ::: computed half from cells which are internal, plus an externally supplied
c ::: tangential derivative centered at location bcl.
c ??? exttd is the externally supplied tangential derivative which is
c ??? centered at bcl from wall.  Since the location of exttd is
c ??? indeterminate (determined by bcl), we arbitrarily index it in the
c ??? cell just outside the wall.
c ??? tander is the calculated tangential derivative which is centered at
c ??? at the wall.  It is indexed as if it were a cell centered quantity
c ??? just outside the wall.  This is anomalous since it is really edge
c ??? centered on the outer wall.
c      
c ::: other notes since previous developers didn't bother to document
c ::: cdir is mnemonic for coordinate direction, i.e. which side this is
c ::: cdir==0 -> left wall
c ::: cdir==2 -> right
c ::: cdir==1 -> bottom
c ::: cdir==3 -> top

c ::: notes about mask
c ::: since we now require other masks, there are 3.
c ::: 1) mask -- the main one covering the cells off the side of valid domain.
c ::: 2) maskp -- the plus end of transverse direction
c ::: 3) maskn -- the minus end of transverse direction

      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(u)
      REAL_T u(DIMV(u),nc)
      integer DIMDEC(den)
      REAL_T den(DIMV(den),nc)
      integer DIMDEC(dent)
      REAL_T dent(DIMV(dent),nc)
      integer DIMDEC(exttd)
      REAL_T exttd(DIMV(exttd),nc)
      integer DIMDEC(tander)
      REAL_T tander(DIMV(tander),nc)
      integer DIMDEC(bcval)
      REAL_T bcval(DIMV(bcval),nc)
      integer DIMDEC(mask)
      integer mask(DIMV(mask))
      integer DIMDEC(maskp)
      integer maskp(DIMV(maskp))
      integer DIMDEC(maskm)
      integer maskm(DIMV(maskm))
      integer bct(nc)
      REAL_T bcl
      REAL_T h(BL_SPACEDIM)
      end
