#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "POROUS_F.H"

c
c***********************************************************************
c*                                                                     *
c*      turn -- general version                                        *
c*      program for generating a three-dimensional, stationary,        *
c*      correlated anisotropic random field                            *
c*                                                                     *
c*      andrew f. b. tompson                                           *
c*      parsons laboratory                                             *
c*      mit                                                            *
c*      cambridge ma 02139                                             *
c*                                                                     *
c*      now at:                                                        *
c*      atmosphere, earth, and energy division, l-208                  *
c*      lawrence livermore national lab                                *
c*      box 808, livermore ca 94551  usa                               *
c*      (925) 422 6348                                                 *
c*      afbt@llnl.gov                                                  *
c*                                                                     *
c*      method: the turning bands algorithm                            *
c*                                                                     *
c*      generate normal field: cond = z(x,y,z): n(0,1)                 *
c*      other fields possible: (1)    z*(x,y,z) = sig*z + f:  n(f,sig) *
c*                             (2)    k(x,y,z) = exp(z*): log normal   *
c*                                                                     *
c*                                                                     *
c*      option (2) implemented in subroutine tlog (if ilog=1) where:   *
c*                 condg = exp(f) = exp(e [ln(k)])                     *
c*                 sig = sqrt(var [ln(k)])                             *
c*                                                                     *
c***********************************************************************
c
      subroutine FORT_TURN(cond,nxc,nyc,nzc,iuc)

#include "TURN_F.h"

      integer nxc,nyc,nzc,iuc
      integer m,n6,n7,n8,n9,n10,n11,ntest,ntestt
      integer inphys,inturn,normout,ntout
      REAL_T  ks(3)
      REAL_T  cond(nxc,nyc,nzc)
c
c.....open log file
c

      open(16,file='turn.log')

      inphys=5
      inturn=7
      normout=8
      ntout=9
c
c.....input
c
      call input(ks,inphys,inturn)
c
      nx = nxc
      ny = nyc
      nz = nzc
      iu = iuc

c$$$      nline  = nx*2
c$$$      delzet = 0.1/10.0
c$$$      nztest = nline*100
c$$$      bigk   = 100
c$$$      dk     = 0.1

      m=int(bigk/dk)
      n1=1+nx*ny*nz
      n2=n1+ndelt+1
      n3=n2+nztest
      n4=n3+m
      n5=n4+m
      n6=n5+m
      n7=n6+nl
      n8=n7+nl
      n9=n8+nl
      n10=n9+nl
      n11=n10+nl
      n12=n11+nl
      n13=n12+nline
      ntestt=n13+nline-1

 1010 call driver(cond,ks,normout,ntout)

      close(16)
      
      end
c
c-----------------------------------------------------------------------
c       driver routine
c-----------------------------------------------------------------------
c
      subroutine driver(cond,ks,normout,ntout)
c
c-----------------------------------------------------------------------
c        b l a n k  c o m m o n  s t o r a g e
c-----------------------------------------------------------------------
c        cond(i,j,k). .generated field                <begins at a(1)>
c        h(i) . . . . .histogram vector               <begins at a(n1)>
c        zline(i) . . .realization of line process    <begins at a(n2)>
c        sk(i). . . . .spectral kernal                <begins at a(n3)>
c        ohmpri(i). . .perturbed wave #               <begins at a(n4)>
c        fik(i) . . . .random angle                   <begins at a(n5)>
c        rx(i). . . . .ensemble x covariance          <begins at a(n6)>
c        rxx(i) . . . .single x covariance            <begins at a(n7)>
c        ry(i). . . . .ensemble y covariance          <begins at a(n8)>
c        ryy(i) . . . .single y covariance            <begins at a(n9)>
c        rz(i). . . . .ensemble z covariance          <begins at a(n10)>
c        rzz(i) . . . .single z covariance            <begins at a(n11)>
c        phi(i) . . . .random phi angle for line i    <begins at a(n12)>
c                         (between 0 and pi)
c        th(i). . . . .random theta angle for line i  <begins at a(n13)>
c                         (between 0 and 2pi)
c-----------------------------------------------------------------------
c        a(ntest)      mass storage vector
c        cmin,cmax     histogram limits
c        condg         kg = exp(f) where f = e [ln(k)]
c                           (use only if ilog = 1)
c        bigk          maximum normalized frequency considered in line
c                           process
c        delzet        distance between points in line process
c        dk            normalized frequency increment for line process
c        dx,dy,dz      nodal spacings
c        ilog          = 1 => convert realization to lognormal
c        iu            seed for random number generator
c        ks(3)         ks(i)=1 => compute statistics in direction i for
c                           each realization
c        n1,n2,...n13  array locations within the mass storage vector a
c        ndelt         histogram limits
c        nfile,tfile   raw data files (unformatted)
c        nl            1+number of lags computed in covariances
c        nline         number of lines chosen
c        nmont         number of monte-carlo runs
c        nx,ny,nz      actual dimensions of cond for x,y,z
c        nzet          number of points on line process
c        nztest        estimate of number of points on line process
c        rm,rssq       random number generator running statistics
c        sig           desired standard deviation of ln(k)
c                           (only if ilog = 1)
c        th,phi        theta-phi line directions
c        theta,fi      single theta-phi direction of line l
c        xl1,xl2,xl3   correlation scales in x,y,z directions
c
c       storage and file assignments
c
c-----------------------------------------------------------------------

#include "TURN_F.h"

      integer normout, ntout
      character*12 nofile
      double precision cond(nx,ny,nz),ks(3),h(ndelt+1)
      double precision rx(nl),ry(nl),rz(nl),rxx(nl)
      double precision ryy(nl),rzz(nl)
      double precision rl, amea, assq,um1,um2,rvar1,rvar2,fu,amm,ass,count
      double precision x1,x2,x3,x4,var,rxr,ryr,rzr,fnru,ru,fnr,randloc,delt
      integer i,j,k,l,ll,m,n,mc,nlevel,istat,nru
c
      nofile='no file     '
      nran=0
      rm=0.
      rssq=0.
      pi=3.141592654
c
c.....zero statistical variables
c
      ndelt1=ndelt+1
      do 100 l=1,ndelt1
      h(l)=0.
100   continue
      delt=(cmax-cmin)/float(ndelt)
      rl=cmin
      do 8 i=1,nl
      rxx(i)=0.
      ryy(i)=0.
      rzz(i)=0.
      rx(i)=0.
      ry(i)=0.
8     rz(i)=0.
      amea=0.
      assq=0.
      um1=0.
      um2=0.
      rvar1=0.
      rvar2=0.
c
c.....start monte carlo loop..............................................
c                                                                        .
c                                                                        .
c
      fu=float(nmont)
      do 1111 mc=1,nmont
      nlevel=0

      call gener(mc,cond)

c.....print normal data ..........(unit=normout)
c
c      if (nfile.ne.nofile) call printout(cond,nfile,normout)
c
c.....print color images
c
c     open(15,file='picture',form='unformatted',
c    +     access='direct',recl=1)
c     write(15,rec=1)((cond(i,j,1),i=1,nx),j=1,ny)
c     close(15)
c
c.....update sample statistics
c
      istat=ks(1)+ks(2)+ks(3)
      if(istat.lt.1) go to 1109
      if(ks(1).eq.1) call statloc(cond,rxx,amm,ass,nl1,1)
      if(ks(2).eq.1) call statloc(cond,ryy,amm,ass,nl2,2)
      if(ks(3).eq.1) call statloc(cond,rzz,amm,ass,nl3,3)
      amea=amea+amm
      assq=assq+ass
      do 22 i=1,nl
      rx(i)=rx(i)+rxx(i)
      ry(i)=ry(i)+ryy(i)
      rz(i)=rz(i)+rzz(i)
   22 continue
c
c.....histogram
c
      call histo(cond,h,mc,delt)
1109  continue
c.....log transform data
c
      if (ilog .eq. 1) call tlog(cond)
c
c.....print transformed data .....(unit=ntout)
c
      if (tfile.ne.nofile) call printout(cond,tfile,ntout)
c
 1111 continue
      close(7)
c                                                                        .
c                                                                        .
c                                                                        .
c.....end monte carlo loop................................................
c
c.....write ensemble statistics
c
      if(istat.lt.1) go to 8888
      x1=float(nmont)
      x3=x1*float(nx*ny*nz)
      x4=x3-1.
      amea=amea/x3
      var=(assq-x3*amea*amea)/x4
      write(16,9150) amea,var
      write(16,9200)
      rxr=0.
      ryr=0.
      rzr=0.
      do 200 i=1,nl
      if(i.gt.nx.or.ks(1).eq.0) go to 811
      nru=(nx-i+1)*ny*nz*nmont
      fnru=float(nru)
      rxr=rx(i)/fnru-amea*amea
      rxr=rxr*x3/x4
811   if(i.gt.ny.or.ks(2).eq.0) go to 821
      nru=nx*(ny-i+1)*nz*nmont
      fnru=float(nru)
      ryr=ry(i)/fnru-amea*amea
      ryr=ryr*x3/x4
821   if(i.gt.nz.or.ks(3).eq.0) go to 831
      nru=nx*ny*(nz-i+1)*nmont
      fnru=float(nru)
      rzr=rz(i)/fnru-amea*amea
      rzr=rzr*x3/x4
831   continue
      write(16,9201) i-1,rxr,ryr,rzr
200   continue
c
c.....write ensemble histogram
c
      write(16,330)
      rl=cmin
      count=0.
      ll=1
      do 2001 l=2,ndelt1
      count=100.*h(l)/x3
      ru=rl+delt
      if(l.eq.ndelt1) ru=cmax
      write(16,870) l-1,rl,ru,count
      rl=ru
2001  continue
c
c.....write uniform random number generator statistics
c
8888  continue
      fnr=float(nran)
      rm=rm/fnr
      rssq=(rssq-fnr*rm*rm)/(fnr-1.)
      write(16,9998) nran,rm,rssq
      write(16,9999)
88    format(' end of monte carlo simulation no. ',i5)
9150  format('  ensemble statistics :',//,
     +   '   global mean        ',e13.4,/
     +   '   global variance    ',e13.4)
9200  format(/,'   covariances',/,
     +  7x,'lag        x dir        y dir        z dir',/)
9201  format(i10,3e13.4)
330   format(//,'  ensemble histogram ',//
     +    '    interval    lower       upper       pct  ',/,
     +    '                limit       limit           ',//)
870   format(i10,2e12.3,f10.2)
901   format(2i6)
9998  format('rng sample size,mean, and variance ',i10,2e12.3)
9999  format(1x,'***normal end of execution***')
      end
c
c-----------------------------------------------------------------------
c     this subroutine generates random field with turning bands method
c-----------------------------------------------------------------------
c
      subroutine gener(mc,cond)

#include "TURN_F.h"

      double precision zline(nztest),cond(nx,ny,nz),th(nline),phi(nline)
      double precision randloc,xx,theta,fi,xli,yli,zli,zetmax,zetmin
      double precision x,y,z,zet,xx1,clambda
      integer i,j,k,l,mc,nzet,index
c
c.....obtain fixed line directions (nline .eq. 15)
c
      if (nline.ne.15) go to 6
      th(1)=0.
      th(2)=0.
      th(3)=1.5708
      th(4)= 1.017222
      th(5)=-1.017222
      th(6)=.3648638
      th(7)=-.3648638
      th(8)=-1.017222
      th(9)=1.017222
      th(10)=1.017222
      th(11)=-1.017222
      th(12)=-1.017222
      th(13)=1.017222
      th(14)=-.3648638
      th(15)=.3648638
c
      phi(1)=1.5708
      phi(2)=0.
      phi(3)=1.5708
      phi(4)=1.256637
      phi(5)=-1.256637
      phi(6)=1.047198
      phi(7)=-1.047198
      phi(8)=.6283185
      phi(9)=-.6283185
      phi(10)=.6283185
      phi(11)=-.6283185
      phi(12)=1.256637
      phi(13)=-1.256637
      phi(14)=1.047198
      phi(15)=-1.047198
      goto 8
c
c.....obtain random line directions (nline .ne. 15)
c
    6 do 7 i=1,nline
      th(i)=randloc(iu)*pi*2.0
      xx=randloc(iu)
      phi(i)=acos(1.-2.*xx)
    7 continue
c
c.....zero the random field array
c
    8 do 9 i=1,nx
      do 9 j=1,ny
      do 9 k=1,nz
      cond(i,j,k) = 0.0
    9 continue
c
c.....start line loop......................................................
c                                                                         .
c.....zeta min and zeta max                                               .
c                                                                         .
      do 3333 l=1,nline
      theta=th(l)
      fi=phi(l)
      xli=cos(theta)*sin(fi)
      yli=sin(theta)*sin(fi)
      zli=cos(fi)
      zetmax=-1.e30
      zetmin=1.e30
      do 111 i=1,2
      x=nx*dx*(i-1)
      do 111 j=1,2
      y=ny*dy*(j-1)
      do 111 k=1,2
      z=nz*dz*(k-1)
      zet=x*xli+y*yli+z*zli
      if(zet.lt.zetmin)zetmin=zet
      if(zet.gt.zetmax)zetmax=zet
111   continue
c
c.....generate random line process
c
      nzet=int((zetmax-zetmin)/delzet)+2

      if (nzet .gt. nztest) then
         print *, "increasing ", nztest, " to ", nzet*2
         nztest = nzet*2
      end if
      
c      if(nzet.le.nztest) go to 112
c      write(16,113) nzet
c      stop

112   call liproc(cond,theta,fi,zetmin,zetmax,zline,nzet)
c
c.....project line process to generate 3-d array
c
      do 1111 i=1,nx
      x=dx*(i-1)
      do 1111 j=1,ny
      y=dy*(j-1)
      do 1111 k=1,nz
      z=dz*(k-1)
      zet=x*xli+y*yli+z*zli
      index=int((zet-zetmin)/delzet+.5)+1
      if (index.le.nzet) goto 1110
      stop
1110  cond(i,j,k)=cond(i,j,k)+zline(index)
1111  continue
3333  continue
c                                                                         .
c                                                                         .
c                                                                         .
c.....end line loop........................................................
c
c.....divide by sqrt(l)
c
      xx1=1./sqrt(float(nline))
      do 4444 i=1,nx
      do 4444 j=1,ny
      do 4444 k=1,nz
      cond(i,j,k)=cond(i,j,k)*xx1
4444  continue
c
c
  113 format (1x,'*****estimate of nzet too small; increase to ',i8,
     +       '*****')
      return
      end
c
c-----------------------------------------------------------------------
c      this subroutine generates a correlated random line process
c      using the methods of rice (1954); shinozuka and jan (1972)
c-----------------------------------------------------------------------
c
      subroutine liproc(cond,theta,fi,zetmin,zetmax,zline,
     +           nzet)
c
#include "TURN_F.h"

      double precision cond(nx,ny,nz),sk(nztest)
      double precision ohmpri(nztest),fik(nztest),zline(nztest)
      double precision um,rvar,x,y,z,clamda,delohm,deldel,randloc
      double precision wk,ss,theta,fi,zet,zetmin,zetmax,sum,arg,fu
      integer i,k,m,nzet	
c
c.....initialization...convert bigk and dk to ohmax and delohm =>
c                      m = # points on line
c
      um=0.
      rvar=0.
      x=cos(theta)*sin(fi)
      y=sin(theta)*sin(fi)
      z=cos(fi)
      clamda=sqrt((xl1*x)**2+(xl2*y)**2+(xl3*z)**2)
      delohm=dk/clamda
      m=int(bigk/dk)
      deldel=delohm/20.
c
c.....calculate random arguments and amplitudes
c.....hole spectrum on line (=ss) used for exponential covariance
c          function in space
c
      do 10 k=1,m
      fik(k)=pi*randloc(iu)*2.
      wk=delohm*(k-0.5)
      ohmpri(k)=deldel*2.*(randloc(iu)-.5)+wk
c
c.....change spectrum here..................................
c                                                          .
c                                                          .
      ss=2.0*xl1*xl2*xl3*wk*wk/pi/
     +    (1.+(xl1*wk*x)**2+(xl2*wk*y)**2+(xl3*wk*z)**2)**2
c                                                          .
c                                                          .
c...........................................................
c
      sk(k)=sqrt(ss*delohm)
10    continue
c
c.....calculate realization
c
      do 1111 i=1,nzet
      zet=zetmin+delzet*float(i-1)
      sum=0.0
      do 2222 k=1,m
      arg=ohmpri(k)*zet+fik(k)
      sum=sum+sk(k)*cos(arg)
2222  continue
      zline(i)=2.0*sum
      um=um+zline(i)
1111  continue
      fu=float(nzet)
      um=um/fu
      return
      end
c
c-----------------------------------------------------------------------
c     function to generate pseudo random numbers, uniform (0,1)
c-----------------------------------------------------------------------
c
      function randloc(iuu)
c
#include "TURN_F.h"
c      common / internal / pi,rm,rssq,n1,n2,n3,n4,n5,n12,n13,nran,ndelt1
c
      double precision randloc
      integer n, m, l,iuu
	
      data m/1048576/
      data l/1027/
      n=l*iuu
      iuu=mod(n,m)
      randloc=float(iuu)/float(m)
      nran=nran+1
      rm=rm+randloc
      rssq=rssq+randloc*randloc
      return
      end
c
c-----------------------------------------------------------------------
c     this subroutine calculates sample statistics of the random field
c-----------------------------------------------------------------------
c
      subroutine statloc(cond,ru,umea,ussq,nli,id)
c
#include "TURN_F.h"

      double precision ru(nl),cond(nx,ny,nz),umea,ussq,ck,fu,um,uv,ut,utu
      double precision uc	
      integer nu(3), i,j,k,nli,id,ntot,nutt,kk,nru
c
c.....initialization
c
       do 2 i=1,nli
       ru(i)=0.0
2      continue
       umea=0.0
       ussq=0.0
c
c.....update mean and sum-squares sums
c
      nu(1)=nx
      nu(2)=ny
      nu(3)=nz
      do 10 i=1,nu(1)
      do 10 j=1,nu(2)
      do 10 k=1,nu(3)
      ck=cond(i,j,k)
      umea=umea+ck
      ussq=ussq+ck*ck
10    continue
      ntot=nu(1)*nu(2)*nu(3)
      fu=float(ntot)
      um=umea/fu
      uv=(ussq-fu*um*um)/(fu-1.)
c
c.....calculate auto covariance
c
      nutt=nu(id)
      do 30 kk=1,nli
      if(kk.gt.nutt) go to 30
      nu(id)=nutt-kk+1
      do 20 i=1,nu(1)
      do 20 j=1,nu(2)
      do 20 k=1,nu(3)
      ut=cond(i,j,k)
      if(id.eq.1) utu=cond(i+kk-1,j,k)
      if(id.eq.2) utu=cond(i,j+kk-1,k)
      if(id.eq.3) utu=cond(i,j,k+kk-1)
      ru(kk)=ru(kk)+ut*utu
20    continue
30    continue
c
c.....write single realization statistics
c
      nu(id)=nutt
      write(16,8874) um,uv
      do 40 kk=1,nli
      nu(id)=nutt-kk+1
      nru=nu(1)*nu(2)*nu(3)
      uc=ru(kk)/float(nru)-um*um
      uc=uc*fu/(fu-1.)
      write(16,8876) kk-1,uc
8874  format(' mean,var',2e13.4)
8876  format(' lag,cov',i10,e13.4)
40    continue
      return
      end
c
c-----------------------------------------------------------------------
c     this subroutine reads all options and turning band inputs
c-----------------------------------------------------------------------
c
      subroutine input(ks,inphys,inturn)
c
#include "TURN_F.h"      
c      common / data1 / dx,dy,dz,xl1,xl2,xl3,condg,sig,nx,ny,nz,idate,
c     +                 idnum
c      common / data2 / delzet,bigk,dk,cmin,cmax,nline,nztest,nmont,
c     +                 ndelt,ilog,iu,nl1,nl2,nl3,nl
c      common / internal / pi,rm,rssq,n1,n2,n3,n4,n5,n12,n13,nran,ndelt1
c      common / file / nfile,tfile
c      character*12 nfile,tfile
      double precision ks(3)
      integer inphys,inturn,idrun
c
c.....open files
c
      open(inphys,file='turn1.dat')
      open(inturn,file='turn2.dat')
c
c.....physical data from turn1.dat........(unit=inphys)
c
      read(inphys,*) idate,idrun
      read(inphys,*) nx,ny,nz
      read(inphys,*) dx,dy,dz
      read(inphys,*) xl1,xl2,xl3
      read(inphys,*) condg,sig
c
c.....turning bands data from turn2.dat...(unit=inturn)
c
      read(inturn,*) nline,delzet,nztest
      read(inturn,*) bigk,dk
      read(inturn,*) nmont,iu,nl
      read(inturn,*) cmin,cmax,ndelt
      read(inturn,*) ks(1),ks(2),ks(3)
      read(inturn,*) ilog
      read(inturn,5) nfile,tfile
c
c.....close files
c
      close(inphys)
      close(inturn)
c
c.....adjust lag
c
      nl1=nl
      nl2=nl
      nl3=nl
      if (((nl-nx)*ks(1)).gt.0)nl1=nx
      if (((nl-ny)*ks(2)).gt.0)nl2=ny
      if (((nl-nz)*ks(3)).gt.0)nl3=nz
c
c.....print data and initial information to log file....unit=6
c
      write(16,10) idate,idnum
      write(16,20) nx,ny,nz,dx,dy,dz,xl1,xl2,xl3
      if(ilog.eq.1) write(16,21) condg,sig
      if(ilog.ne.1) write(16,22)
      write(16,30) nline,delzet,bigk,dk,nmont,iu,nl1,cmin,cmax,ndelt,
     +            nfile,tfile
    5 format(a12/a12)
   10 format ('1'/////,' ',70(1h*)//' turning bands program'//' ',
     +       70(1h*)//////////' idate = ',i10///' idnum = ',i10/'1')
   20 format (//' ',70(1h-)//' physical data'//' ',70(1h-)//
     +       ' number of points. . . x-direction...',i10/
     +       '                       y-direction...',i10/
     +       '                       z-direction...',i10//
     +       ' nodal spacing . . . . x-direction...',f10.4/
     +       '                       y-direction...',f10.4/
     +       '                       z-direction...',f10.4//
     +       ' correlation length. . x-direction...',f10.4/
     +       '                       y-direction...',f10.4/
     +       '                       z-direction...',f10.4//)
   21 format (//
     +       ' kg = exp (mean of ln(k) field.......',f10.4/
     +       ' standard deviation of ln(k) field...',f10.4//)
   22 format (//
     +       ' normal (0,1) field being generated ')
   30 format ('1'//' ',70(1h-)//' turning bands data'//' ',70(1h-)//
     +       ' number of lines ....................',i10/
     +       '     (chosen in random directions if '/
     +       '      different from 15)'/
     +       ' nodal spacing along each line.......',f10.4/
     +       ' maximum normalized spectral freq....',f10.4/
     +       ' normalized spectral freq increment..',f10.4/
     +       ' number of monte carlo runs..........',i10/
     +       ' starting seed for random numbers....',i10/
     +       ' maximum computed lag................',i10/
     +       ' lower histogram limit...............',f10.4/
     +       ' upper histogram limit...............',f10.4/
     +       ' number of histogram increments......',i10//
     +       ' file to which unformatted normal'/
     +       '      data routed ...................',a12/
     +       ' file to which unformatted transformed'/
     +       '      data routed ...................',a12//'1')
      return
      end
c
c-----------------------------------------------------------------------
c       this subroutine transforms the generated normal field
c       to a lognormal field with specified mean and variance
c-----------------------------------------------------------------------
c
      subroutine tlog(cond)
c
#include "TURN_F.h"

      double precision cond(nx,ny,nz)
      integer i,j,k
c
c.....note:      condg = exp(f) where f = e [ln(k)]
c                sig   = sqrt( var [ln(k)] )
c
      do 100 i=1,nx
      do 100 j=1,ny
      do 100 k=1,nz
      cond(i,j,k)=condg*exp(sig*cond(i,j,k))
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c       this subroutine updates histogram data
c-----------------------------------------------------------------------
c
      subroutine histo(cond,h,mc,delt)
c
#include "TURN_F.h"
c      common / data1 / dx,dy,dz,xl1,xl2,xl3,condg,sig,nx,ny,nz,idate,
c     +                 idnum
c      common / data2 / delzet,bigk,dk,cmin,cmax,nline,nztest,nmont,
c     +                 ndelt,ilog,iu,nl1,nl2,nl3,nl
c      common / internal / pi,rm,rssq,n1,n2,n3,n4,n5,n12,n13,nran,ndelt1
      double precision cond(nx,ny,nz),h(ndelt+1),rl,ru,htest,delt
      integer i,j,k,l,mc
c
c.....update histogram counters
c
444   continue
      rl=cmin
      do 200 l=1,ndelt
      ru=rl+delt
      if(l.eq.ndelt) ru=cmax
      do 210 i=1,nx
      do 210 j=1,ny
      do 210 k=1,nz
      htest=cond(i,j,k)
      if(htest.gt.rl.and.htest.le.ru) h(l+1)=h(l+1)+1.
210   continue
      rl=ru
200   continue
      return
      end
c
c
c-----------------------------------------------------------------------
c       this subroutine prints out the normal or transformed data
c       to unformatted data file(s)
c-----------------------------------------------------------------------
c
      subroutine printout(cond,fileout,iunit)
c
#include "TURN_F.h"
c      common / data1 / dx,dy,dz,xl1,xl2,xl3,condg,sig,nx,ny,nz,idate,
c     +                 idnum
      character*12 fileout
      double precision cond(nx,ny,nz)
      integer iunit,i,j,k
c
c      data linp5/1/
c      data lfi5/5/
c      data ktym/0/
c      data tym/0./
c      data grav1,grav2,grav3/0.,0.,1./
c
c.....open file: 'fileout' as unit: 'iunit'
c
c     open(iunit,file=fileout,access='sequential',
c    +     form='unformatted',status='new')
c
c     write(iunit)idnum,idate,lfi5,linp5,ktym,tym,
c    +                 nx,ny,nz,dx,dy,dz,grav1,grav2,grav3,
c    +                 condg,sig,xl1,xl2,xl3
c     write(iunit) (((cond(i,j,k),i=1,nx),j=1,ny),k=1,nz)
c
c      open(iunit,file=fileout,status='new')
       open(iunit,file=fileout)
       rewind iunit
       write(iunit,*) (((cond(i,j,k),i=1,nx),j=1,ny),k=1,nz)
       close(iunit)

      return

      end

