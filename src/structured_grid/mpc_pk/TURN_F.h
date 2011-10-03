      double precision  dx,dy,dz,xl1,xl2,xl3,condg,sig
      double precision  delzet,bigk,dk,cmin,cmax
      double precision  pi,rm,rssq
      integer nx,ny,nz,idate,idnum
      integer nline,nztest,nmont,ndelt,ilog,iu,nl1,nl2,nl3,nl
      integer n1,n2,n3,n4,n5,n12,n13,nk,nran,ndelt1
      character*12 nfile,tfile
      common / data1 / dx,dy,dz,xl1,xl2,xl3,condg,sig,nx,ny,nz,idate,
     &                 idnum
      common / data2 / delzet,bigk,dk,cmin,cmax,nline,nztest,nmont,
     &                 ndelt,ilog,iu,nl1,nl2,nl3,nl,nk
      common / internal / pi,rm,rssq,n1,n2,n3,n4,n5,n12,n13,nran,ndelt1
      common / file / nfile,tfile
