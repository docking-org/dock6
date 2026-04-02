c ====================================================================
c                             READGRID()
c  Function to read in chemgrid data files: *.bmp, *.vdw, *.esp, *.abs
c
c  --reads in grids for force field type scoring.          ECM  4/91
c  modified by Todd Ewing 2/93
c
c --------------------------------------------------------------------

#ifdef BUILD_FOR_RDSOL
      integer function readgrid(grdfil, iflag37, iflagqnifft, ngrd, 
     &        grddiv, grdpts,
     &        offset, dnsize, aval, bval, esval, dslb, dslx, bump,
     &        ndgrd, pscal, dgrd, pmid, phi)
#else
      integer function readgrid(grdfil, iflag37, iflagqnifft, ngrd, 
     &        grddiv, grdpts, 
     &        offset, dnsize, aval, bval, esval,             bump, 
     &        ndgrd, pscal, dgrd, pmid, phi)
#endif

#include "max.h"
#include "phimap.h"

      integer ngrd, ndgrd, igrid
      integer iflag37 !boolean to use grid generated with DOCK 3.7
      logical flag37 !boolean to use grid generated with DOCK 3.7
      integer iflagqnifft !boolean to use grid generated with qnifft (true) or delphi (false)
      logical flagqnifft !boolean to use grid generated with qnifft (true) or delphi (false)
      integer dnsize
      integer tmpint1, tmpint2
      real tmpfloat1(274625), tmpfloat2(274625)
      real tmp2
c  ngrd--number of grid points
      real aval(1000000000), bval(1000000000), esval(1000000000),  !  maxsize
     & phi(1000000000)
#ifdef BUILD_FOR_RDSOL
      real dslb(1000000000), dslx(1000000000)
#endif
      character*1 bump(1000000000)
c  aval(), bval(), esval(), bump()--values stored "at" grid points
      real grddiv, pscal
      real xstart, xend, ystart, yend, zstart, zend
c  grddiv--spacing of grid points in angstroms
      real offset(3), pmid(3), temp(6)
c  offset()--box xmin, ymin, zmin in angstroms
      integer grdpts(3), dgrd(3)
c  grdpts()--number of grid points along box dimensions (x,y,z)
c
      character*80 grdfil,append,tmpfil
      character*17 header
      character*80 tmpval
      character*20 toplabel
      character*10 head 
      character*60 title
      character*16 botlabel
      integer i, j, k, l, namend, unitno, gres
      real x,y,z
      integer il != 1
      integer idiff != 0
      integer maxsize != 1000000000
      !flag37 = LOGICAL(iflag37)
      !flag37 = iflag37
      if (iflag37.eq.1) then
          flag37=.true.
      else if (iflag37.eq.0) then
          flag37=.false.
      else
          write(*,*) "Error... iflag37 is not 0 or 1. "
          stop
      endif
c
      if (iflagqnifft.eq.1) then
          flagqnifft=.true.
      else if (iflagqnifft.eq.0) then
          flagqnifft=.false.
      else
          write(*,*) "Error... iflag37 is not 0 or 1. "
          stop
      endif

c
      nsize = dnsize
      allocate(phimap(dnsize,dnsize,dnsize), stat=istat) 
      il = 1
      idiff = 0
      maxsize = 1000000000
      unitno = 2
c     maxgrid = ngrd
c
      namend=80
      do 100 i=2,80
	if (grdfil(i:i) .eq. ' ') then
	  namend=i-1
	  go to 105
        endif
  100 continue
  105 continue
c
c   1 format (A17)
c   2 format (4F8.3, 3I4)
c   3 format (80A1)
c      character (len=17), parameter :: vdw1 = "(A17)"
c      character (len=45), parameter :: vdw2 = "(4F8.3, 3I4)"

c     if (flag37) then
c       open (unit=unitno, file=grdfil(1:namend)//'.bmp', status='old',
c    &  form='unformatted')
c     else 
c       open (unit=unitno, file=grdfil(1:namend)//'.bmp', status='old')
c     endif
      open (unit=unitno, file=grdfil(1:namend)//'.bmp', status='old')

      read (unitno, "(A17)") header
      read (unitno,"(4F8.3, 3I4)") grddiv, (offset(i), i=1,3), 
     & (grdpts(i), i=1,3)
      ngrd=grdpts(1)*grdpts(2)*grdpts(3)
      if (maxsize.lt.ngrd) then
         write(*,*), "maxsize < ngrd",maxsize,".lt.",ngrd
         Stop
      endif
      if(ngrd.gt.maxgrd)then
        write(6,*)'too many points in force field grid.'
        write(6,*)'rerun with grid_points at least ', ngrd
        write(6,*)'program stops'
        stop
      endif
c     flush(6)
c      read (unitno) (bump(i), i=1, ngrd)
c      read (unitno,"(80A1)") (bump(i), i=1, ngrd-1)
c     do j = 1,int(floor(float(ngrd)/80)) ! read it in line by line
      do j = 1,int(floor(float(ngrd)/80.0)) ! read it in line by line
c        write(*,*), j, il,il+79,ngrd
         read (unitno,"(80A1)") (bump(i), i=il, il+79)
         !read (unitno,"(A1)") temp ! read newline
c        write(*,*) (bump(i), i=il, il+79)
         !il = il+79
         il = il+80
      end do
      !il=il+79
      if ( (j-1) * 80 < ngrd) then ! read in the remainder
         idiff = ngrd - il
c        write(*,*), "I AM HERE", il, j, j*80.0, ngrd, idiff
         read (unitno,"(80A1)") (bump(i), i=il, ngrd-1)
c        write(*,*) (bump(i), i=il, ngrd-1)
      endif
      close (unitno)
      close (14)
 
      tmpfil = grdfil(1:namend)//'.vdw'
      if (flag37) then
        open (unit=unitno, file=tmpfil, status='old',
     &   form='unformatted', action='read',
     &   CONVERT='SWAP')
      else
        open (unit=unitno, file=tmpfil, status='old',
     &  form='unformatted', action='read')

      endif
c     open (unit=unitno, file=tmpfil, status='old',
c    & form='unformatted', action='read')
c    & form='unformatted', action='read',
c    & CONVERT='SWAP')
      read (unitno) (aval(i), i=1, ngrd)
      read (unitno) (bval(i), i=1, ngrd)
      close (unitno)
c for debuging // TEB 2021/10/09
c     do i = 1,ngrd
c         if (aval(i)-bval(i) < -1.0) then
c            write(*,*) "aval = ", aval(i), "bval = ", bval(i)
c         endif
c     end do
 200  continue 

      tmpfil = grdfil(1:namend)//'.esp'
      if (flag37) then
        open (unit=unitno, file=tmpfil, status='old', 
     &   form='unformatted', action='read',
     &   CONVERT='SWAP')
      else
        open (unit=unitno, file=tmpfil, status='old', 
     &   form='unformatted', action='read')
      endif
c     open (unit=unitno, file=tmpfil, status='old', 
c    &form='unformatted')
c    & form='unformatted', action='read')
c    & form='unformatted', action='read',
c    & CONVERT='SWAP')
      read (unitno) (esval(i), i=1, ngrd)
      close (unitno)

#ifdef BUILD_FOR_RDSOL
      tmpfil = grdfil(1:namend)//'.dsl'
      open (unit=unitno, file=tmpfil, status='old',
     &form='unformatted')
      read (unitno) (dslb(i), i=1, ngrd)
      read (unitno) (dslx(i), i=1, ngrd)
      close (unitno)
 300  continue 
#endif
 
c     tmpfil = 'rec+sph.phi'
      tmpfil = grdfil(1:namend)//'.phi'
c     open (unit=unitno, file=tmpfil, status='old',
c    &form='unformatted',action='read',CONVERT='SWAP')
c    &form='unformatted',action='read')
      !read(unitno)toplabel,head,title,phimap,temp,botlabel,scale,oldmid,
      !read(unitno)toplabel,head,title,temp,phimap,botlabel,
      
      if (flagqnifft) then ! use qnifft grid
         open (unit=unitno, file=tmpfil, status='old',
     &   form='unformatted',action='read',CONVERT='SWAP')
         !read(unitno)toplabel,head,title
         read(unitno)toplabel
         read(unitno)head,title
         !write(*,*) 'head: ', head
         !write(*,*) 'title: ', title
         !flush(6)
         read(unitno)phimap
         read(unitno)botlabel
         read(unitno)scale,oldmid
         !write(*,*) toplabel,head,title,botlabel,scale,oldmid,igrid
         !flush(6)
      else ! use DelPhi grid
c        open (unit=unitno, file=tmpfil, status='old',
c    &   form='unformatted',action='read',CONVERT='SWAP')
c    &   form='unformatted',action='read')
         open (unit=unitno, file=tmpfil, status='old',
     &   form='unformatted')
         read(unitno)toplabel
         read(unitno)head,title
         read(unitno)phimap
         read(unitno)botlabel
         read(unitno)scale,oldmid,igrid
         !write(*,*) "toplabel,head,title =", toplabel,head,title
         !write(*,*) "botlabel =", botlabel 
         !write(*,*) "scale,oldmid,igrid,temp =", scale,oldmid,igrid,temp
         !write(*,*) "oldmid =", oldmid(1), oldmid(2), oldmid(3)
         !flush(6)
      endif
      close(unitno)

!      write(*,*) oldmid(1), oldmid(2), oldmid(3)
!      pmid(1) = oldmid(1) - ((igrid-1)/2)*scale
!      pmid(2) = oldmid(2) - ((igrid-1)/2)*scale
!      pmid(3) = oldmid(3) - ((igrid-1)/2)*scale

      pmid(1) = oldmid(1)
      pmid(2) = oldmid(2)
      pmid(3) = oldmid(3)
      pscal = scale

      l = 0
!     print *,'writing phimap'
      do i = 1,dnsize
        do j = 1, dnsize
          do k = 1, dnsize
           l = l+1
           phi(l) = phimap(i,j,k)
           !write(*,*) l, phi(l)
           !stop
          end do
        end do
      end do
      dgrd(1) = dnsize
      dgrd(2) = dnsize
      dgrd(3) = dnsize
      ndgrd = dgrd(1)*dgrd(2)*dgrd(3)

!     crap testing code
!      open(85,file='fgrd.in',status='old')
!       do i = 1,274625
!         read(85,*)tmpfloat1(i)
!       end do
!       close(85)
!
!      open(86,file='cgrd.out',status='old')
!      do j = 1,274625
!         read(86,*)tmpfloat2(j)
!      end do
!      close(86)
!
!      do i=1,274624
!       write(*,'(3f17.8)')tmpfloat1(i), tmpfloat2(i), 
!    &  tmpfloat1(i)-tmpfloat2(i)
!      end do 

!     Write out grids in PDB format

!     X = offset(1)
!     Y = offset(2)
!     Z = offset(3)
!     l = 0
!     gres = 0
!     tmp2 = 1.0

!      open(45,file='grd-aval.pdb',status='new')

!     do i = 1, nsize
!       if(i > 1) Z = Z + grddiv
!       do j = 1, nsize
!         if(j > 1) Y = Y + grddiv
!         do k = 1, nsize
!           l = l + 1
!           gres = gres + 1
!           !if(phi(l).lt.-0.5) then
!           write(45,'(a6,i5,2x,a2,2x,a3,i6,4x,3f8.3,2f6.2)')'HETATM',l,
!    &'O2','WAT',gres,X,Y,Z,tmp2,aval(l)
!	      write(45,'(a3)')'TER'
!           !endif
!           X = X + grddiv
!         end do
!         X = offset(1)
!       end do
!       X = offset(1)
!      Y = offset(2)
!     end do

!     close(45)

      
      return
      end
