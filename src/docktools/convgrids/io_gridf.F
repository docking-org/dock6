c ====================================================================
c                             READGRID()
c  Function to read in chemgrid data files: *.bmp, *.vdw, *.esp, *.abs
c
c  --reads in grids for force field type scoring.          ECM  4/91
c  modified by Todd Ewing 2/93
c
c --------------------------------------------------------------------

#ifdef BUILD_FOR_RDSOL
      integer function readgrid(grdfil, ngrd, grddiv, grdpts,
     &        offset, aval, bval, esval, dslb, dslx, bump, ndgrd,
     &        pscal, dgrd, pmid, phi)
#else
      integer function readgrid(grdfil, ngrd, grddiv, grdpts, 
     &        offset, aval, bval, esval, bump, ndgrd, pscal,
     &        dgrd, pmid, phi)
#endif

#include "max.h"
#include "phimap.h"

      integer ngrd, ndgrd, igrid
      integer tmpint1, tmpint2
      real tmpfloat1(274625), tmpfloat2(274625)
      real tmp2
c  ngrd--number of grid points
      real aval(*), bval(*), esval(*), phi(*)
#ifdef BUILD_FOR_RDSOL
      real dslb(*), dslx(*)
#endif
      character*1 bump(*)
c  aval(), bval(), esval(), bump()--values stored "at" grid points
      real grddiv, pscal
      real xstart, xend, ystart, yend, zstart, zend
c  grddiv--spacing of grid points in angstroms
      real offset(3), pmid(3)
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
c
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
    1 format (A17)
    2 format (4F8.3, 3I4)
    3 format (80A1)
      open (unit=unitno, file=grdfil(1:namend)//'.bmp', status='old',
     &form='unformatted')
      read (unitno) header
      read (unitno) grddiv, (offset(i), i=1,3), (grdpts(i), i=1,3)
      ngrd=grdpts(1)*grdpts(2)*grdpts(3)
      write(6,*) ngrd, grddiv, (grdpts(i), i=1,3), (offset(i), i=1,3)
      write(6,*) 'maxgrid ', maxgrd
      if(ngrd.gt.maxgrd)then
        write(6,*)'too many points in force field grid.'
        write(6,*)'rerun with grid_points at least ', ngrd
        write(6,*)'program stops'
        stop
      endif
      read (unitno) (bump(i), i=1, ngrd)
      close (unitno)
      close (14)

  
      tmpfil = grdfil(1:namend)//'.vdw'
      open (unit=unitno, file=tmpfil, status='old',
     &form='unformatted')
      read (unitno) (aval(i), i=1, ngrd)
      read (unitno) (bval(i), i=1, ngrd)
      close (unitno)
 200  continue 

      tmpfil = grdfil(1:namend)//'.esp'
      open (unit=unitno, file=tmpfil, status='old', 
     &form='unformatted')
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
 
      tmpfil = 'rec+sph.phi'
      open (unit=unitno, file=tmpfil, status='old',
     &form='unformatted')
      read(unitno)toplabel
      read(unitno)head,title
      read(unitno)phimap
      read(unitno)botlabel
      read(unitno)scale,oldmid,igrid
      close(unitno)

!      write(*,*) oldmid(1), oldmid(2), oldmid(3)
!     pmid(1) = oldmid(1) - ((igrid-1)/2)*scale
!	pmid(2) = oldmid(2) - ((igrid-1)/2)*scale
!	pmid(3) = oldmid(3) - ((igrid-1)/2)*scale
      pmid(1) = oldmid(1)
      pmid(2) = oldmid(2)
      pmid(3) = oldmid(3)
      pscal = scale

      l = 0
!     print *,'writing phimap'
      do i = 1,nsize
        do j = 1, nsize
          do k = 1, nsize
           l = l+1
           phi(l) = phimap(i,j,k)
           !write(*,*) l, phi(l)
          end do
        end do
      end do
      dgrd(1) = nsize
      dgrd(2) = nsize
      dgrd(3) = nsize
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
