c ====================================================================
c                             READGRID()
c  Function to read in chemgrid data files: *.bmp, *.vdw, *.esp, *.abs
c
c  --reads in grids for force field type scoring.          ECM  4/91
c  modified by Todd Ewing 2/93
c
c --------------------------------------------------------------------

      integer function readgrid(grdfil, ngrd, grddiv, grdpts, 
     &        offset, aval, bval, esval, bump)

      integer ngrd
c  ngrd--number of grid points
      real aval(*), bval(*), esval(*)
      character*1 bump(*)
c  aval(), bval(), esval(), bump()--values stored "at" grid points
      real grddiv
c  grddiv--spacing of grid points in angstroms
      real offset(3)
c  offset()--box xmin, ymin, zmin in angstroms
      integer grdpts(3)
c  grdpts()--number of grid points along box dimensions (x,y,z)
c
      character*80 grdfil,append
      character*17 header
      integer i, namend, unitno
      real maxgrd
c
      unitno = 2
      maxgrd = ngrd
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
      open (unit=unitno, file=grdfil(1:namend)//'.bmp', status='old')
      read (unitno, 1) header
      read (unitno, 2) grddiv, (offset(i), i=1,3), (grdpts(i), i=1,3)
      ngrd=grdpts(1)*grdpts(2)*grdpts(3)
      if(ngrd.gt.maxgrd)then
        write(6,*)'too many points in force field grid.'
        write(6,*)'rerun with grid_points at least ', ngrd
        write(6,*)'program stops'
        stop
      endif
      read (unitno, 3) (bump(i), i=1, ngrd)
      close (unitno)

      open (unit=unitno, file=grdfil(1:namend)//'.vdw', status='old',
     &form='unformatted')
      read (unitno) (aval(i), i=1, ngrd)
      read (unitno) (bval(i), i=1, ngrd)
      close (unitno)

      open (unit=unitno, file=grdfil(1:namend)//'.esp', status='old',
     &form='unformatted')
      read (unitno) (esval(i), i=1, ngrd)
      close (unitno)

      return
      end


