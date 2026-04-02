c----------------------------------------------------------------------
c*********************** modified ************************************
      subroutine grdout(grdfil, unitno, npts, grddiv, grdpts, offset,
     &grdcut, grdcuto, natm_total, estype)
c************************* ended *************************************
c
c  --called from CHEMGRID
c  --writes out grids; makes a formatted "bump" file and unformatted
c    van der Waals and electrostatics files
c                                                  ECMeng    4/91
c----------------------------------------------------------------------
c zou, 7/96
      include 'nchemgrid.h'

c
      character*80 grdfil

c zou, 7/96
      integer i, k, namend, unitno,n

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
    4 format (3I5, 2F9.3, F11.3)
      open (unit=unitno, file=grdfil(1:namend)//'.bmp', 
     &status='unknown')
      write (unitno, 1) 'bump map         '
      write (unitno, 2) grddiv, (offset(i), i=1,3), (grdpts(i), i=1,3)

c zou, 7/96
      write (unitno, 4) grdcut, grdcuto, natm_total, cutsq, f_scale,
     *                  solv_rec
      write (unitno, 3) (bump(i), i=1, npts)
      close (unitno)

      open (unit=unitno, file=grdfil(1:namend)//'.vdw',
     &status='unknown', form='unformatted')
      write (unitno) (aval(i), i=1, npts)
      write (unitno) (bval(i), i=1, npts)
      close (unitno)

      if (estype.ne.2) then
      open (unit=unitno, file=grdfil(1:namend)//'.esp',
     &status='unknown', form='unformatted')
      write (unitno) (esval(i), i=1, npts)
      close (unitno)
      else


c zou, 7/96
      open (unit=unitno, file=grdfil(1:namend)//'.sol',
     &status='unknown', form='unformatted')
      write (unitno) (mark_occup(i), i=1, npts)
      close (unitno)

c srb 7/2009, from unformatted to formatted for 64 bit platforms
      open (unit=unitno, file=grdfil(1:namend)//'.rec',
     &status='unknown' )
   66 format (7F16.8)
      write (unitno, 66) ((crd_rec(i,k), k=1,3), rcrg_rec(i),
     & rvdweff(i), inv_a(i), descreen(i), i=1, natm_total)

c zou, 7/99
      close (unitno)
      end if
c

      return
      end
c----------------------------------------------------------------------
