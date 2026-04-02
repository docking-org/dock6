c----------------------------------------------------------------------
      subroutine grdout(grdfil, unitno, npts, grddiv, grdpts, offset)
c
c  --called from CHEMGRID
c  --writes out grids; makes a formatted "bump" file and unformatted
c    van der Waals and electrostatics files
c    explicit an bulk desolvation grids for the receptor  kxr 12/05
c                                                  ECMeng    4/91
c----------------------------------------------------------------------
      include 'chemgrid.h'
c
      character*80 grdfil
      integer i, namend, unitno
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
      open (unit=unitno, file=grdfil(1:namend)//'.bmp', status='new', 
     &form='unformatted')
      write (unitno) 'bump map         '
      write (unitno) grddiv, (offset(i), i=1,3), (grdpts(i), i=1,3)
      write (unitno) (bump(i), i=1, npts)
      close (unitno)
      open (unit=unitno, file=grdfil(1:namend)//'.vdw', status='new',
     &form='unformatted')
      write (unitno) (aval(i), i=1, npts)
      write (unitno) (bval(i), i=1, npts)
      close (unitno)
      open (unit=unitno, file=grdfil(1:namend)//'.esp', status='new',
     &form='unformatted')
      write (unitno) (esval(i), i=1, npts)
      close (unitno)
      if(idsol .eq. 1) then 
c       write out bulk and explicit desolvation grids
        open (unit=unitno, file=grdfil(1:namend)//'.dsl', status='new',
     &  form='unformatted') 
        write (unitno) (dsol_bl(i), i=1, npts)
        write (unitno) (dsol_ex(i), i=1, npts)
        close (unitno)
      else if(idsol .eq. 2) then 
c       write out delphi desolvation grid
        open (unit=unitno, file=grdfil(1:namend)//'.dsl', status='new',
     &  form='unformatted')
        write (unitno) (dsol_bl(i), i=1, npts)
      else if(idsol .eq. 4) then 
c       write out gaussian scoring grid
        open (unit=unitno, file=grdfil(1:namend)//'.gsf', status='new',
     &  form='unformatted')
        write (unitno) (volp(i), i=1,npts)
        write (unitno) (ndf(i), i=1,npts)
        close(unitno)
      endif
c
      return
      end
c----------------------------------------------------------------------
