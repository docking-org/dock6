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
c**********modified********
      include 'nchemgrid.h'
c***********ended**********
c
      character*80 grdfil
c**********modified**************
      integer i, namend, unitno, n, k
c***********ended****************
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
    4 format (3I5, 3F9.3, I3)
      open (unit=unitno, file=grdfil(1:namend)//'.bmp', 
     &status='unknown')
      write (unitno, 1) 'bump map         '
      write (unitno, 2) grddiv, (offset(i), i=1,3), (grdpts(i), i=1,3)
      write (unitno, 4) grdcut, grdcuto, natm_total, cutsq, f_scale,
     *                  solv_rec, nvtyp
      write (unitno, 3) (bump(i), i=1, npts)
      close (unitno)

      if (estype.ne.2) then
      open (unit=unitno, file=grdfil(1:namend)//'.vdw', 
     &status='unknown', form='unformatted')
      write (unitno) (aval(i), i=1, npts)
      write (unitno) (bval(i), i=1, npts)
      close (unitno)
      open (unit=unitno, file=grdfil(1:namend)//'.esp',
     &status='unknown', form='unformatted')
      write (unitno) (esval(i), i=1, npts)
      close (unitno)
      else
c********************** added ***************************************
c      open (unit=unitno, file=grdfil(1:namend)//'.sol',
c     &status='unknown', form='unformatted')
c      write (unitno) (mark_occup(i), i=1, npts)
c      close (unitno)
c      open (unit=unitno, file=grdfil(1:namend)//'.rec',
c     &status='unknown', form='unformatted')
ccc   write (unitno) ((rcrd_rec(i,n), n=1, 3), i=1, natm_total)
c      write (unitno) (crd_rec(i,1), i=1, natm_total)
c      write (unitno) (crd_rec(i,2), i=1, natm_total)
c      write (unitno) (crd_rec(i,3), i=1, natm_total)
c      write (unitno) (rcrg_rec(i), i=1, natm_total)
c      write (unitno) (rvdweff(i), i=1, natm_total)
c      write (unitno) (sum_rec(i), i=1, natm_total)

c srb 7/2009, from unformatted to formatted for 64 bit platforms
      open (unit=unitno, file=grdfil(1:namend)//'.sas',
     &status='unknown' )
   66 format (4F16.8,1I12)
      write (unitno, 66) ((crd_rec(i,k), k=1,3), rvdweff(i),vdwn_rec(i),
     &i=1, natm_total)
   67 format (2I12)
      write (unitno, 67) (nsphgrid(i), nhp(i), i=1, nvtyp)
      close (unitno)

    5 format (3F9.4)
    6 format (6F10.6)
    7 format (40I2)
    8 format (3I8)
      open (unit=unitno, file=grdfil(1:namend)//'.sasmark', 
     &status='unknown')
      write (unitno, 8) nsas, nsas_hp, nsas_pol
      write (unitno, 5) r_probe, spacing, r2_cutoff
      write (unitno, 6) (((sphgrid_crd(n,i,k), k=1,3), 
     &i=1, nsphgrid(n)), n=1, nvtyp)
      write (unitno, 7) ((mark_sas(i,n), n=1, nsphgrid(vdwn_rec(i))), 
     &i=1, natm_total)
      close (unitno)

      end if
c********************** ended ***************************************
c
      return
      end
c----------------------------------------------------------------------
