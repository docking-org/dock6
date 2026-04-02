c------------------------------------------------------------------------
c                          PROGRAM SOLVMAP
c
c       Copyright (C) 2003 Regents of the University of California
c                         All Rights Reserved.
c
c  This program creates the grid for contact scoring.
c  The core of the program is Brian's distmap.f, plus....
c  changes by DLB:
c   - new algorithm for increased speed.
c   - grid points assigned a negative score only get a value of -128.
c     (in old version, -128 is modified by the number of "good" contacts)
c  scoring:
c   - grid point gets -128 if it is within rad_o(n,c,s,p,q)+probe angstroms 
c     of an approp. atom, else
c   - grid point gets a score of 1 for each atom within dislim angstroms,
c   - else 0.
c   - to make compatible with distmap, added routine bdump to write out
c   the grid in the same format as in dockv1.1
c  changes by ECM:
c   - program ignores hydrogens and deuteriums, if any, in receptor file
c   - unnecessary variables and passing of variables removed
c   - pdb file of box outlining grid, 'distmap.box', written out
c   - user may specify a limiting box (allows use of entire receptor file
c     as input, but inclusion of only the boxed region in the calculation)
c  changes by BKS (2003, who can believe it?)
c   - program now calculates an displaced volume (by solute = protein)
c   - based partial desolvation energy.  Based on calculation of GB
c   - dielectric
c------------------------------------------------------------------------

        include 'dist.h'

       character*80 pdbnam, scoren, boxlin
c  pdbnam: name of pdb file containing site to be scored.
c  scoren: name of output file to write the grid.
c  boxlin: name of pdb file of limiting box (optional); allows one to use
c    the total receptor file without forcing the grid to enclose it;
c    only the receptor atoms within the limiting box will be included in
c    the calculation.
c  polcon: bad contact limit to a polar receptor atom.
c  ccon: bad contact limit to a non-polar receptor atom.
c  polcon2 is polcon squared, ccon2 is ccon squared.
c  discut: cutoff distance for a "good" contact 
c  discut2: discut squared
       integer grid(0:maxdix,0:maxdiy,0:maxdiz)
c  grid: array covering the site of pdbnam containing the scores.
       integer perang
c  perang: number of grid points per angstrom
       real invgrid
c  invgrid: inverse of perang
       integer n(3), il, ierr
c  n: number of grid points in x,y,z directions.
c  il: number of atoms 
c  ierr: whether there is an error opening the indicated limiting box
c        file (0=no error, 1=error)
       real ca(2,3),pdbcrd(maxlin,3)
c  ca: coordinates of the maximum (ca(1,*)) and minimum (ca(2,*))
c      points of the grid.
c  pdbcrd: pdb coordinates of the atoms.
       character*4 atypr(maxlin)
c  atypr: atom type 
       logical limit
c  limit: whether or not the user has specified a limiting box
c--------------------------------------------------------------------------
c  open input and output files
c  these will be indist.dat and outdist.dat on vms machines
       open(unit=5,file='INSOLV',status='old')
       open(unit=6,file='OUTSOLV',status='new')
c
       write(6,'(t5,a,$)') 'name of input pdb file:'
       read(5,10) pdbnam
 10    format(a80)
       write(6,*) pdbnam
       write(6,'(t5,a,$)') 'name of output contact grid file:'
       read(5,10) scoren
       write(6,*) scoren
   
c  read input parameters from INSOLV.
       write(6,'(t5,a)') 'atomic radii:'
       write(6,'(t10,a,$)') '1. for O N C S P Other [1.4 1.3 1.7 2.2'
       write(6,'(a)') '2.2 1.8]?'
       read(5,*) rad_o, rad_n, rad_c, rad_s, rad_p, rad_q
       if (rad_o .le. 0.0) rad_o = 1.4
       if (rad_n .le. 0.0) rad_n = 1.3
       if (rad_c .le. 0.0) rad_c = 1.7
       if (rad_s .le. 0.0) rad_s = 2.2
       if (rad_p .le. 0.0) rad_p = 2.2
       if (rad_q .le. 0.0) rad_q = 1.8
       write(6,*) rad_o, rad_n, rad_c, rad_s, rad_p, rad_q
       write(6,'(t10,a,$)') '2. for probe (exclusion) [1.4]?'
       read(5,*) probe
       if (probe .lt. 0.0) probe = 1.4
       write(6,*) probe
       write(6,'(t5,a,$)') 'number of grid points per angstrom?(1-4)'
       read(5,*) perang
       if (perang .le. 0) perang = 2
       if (perang .gt. 10) perang = 2
       write(6,*)perang
       limit=.false.
       boxlin(1:10)='          '
       read(5,'(a80)',end=21)boxlin
       if (boxlin(1:10).eq.'          ') go to 21
       write(6,*)'pdb file of limiting box:',boxlin
       limit=.true.
 21    continue
       close(5)

c  square for speed in the distance calculations. 
       rad_o = (rad_o+probe)**2
       rad_n = (rad_n+probe)**2
       rad_c = (rad_c+probe)**2
       rad_s = (rad_s+probe)**2
       rad_p = (rad_p+probe)**2
       rad_q = (rad_q+probe)**2
       probe2 = probe**2

       open(14,file=pdbnam,status='old')
       open(16,file=scoren,status='new')

       invgrid = 1.0/real(perang)
       write(6,*) 'grid resolution: ',invgrid
 
       if (limit) then
         ierr=0
         call readbox(15,boxlin,ca,ierr)
         if (ierr.eq.1) then
           limit=.false.
           write(6,'("WARNING -> ",A)') 'error opening limiting box'
           write(6,'("WARNING -> ",A)') 
     &		'   attempting to use entire input pdb file'
           go to 22
         endif
         call readpdb(14,ca,maxlin,pdbcrd,il,atypr)
       endif

 22    continue
       if (.not.limit) then
         call highlow(14,ca,maxlin,pdbcrd,il,atypr)
       endif

c  make pdb-format file of box outlining the grid
       call mkbox(15,ca)
c  calculate number of grid points per dimension
       call gridinfo(perang,n,maxdix,maxdiy,maxdiz,ca)

c  determine score of each grid point
        call score(grid,ca,pdbcrd,il,atypr,n,
     &  perang,invgrid,discut2,discut,polcon2)

         call flush(6)
c  dump grid in format of distmap from dockv1.1
       call bdump(n,perang,ca)

       write(6,*)'successful completion'
       close(6)
       close(16)
       close(14)
       stop
       end
c------------------------------------------------------------------------
c   subroutine highlow
c  -finds the highest and lowest x,y, and z values present in the
c   receptor atom file.  These determine the dimensions of the grid.
c  -stores the coordinates of all atoms in
c   the array pdbcrd and their atom types in atypr.
c------------------------------------------------------------------------
       subroutine highlow(lun,ca,maxlin,pdbcrd,il,atypr)

       integer lun,il,maxlin,i,j
       character*4 atypr(*)
       character*80 string
       real ca(2,3), pdbcrd(maxlin,3)

c  initialize mins to high value, maxs to low value.
       do 10 i = 1,3
         ca(1,i) = -1000.
         ca(2,i) =  1000.
 10    continue

       il = 0
 15    read(lun,'(a80)',end=150) string
       if (string(1:4).ne.'ATOM' .and. string(1:4).ne.'HETA') go to 15
c  ignore hydrogens: position 13 or 14 in AMBER PDB format, position
c  14 in standard PDB format; this will not ignore heavy (non-H,D) atoms
c  in PDB files containing standard atom names (ECM).
       if  (string(13:13).eq.'H'.or.string(13:13).eq.'D'.or.
     & string (14:14).eq.'H'.or. string(14:14).eq.'D') go to 15
       il = il + 1
       if (il .gt. maxlin) then
         write(6,*) 'too many atoms - limit is ', maxlin
       endif
       read(string,80,end=150) atypr(il),(pdbcrd(il,i),i=1,3) 
 80    format(12x,a4,14x,3f8.3)
       do 20 i=1,3 
         if (pdbcrd(il,i) .gt. ca(1,i)) ca(1,i) = pdbcrd(il,i)
         if (pdbcrd(il,i) .lt. ca(2,i)) ca(2,i) = pdbcrd(il,i)
 20    continue 
       go to 15
 150   continue
       write(6,*)'maximum and minimum x y z:'
       do 200 i=1,2
         write(6,*)(ca(i,j),j=1,3)
200    continue
       write(6,*) 'number of nonhydrogen atoms = ',il
       return
       end
c------------------------------------------------------------------------
c   subroutine gridinfo
c  --calculates number of grid points in each dimension.
c------------------------------------------------------------------------
        subroutine gridinfo(perang,n,maxdix,maxdiy,maxdiz,ca)

        real ca(2,3)
        integer perang,maxdix,maxdiy,maxdiz
        integer n(3)

        n(1) = (nint(ca(1,1)*perang)) - (nint(ca(2,1)*perang))
        n(2) = (nint(ca(1,2)*perang)) - (nint(ca(2,2)*perang))
        n(3) = (nint(ca(1,3)*perang)) - (nint(ca(2,3)*perang))
        write(6,*) 'grid x, y, z dimensions:', n(1),n(2),n(3)
        if (n(1).gt.maxdix .or. n(2).gt.maxdiy .or. n(3).gt.maxdiz)then
          write(6,*)'too many grid points - use fewer receptor atoms'
          write(6,*)n(1),maxdix,n(2),maxdiy,n(3),maxdiz
          stop
        endif
        return
        end
c------------------------------------------------------------------------
c   subroutine score
c------------------------------------------------------------------------
        subroutine score(grid,ca,pdbcrd,il,atypr,n,perang,invgrid,
     &  discut2,discut,polcon2)

        include 'dist.h'

        integer i,j,k,il,ii,jj,kk
        real ca(2,3),pdbcrd(maxlin,3)
        real dist2,invgrid
        real polcon2,distlim2
        real discut,discut2
        real ptpdb(3)
        integer grid(0:maxdix,0:maxdiy,0:maxdiz)
        integer n(3),min(3),max(3),l,perang
        character*4 atypr(maxlin)

c  initialize array grid to all 0's.
        do 120 i = 0,n(1)
        do 110 j = 0,n(2)
        do 100 k = 0,n(3)
          grid(i,j,k) = 0 
          solgrid(i,j,k) = 0. 
 100    continue
 110    continue
 120    continue

        dV = invgrid**3/(4*3.14159265) 
        dlim = (1.4)**2
        write(6,*) 'dV is', dV
        call flush (6)

        do 500 l = 1,il
c  find min,max possible grid points for x,y,z
c  any grid point within discut angstroms of an atom could potentially
c  affect the score.
c  note that the larger coordinate (max) has the smaller grid point
c  number.
c  excluded volume distance is determined by the atom type.
          if(atypr(l)(2:2).eq.'O') then
            discut = rad_o
          elseif(atypr(l)(2:2).eq.'N') then
            discut = rad_n
          elseif(atypr(l)(2:2).eq.'C') then
            discut = rad_c
          elseif(atypr(l)(2:2).eq.'S') then
            discut = rad_s
          elseif(atypr(l)(2:2).eq.'P') then
            discut = rad_p
          else 
            discut = rad_q
          endif
          max(1)=nint((ca(1,1)-(pdbcrd(l,1)+sqrt(discut)))*perang)
          max(2)=nint((ca(1,2)-(pdbcrd(l,2)+sqrt(discut)))*perang)
          max(3)=nint((ca(1,3)-(pdbcrd(l,3)+sqrt(discut)))*perang)
          min(1)=nint((ca(1,1)-(pdbcrd(l,1)-sqrt(discut)))*perang)
          min(2)=nint((ca(1,2)-(pdbcrd(l,2)-sqrt(discut)))*perang)
          min(3)=nint((ca(1,3)-(pdbcrd(l,3)-sqrt(discut)))*perang)

c  check bounds.
          do 130 i = 1,3
            if (max(i) .lt. 0 ) max(i) = 0
            if (min(i) .gt. n(i) ) min(i) = n(i)
 130      continue

c  loop over all grid points in the cube between min and max.
c  for those that come too close to protein, set score to -127.
c  can skip points that already have a negative score.
          do 320 i = max(1), min(1)
          do 310 j = max(2), min(2)
          do 300 k = max(3), min(3)
            if ( grid(i,j,k) .ge. 0) then
               ptpdb(1) = ca(1,1) - i*invgrid 
               ptpdb(2) = ca(1,2) - j*invgrid 
               ptpdb(3) = ca(1,3) - k*invgrid 
               dist2 = (pdbcrd(l,1) - ptpdb(1))**2 +
     &              (pdbcrd(l,2) - ptpdb(2))**2 +
     &              (pdbcrd(l,3) - ptpdb(3))**2 
               if (dist2 .le. discut) then
                 grid(i,j,k) = -127
               endif
            endif
  300     continue
  310     continue
  320     continue
  500    continue
         write(6,*) 'end contact grid calc.'
         call flush (6)

          do 560 i = 0, n(1)
          write(6,*) i,'/',n(1),' the way there!'
          call flush (6)
          do 550 j = 0, n(2)
          do 540 k = 0, n(3)
          if(grid(i,j,k) .ge. 0) goto 540
             max(1)=i-(10*perang)
             max(2)=j-(10*perang)
             max(3)=k-(10*perang)
             min(1)=i+(10*perang)
             min(2)=j+(10*perang)
             min(3)=k+(10*perang)

c            check bounds.
             do 505 ii = 1,3
               if (max(ii) .lt. 0 ) max(ii) = 0
               if (min(ii) .gt. n(ii) ) min(ii) = n(ii)
  505        continue

             do 530 ii = max(1), min(1)
             do 520 jj = max(2), min(2)
             do 510 kk = max(3), min(3)
                dist2 = (real(i-ii)*invgrid)**2 + 
     &          (real(j-jj)*invgrid)**2 + (real(k-kk)*invgrid)**2
c      only count grid points outside of the vdw of the probe atom.
*               if (dist2.le.dlim) then
                if (dist2.le.dlim) dist2 = dlim
                solgrid(ii,jj,kk) = dV/(dist2**2) + 
     &                                 solgrid(ii,jj,kk)
  510        continue
  520        continue
  530        continue
  540     continue
  550     continue
  560     continue

        return
        end
c------------------------------------------------------------------------
c   subroutine bdump
c  -writes out grid in distmap format for compatibility with output by 
c   distmap.f.
c------------------------------------------------------------------------
        subroutine bdump(n,perang,ca)

        include 'dist.h'

        integer nlist,i,j,k,perang,neg, kk
        integer n(3)
        real ca(2,3)
        logical wlist
        data nlist / 1 /

        write(16,740)n(1),n(2),n(3),perang,nint(ca(1,1)*perang),
     &         nint(ca(1,2)*perang),nint(ca(1,3)*perang)
 740    format(7(i4,2x))
         do 550 i = 0,n(1)
           do 350 j = 0,n(2)
             do 150 k = 0,n(3),13
               if(k+12.le.n(3))then
                  write(16,761) (solgrid(i,j,kk),kk=k,k+12)
               else
                  write(16,762) (solgrid(i,j,kk),kk=k,n(3)-1)
                  write(16,'(f6.3)') solgrid(i,j,n(3))
               endif
 761              format(13f6.3)    
 762              format(f6.3,$)    
150           continue
350       continue     
550     continue      
        nlist=1
        return
        end
c------------------------------------------------------------------------
      subroutine mkbox(unitno, ca)
c  --makes a pdb format box to show the size and location of the grid.
c------------------------------------------------------------------------
      real boxdim(3), cent(3), ca(2,3)
      integer unitno, i
c
      do 10 i=1,3
        cent(i) = (ca(1,i) + ca(2,i))/2.0
        boxdim(i) = ca(1,i) - ca(2,i)
   10 continue
c
      open (unit=unitno, file='distmap.box', status='unknown')
      write (unitno, '(A24)') 'HEADER    CORNERS OF BOX'
      write (unitno, 1) 'REMARK    CENTER (X Y Z) ', (cent(i), i=1,3)
    1 format (A25, 3F8.3)
      write (unitno, 2) 'REMARK    DIMENSIONS (X Y Z) ',
     &(boxdim(i), i=1,3)
    2 format (A29, 3F8.3)
      write (unitno, 3) 'ATOM', 1, 'DUA', 'BOX', 1,
     &ca(2,1), ca(2,2), ca(2,3)
      write (unitno, 3) 'ATOM', 2, 'DUB', 'BOX', 1,
     &ca(1,1), ca(2,2), ca(2,3)
      write (unitno, 3) 'ATOM', 3, 'DUC', 'BOX', 1,
     &ca(1,1), ca(2,2), ca(1,3)
      write (unitno, 3) 'ATOM', 4, 'DUD', 'BOX', 1,
     &ca(2,1), ca(2,2), ca(1,3)
      write (unitno, 3) 'ATOM', 5, 'DUE', 'BOX', 1,
     &ca(2,1), ca(1,2), ca(2,3)
      write (unitno, 3) 'ATOM', 6, 'DUF', 'BOX', 1,
     &ca(1,1), ca(1,2), ca(2,3)
      write (unitno, 3) 'ATOM', 7, 'DUG', 'BOX', 1,
     &ca(1,1), ca(1,2), ca(1,3)
      write (unitno, 3) 'ATOM', 8, 'DUH', 'BOX', 1,
     &ca(2,1), ca(1,2), ca(1,3)
    3 format (A4, 6x, I1, 2x, A3, 1x, A3, 5x, I1, 4x, 3F8.3)
      write (unitno, 4) 'CONECT', 1, 2, 4, 5
      write (unitno, 4) 'CONECT', 2, 1, 3, 6
      write (unitno, 4) 'CONECT', 3, 2, 4, 7
      write (unitno, 4) 'CONECT', 4, 1, 3, 8
      write (unitno, 4) 'CONECT', 5, 1, 6, 8
      write (unitno, 4) 'CONECT', 6, 2, 5, 7
      write (unitno, 4) 'CONECT', 7, 3, 6, 8
      write (unitno, 4) 'CONECT', 8, 4, 5, 7
    4 format (A6, 4I5)
      close (unitno)
      return
      end
c------------------------------------------------------------------------
      subroutine readbox(unitno,boxlin,ca,ierr)
c  --reads in a pdb file of a limiting box, as written using the program
c    showbox or the program chemgrid.
c------------------------------------------------------------------------
      real ca(2,3)
      integer unitno, i, ierr
      character*80 boxlin,line
c
      ierr=0
      open (unit=unitno, file=boxlin, status='old',err=80)
   10 read (unitno, 1000) line
 1000 format (a80)
      if (line(1:11).eq.'ATOM      1') then
        read (line, '(30x,3f8.3)') (ca(2,i), i=1,3)
        go to 10
      else if (line(1:11).eq.'ATOM      7') then
        read (line, '(30x,3f8.3)') (ca(1,i), i=1,3)
        go to 90
      else
        go to 10
      endif
   80 continue
      ierr=1
   90 continue
      return
      end
c------------------------------------------------------------------------
       subroutine readpdb(lun,ca,maxlin,pdbcrd,il,atypr)
c  --reads in the receptor pdb file; ignores hydrogens and atoms outside
c    of the limiting box.
c  --stores the coordinates in array pdbcrd and atom types in atypr.
c------------------------------------------------------------------------
       integer lun,il,maxlin,i,j
       character*4 atypr(*),atype
       character*80 string
       real ca(2,3), pdbcrd(maxlin,3),crd(3)

       il = 0
 15    read(lun,'(a80)',end=150) string
       if (string(1:4).ne.'ATOM' .and. string(1:4).ne.'HETA') go to 15
c  ignore hydrogens: position 13 or 14 in AMBER PDB format, position
c  14 in standard PDB format; this will not ignore heavy (non-H,D) atoms
c  in PDB files containing standard atom names (ECM).
       if  (string(13:13).eq.'H'.or.string(13:13).eq.'D'.or.
     & string (14:14).eq.'H'.or. string(14:14).eq.'D') go to 15
       read(string,80,end=150) atype,(crd(i),i=1,3) 

cbks   read all pdb atoms in to use in desolvation calculation.
c  ignore atoms outside the limiting box
*      do 18 i=1,3
*        if (crd(i).gt.ca(1,i)) go to 15
*        if (crd(i).lt.ca(2,i)) go to 15
*18    continue
       il = il + 1
       if (il .gt. maxlin) then
         write(6,*) 'too many atoms - limit is ', maxlin
       endif
       read(string,80,end=150) atypr(il),(pdbcrd(il,i),i=1,3) 
 80    format(12x,a4,14x,3f8.3)
       go to 15
 150   continue
       write(6,*)'maximum and minimum x y z:'
       do 200 i=1,2
         write(6,*)(ca(i,j),j=1,3)
200    continue
       write(6,*) 'number of nonhydrogen atoms = ',il
       return
       end
c------------------------------------------------------------------------
