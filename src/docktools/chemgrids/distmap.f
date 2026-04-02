c------------------------------------------------------------------------
c                          PROGRAM DISTMAP
c
c       Copyright (C) 1991 Regents of the University of California
c                         All Rights Reserved.
c
c  This program creates the grid for contact scoring.
c  The core of the program is Brian's distmap.f, plus....
c  changes by DLB:
c   - new algorithm for increased speed.
c   - grid points assigned a negative score only get a value of -128.
c     (in old version, -128 is modified by the number of "good" contacts)
c  scoring:
c   - grid point gets -128 if it is within polcon angstroms of a polar
c     atom or ccon angstroms of a nonpolar atom, else
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
c------------------------------------------------------------------------
       integer maxlin,maxdix,maxdiy,maxdiz
       parameter (maxlin=5000)
       parameter (maxdix=180,maxdiy=180,maxdiz=180)
c  maxlin: maximum number of atoms in site 
c  maxdix,maxdiy,maxdiz: maximum number of grid points per x/y/z dimension
       character*80 pdbnam, scoren, boxlin
c  pdbnam: name of pdb file containing site to be scored.
c  scoren: name of output file to write the grid.
c  boxlin: name of pdb file of limiting box (optional); allows one to use
c    the total receptor file without forcing the grid to enclose it;
c    only the receptor atoms within the limiting box will be included in
c    the calculation.
       real polcon, ccon, polcon2, ccon2
       real discut,discut2
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
       open(unit=5,file='INDIST',status='old')
       open(unit=6,file='OUTDIST',status='new')
c
       write(6,'(t5,a,$)') 'name of input pdb file:'
       read(5,10) pdbnam
 10    format(a80)
       write(6,*) pdbnam
       write(6,'(t5,a,$)') 'name of output contact grid file:'
       read(5,10) scoren
       write(6,*) scoren
       WRITE(6,*)'BLEEEEEEEEEEEE'
   
c  read input parameters from INDIST.
       write(6,'(t5,a)') 'bad contact distance:'
       write(6,'(t10,a,$)') '1. for polar atoms [2.4]?'
       read(5,*) polcon
       if (polcon .lt. 0.0) polcon = 2.4
       write(6,*)polcon
       write(6,'(t10,a,$)') '2. for non-polar atoms [2.8]?'
       read(5,*) ccon
       if (ccon .lt. 0.0) ccon = 2.8
       write(6,*) ccon
       write(6,'(t5,a,$)') 'scoring distance cutoff [4.5]?'
       read(5,*) discut
       if (discut .lt. 0.0) discut = 4.5
       write(6,*) discut
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
       polcon2=polcon**2
       ccon2 =ccon**2
       discut2=discut**2

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
        call score(grid,ca,pdbcrd,il,atypr,n,maxdix,maxdiy,maxdiz,
     &  maxlin,perang,invgrid,discut2,discut,polcon2,ccon2)

c  dump grid in format of distmap from dockv1.1
       call bdump(grid,n,perang,ca,maxdix,maxdiy,maxdiz)

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
        subroutine score(grid,ca,pdbcrd,il,atypr,n,
     &  maxdix,maxdiy,maxdiz,maxlin,perang,invgrid,discut2,discut,
     &  polcon2,ccon2)

        integer i,j,k,il,maxdix,maxdiy,maxdiz,maxlin
        real ca(2,3),pdbcrd(maxlin,3)
        real dist2,invgrid
        real polcon2,ccon2,distlim2
        real discut,discut2
        real ptpdb(3)
        integer grid(0:maxdix,0:maxdiy,0:maxdiz)
        integer n(3),min(3),max(3),l,perang
        character*4 atypr(maxlin)

c  initialize array grid to all 0's.
        do 120 i = 1,n(1)
        do 110 j = 1,n(2)
        do 100 k = 1,n(3)
          grid(i,j,k) = 0 
 100    continue
 110    continue
 120    continue

        do 1000 l = 1,il
c  find min,max possible grid points for x,y,z
c  any grid point within discut angstroms of an atom could potentially
c  affect the score.
c  note that the larger coordinate (max) has the smaller grid point
c  number.
          max(1) = nint((ca(1,1)- (pdbcrd(l,1)+discut) )*perang)
          max(2) = nint((ca(1,2)- (pdbcrd(l,2)+discut) )*perang)
          max(3) = nint((ca(1,3)- (pdbcrd(l,3)+discut) )*perang)
          min(1) = nint((ca(1,1)- (pdbcrd(l,1)-discut) )*perang)
          min(2) = nint((ca(1,2)- (pdbcrd(l,2)-discut) )*perang)
          min(3) = nint((ca(1,3)- (pdbcrd(l,3)-discut) )*perang)

c  close contact distance is determined by the atom type.
          if(atypr(l)(2:2).eq.'N'.or.atypr(l)(2:2).eq.'O')then
            distlim2 = polcon2
          else
            distlim2 = ccon2
          endif

c  check bounds.
          do 130 i = 1,3
            if (max(i) .lt. 0 ) max(i) = 0
            if (min(i) .gt. n(i) ) min(i) = n(i)
 130      continue

c  loop over all grid points in the square between min and max.
c  for those within distance limits, add 1 to the score.
c  can skip points that already have a negative score.
          do 320 i = max(1),min(1)
          do 310 j = max(2),min(2)
          do 300 k = max(3),min(3)
            if ( grid(i,j,k) .ge. 0) then
              ptpdb(1) = ca(1,1) - i*invgrid 
              ptpdb(2) = ca(1,2) - j*invgrid 
              ptpdb(3) = ca(1,3) - k*invgrid 
              dist2 = (pdbcrd(l,1) - ptpdb(1))**2 +
     &                (pdbcrd(l,2) - ptpdb(2))**2 +
     &                (pdbcrd(l,3) - ptpdb(3))**2 
              if (dist2 .le. distlim2) then
                grid(i,j,k) = -127
              elseif (dist2 .le. discut2) then
                grid(i,j,k) = grid(i,j,k) + 1
              endif
            endif
  300     continue
  310     continue
  320     continue
 1000    continue
        return
        end
c------------------------------------------------------------------------
c   subroutine bdump
c  -writes out grid in distmap format for compatibility with output by 
c   distmap.f.
c------------------------------------------------------------------------
        subroutine bdump(grid,n,perang,ca,maxdix,maxdiy,maxdiz)

        integer nlist,i,j,k,perang,maxdix,maxdiy,maxdiz,neg
        integer grid(0:maxdix,0:maxdiy,0:maxdiz),n(3)
        real ca(2,3)
        logical wlist
        data nlist / 1 /

        neg = 128
        write(16,740)n(1),n(2),n(3),perang,nint(ca(1,1)*perang),
     &         nint(ca(1,2)*perang),nint(ca(1,3)*perang)
 740    format(7(i4,2x))
         do 550 i = 0,n(1)
           do 350 j = 0,n(2)
             do 150 k = 0,n(3)-1
               if(grid(i,j,k).eq.grid(i,j,k+1))then
                  nlist=nlist+1
c  nlist counts the number of consecutive
c  elements of score which have the same value.
                  wlist=.true.
                else
                  if (grid(i,j,k) .ge. 0) then
                    write(16,760)nlist,grid(i,j,k)
                  else
                    write(16,760)nlist,neg
                  endif
 760              format(i3,i3)    
                  nlist=1
                  wlist=.false.
                endif
150           continue
c  logical wlist is used so as not to
c  write the same nlist, score twice.
              if(wlist)then
                  if (grid(i,j,k) .ge. 0) then
                    write(16,760)nlist,grid(i,j,k)
                  else
                    write(16,760)nlist,neg
                  endif
             elseif(grid(i,j,k-1).ne.grid(i,j,k))then
                  if (grid(i,j,k) .ge. 0) then
                    write(16,760)nlist,grid(i,j,k)
                  else
                    write(16,760)nlist,neg
                  endif
              endif                            
              nlist=1
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
c  ignore atoms outside the limiting box
       do 18 i=1,3
         if (crd(i).gt.ca(1,i)) go to 15
         if (crd(i).lt.ca(2,i)) go to 15
 18    continue
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
