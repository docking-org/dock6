c---------------------------------------------------------------------
c                          PROGRAM CHEMGRID
c
c       Copyright (C) 1991 Regents of the University of California
c                         All Rights Reserved.
c
c  --calculates the grids used in DOCK 3.0 for force field type scoring,
c  given a receptor pdb file with hydrogens on the polar atoms, and the
c  appropriate parameter tables.  4 grids are created and stored in 3
c  files:  a bump grid/file, an electrostatic potential grid/file and
c  van der Waals repulsion and attraction grids stored in the same 
c  file.
c  ECM  9/93  --version to read input box from SHOWBOX; no output box
c
c  --informative output files:
c    OUTCHEM    restatement of input parameters; messages pertaining
c               to calculation of the grids
c    OUTPARM    messages pertaining to parameterization of receptor
c               atoms; net charge on the receptor molecule (including
c               any ions or waters in the receptor pdb file); created
c               in subroutine parmrec
c    PDBPARM    shows which parameters have been associated with each
c               atom in the receptor pdb file; created in subroutine
c               parmrec, and used in subroutines dconst (constant
c               dielectric) and ddist (distance-dependent dielectric)
c---------------------------------------------------------------------
c zou, 7/96 
      include 'nchemgrid.h'
c
      character*80 vdwfil
      real scrd(3), sumcrd(3), com(3)
c  scrd()--coordinates of a sphere center
c  sumcrd()--sum of sphere center coordinates so far
c  com()--grid center coordinates; user input or sphere cluster center
c   of mass
c zou, 8/98
      character*80 recfil, sphfil, boxfil, grdfil, cavfil
c  recfil--pdb-format receptor file (input file)
c  sphfil--file containing one or more sphere clusters, read when the
c   user wishes to center the grids on a sphere cluster center of mass
c   (input file)
c  boxfil--pdb-format file for displaying the grid boundaries (output
c   file)
c  grdfil--prefix name for grid files (output)
      character*80 table, dumlin
c  table--table containing receptor atom parameters
      character*1 ctrtyp
c  ctrtyp--character flag for grid center option; 'u' or 'U' for user
c   input, otherwise a sphere cluster center of mass
      integer i, j, n
c  variables for reading sphere coordinates:
      integer cntemp, nstemp, clnum, nsph
      logical done
c  done--whether or not sphere cluster center of mass has been calculated
c zou, 8/96 
      real atime,btime
c zou, 7/98 
      real epsilon_solvent, epsilon_cavity
c
      open (unit=1, file='INCHEM', status='old')
      open (unit=2, file='OUTCHEM', status='unknown')
c

c zou, 7/98 
      write (2, *) '***From SCHEMGRID, Version 4.0 ***'
      write (2, *) '__________________________________'
      write (2, *) ' '


      read (1, 1000) recfil
c zou, 7/98
 1000 format (A50)
      write (2, *) 'receptor pdb file:'
      write (2, 1000) recfil
      read (1, 1000) cavfil
      write (2, *) 'cavity pdb file:'
      write (2, 1000) cavfil
      read (1, 1000) table
      write (2, *) 'receptor parameters will be read from:'
      write (2, 1000) table
      read (1, 1000) vdwfil
      write (2, *) 'van der Waals parameter file:'
      write (2, 1000) vdwfil
c  
      call parmrec(recfil, cavfil, table, vdwfil, 2)
c  
      read (1, 1000) boxfil
      write (2, *) 'input box file defining grid location:'
      write (2, 1000) boxfil
c
      open (unit=3, file=boxfil, status='old')
      read (3, 1000) dumlin
      read (3, 1001) (com(i), i=1,3)
 1001 format (25x, 3f8.3)
      read (3, 1002) (boxdim(i), i=1,3)
 1002 format (29x, 3f8.3)
      close (3)
c
      write (2, '(A)') 'box center coordinates [x y z]:'
      write (2, '(3G13.7)') (com(i), i=1,3)
      write (2, '(A,G13.7)') 'box x-dimension = ', boxdim(1)
      write (2, '(A,G13.7)') 'box y-dimension = ', boxdim(2)
      write (2, '(A,G13.7)') 'box z-dimension = ', boxdim(3)
c
c  --set offset to xmin, ymin, zmin of box
c
      do 65 i=1,3
        offset(i)=com(i) - boxdim(i)/2.0
   65 continue
c
      read (1, *) grddiv
      write (2, *) 'grid spacing in angstroms'
      write (2, *) grddiv
c
c  --convert box dimensions to grid units, rounding upwards
c  --note that points per side .ne. side length in grid units,
c    because lowest indices are (1,1,1) and not (0,0,0)
c
      npts=1
      do 70 i=1,3
        grddim(i)=int(boxdim(i)/grddiv + 1.0)
        grdpts(i)=grddim(i) + 1
        npts=npts*grdpts(i)
   70 continue
      if (npts .gt. maxpts) then
        write (2, *) 'maximum number of grid points exceeded--'
        write (2, *) 'decrease box size, increase grid spacing, or'
        write (2, *) 'increase parameter maxpts'
        write (2, *) 'program stops'
        stop
      endif
      write (2, *) 'grid points per side [x y z]:'
      write (2, *) (grdpts(i), i=1,3)
      write (2, *) 'total number of grid points = ', npts

      read (1, *) estype
c zou, 7/96
c 0: constant dielectric; 1: distance-dependent dielectric;
c 2: Generalized Born polarization
      if (estype .eq. 0) then
        write (2, *) 'a constant dielectric will be used'
      else if (estype .eq. 1) then
        write (2, *) 'a distance-dependent dielectric will be used'
      else
        write (2, *) 'Generalized Born polarization will be used'
      endif

      read (1, *) esfact
      write (2, '(A31, A15, F6.2)') 'the dielectric function will be',
     &' multiplied by ', esfact
   75 continue

      read (1, *) cutoff, cutoffb
      write (2, *) 'cutoff distance for energy calculations:'
      write (2, *) cutoff
      write (2, *) 'cutoff distance for receptors outside the box:'
      write (2, *) cutoffb
      cutsq=cutoff*cutoff
      grdcut=int(cutoff/grddiv + 1.0)
      grdcutb=int(cutoffb/grddiv + 1.0)

c      read (1, *) cutoffo
      cutoffo = 3.0
      write (2, *) 'cutoff distance for overlap marking:'
      write (2, *) cutoffo
      grdcuto=int(cutoffo/grddiv + 1.0)
c
c  --convert cutoff to grid units, rounding up (only add 1 rather
c    than 2, because differences in indices rather than the 
c    absolute indices are required)
c


c zou, 10/96 
c      read (1, *) f_scale
      f_scale = 4.0
      write (2, *) 'scale factor in exponent term of f_gb:'
      write (2, *) f_scale


c zou, 7/98 
      read (1, *) epsilon_solvent, epsilon_cavity
      cav_scale = (1.0 / epsilon_cavity - 1.0 / epsilon_solvent)
     +            /(1.0 -1.0 / epsilon_solvent)
      write (2, *) 'dielectric contants of solvent and cavity:'
      write (2, *) epsilon_solvent, epsilon_cavity
      write (2, *) 'cav_scale:'
      write (2, *) cav_scale


      read (1, *) pcon, ccon
      write (2, *) 'distances defining bumps with receptor atoms:'
      write (2, '(A21, F5.2)') 'receptor polar atoms ', pcon
      write (2, '(A22, F5.2)') 'receptor carbon atoms ', ccon
      pconsq=pcon*pcon
      cconsq=ccon*ccon

      read (1, 1000) grdfil
      write (2, *) 'output grid prefix name:'
      write (2, 1000) grdfil

c zou, 7/99
      read (1, *) mark_pair
      if (mark_pair.eq.1) then
        write (2, *) 'Use the pairwise formula to compute Born radii'
      else
        write (2, *) 'Use the original formula to compute Born radii'
      end if


      close (1)
c
c  --initialize grid
c
      do 80 n=1, npts
        aval(n)=0.0
        bval(n)=0.0
        esval(n)=0.0
        bump(n)='F'
c zou, 7/96
        mark_occup(n)=0
   80 continue
c

c zou, 7/96
      call second(atime)

      write (2, *) '***elapsed time before grid computation: ',
     &              atime

      if (estype .eq. 0) then
        call dconst(3, grdcut, grddiv, grdpts, esfact, offset) 
      else if (estype .eq. 1) then
        call ddist(3, grdcut, grddiv, grdpts, esfact, offset)
      else
        call polariz(3, grdcut, grddiv, grdpts, offset, 
     &grdcuto, natm_total)
      endif

      call second(btime)
      write (2, *) '***elapsed time on grid computation: ',
     &              btime-atime
c
      call grdout(grdfil, 3, npts, grddiv, grdpts, offset, grdcut,
     &grdcuto, natm_total, estype)
      call second(btime)
      write (2, *) '***total elapsed time: ', btime

c
      close (2)
      end
c---------------------------------------------------------------------
