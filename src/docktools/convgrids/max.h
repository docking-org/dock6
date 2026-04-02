c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c                      distribution version
c  --header file containing all nonlocal parameters used in dock3
c                                                      ECM   3/92
c-----------------------------------------------------------------------
c
c GRP = unique levels and sublevels of molecule
c MATM = # of atoms in molecule, EATM = # of atoms in ensemble
      integer MAXMGRP, MAXEGRP, MAXMATM
      integer MAXEATM, MAXMBRANCH, MAXOR, MAXMOL, MAXFBRANCH
      parameter (MAXMGRP = 1000, MAXEGRP = 25000)
      parameter (MAXMATM = 4100, MAXEATM = 750000)
      parameter (MAXMBRANCH = 101, MAXFBRANCH = 1000)
      parameter (MAXMOL = 100, MAXOR = 100000)
c SAR_PICK = which orientation to write for SAR DOCK
c introduced by jji in summer 2002. currently unused.
c see save_branch_ptrs.f
      integer SAR_PICK
      parameter (SAR_PICK = 6262)
c
      integer numconfs
c      integer maxlig
c      parameter (maxlig=100)
c maxlig was 200
c  maxlig--maximum number of ligand atoms
      integer maxlst
      parameter (maxlst=10000)
c  maxlst--maximum number of ligands saved (per list)
      integer maxgrd
      !parameter (maxgrd=2500000)
      parameter (maxgrd=10000000)
c  maxgrd--maximum number of points in grids for force field scoring
      integer maxtyv
      parameter (maxtyv=50)
c  maxtyv--maximum number of entries in van der Waals parameter file
      integer maxclu
      parameter (maxclu=10)
c  maxclu--maximum number of sphere clusters that can be docked to in a
c   given run
      integer maxpts,maxwid,maxnod
      parameter (maxpts=440,maxwid=70,maxnod=60)
c  maxpts--maximum number of ligand or receptor centers
c  maxwid--maximum number of items in a bin (number of nodes + 2)
c  maxnod--maximum number of bins per ligand or receptor center
      integer nsize
      !parameter (nsize=65)
      !parameter (nsize=100)
      !parameter (nsize=193)
      !parameter (nsize=194)
c  nsize--points per side of DelPhi electrostatic potential grid
      integer maxdix,maxdiy,maxdiz
      parameter (maxdix=200, maxdiy=200, maxdiz=200)
c  maxdim--maximum points per side of grid for contact scoring
      integer maxctr
      parameter (maxctr=12)
c  maxctr--maximum number of ligand/receptor pairs that the program will
c   attempt to find; maximum value of nodlim
      real inithi, initlo
      parameter (inithi = 1.0E+10, initlo = -127.0)
c  inithi--initialization value for electrostatic and ff scores
c  initlo--initialization value for contact scores
c-----------------------------------------------------------------------
