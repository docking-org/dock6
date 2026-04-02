c-------------------------------------------------------------------------
c                        distribution version
c      header file for CHEMGRID                             ECMeng  4/91
c-------------------------------------------------------------------------
      integer maxpts
c      parameter (maxpts=16700000)
      parameter (maxpts=5000000)
c  maxpts--maximum number of grid points
      integer npts
c  npts--number of grid points
      real aval(maxpts), bval(maxpts), esval(maxpts)
      character*1 bump(maxpts)
c  aval(), bval(), esval(), bump()--values stored "at" grid points
c***********************added*******************
      integer mark_occup(maxpts)
c  values stored "at" grid points for solvation
      integer natm_total
c total # of receptor atoms
c************************end********************
      real rsra, rsrb, rcrg, rcrd(3)
      integer natm, vdwn
c  rsra, rsrb, rcrg, rcrd(), natm, vdwn--values for current receptor atom
c***************************************************************************
c add in original, offset and effective vdw radius for current receptor atom
      real rvdw, vdwoff, rvdweff0
      integer maxatm
      parameter (maxatm=6000)
      parameter (vdwoff=-0.09)
      real rvdweff(maxatm), rcrg_rec(maxatm)
      real crd_rec(maxatm,3)
      real sum_rec(maxatm), inv_a(maxatm)
      real cutoffo, cutoffb 
c  cutoffo, cutoffb --cutoff distance for overlap marking and margin outside box
      integer grdcuto, grdcutb
c  grdcuto and grdcutb--cutoffo and cutoffb in grid units
      real solv_rec
c  solv_rec--solvation parameter of receptor
      real f_scale
c  f_scale--scale factor in the exponent term of f_gb (originally set to be 4.0)
      real r_probe
c  r_probe--probe radius
      integer maxnear, maxgrd, maxtyv
      parameter (maxnear=100, maxgrd=5000, maxtyv=30)
c  maxnear--maximum number of neighboring atoms
c  maxgrd--maximum number of spherical grid pts on each atom 
c  maxtyv--maximum number of vdw type
      integer vdwn_rec(maxatm)
c  vdwn_rec()--mark of vdw type of receptor atoms
      integer mark_sas(maxatm, maxgrd)
c  mark_sas()--mark of grid pts. on SAS 
      real sphgrid_crd(maxtyv,maxgrd,3)
c  sphgrid_crd()--coord. of grid pts. on each atom
      integer nsphgrid(maxtyv)
c  nsphgrid()--number of spherical grid pts on each atom
      integer nhp(maxtyv)
c  nhp()--type of polarity of each VDW type
      real r_cutoff, r2_cutoff
c  r_cutoff--cutoff of inter-atom distance for SAS evaluation
c  r2_cutoff--square of r_cut
      real spacing
c  spacing--square root of area per grid point
      integer nsas, nsas_hp, nsas_pol
c  nsas's--total # of grid pts on SAS (of different polarity types)
      integer nvtyp
c  nvtyp--number of entries in 'vdw.parms' with nonzero radius
c***************************************************************************
      integer nearpt(3)
c  nearpt()--3D indices of grid point closest to current receptor atom
      real gcrd(3)
c  gcrd()--coordinates in angstroms of current grid point
      real grddiv
c  grddiv--spacing of grid points in angstroms
      real boxdim(3)
c  boxdim()--box dimensions in angstroms (x,y,z)
      real offset(3)
c  offset()--box xmin, ymin, zmin in angstroms
      integer grddim(3)
c  grddim()--box dimensions in grid units (x,y,z)
      integer grdpts(3)
c  grdpts()--number of grid points along box dimensions (x,y,z)
c    NOTE: grdpts(i)=griddim(i) + 1   (lowest indices are (1,1,1))
      integer estype
c  estype--type of electrostatic calculation desired:
c    0 = use constant dielectric function
c    1 = use distance-dependent dielectric function
c    2 = use previously generated (DelPhi) electrostatic potential map
      real esfact
c  esfact--factor to multiply dielectric by when estype = 0 or estype = 1;
c    not read or used when estype = 2
c    examples:  D = 1    estype = 0, esfact = 1
c               D = 4    estype = 0, esfact = 4
c               D = r    estype = 1, esfact = 1
c               D = 4r   estype = 1, esfact = 4
      real cutoff, cutsq, pcon, ccon, pconsq, cconsq
c  cutoff--cutoff distance for energy calculations
c  cutsq--cutoff distance squared
c  pcon (ccon)--distance defining a bump with a polar atom (a carbon)
c    of the receptor
c  pconsq, cconsq--the squares of pcon and ccon, respectively
      integer grdcut
c  grdcut--cutoff, in grid units
      real dist2
c  dist2--function to calculate distance squared
      integer indx1
c  indx1--function to convert the 3-dimensional (virtual) indices of a
c    grid point to the actual index in a 1-dimensional array

c
      common
     &/rmaps/ aval, bval, esval
c********************* added *********************************
     &/smaps/ rcrg_rec, rvdweff, sum_rec, crd_rec, mark_occup
     &/vals/ cutsq, pconsq, cconsq, grdcutb, f_scale, solv_rec
     &/solv/ r_probe, sphgrid_crd, nsphgrid, nhp, nvtyp
     &/solv1/ spacing, r2_cutoff, vdwn_rec, mark_sas
     &/solv2/ nsas, nsas_hp, nsas_pol
c********************* ended *********************************
     &/cmap/ bump
     &/tmap/ boxdim
c-------------------------------------------------------------------------
