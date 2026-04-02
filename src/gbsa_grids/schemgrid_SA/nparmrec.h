c----------------------------------------------------------------
c                  distribution version
c     header file for subroutine parmrec          ECMeng   4/91
c----------------------------------------------------------------
c Modified by Zou, 11/96
      integer maxtyp, nptyp
      parameter (maxtyp=1000)
c  maxtyp--maximum number of entries in 'prot.table' or 'na.table'
c  nptyp--number of entries in 'prot.table' or 'na.table' so far
      integer inum(maxtyp), ilink(maxtyp)
c  inum()--id numbers in hash table
c  ilink()--links for hash table
      character*1 chain(maxtyp), chn, schn
      character*3 res(maxtyp), resid, sresid
      character*4 atm(maxtyp), resnum(maxtyp), atom, resno, sresno
      real crg(maxtyp)
      integer vdwtyp(maxtyp)
c  vdwtyp()--integer vdw type indicators
      integer maxtyv
c      parameter (maxtyv=50)
      parameter (maxtyv=30)
c  maxtyv--maximum number of entries in 'vdw.parms'
      integer nvtyp
c  nvtyp--number of entries in 'vdw.parms' so far
c************************************************
c add in vdw(maxtyv) (VDW radius of each atom) and polarity type nhp(maxtyv)
      real sra(maxtyv), srb(maxtyv), vdw(maxtyv)
      integer nhp(maxtyv)
      integer maxgrd
      parameter (maxgrd=5000)
      real r_probe, sphgrid_crd(maxtyv,maxgrd,3)
      integer nsphgrid(maxtyv)
      real sgcrd(maxgrd,3)
c cardial coordinates in angstroms of grid points relative to atomic center
c************************************************
c  sra(), srb()--vdw parameters, sqrt of A and sqrt of B
      integer maxatm
      parameter (maxatm=30000)
c  maxatm--maximum number of receptor atoms
      logical found
      character*80 line
      real crgtot
c
      common
     &/chem_link/ inum, ilink
     &/name/ atm, res, resnum, chain
c     &/value/ crg, vdwtyp, sra, srb, vdw
c**************add*******************************
     &/solv/ r_probe, sphgrid_crd, nsphgrid, nhp, nvtyp
c************************************************
c----------------------------------------------------------------
