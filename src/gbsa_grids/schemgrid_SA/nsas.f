c**************************** NEW *************************************
      subroutine sas(unitno, grddiv, grdpts, offset, grdcuto, 
     &natm_total)
c
c  --called from CHEMGRID
c  --preparation for hydrophobic energy calculation
c----------------------------------------------------------------------
c
      include 'nchemgrid.h'
c
      integer unitno, i, j, k, n, natm_all, natm1
      integer ilo, ihi(3)
      real rcrd0(3), rcrd1(3), r2, rij2
      integer mark_in(maxatm), vdwnG(maxatm)
      integer maxneighb
      parameter (maxneighb=300)
      integer neighbor_num(maxatm), neighbors(maxatm, maxneighb)
      real reff2G(maxatm), crd_recG(maxatm,3)
      integer nsum
c

c
c  --open parameterized receptor file (from subroutine parmrec)
c
      open (unit=unitno, file='PDBPARM', status='old')
c
c initialize variables
      natm_total=0      ! total # of R atoms to be saved
      natm_all=0        ! total # of R atoms in consideration
c         print*, 'maxatm=',maxatm,'maxneighb=',maxneighb
      do n=1,maxatm
        neighbor_num(n)=0
        do i=1,maxneighb
        neighbors(n,i)=0
        end do
      end do

c factor 2 is to account for all atoms whose SAS make contact
c with SAS of atoms in the grid box
      ilo=1-grdcuto*2
      do i=1,3
      ihi(i)=grdpts(i)+grdcuto*2
      end do


c add in rvdw (vdw radius)
  100 read (unitno, 1006, end=500) natm, vdwn, rsra, rsrb, rvdw,
     &rcrg, (rcrd(i), i=1,3)
 1006 format (2I5, 2(1x, F8.2), 1x, F5.3, 1x, F8.3, 1x, 3F8.3)
      if (vdwn .le. 0) go to 100    ! ignore atoms with 0 size

c
c  --subtract offset from receptor atom coordinates, find the 3D indices
c    of the nearest grid point (adding 1 because the lowest indices
c    are (1,1,1) rather than (0,0,0)); ignore receptor atoms farther
c    from the grid than the cutoff distance
c
      do 110 i=1,3
        rcrd(i)=rcrd(i) - offset(i)
        nearpt(i)=nint(rcrd(i)/grddiv) + 1
        if (nearpt(i) .gt. ihi(i)) go to 100 
        if (nearpt(i) .lt. ilo) go to 100 
  110 continue

      natm_all=natm_all+1    ! count # of R atoms in consideration
      r2=(rvdw+r_probe)**2   
      reff2G(natm_all)=r2    ! square of effective radius
      vdwnG(natm_all)=vdwn   ! vdw type
      do i=1,3
         crd_recG(natm_all,i)=rcrd(i)  ! grid coordinates
      end do


      if (nearpt(1) .le. grdpts(1) .and. nearpt(1) .gt. 0) then
        if (nearpt(2) .le. grdpts(2) .and. nearpt(2) .gt. 0) then
          if (nearpt(3) .le. grdpts(3) .and. nearpt(3) .gt. 0) then
            natm_total=natm_total+1    ! count # of R atoms in the box
            mark_in(natm_all)=1        ! mark for R atoms in the box
            rvdweff(natm_total)=r2     ! square of effective radius
            vdwn_rec(natm_total)=vdwn  ! vdw type
            do i=1,3
              crd_rec(natm_total,i)=rcrd(i)  ! grid coordinates
            end do
          end if
        end if
      end if

  200     continue
  300   continue
  400 continue
      go to 100
  500 continue

      write(2,*) 'total # of R atoms in the box = ', natm_total
      write(2,*) 'total # of R atoms in consideration = ', natm_all

c
c -- Reduce computational time by singling out neighbors of each atom;
c    SAS of each atom will be affected only by these neighbors.
c

      do i=1,natm_all-1
        do k=1,3
          rcrd(k) = crd_recG(i,k)
        end do
        do j=i+1,natm_all
          do k=1,3
            rcrd1(k) = crd_recG(j,k)
          end do
          rij2 = dist2(rcrd,rcrd1)
          if (rij2 .lt. r2_cutoff) then   ! screen for neighbors
            neighbor_num(i)=neighbor_num(i)+1
            neighbor_num(j)=neighbor_num(j)+1
            neighbors(i,neighbor_num(i))=j
            neighbors(j,neighbor_num(j))=i
          end if
        end do
      end do


c count # of grid pts on SAS

      nsas=0
      nsas_hp=0
      nsas_pol=0
      natm_total=0

      do natm=1,natm_all
        if (mark_in(natm).eq.1) then           ! atom inside box  1111
          nsum=0
          natm_total=natm_total+1  
          vdwn=vdwnG(natm)
          do k=1,3
            rcrd0(k) = crd_recG(natm,k)
          end do

          do 1 n=1,nsphgrid(vdwn)                ! 2222
            mark_sas(natm,n)=0                   ! initialization
            do k=1,3
              rcrd(k) = rcrd0(k) + sphgrid_crd(vdwn,n,k)
            end do

            do natm1=1,neighbor_num(natm)      ! 3333
              i=neighbors(natm,natm1)          ! can be in/out of box
              do k=1,3
                rcrd1(k) = crd_recG(i,k)
              end do

              rij2 = dist2(rcrd,rcrd1)
              if (rij2 .lt. reff2G(i)) goto 1  ! grid point not accessible
            end do                          ! 3333

            mark_sas(natm_total,n)=1        ! accessible to solvent
            nsas=nsas+1                     ! count # of grid pts on SAS
            nsum=nsum+1                     ! count # of grid pts on SAS
1         continue                          ! 2222

          if (nhp(vdwn) .eq. 1) then        ! atom is nonpolar
             nsas_hp=nsas_hp+nsum
          else                              ! atom is polar
             nsas_pol=nsas_pol+nsum
          end if

        end if                              ! 1111
      end do

      close (unitno)
      return
      end
c----------------------------------------------------------------------
