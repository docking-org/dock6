c**************************** NEW *************************************
      subroutine polariz(unitno, grdcut, grddiv, grdpts, offset,
     &grdcuto, natm_total)
c
c  --called from CHEMGRID
c  --preparation for polarization energy calculation
c----------------------------------------------------------------------
c
      include 'nchemgrid.h'
c
      real mincon, minsq
      parameter (mincon=0.0001)
      integer unitno, i, j, k, n, natm_all
      integer ilo, ilo1(3), ihi(3), ihi1(3)
      real rcrd1(3), r2, rvdweff2, tem, sum, inv_aij2, dij, f_gb
      integer mark_Br(maxatm)
      integer max1, max2
c      parameter (max1=-80, max2=170)
      parameter (max1=-50, max2=204)
      integer mark_occupG(max1:max2,max1:max2,max1:max2)
c

      minsq=mincon*mincon
c
c  --open parameterized receptor file (from subroutine parmrec)
c
      open (unit=unitno, file='PDBPARM', status='old')
c
c initialize variables
      natm_all=0        ! total # of R atoms in consideration
      natm_total=0      ! total # of R atoms to be saved

      do i=1,maxatm
        mark_Br(i)=0
      end do

      ilo=1-grdcutb
      do i=1,3
      ihi(i)=grdpts(i)+grdcutb
      end do



      if (ilo.lt.max1) then
          write(2,*) 'ilo = ', ilo, 'max1 =', max1
          write (2, *) 'Lower max1, program stops'
          stop
      end if
      do i=1,3
        if (ihi(i).gt.max2) then
          write(2,*) 'ihi(',i,') = ', ihi(i), 'max2 =', max2
          write (2, *) 'Increase max2, program stops'
          stop
        end if
      end do

      do 1 i=ilo, ihi(1)
      do 1 j=ilo, ihi(2)
      do 1 k=ilo, ihi(3)
        mark_occupG(i,j,k)=0
    1 continue



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

      natm_all=natm_all+1        ! count # of R atoms in consideration
      rvdweff0=rvdw+vdwoff  ! effective vdw radius in grid unit
      rvdweff2=rvdweff0*rvdweff0

      if (nearpt(1) .le. grdpts(1) .and. nearpt(1) .gt. 0) then
        if (nearpt(2) .le. grdpts(2) .and. nearpt(2) .gt. 0) then
          if (nearpt(3) .le. grdpts(3) .and. nearpt(3) .gt. 0) then
c      if (rcrd(1) .le. boxdim(1) .and. rcrd(1) .ge. 0) then
c        if (rcrd(2) .le. boxdim(2) .and. rcrd(2) .ge. 0) then
c          if (rcrd(3) .le. boxdim(3) .and. rcrd(3) .ge. 0) then
            natm_total=natm_total+1    ! count # of R atoms in the box
            mark_Br(natm_all)=1        ! symbol to compute Born radius
            rcrg_rec(natm_total)=rcrg     ! charge 
            rvdweff(natm_total)=rvdweff0  ! effective vdw radius
            do i=1,3
              crd_rec(natm_total,i)=rcrd(i)  ! grid coordinates
            end do
          end if
        end if
      end if



c
c --loop through grid points within the cutoff cube (not sphere) of
c   the current receptor atom, but only increment values if the grid
c   point is within the cutoff sphere for the atom
c

      do i=1,3
      ilo1(i)=nearpt(i)-grdcuto
      ihi1(i)=nearpt(i)+grdcuto
      end do

      do 400 i=max(ilo,ilo1(1)), min(ihi(1),ihi1(1))
        gcrd(1)=float(i-1)*grddiv
        do 300 j=max(ilo,ilo1(2)), min(ihi(2),ihi1(2))
          gcrd(2)=float(j-1)*grddiv
          do 200 k=max(ilo,ilo1(3)), min(ihi(3),ihi1(3))
            gcrd(3)=float(k-1)*grddiv
            r2 = dist2(rcrd,gcrd)
c
c add in marks for gridpoints overlapped with the current receptor atom
c
            if (r2 .le. rvdweff2) mark_occupG(i,j,k)=1


  200     continue
  300   continue
  400 continue
      go to 100
  500 continue

c**********************************************************************
c add in calculation of sum of (1/r_ij^4) over j (i,j: receptor atoms)
c compute Born radius (alpha)

      write(2,*) 'total # of R atoms in the box = ', natm_total
      write(2,*) 'total # of R atoms for evaluation of Born radius = ', 
     &           natm_all

      tem = grddiv**3/4.0/3.141592654
c      tem = 1.0/(4.0*3.141592654)
c inv_a is in grid unit


      do 700 natm=1,natm_total
        rvdweff0=rvdweff(natm)
        rvdweff2=rvdweff0*rvdweff0
        do i=1,3
          rcrd(i)=crd_rec(natm,i)
          nearpt(i)=nint(rcrd(i)/grddiv) + 1
          ilo1(i)=nearpt(i)-grdcut
          ihi1(i)=nearpt(i)+grdcut
        end do
        sum=0.0
        do 600 i=max(ilo,ilo1(1)), min(ihi(1),ihi1(1))
          gcrd(1)=float(i-1)*grddiv
          do 600 j=max(ilo,ilo1(2)), min(ihi(2),ihi1(2))
            gcrd(2)=float(j-1)*grddiv
            do 600 k=max(ilo,ilo1(3)), min(ihi(3),ihi1(3))
              gcrd(3)=float(k-1)*grddiv
              if (mark_occupG(i,j,k).eq.1) then  ! grid pts overlap with atoms
                r2 = dist2(rcrd,gcrd)
                if (r2 .gt. rvdweff2) sum=sum+1.0/(r2*r2) 
c                       no overlap with current atom enveloped by a probe atom
              end if
  600   continue
        sum_rec(natm)=sum
        inv_a(natm) = 1.0/rvdweff0 - sum*tem ! inverse of Born radius
c          write(100,*) natm, rvdweff0, 1./inv_a(natm)
  700 continue



c compute solvation component -- polarization energy term of pure receptor

      solv_rec = 0.0

      do 800 i=1,natm_total
      do 800 j=1,natm_total

        do k=1,3
          rcrd(k) = crd_rec(i,k)
          rcrd1(k) = crd_rec(j,k)
        end do
        r2 = dist2(rcrd,rcrd1)
        if (r2 .gt. cutsq) goto 800
        inv_aij2 = inv_a(i) * inv_a(j)
c        dij = r2 * inv_aij2 / 4.0
        dij = r2 * inv_aij2 / f_scale
        f_gb= sqrt(r2 + exp(-dij)/inv_aij2)
        solv_rec = solv_rec + rcrg_rec(i)*rcrg_rec(j)/f_gb
  800 continue

c      tem = 166.0 * (1.0-1.0/78.3)
c      solv_rec = solv_rec * tem

c in real unit, should be 166.0 * (1.0-1.0/78.3) / grddiv * solv_rec
      write (2, '(A,G13.7)') 'solv_rec = ', solv_rec



c in real unit, should be 332.0 / grddiv * esval

      do 900 natm=1,natm_total
        do i=1,3
          rcrd(i)=crd_rec(natm,i)
          nearpt(i)=nint(rcrd(i)/grddiv) + 1
          ilo1(i)=nearpt(i)-grdcut
          ihi1(i)=nearpt(i)+grdcut
        end do
        rcrg=rcrg_rec(natm)
        do i=max(1,ilo1(1)), min(grdpts(1),ihi1(1))
          gcrd(1)=float(i-1)*grddiv
          do j=max(1,ilo1(2)), min(grdpts(2),ihi1(2))
            gcrd(2)=float(j-1)*grddiv
            do k=max(1,ilo1(3)), min(grdpts(3),ihi1(3))
              gcrd(3)=float(k-1)*grddiv
              n = indx1(i,j,k,grdpts)
              if (mark_occupG(i,j,k).eq.1) mark_occup(n)=1
              r2 = dist2(rcrd,gcrd)
              if (r2 .gt. cutsq) goto 2
              if (r2 .lt. minsq) then
                bump(n)='X'
                r2 = minsq
              else if(((r2 .lt. cconsq .and. vdwn .le. 5) .or. (r2 .lt.
     &        pconsq .and. vdwn .ge. 8)) .and. bump(n) .eq. 'F') then
                bump(n)='T'
              endif
              esval(n)=esval(n) + rcrg/sqrt(real(r2))
    2         continue
            end do
          end do
        end do
  900 continue


      close (unitno)
      return
      end
c----------------------------------------------------------------------
