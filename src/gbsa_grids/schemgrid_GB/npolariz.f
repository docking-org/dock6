c zou, 7/96
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
      real rcrd1(3), r2, r6, rvdweff2, tem, sum, inv_aij2, dij, f_gb
      integer vdwn_rec(maxatm)
      real sumc_rec(maxatm) 
      real rsra_rec(maxatm), rsrb_rec(maxatm)
      integer max1, max2
c      parameter (max1=-80, max2=170)
      parameter (max1=-50, max2=204)
      integer mark_occupG(max1:max2,max1:max2,max1:max2)
c zou, 7/99
      real xi,yi,zi,ri,si,si2, screen
      real rj,sj,sj2,rij,RLij,RUij,RLij2,RUij2,tmp

      real atime, btime
cliuhy debug*********************************
      integer kkk
      integer natm_rec(maxatm)
cliuhy add***********************
      real crd_rec_aux(maxatm,3)
      real rvdweff_aux(maxatm),descreen_aux(maxatm)
      real r_dummy,descreen_dummy
      parameter(r_dummy=1.85,descreen_dummy=0.8)
      integer natm_cavity
      real INVA_NEG
  
      INVA_NEG=0.1
c**************************************
c
      minsq=mincon*mincon
c
c  --open parameterized receptor file (from subroutine parmrec)
c
      open (unit=unitno, file='PDBPARM', status='old')
      open (unit=15, file='PDBCAV', status='old')
c
c initialize variables
      natm_all=0        ! total # of R atoms in consideration
      natm_total=0      ! total # of R atoms to be saved

      do i=1,maxatm
        sumc_rec(i)=0
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

       call second(atime)
c      write (6, *) '***elapsed time before pair: ',atime
     

c read in receptor information
c add in rvdw (vdw radius)
c zou, 7/99
c add in screen (descreening parameter)
  100 read (unitno, 1006, end=500) natm, vdwn, rsra, rsrb, rvdw, screen,
     &rcrg, (rcrd(i), i=1,3)
 1006 format (2I5, 2(1x, F8.2), 1x, F5.3, 1x, F5.3, 1x, F8.3, 1x, 3F8.3)
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
cliuhy add****************
      descreen_aux(natm_all)=screen
      rvdweff_aux(natm_all)=rvdweff0
      do i=1,3
        crd_rec_aux(natm_all,i)=rcrd(i)
      end do
c************************

      if (nearpt(1) .le. grdpts(1) .and. nearpt(1) .gt. 0) then
        if (nearpt(2) .le. grdpts(2) .and. nearpt(2) .gt. 0) then
          if (nearpt(3) .le. grdpts(3) .and. nearpt(3) .gt. 0) then
c      if (rcrd(1) .le. boxdim(1) .and. rcrd(1) .ge. 0) then
c        if (rcrd(2) .le. boxdim(2) .and. rcrd(2) .ge. 0) then
c          if (rcrd(3) .le. boxdim(3) .and. rcrd(3) .ge. 0) then
            natm_total=natm_total+1    ! count # of R atoms in the box
            rcrg_rec(natm_total)=rcrg     ! charge 
            rvdweff(natm_total)=rvdweff0  ! effective vdw radius
c zou, 7/99
            descreen(natm_total)=screen   ! descreening parameter
            vdwn_rec(natm_total)=vdwn     ! vdw type
            rsra_rec(natm_total)=rsra     ! vdw parameter
            rsrb_rec(natm_total)=rsrb     ! vdw parameter
            do i=1,3
              crd_rec(natm_total,i)=rcrd(i)  ! grid coordinates
            end do
cliuhy debug********************
           natm_rec(natm_total)=natm
c************************
          end if
        end if
      end if



c
c --loop through grid points within the cutoff cube (not sphere) of
c   the current receptor atom, but only increment values if the grid
c   point is within the cutoff sphere for the atom
c
cliuhy add******************
c pairwise calculation do not need the grid

      if (mark_pair.ne.1) then
c*******************************

      do i=1,3
      ilo1(i)=max(nearpt(i)-grdcuto, ilo)
      ihi1(i)=min(nearpt(i)+grdcuto, ihi(i))
      end do

      do 400 i=ilo1(1), ihi1(1)
        gcrd(1)=float(i-1)*grddiv
        do 300 j=ilo1(2), ihi1(2)
          gcrd(2)=float(j-1)*grddiv
          do 200 k=ilo1(3), ihi1(3)
            gcrd(3)=float(k-1)*grddiv
            if (mark_occupG(i,j,k).eq.1) goto 200
            r2 = dist2(rcrd,gcrd)
c
c add in marks for gridpoints overlapped with the current receptor atom
c
            if (r2 .le. rvdweff2) mark_occupG(i,j,k)=1


  200     continue
  300   continue
  400 continue
cliuhy add*************
      end if
c*********************
      go to 100
  500 continue

c**********************************************************************
c add in calculation of sum of (1/r_ij^4) over j (i,j: receptor atoms)
c compute Born radius (alpha)

      write(2,*) 'total # of R atoms in the box = ', natm_total
      write(2,*) 'total # of R atoms for evaluation of Born radius = ', 
     &           natm_all

cliuhy debug*********************
c add to calculate the speed of pairwise
       call second(atime)
c*************************

cliuhy add*****************************
c open an error file to record the "supposed" embeded atoms
      open(24,file='NEG_INVA')
c zou, 7/99
      if (mark_pair.eq.1) then      ! pairwise calculations  1111111

cliuhy add**********************
        do i=1,natm_total
          inv_a(i)=1.0/(rvdweff(i)-vdwoff)
cc         descreen(i)=0.5
        end do

        do i=1,natm_total
          xi=crd_rec(i,1)
          yi=crd_rec(i,2)
          zi=crd_rec(i,3)
          ri=rvdweff(i)-vdwoff
          do j=1,natm_all
cliuhy add*************************
c need not the follow line for two reasons:
c 1. j and i are not in the same sequence (box flag), "i" is in grid box
c     "j" is in auxiliary box
c 2. if they are same atom, rij=0. and sj<ri (screen<1.0)
c
c           if (j.eq.i) goto 876
c*********************************
            rj=rvdweff_aux(j)-vdwoff
            sj=descreen_aux(j)*rj
            sj2=sj*sj
            rij=sqrt((crd_rec_aux(j,1)-xi)**2+(crd_rec_aux(j,2)-yi)**2
     &               +(crd_rec_aux(j,3)-zi)**2)
            if (ri .lt. rij+sj) then
              RLij=1./max(ri,rij-sj)
              RUij=1./(rij+sj)
              RLij2=RLij*RLij
              RUij2=RUij*RUij
              tmp=RLij-RUij+0.25*(sj2/rij-rij)*(RLij2-RUij2)
     &            -0.5/rij*alog(RLij/RUij)
            end if
            if (sj .le. ri-rij) then
              tmp=0.
            end if
cliuhy debug******************
c      if (i.eq.427) then
c         print*, i,j,ri,sj,rij,'******',rj
c         print*, i,natm_rec(i)
c         print*, natm_rec(427),natm_rec(436),natm_rec(589)
c         stop
c      end if

             if (tmp.lt.-1.E-6) print*, i,j,sj,rij,RLij,RUij,0.5*tmp
cliuhy debug***********************
c     do kkk=1,10
c       write(*,*) descreen(kkk)
c     end do
ccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc
c             if (tmp.lt.0.) print*, i,j,sj,rij,RLij,RUij,0.5*tmp
c             if (tmp.lt.0.) print*, ri,sj,rij,'******',rj
c             if (tmp.lt.0) then
c               print*, RLij-RUij,0.25*(sj2/rij-rij)*(RLij2-RUij2),
c    &            -0.5/rij*alog(RLij/RUij)
c               print*, 'RLij2,RUij2,RLij2-RUij2,sj2/rij,sj2/rij-rij'
c               print*, RLij2,RUij2,RLij2-RUij2,sj2/rij,sj2/rij-rij
c               print*, RLij/RUij,alog(RLij/RUij)
c             end if
cccccccccccccccccccccccccccccccccc
c************************************

              inv_a(i)=inv_a(i)-0.5*tmp
876       continue
          end do
cliuhy add*************************
c          pairwise
c if inv_a(i) is negative, that means it is embeded deeply
c and we just consider inv_a(i) = 0.
c meanwhile, we write out a file to record these atoms
        if (inv_a(i).lt.0) then
          write(24,2424) i,natm_rec(i),vdwn_rec(i),ri
          write(*,2424) i,natm_rec(i),vdwn_rec(i),ri
          write(*,*) xi+offset(1),yi+offset(2),zi+offset(3)
          write(*,*) rcrg_rec(i)
          inv_a(i)=INVA_NEG
        end if
2424      format(3i5,f8.3)
        end do
       call second(btime)
      write (6, *) '***elapsed time AFTER pair: ', btime-atime


      else                      ! original calculations  1111111

      tem = grddiv**3/4.0/3.141592654
c      tem = 1.0/(4.0*3.141592654)
c inv_a is in grid unit


      do 700 natm=1,natm_total
        rvdweff0=rvdweff(natm)
        rvdweff2=rvdweff0*rvdweff0
        do i=1,3
          rcrd(i)=crd_rec(natm,i)
          nearpt(i)=nint(rcrd(i)/grddiv) + 1
          ilo1(i)=max(nearpt(i)-grdcut, ilo)
          ihi1(i)=min(nearpt(i)+grdcut, ihi(i))
        end do
        sum=0.0
        do 600 i=ilo1(1), ihi1(1)
          gcrd(1)=float(i-1)*grddiv
          do 600 j=ilo1(2), ihi1(2)
            gcrd(2)=float(j-1)*grddiv
            do 600 k=ilo1(3), ihi1(3)
              gcrd(3)=float(k-1)*grddiv
              if (mark_occupG(i,j,k).eq.1) then  ! grid pts overlap with atoms
                r2 = dist2(rcrd,gcrd)
                if (r2 .gt. rvdweff2) sum=sum+1.0/(r2*r2) 
c                       no overlap with current atom enveloped by a probe atom
              end if
  600   continue
        inv_a(natm)=1.0/rvdweff0 - sum*tem     ! inverse of Born radius
cliuhy add*************************
c          grid-base
c if inv_a(i) is negative, that means it is embeded deeply
c and we just consider inv_a(i) = 0.
c meanwhile, we write out a file to record these atoms
        if (inv_a(natm).lt.0) then
          write(24,2424) natm,natm_rec(natm),
     &                   vdwn_rec(natm),rvdweff0-vdwoff
          inv_a(natm)=0.001
        end if
  700 continue

       call second(btime)
      write (6, *) '***elapsed time BEFORE pair: ', btime-atime
      end if                    ! 1111111




c in real unit, should be 332.0 / grddiv * esval

      do 900 natm=1,natm_total
        do i=1,3
          rcrd(i)=crd_rec(natm,i)
          nearpt(i)=nint(rcrd(i)/grddiv) + 1
          ilo1(i)=max(nearpt(i)-grdcut, 1)
          ihi1(i)=min(nearpt(i)+grdcut, grdpts(i))
        end do
c        rcrg=rcrg_rec(natm)
        vdwn=vdwn_rec(natm)
        rsra=rsra_rec(natm)
        rsrb=rsrb_rec(natm)
        do i=ilo1(1), ihi1(1)
          gcrd(1)=float(i-1)*grddiv
          do j=ilo1(2), ihi1(2)
            gcrd(2)=float(j-1)*grddiv
            do k=ilo1(3), ihi1(3)
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
              r6 = r2*r2*r2
              aval(n)=aval(n) + rsra/(r6*r6)
              bval(n)=bval(n) + rsrb/r6
    2         continue
            end do
          end do
        end do
  900 continue

       call second(btime)

c mark cavity
cliuhy add*****************************
      natm_cavity=0
c*************************************

   99 read (15, 3001, end=499) natm, (rcrd(i), i=1,3)
 3001 format (I5, 1x, 3F8.3)
      do i=1,3
        rcrd(i)=rcrd(i) - offset(i)
        nearpt(i)=nint(rcrd(i)/grddiv) + 1
      end do
cliuhy add********************************
c record the dummy cavity atoms
      natm_cavity=natm_cavity+1
      do i=1,3
        crd_rec_aux(natm_cavity,i)=rcrd(i)
      end do
c*******************************************

cliuhy add******************
c pairwise do not need the grid
      if (mark_pair.ne.1) then
c****************************
      do i=1,3
      ilo1(i)=nearpt(i)-grdcuto
      ihi1(i)=nearpt(i)+grdcuto
      end do

      do i=ilo1(1), ihi1(1)
        gcrd(1)=float(i-1)*grddiv
        do j=ilo1(2), ihi1(2)
          gcrd(2)=float(j-1)*grddiv
          do k=ilo1(3), ihi1(3)
            gcrd(3)=float(k-1)*grddiv
            r2 = dist2(rcrd,gcrd)
c
c add in marks for gridpoints overlapped with the dummy atom in the cavity
c
            if (r2 .le. rdum2) then
              mark_occupG(i,j,k)=2
            end if
      end do
      end do
      end do

cliuhy add********************
      end if
c******************************
      goto 99
  499 continue
      close (15)

c effect of cavity on Born radii of receptor atoms


cliuhy add*************************
      if (mark_pair.eq.1) then      ! pairwise calculations  1111111

        do i=1,natm_total
          xi=crd_rec(i,1)
          yi=crd_rec(i,2)
          zi=crd_rec(i,3)
          ri=rvdweff(i)-vdwoff
          sum=0.
          do j=1,natm_cavity
cliuhy add*************************
c need not the follow line for two reasons:
c 1. j and i are not in the same sequence (box flag), "i" is in grid box
c     "j" is in auxiliary box
c 2. if they are same atom, rij=0. and sj<ri (screen<1.0)
c
c           if (j.eq.i) goto 876
c*********************************
            rj=r_dummy
            sj=rj*descreen_dummy
            sj2=sj*sj
            rij=sqrt((crd_rec_aux(j,1)-xi)**2+(crd_rec_aux(j,2)-yi)**2
     &               +(crd_rec_aux(j,3)-zi)**2)
            if (ri .lt. rij+sj) then
              RLij=1./max(ri,rij-sj)
              RUij=1./(rij+sj)
              RLij2=RLij*RLij
              RUij2=RUij*RUij
              tmp=RLij-RUij+0.25*(sj2/rij-rij)*(RLij2-RUij2)
     &            -0.5/rij*alog(RLij/RUij)
            end if
            if (sj .le. ri-rij) then
              tmp=0.
            end if

            if (tmp.lt.-1.E-6) print*, i,j,sj,rij,RLij,RUij,0.5*tmp
            sum=sum+tmp
          end do
          sumc_rec(i)=sum
cliuhy debug**********************
c         print*,i,sum,cav_scale
c**********************************

          inv_a(i) = inv_a(i) - 0.5*sum*cav_scale ! inv of Born radius

        end do

      else           !original calculation

      do 701 natm=1,natm_total
        rvdweff0=rvdweff(natm)
        rvdweff2=rvdweff0*rvdweff0
        do i=1,3
          rcrd(i)=crd_rec(natm,i)
          nearpt(i)=nint(rcrd(i)/grddiv) + 1
          ilo1(i)=max(nearpt(i)-grdcut, ilo)
          ihi1(i)=min(nearpt(i)+grdcut, ihi(i))
        end do
        sum=0.0
        do 601 i=ilo1(1), ihi1(1)
          gcrd(1)=float(i-1)*grddiv
          do 601 j=ilo1(2), ihi1(2)
            gcrd(2)=float(j-1)*grddiv
            do 601 k=ilo1(3), ihi1(3)
              gcrd(3)=float(k-1)*grddiv
              n = indx1(i,j,k,grdpts)
              if (mark_occupG(i,j,k).eq.2 .and. mark_occup(n).eq.0) then  ! grid pts occupied by cavity but not by receptor atoms
                r2 = dist2(rcrd,gcrd)
                if (r2 .gt. rvdweff2) sum=sum+1.0/(r2*r2) 
c                       no overlap with current atom
              end if
  601   continue
        sumc_rec(natm)=sum
        inv_a(natm) = inv_a(natm) - sum*tem*cav_scale ! inverse of Born radius
  701 continue
      
      end if


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


      solv_rec = - 166.0 * (1.0-1.0/78.3) * solv_rec
      write (2, '(A,G13.7)') 'solv_rec = ', solv_rec

cliuhy debug*************************
c just written out the inv_a for compare

      open(22,file='inva')
      do natm=1,natm_total
        write(22,*) natm,inv_a(natm)
      end do
      close(22)

c*****************************

      do natm=1,natm_total
        inv_a(natm) = inv_a(natm) - sumc_rec(natm)*tem*(1.0-cav_scale) ! cavity correction
      end do



        do i=1, grdpts(1)
        do j=1, grdpts(2)
        do k=1, grdpts(3)
          n = indx1(i,j,k,grdpts)
          if (mark_occupG(i,j,k).eq.2) mark_occup(n)=1
        end do
        end do
        end do


c        do i=1,natm_total
c          write (100,*) i, inv_a(i)
c        end do


      close (unitno)
      close (15)
      return
      end
c----------------------------------------------------------------------
