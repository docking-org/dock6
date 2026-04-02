c----------------------------------------------------------------------
c                       PROGRAM SHOWBOX            Elaine Meng   9/91
c
c     This program makes a pdb file for a box after interactively
c     requesting box dimensions and center.  Alternatively, the
c     center of mass of a user-specified sphere cluster can be used
c     as the box center.  A final option is to enclose all of the
c     spheres in a cluster, with a user-specified extra "cushion" of
c     space.
c----------------------------------------------------------------------
c
      real crd(3), crdsum(3), com(3), mincrd(3), maxcrd(3)
      real margin
      real boxdim(3)
      character*80 sphfil, boxfil
      character*80 line
      character*1 ctrtyp, yn
      integer i, j
      integer cntemp, nstemp, clnum, nsph
      logical auto
c
      write (6,*)
     &'automatically construct box to enclose spheres [Y/N]?'
      read (5, '(a1)') yn
      if (yn.eq.'Y'.or.yn.eq.'y') then
        auto=.true.
        write (6,*) 'extra margin to also be enclosed (angstroms)?'
        write (6,*) '  (this will be added in all 6 directions)'
        read (5,*) margin
      else
        auto=.false.
        write (6,*) 'box center user-defined [U] or sphere cluster '
        write (6, *) 'center of mass [S]?'
        read (5, '(A1)') ctrtyp
        if (ctrtyp .eq. 'U' .or. ctrtyp .eq. 'u') then
          write (6, *) 'enter box center coordinates [x y z] '
          read (5, *) (com(i), i=1,3)
          go to 600
        endif 
      endif
      write (6, *) 'sphere file-'
      read (5, 100) sphfil
  100 format (A80)
      write (6, *) 'cluster number-'
      read (5, *) clnum 
c
c  --initialize coordinate sums for calculating center of mass
c    and min/max coords for automatic option
c
      do 110 i=1,3
        crdsum(i)=0
  110 continue
      do 120 i=1,3
        mincrd(i)=10000.0
        maxcrd(i)=-10000.0
  120 continue
c
c  --open the sphere file; read cluster and calculate center of mass
c    or min/max coords for automatic option
c
        open (unit=1, file=sphfil, status='old')
c skip header records of new-format sphere cluster file
130     read(1,'(a80)') line
        if (line(1:7).ne.'cluster') goto 130
        read (line, 161, end=500) cntemp, nstemp
  161   format (8x, I5, 32x, I5)
        if (cntemp .eq. clnum) then
          nsph=nstemp
          do 300 i=1,nsph
            read (1,172) (crd(j), j=1,3)
  172       format (5x, 3F10.5)
            if (auto) then
              do 212 j=1,3
                if (crd(j).lt.mincrd(j)) mincrd(j)=crd(j)
                if (crd(j).gt.maxcrd(j)) maxcrd(j)=crd(j)
  212         continue
            else
              do 220 j=1,3
                crdsum(j)=crdsum(j) + crd(j)
  220         continue
            endif
  300     continue
          close (1)
          if (auto) then
            do 320 i=1,3
              mincrd(i)=mincrd(i) - margin
              maxcrd(i)=maxcrd(i) + margin
              com(i)=(maxcrd(i) + mincrd(i))/2.0
              boxdim(i)=maxcrd(i) - mincrd(i)
  320       continue
            go to 700
          else
            do 340 i=1,3
              com(i)=crdsum(i)/real(nsph)
  340       continue
            go to 600
          endif
        else
          goto 130
        endif
  500   continue
        close (1)
        write (6,*) 'sphere cluster not found'
        stop
  600 continue
      write (6, *) 'enter box dimensions [x y z] '
      read (5, *) (boxdim(i), i=1,3)
  700 continue
      write (6, *) 'output filename?'
      read (5, 100) boxfil
      open (unit=2, file=boxfil, status='unknown')
      write (2, '(A24)') 'HEADER    CORNERS OF BOX'
      write (2, '(A25, 3F8.3)') 'REMARK    CENTER (X Y Z) ',
     &(com(i), i=1,3)
      write (2, '(A29, 3F8.3)') 'REMARK    DIMENSIONS (X Y Z) ',
     &(boxdim(i), i=1,3)
      write (2, 1003) 'ATOM', 1, 'DUA', 'BOX', 1,
     &(com(1) - boxdim(1)/2.0), (com(2) - boxdim(2)/2.0),
     &(com(3) - boxdim(3)/2.0)
      write (2, 1003) 'ATOM', 2, 'DUB', 'BOX', 1,
     &(com(1) + boxdim(1)/2.0), (com(2) - boxdim(2)/2.0),
     &(com(3) - boxdim(3)/2.0)
      write (2, 1003) 'ATOM', 3, 'DUC', 'BOX', 1,
     &(com(1) + boxdim(1)/2.0), (com(2) - boxdim(2)/2.0),
     &(com(3) + boxdim(3)/2.0)
      write (2, 1003) 'ATOM', 4, 'DUD', 'BOX', 1,
     &(com(1) - boxdim(1)/2.0), (com(2) - boxdim(2)/2.0),
     &(com(3) + boxdim(3)/2.0)
      write (2, 1003) 'ATOM', 5, 'DUE', 'BOX', 1,
     &(com(1) - boxdim(1)/2.0), (com(2) + boxdim(2)/2.0),
     &(com(3) - boxdim(3)/2.0)
      write (2, 1003) 'ATOM', 6, 'DUF', 'BOX', 1,
     &(com(1) + boxdim(1)/2.0), (com(2) + boxdim(2)/2.0),
     &(com(3) - boxdim(3)/2.0)
      write (2, 1003) 'ATOM', 7, 'DUG', 'BOX', 1,
     &(com(1) + boxdim(1)/2.0), (com(2) + boxdim(2)/2.0),
     &(com(3) + boxdim(3)/2.0)
      write (2, 1003) 'ATOM', 8, 'DUH', 'BOX', 1,
     &(com(1) - boxdim(1)/2.0), (com(2) + boxdim(2)/2.0),
     &(com(3) + boxdim(3)/2.0)
 1003 format (A4, 6x, I1, 2x, A3, 1x, A3, 5x, I1, 4x, 3F8.3)
      write (2, 1004) 'CONECT', 1, 2, 4, 5
      write (2, 1004) 'CONECT', 2, 1, 3, 6
      write (2, 1004) 'CONECT', 3, 2, 4, 7
      write (2, 1004) 'CONECT', 4, 1, 3, 8
      write (2, 1004) 'CONECT', 5, 1, 6, 8
      write (2, 1004) 'CONECT', 6, 2, 5, 7
      write (2, 1004) 'CONECT', 7, 3, 6, 8
      write (2, 1004) 'CONECT', 8, 4, 5, 7
 1004 format (A6, 4I5)
      close (2)
      write (6, *) 'finished making pdb-format box file'
      end
c----------------------------------------------------------------------
