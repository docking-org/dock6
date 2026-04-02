      integer function orient_gk
     & (n,x0,x1,cg1_corr,rot,trans,refl)
c       rotates coords in x1 onto x0
c     notes/cmo 
c       golub-kabsch algorithm, with checks by cmo
c       returns early if discarding mirror'd ligand
c
c     implicit none
c
c
c     variables --
      integer n
c        n:  number of points to be matched.
      integer i,j,k
c        i:  do loop index.
c        j:  do loop index.
c        k:  do loop index.
	integer job,info  ! for dsvdc
      real an
c        an:  real number equal to the inverse of the number
c             of points to be matched.
c     arrays --
      real x0(3,*)
c        x0:  stationary coordinates.
      real x1(3,*)
c        x1:  coordinates to be least squares fit to x0.
      real*8 aa(3,3) 
c        aa:  correlation matrix used in calculation of rotation
c             matrix.
      real rot(3,3)
c        rot:  the rotation matrix.
      real cg1(3),cg1_corr(3),cg0(3)
c        trans: translation vector.
      real trans(3)
c        cg1:  center of mass of x1.
c        cg1_corr:  correction to the center of mass of x1.
c        cg0:  center of mass of x0.
      integer refl
c        refl: flag for whether a reflection is ok
c for new calc
	real*8 m,p,q,norm,temp,temp1,temp2,temp0,det
     x    ,phi,cphi,sphi
	real*8 lam(3)   ! eigval
	real*8 c(3,3),v(3,3),w(3,3),test(3,3)
     x    ,a1,a2,a3
c
	logical fh_orient
c
c
c       center of gravity calcn
      an=1./float(n)
      do 10 i=1,3
         cg0(i)=0.
         cg1(i)=0.
   10 continue
      do 12 i=1,n
         do 11 j=1,3
            cg0(j)=cg0(j)+x0(j,i)
            cg1(j)=cg1(j)+x1(j,i)
   11    continue
   12 continue
      do 14 i=1,3
         cg0(i)=cg0(i)*an
         cg1(i)=cg1(i)*an
   14 continue
c
      do 21 i = 1,3
         do 20 j = 1,3
            aa(i,j) = 0.0
   20    continue
   21 continue
c       calcn of correlation matrix
      do 32 k=1,n
         aa(1,1)=aa(1,1)+(x1(1,k)-cg1(1))*(x0(1,k)-cg0(1))
         aa(2,1)=aa(2,1)+(x1(2,k)-cg1(2))*(x0(1,k)-cg0(1))
	 aa(3,1)=aa(3,1)+(x1(3,k)-cg1(3))*(x0(1,k)-cg0(1))
         aa(1,2)=aa(1,2)+(x1(1,k)-cg1(1))*(x0(2,k)-cg0(2))
         aa(2,2)=aa(2,2)+(x1(2,k)-cg1(2))*(x0(2,k)-cg0(2))
         aa(3,2)=aa(3,2)+(x1(3,k)-cg1(3))*(x0(2,k)-cg0(2))
         aa(1,3)=aa(1,3)+(x1(1,k)-cg1(1))*(x0(3,k)-cg0(3))
         aa(2,3)=aa(2,3)+(x1(2,k)-cg1(2))*(x0(3,k)-cg0(3))
         aa(3,3)=aa(3,3)+(x1(3,k)-cg1(3))*(x0(3,k)-cg0(3))
   32 continue
c
	det = aa(1,1)*aa(2,2)*aa(3,3) + aa(2,1)*aa(3,2)*aa(1,3)
     x   + aa(3,1)*aa(2,3)*aa(1,2) - aa(1,3)*aa(2,2)*aa(3,1)
     x   - aa(2,3)*aa(3,2)*aa(1,1) - aa(3,3)*aa(1,2)*aa(2,1)
c	print *, "det (aa) = ", det
c
c
	do i=1,3
	do j=1,3
		c(i,j)=0.d0
	enddo
	enddo
c
	do k=1,3
		c(1,1) = c(1,1) + aa(k,1)*aa(k,1)
		c(1,2) = c(1,2) + aa(k,1)*aa(k,2)
		c(2,2) = c(2,2) + aa(k,2)*aa(k,2)
		c(1,3) = c(1,3) + aa(k,1)*aa(k,3)
		c(2,3) = c(2,3) + aa(k,2)*aa(k,3)
		c(3,3) = c(3,3) + aa(k,3)*aa(k,3)
	enddo
c eigenvalues
	m = ( c(1,1) + c(2,2) + c(3,3) )/3.d0
	q = ( (c(1,1)-m)*(c(2,2)-m)*(c(3,3)-m) ) +  c(1,2)*c(2,3)*c(1,3) 
     x	 + c(1,3)*c(2,3)*c(1,2) - ( (c(2,2)-m)*c(1,3)*c(1,3) ) 
     x   - ( (c(1,1)-m)*c(2,3)*c(2,3) ) - ( (c(3,3)-m)*c(1,2)*c(1,2) ) 
	q = q/2.d0
	p = ( (c(1,1)-m)*(c(1,1)-m) ) + ( (c(2,2)-m)*(c(2,2)-m) ) +
     x 		( (c(3,3)-m)*(c(3,3)-m) )
	p = p+2.d0*c(1,2)*c(1,2)+2.d0*c(1,3)*c(1,3)+2.d0*c(2,3)*c(2,3)
	p = p/6.d0

	fh_orient = .FALSE.
	if ((p*p*p - q*q).lt.0.d0) then
	    fh_orient = .TRUE.
	    goto 50
	endif

	temp = (dsqrt(p*p*p - q*q) )/q
	phi = ( atan(temp) )/3.d0
	if (temp.lt.0.d0) phi = 3.14159265358979d0 + phi

c  Check to make sure p not close to zero (100 * machine eps)
	if (phi.lt.1.75d-5) fh_orient = .TRUE.

50	if (fh_orient) then
		call orient_fh(aa,rot)
		goto 100
	endif
c
	p = dsqrt(p)
	cphi = dcos(phi)
	sphi = dsin(phi)

	lam(1) = m + 2.d0*p*cphi
	lam(2) = m - p*(cphi - dsqrt(3.d0)*sphi)
	temp0 = dsqrt(p)
	lam(3) = m - p*(cphi + dsqrt(3.d0)*sphi)
c
	if (temp.lt.0.d0) then
		temp0 = lam(1)
		lam(1) = lam(2)
		lam(2) = lam(3)
		lam(3) = temp0
	endif	
c
c eigvectors
c23456
	do i=1,3
	   temp1 = ( -c(1,3)*c(1,2) + (c(1,1)-lam(i))*c(2,3) ) /
     x		( c(1,2)*c(1,2) - (c(1,1)-lam(i))*(c(2,2)-lam(i)) )
	   temp2 = (-c(1,3)-c(1,2)*temp1)/(c(1,1)-lam(i))
	   norm = 1.d0/dsqrt( temp1*temp1 + temp2*temp2 + 1.d0)
	   v(1,i) = temp2*norm
	   v(2,i) = temp1*norm
	   v(3,i) = norm
	enddo 
c aa*v = w = u s
	do i=1,3 
	   do j=1,3
		temp = 0.
		do k=1,3
		   temp = temp + aa(i,k)*v(k,j)
		enddo
		w(i,j) = temp
	   enddo
	enddo
c
	do j=1,3
	   norm = w(1,j)*w(1,j)+w(2,j)*w(2,j)+w(3,j)*w(3,j)
	   norm = 1.d0/dsqrt(norm)
  	   w(1,j) = w(1,j)*norm
  	   w(2,j) = w(2,j)*norm
  	   w(3,j) = w(3,j)*norm
	enddo
        if ((det.lt.0.1d-2).and.(det.gt.-.1d-2)) then
           a1 = 1.d0 - (w(1,1)+w(2,1)+w(3,1))*w(1,1)
     x           - (w(1,2)+w(2,2)+w(3,2))*w(1,2)
           a2 = 1.d0 - (w(1,1)+w(2,1)+w(3,1))*w(2,1)
     x           - (w(1,2)+w(2,2)+w(3,2))*w(2,2)
           a3 = 1.d0 - (w(1,1)+w(2,1)+w(3,1))*w(3,1)
     x           - (w(1,2)+w(2,2)+w(3,2))*w(3,2)
           norm = a1*a1 + a2*a2 + a3*a3
           norm = dsqrt(norm)
           w(1,3) = a1/norm
           w(2,3) = a2/norm
           w(3,3) = a3/norm
	endif
c
c calculate rotation matrix
c
80	do i=1,3
	   do j=1,3
		temp = 0.
		do k=1,3
		   temp = temp + w(i,k)*v(j,k)
		enddo
		rot(j,i)=temp   ! to make consistent w orient
	   enddo
	enddo
c
c  Check to make sure transformation will not invert coordinates
c
	det = rot(1,1)*rot(2,2)*rot(3,3) + rot(2,1)*rot(3,2)*rot(1,3)
     x  	+ rot(3,1)*rot(2,3)*rot(1,2) - rot(1,3)*rot(2,2)*rot(3,1)
     x  	- rot(2,3)*rot(3,2)*rot(1,1) - rot(3,3)*rot(1,2)*rot(2,1)

	if (((det.lt.0.).and.((refl.eq.0).or.(refl.eq.-1))).or.
     x		((det.ge.0.).and.(refl.eq.1))) then
                w(1,3) = -w(1,3)
                w(2,3) = -w(2,3)
                w(3,3) = -w(3,3)
c
c recalculate rotation matrix
c
		do i=1,3
		   do j=1,3
			temp = 0.
			do k=1,3
			   temp = temp + w(i,k)*v(j,k)
			enddo
			rot(j,i)=temp   ! to make consistent w orient
		   enddo
		enddo
	endif
c
100	continue

	call get_translation (trans, rot, cg0, cg1, cg1_corr)
c
	orient_gk = 1
	return
	end
c
c//////////////////////////////////////////////////////////////////
      subroutine orient_fh(aa,rot)
c       rotates coords in x1 to give least-
c       squares fit to coords in x0
c       uses hermans and ferro program.
c
c
c     implicit none
c
      real tol
      parameter (tol = 0.001)
c        tol:  tolerance level check for whether iteration
c              can be stopped.
c
c     variables --
      integer i,j,k
c        i:  do loop index.
c        j:  do loop index.
c        k:  do loop index.
      integer ict
c        ict:  counter for number of iterations.  must do at
c              least 3 iterations.
      integer ix,iy,iz
c        ix:  pointer used in iterative least squares.  may
c             have the value 1, 2, or 3.  the value changes on
c             each iteration.
c        iy:  pointer used in iterative least squares.  may
c             have the value 1, 2, or 3.  the value changes on
c             each iteration.
c        iz:  pointer used in iterative least squares.  may
c             have the value 1, 2, or 3.  the value changes on
c             each iteration.
      integer iflag
c        iflag:  indicator variable for whether iterations can be
c                stopped.  must do at least 3 iterations and iflag
c                must equal 0 to be finished.
c             of points to be matched.
      real sig,gam,sg,bb,cc
c        sig:  temporary variable used in calculation of rotation
c              matrix.
c        gam:  temporary variable used in calculation of rotation
c              matrix.
c        sg:  temporary variable used in calculation of rotation
c              matrix.
c        bb:  temporary variable used in calculation of rotation
c              matrix.
c        cc:  temporary variable used in calculation of rotation
c              matrix.
c     arrays --
      real*8 aa(3,3)
c        aa:  correlation matrix used in calculation of rotation
c             matrix.
      real rot(3,3)
c        rot:  the rotation matrix.
c
      do 41 i=1,3
         do 40 j=1,3
            rot(i,j)=0.0
   40    continue
         rot(i,i)=1.0
   41 continue
c       from here to 70, iterative rotation scheme
      ict=0
      goto 51
   50 continue
      ix=ix+1
      if(ix.lt.4)goto 52
      if(iflag.eq.0)goto 70
   51 iflag=0
      ix=1
   52 continue
      ict=ict+1
      if(ict.gt.100)goto 70   ! changed
      iy=ix+1
      if(iy.eq.4) iy=1
      iz=6-ix-iy
      sig=aa(iz,iy)-aa(iy,iz)
      gam=aa(iy,iy)+aa(iz,iz)
      sg=sqrt(sig*sig+gam*gam)
      if(sg.eq.0.)goto 50
      sg=1./sg
      if(abs(sig).le.tol*abs(gam))goto 50
      do 60 k=1,3
         bb=gam*aa(iy,k)+sig*aa(iz,k)
         cc=gam*aa(iz,k)-sig*aa(iy,k)
         aa(iy,k)=bb*sg
         aa(iz,k)=cc*sg
         bb=gam*rot(iy,k)+sig*rot(iz,k)
         cc=gam*rot(iz,k)-sig*rot(iy,k)
         rot(iy,k)=bb*sg
         rot(iz,k)=cc*sg
   60 continue
      iflag=1
      goto 50
   70 continue
      return
      end
c
c-----------------------------------------------------------------------
      subroutine get_translation (trans, rot, cg0, cg1, cg1_corr)
c
      integer i
      real trans(3), rot(3,3), cg0(3), cg1(3), cg1_corr(3)
c     Construct translation vector
      do 110 i=1,3
        trans(i) = rot(i,1) * (cg1_corr(1) - cg1(1))
     &           + rot(i,2) * (cg1_corr(2) - cg1(2))
     &           + rot(i,3) * (cg1_corr(3) - cg1(3))
     &           + cg0(i) - cg1_corr(i)
110	continue
c
	return
	end
c-----------------------------------------------------------------------
      subroutine transform (natl,cl,coml,rot,trans,xatm)
c
c     this subroutine rotates and translates coordinates
c     according to a rotation matrix and translation vector
c
c     implicit none
c
c     variables --
      integer i, j
c        i, j:  do loop index.
      integer natl
c        natl:  number of atoms in ligand.
c
c     arrays --
      real cl(3,*)
      real  coml(3)
c        cl:  ligand atomic coordinates.
c        coml:  ligand center of mass.
      real xatm(3,*)
c        xatm:  rotated and translated ligand atomic coordinates.
      real rot(3,3)
c        rot:  rotation matrix determined by orient for this match.
      real  trans(3)
c        trans:  translation vector
c
c     apply rotation matrix calculated by orient
      do 30 i=1,natl
         xatm(1,i) = rot(1,1) * (cl(1,i) - coml(1))
     &             + rot(1,2) * (cl(2,i) - coml(2))
     &             + rot(1,3) * (cl(3,i) - coml(3))
     &             + coml(1) + trans(1)
         xatm(2,i) = rot(2,1) * (cl(1,i) - coml(1))
     &             + rot(2,2) * (cl(2,i) - coml(2))
     &             + rot(2,3) * (cl(3,i) - coml(3))
     &             + coml(2) + trans(2)
         xatm(3,i) = rot(3,1) * (cl(1,i) - coml(1))
     &             + rot(3,2) * (cl(2,i) - coml(2))
     &             + rot(3,3) * (cl(3,i) - coml(3))
     &             + coml(3) + trans(3)
   30 continue
      return
      end

c-----------------------------------------------------------------------
      subroutine transform_atom (atom,rot,orig)
c
c     this subroutine translates the ligand the center of mass of the
c     receptor site and rotates it according to the rotation matrcounteror ix
c     calculated by orient.
c
      implicit none
c
c     arrays --
      real atom(3)
      real atom_tmp(3)
c        atom:  atom coordinates.
c        atom_tmp:  rotated and translated atom coordinates.
      real rot(3,3)
c        rot:  rotation matrix
      real  orig(3)
c        orig:  point of rotation.
c
c     apply rotation matrix calculated by orient
         atom_tmp(1)=orig(1)+rot(1,1)*(atom(1)-orig(1))
     &                    +rot(1,2)*(atom(2)-orig(2))
     &                    +rot(1,3)*(atom(3)-orig(3))
         atom_tmp(2)=orig(2)+rot(2,1)*(atom(1)-orig(1))
     &                    +rot(2,2)*(atom(2)-orig(2))
     &                    +rot(2,3)*(atom(3)-orig(3))
         atom_tmp(3)=orig(3)+rot(3,1)*(atom(1)-orig(1))
     &                    +rot(3,2)*(atom(2)-orig(2))
     &                    +rot(3,3)*(atom(3)-orig(3))
	atom(1) = atom_tmp(1)
	atom(2) = atom_tmp(2)
	atom(3) = atom_tmp(3)
      return
      end
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C//////////////////////////////////////////////////////////////////////////////

        SUBROUTINE get_matrix_from_angles (rot_matrix, angle)

C       Computes rotation matrix given Euler angles.

        IMPLICIT NONE

        REAL            angle(3), rot_matrix(3,3)

        rot_matrix(1,1)= COS(angle(3))*COS(angle(1))-
     +                   COS(angle(2))*SIN(angle(1))*SIN(angle(3))
        rot_matrix(1,2)=-COS(angle(3))*SIN(angle(1))-
     +                   COS(angle(2))*COS(angle(1))*SIN(angle(3))
        rot_matrix(1,3)= SIN(angle(3))*SIN(angle(2))

        rot_matrix(2,1)= SIN(angle(3))*COS(angle(1))+
     +                   COS(angle(2))*SIN(angle(1))*COS(angle(3))
        rot_matrix(2,2)=-SIN(angle(3))*SIN(angle(1))+
     +                   COS(angle(2))*COS(angle(1))*COS(angle(3))
        rot_matrix(2,3)=-COS(angle(3))*SIN(angle(2))

        rot_matrix(3,1)= SIN(angle(2))*SIN(angle(1))
        rot_matrix(3,2)= SIN(angle(2))*COS(angle(1))
        rot_matrix(3,3)= COS(angle(2))

        RETURN
        END
C//////////////////////////////////////////////////////////////////////////////
