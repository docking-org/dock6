c................NEW.................................
      subroutine sphgrid(radius, ngrid, sgcrd)
c
c  --called from CHEMGRID
c  --Assign grids on the spherical surface of the atom
c----------------------------------------------------------------------
c

      include 'nchemgrid.h'

      real radius
c effective atomic radius (r_vdw + r_probe)
      integer ngrid
c total number of grid points on the atomic spherical surface 
      real sgcrd(maxgrd,3)
c cardial coordinates in angstroms of grid points relative to atomic center
      real PI, PI2
      parameter (PI=3.141592654, PI2=PI*2)
      integer i, j, ntheta, nphi
      real theta, delta_theta, phi, delta_phi, r, z
c
      delta_theta=spacing/radius 

      sgcrd(1,1)=0.               ! corresp. to top pt (theta=phi=0)
      sgcrd(1,2)=0.
      sgcrd(1,3)=radius
      ntheta=PI/delta_theta       ! theta<=PI

      ngrid=1
      theta=delta_theta
      do i=1,ntheta
        z=radius*cos(theta)
        r=radius*sin(theta)
        delta_phi=spacing/r
        nphi=nint(PI2/delta_phi+0.499)  ! phi<2PI 

        phi=0.0
        do j=1,nphi
          ngrid=ngrid+1
          sgcrd(ngrid,1)=r*cos(phi)
          sgcrd(ngrid,2)=r*sin(phi)
          sgcrd(ngrid,3)=z
          phi=phi+delta_phi
        end do
        theta=theta+delta_theta
      end do

      if (ngrid.gt.maxgrd) pause 'Maximum spherical grid pts exceeded.
     & Increase maxgrd or reduce spacing.'

      return
      end
c----------------------------------------------------------------------
