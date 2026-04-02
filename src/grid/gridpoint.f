
c---------------------------------------------------------------------
      subroutine get_index(crd0,offset,grddiv,grdpts,out,index,ccrd)
c
c  --does trilinear interpolations, using the values at the 8
c  vertices of the tiny grid cube enclosing each point of interest;
c  returns the indices and cubic coordinates for use when retreiving
c  the grid value using get_value.
c  Interpolation algorithm taken from Honig et al.'s DelPhi program.
c  ECM   11/91
c  Modified by Todd Ewing to only compute indices
C	
C      W = A1.XYZ + A2.XY + A3.XZ + A4.YZ + A5.X +
C		A6.Y + A7.Z + A8
C	
C      where Ai coefficients are linear combinations of Wi at
C      the cube corners
c
c---------------------------------------------------------------------
c
c
      real offset(3), grddiv
      integer grdpts(3)
      integer out
      real crd(3), crd0(3), ccrd(3)
      integer low(3)
      integer i, index(8)
      integer indx1
c
      out = 0
c
c  --find lower left bottom grid point, return 0.0 if coordinates are
c    outside the grid; find nearest grid point (1 added in both cases
c    since lowest indices are (1,1,1) rather than (0,0,0))
c
      do 20 i=1,3
        crd(i)=crd0(i) - offset(i)
        low(i)=int(crd(i)/grddiv) + 1
        if (low(i) .lt. 1  .or. low(i) .ge. grdpts(i) .or.
     &  crd(i) .lt. 0.0) then
          out=1
          go to 60
        endif
   20 continue
c
c  --calculate cube coordinates of point (between 0 and 1)
c
      do 40 i=1,3
        ccrd(i)=crd(i)/grddiv - float(int(crd(i)/grddiv))
   40 continue
c
c  --get 1-dimensional indices for the points of interest
c
      index(8) = indx1(low(1),low(2),low(3),grdpts)
      index(7) = indx1(low(1),low(2),low(3)+1,grdpts)
      index(6) = indx1(low(1),low(2)+1,low(3),grdpts)
      index(5) = indx1(low(1)+1,low(2),low(3),grdpts)
      index(4) = indx1(low(1),low(2)+1,low(3)+1,grdpts)
      index(3) = indx1(low(1)+1,low(2),low(3)+1,grdpts)
      index(2) = indx1(low(1)+1,low(2)+1,low(3),grdpts)
      index(1) = indx1(low(1)+1,low(2)+1,low(3)+1,grdpts)
c
   60 continue
      return
      end


c---------------------------------------------------------------------
      real function get_value(grid,ccrd,index)
c
c  Uses indices and cubic coordinates to retrieve an interpolated
c  value from a energy grid.
c  written by Todd Ewing based on code by Elaine Meng
c
c---------------------------------------------------------------------
c
c
      real grid(*)
      real ccrd(3)
      integer index(8)
      real a1, a2, a3, a4, a5, a6, a7, a8
c
c  --calculate coefficients of trilinear function for aval
c
      a8 = grid(index(8))
      a7 = grid(index(7)) - a8
      a6 = grid(index(6)) - a8
      a5 = grid(index(5)) - a8
      a4 = grid(index(4)) - a8 - a7 - a6
      a3 = grid(index(3)) - a8 - a7 - a5
      a2 = grid(index(2)) - a8 - a6 - a5
      a1 = grid(index(1)) - a8 - a7 - a6 - a5 - a4 - a3 - a2
c
c  --determine interpolated grid value
c
      get_value = a1*ccrd(1)*ccrd(2)*ccrd(3) + a2*ccrd(1)*ccrd(2) +
     &a3*ccrd(1)*ccrd(3) + a4*ccrd(2)*ccrd(3) + a5*ccrd(1) +
     &a6*ccrd(2) + a7*ccrd(3) + a8

      return
      end


c
c-----------------------------------------------------------------------
      integer function indx1(i,j,k,grdpts)
c
c  --converts the 3-dimensional (virtual) indices of a grid point to the
c    actual index in a 1-dimensional array
c
      integer i, j, k
      integer grdpts(3)
c
      indx1 = grdpts(1)*grdpts(2)*(k-1) + grdpts(1)*(j-1) + i
      return
      end
c-----------------------------------------------------------------------

