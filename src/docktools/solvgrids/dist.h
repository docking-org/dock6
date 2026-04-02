c--------------------------------------------------------------------
c         header file for solmap.  BKS 11/2002
c--------------------------------------------------------------------
        integer maxlin,maxdix,maxdiy,maxdiz
        parameter (maxlin=10000)
        parameter (maxdix=150,maxdiy=150,maxdiz=150)
c  maxlin: maximum number of atoms in site
c  maxdix,maxdiy,maxdiz: maximum number of grid points per x/y/z dimension
        real dlim
c        dlim: minimum distance away from i'th grid to be counted in
c           xcluded volume calculation.  
        real    rad_o, rad_c, rad_s, rad_p, rad_n, rad_q, probe
c        rad_o,c,s,p,q: atom radii.  q is for everything undefined.
c        probe: the radius of the solvent probe 
        real    dV, probe2
c        dv  : the volume/4*pi of a grid element
c        probe2: probe^2
        real    solgrid(0:maxdix,0:maxdiy,0:maxdiz)
c        solgrid: grid of effective dielectric, based on displace volume.

        common /sol/ solgrid
        common /rads/ rad_o,rad_c,rad_s,rad_p,rad_n,rad_q,probe
