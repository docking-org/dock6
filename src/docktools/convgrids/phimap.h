c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
c
 
C             PHIMAP MAP MODULE

c DelPhi map

c also see chemscore.h and mscore.h


c integer/angstrom conversions:

c     scale: from delphi output.
      real scale

c oldmid (center of DelPhi phimap)
      real oldmid(3)

      character*80 phifil

c DelPhi
      real  phimap(NSIZE, NSIZE, NSIZE)

      common /phiv/scale, oldmid 

      common /phic/phifil

      common /phia/phimap

 
