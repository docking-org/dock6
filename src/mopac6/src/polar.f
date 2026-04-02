      SUBROUTINE POLAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C**********************************************************************
C
C   POLAR SETS UP THE CALCULATION OF THE MOLECULAR ELECTRIC RESPONSE
C   PROPERTIES BY FFHPOL.
C
C**********************************************************************
      CHARACTER*2 ELEMNT
      COMMON /POLVOL/ POLVOL(107)
      COMMON /KEYWRD/ KEYWRD
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR),IDUMY,XPARAM(MAXPAR)
      COMMON /SCFTYP/ EMIN,LIMSCF
      COMMON /TIME  / TIME0
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /CORE  / CORE(107)
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),
     1                    LOCDEP(MAXPAR)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /LAST  / LAST
      COMMON /COORD / COORD(3,NUMATM)
      COMMON /EULER / TVEC(3,3),IDTVEC
      COMMON /SCRACH/ RXYZ(MPACK),XDUMY(MAXPAR**2-MPACK)
      DIMENSION GRAD(MAXPAR),ROTVEC(3,3),DIPVEC(3),
     1          TEMPV(3,3)
      CHARACTER  KEYWRD*241
      LOGICAL LET, LIMSCF
      LIMSCF=.FALSE.
      LET=(INDEX(KEYWRD,'LET').NE.0)
      WRITE (6,10)
   10 FORMAT ('1',20('*'),' FINITE-FIELD POLARIZABILITIES ',
     1        20('*'),//,
     2      '    THE F-F METHOD IS PERFORMED USING BOTH AN ENERGY',/,
     3      '    AND DIPOLE MOMENT EXPANSION.  THESE RESULTS ARE',/,
     4      '    LISTED BELOW AS "E4" AND "DIP", RESPECTIVELY.')
      CALL GMETRY(GEO,COORD)
C
C  Orient the molecule with the moments of inertia.
C  This is done to ensure a unique, reproduceable set of directions.
C  If LET is specified, the input orientation will be used.
C
      IF (.NOT.LET) THEN
         MASS = 1.0D00
         CALL AXIS (COORD,NUMAT,A,B,C,SUMW,MASS,ROTVEC)
         WRITE (6,20)
   20    FORMAT (/' ROTATION MATRIX FOR ORIENTATION OF MOLECULE:'/)
         DO 40 I = 1,3
            WRITE (6,30) (ROTVEC(I,J),J=1,3)
   30       FORMAT (5X,3F12.6)
   40    CONTINUE
C
C  ROTATE ATOMS
C
         DO 70 I = 1,NUMAT
            DO 60 J = 1,3
               SUM = 0.0D00
               DO 50 K = 1,3
                  SUM = SUM + COORD(K,I)*ROTVEC(K,J)
   50          CONTINUE
               GEO(J,I) = SUM
   60       CONTINUE
   70    CONTINUE
         DO 90 I = 1,NUMAT
            DO 80 J = 1,3
               COORD(J,I) = GEO(J,I)
   80       CONTINUE
   90    CONTINUE
         WRITE(6,'(//10X,''CARTESIAN COORDINATES '',/)')
         WRITE(6,'(4X,''NO.'',7X,''ATOM'',9X,''X'',
     1  9X,''Y'',9X,''Z'',/)')
         L=0
         DO 100 I=1,NUMAT
            IF(NAT(I).EQ.99.OR.NAT(I).EQ.107) GOTO 100
            L=L+1
            WRITE(6,'(I6,8X,A2,4X,3F10.4)')
     1           L,ELEMNT(NAT(I)),(COORD(J,L),J=1,3)
  100    CONTINUE
C
C  IF POLYMER, ROTATE TVEC
C  (BEWARE:  THE POLYMER SECTIONS MAY NOT WORK YET)
C
         IF (IDTVEC.GT.0) THEN
            DO 130 I = 1,IDTVEC
               DO 120 J = 1,3
                  SUM = 0.0D00
                  DO 110 K = 1,3
                     SUM = SUM + TVEC(K,I)*ROTVEC(K,J)
  110             CONTINUE
                  TEMPV(J,I) = SUM
  120          CONTINUE
  130       CONTINUE
            DO 150 I = 1,3
               DO 140 J = 1,IDTVEC
                  TVEC(I,J) = TEMPV(I,J)
  140          CONTINUE
  150       CONTINUE
            WRITE (6,160) ((TVEC(J,I),J=1,3),I=1,IDTVEC)
  160       FORMAT (/' NEW TRANSLATION VECTOR:'/,
     1           ' ',3(3F15.5))
         ENDIF
      ENDIF
C
      LAST=1
      NA(1)=99
C
C  SET UP THE VARIABLES IN XPARAM AND LOC, THESE ARE IN CARTESIAN
C  COORDINATES.
C
      NDEP=0
      NUMAT=0
      SUMX=0.D0
      SUMY=0.D0
      SUMZ=0.D0
      DO 180 I=1,NATOMS
         IF((LABELS(I).NE.99).AND.(LABELS(I).NE.107)) THEN
            NUMAT=NUMAT+1
            LABELS(NUMAT)=LABELS(I)
            SUMX=SUMX+COORD(1,NUMAT)
            SUMY=SUMY+COORD(2,NUMAT)
            SUMZ=SUMZ+COORD(3,NUMAT)
            DO 170 J=1,3
  170       GEO(J,NUMAT)=COORD(J,NUMAT)
         ENDIF
  180 CONTINUE
      SUMX=SUMX/NUMAT
      SUMY=SUMY/NUMAT
      SUMZ=SUMZ/NUMAT
      SUMMAX=0.D0
      ATPOL=0.D0
      DO 190 I=1,NUMAT
         IF (LABELS(I).NE.107) THEN
            ATPOL=ATPOL+POLVOL(LABELS(I))
         ENDIF
         GEO(1,I)=GEO(1,I)-SUMX
         IF(SUMMAX.LT.ABS(GEO(1,I))) SUMMAX=ABS(GEO(1,I))
         GEO(2,I)=GEO(2,I)-SUMY
         IF(SUMMAX.LT.ABS(GEO(2,I))) SUMMAX=ABS(GEO(2,I))
         GEO(3,I)=GEO(3,I)-SUMZ
         IF(SUMMAX.LT.ABS(GEO(3,I))) SUMMAX=ABS(GEO(3,I))
  190 CONTINUE
C
      NVAR=0
      NATOMS = NUMAT
      CALL COMPFG(GEO, .TRUE., HEAT0, .TRUE., GRAD, .FALSE.)
      WRITE (6,200) HEAT0
  200 FORMAT (//' ENERGY OF "REORIENTED" SYSTEM WITHOUT FIELD:',
     1        F20.10)
C
      CALL FFHPOL (HEAT0,ATPOL,DIPVEC)
C
      RETURN
      END
