#!/bin/csh -f
# input:
# $1   receptor coordinates file;
# $2   grid size;
# $3   # of iterations for the first run;
# $4   # of iterations for the second run;
# $5   # of iterations for the third run;
# output:
# delphi map: $1.phi
# log file: delphi.out

set RECCRG = rec+sph.crg
set GRID = 65
set n1 = 50
set n2 = 120
set n3 = 180

set DELPHI = $DELPHI_HOME/delphi
set VDW = vdw.siz
set CRG = amb.crg.oxt
set PRM = genric.prm

#    three step focussing calc.
     if (-e ARCDAT) /bin/rm ARCDAT
     cp $VDW  fort.11
     cp $CRG  fort.12
     cp $RECCRG fort.13
     echo "gsize=$GRID, perfil=20, bndcon=2, linit=$n1\n" >! fort.10
     cat $PRM >> fort.10
     time $DELPHI
     if ($status != 0) exit(100)
     echo "gsize=$GRID, perfil=60, bndcon=3, linit=$n2\n" >! fort.10
     cat $PRM >> fort.10
     mv fort.14 fort.18
     time $DELPHI
     if ($status != 0) exit(100)
     echo "gsize=$GRID, perfil=90, bndcon=3, linit=$n3\n" >! fort.10
     cat $PRM >> fort.10
     mv fort.14 fort.18
     time $DELPHI
     if ($status != 0) exit(100)
/bin/rm fort.1[0-3] fort.18 ARCDAT
mv fort.14 $RECCRG:r.phi
