#!/bin/csh -fe
# doamsol.csh

# modified by Trent E Balius Nov, 2013
# modified by T.B. Adler, Shoichet group
# modified by Teague Sterling, Sept. 2014 (Fix missing file when antechamber is missing)

set DOCK6BASE = ${DOCK_HOME}

#set VERBOSE = ""
set mol2file = $1
set amsoltlimit = 1m

# obtain the name of the protonated molecule from the mol2-file name:
set mol2_file_FullName = "$1" #e.g. ZINC00000.0.mol2
set constituents = `echo $mol2_file_FullName:q | sed 's/\./ /g'`
set ProtonatedMoleculeName = "$constituents[1].$constituents[2]"
if $?VERBOSE then
echo "ProtonatedMoleculeName: $ProtonatedMoleculeName"
endif
set MoleculeName = $ProtonatedMoleculeName


# obtain the name of the protonated molecule from the mol2-file name:
set mol2_file_FullName = "$1" #e.g. ZINC00000.0.mol2
set constituents = `echo $mol2_file_FullName:q | sed 's/\./ /g'`
set ProtonatedMoleculeName = "$constituents[1].$constituents[2]"
if $?VERBOSE then
echo "ProtonatedMoleculeName: $ProtonatedMoleculeName"
endif
set MoleculeName = $ProtonatedMoleculeName


if -e temp.mol2 then
    echo "warning: temp.mol2. Rewriting link."
endif

if $?VERBOSE then
echo " Preparing AMSOL7.1 input for $mol2file (first: transformation to ZmatMOPAC by openbabel)"
#echo " If there is any trouble, make sure that your DOCK6BASE and OBABELBASE"
echo " If there is any trouble, make sure that your DOCK_HOME and OBABELBASE"
echo " is correctly set in ~/.cshrc or ~/.bashrc"
endif

ln -svfn $mol2file temp.mol2
#ln -svfn $mol2file:r.smi temp.smi
if ! $?AMSOLEXE then
    set AMSOLEXE = $DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/amsol7.1
endif
if $?VERBOSE then
echo "AMSOLEXE is $AMSOLEXE ."
endif

if ! $?OBABELEXE then
   set OBABELEXE=$OBABELBASE/obabel
   if ! -e $OBABELEXE then
       set OBABELEXE=$OBABELBASE/bin/obabel
   endif
endif
if $?VERBOSE then
echo " OBABELEXE is ${OBABELEXE} ."
endif
if ! -e ${OBABELEXE} then
    echo "Couldn't fine OBABLE at ${OBABELEXE}"
    exit -1
endif

# obtain the entire current path:
set current_path = `pwd`
if $?VERBOSE then
echo "current_path ::: $current_path"
endif

# omega does not always assign proper SYBYL atom types (e.g., just S, instead of S.O2 for example): 
# This can lead to obabel warnings during conversion from mol2 to ZmatMOPAC format:
# Therefore, if antechamber is available, we use it to avoid these warnings and we trust antechamber's ability
# to assign correct atom types.

set TEMP_FILE = "${current_path}/temp.mol2"
set TEMP2_FILE = "${current_path}/temp.2.mol2"
if( `where  antechamber || echo ''` == "" ) then
      echo "antechamber (part of ambertools (downloadable for free!)) is not available on your computer";
      echo "obabel might write out a warning since atom types cannot be translated/interpreted correctly.";
      echo "The obabel warnings can be confidently disregarded. They don't affect the docking.";
else
      if $?VERBOSE then
          echo "antechamber is used to create reliable SYBYL atom types in temp.mol2 file: --> temp_AtomTypesFixed.mol2";
      endif
      antechamber -i ${TEMP_FILE} -fi mol2 -o "${current_path}/temp_AtomTypesFixed_bad_coord.mol2" -fo mol2 -at sybyl;
      antechamber -i ${TEMP_FILE} -fi mol2 -a "${current_path}/temp_AtomTypesFixed_bad_coord.mol2" -fa mol2 -ao type -o "${current_path}/temp_AtomTypesFixed.mol2" -fo mol2 -j 0;
      set TEMP_FILE = "${current_path}/temp_AtomTypesFixed.mol2"
endif

if $?VERBOSE then
echo "obabel -i mol2 temp_AtomTypesFixed.mol2 -o mopin -O temp.ZmatMOPAC"
endif
${OBABELEXE} -i mol2 ${TEMP_FILE} -o mopin -O ${current_path}/temp.ZmatMOPAC

# prepare the use of python:
if ! $?PYTHONPATH then
    setenv PYTHONPATH $DOCK6BASE/template_pipeline/hdb_lig_gen/common
else
    setenv PYTHONPATH "$DOCK6BASE/template_pipeline/hdb_lig_gen/common:$PYTHONPATH"
endif

# create AMSOL7.1 input files (SM5.42R calculations in water and hexadecane solvents) using a Z-matrix in MOPAC style:
#
python $DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/make_amsol71_input_use_mol2_charges.py ${current_path}/temp.ZmatMOPAC ${MoleculeName}
#python $DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/make_amsol71_input_use_mol2_charges.py ${current_path}/temp.ZmatMOPAC ${MoleculeName}
#python $DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/make_amsol71_input.py ${current_path}/temp.ZmatMOPAC ${MoleculeName}
#python $DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/make_amsol71_input.py3.py ${current_path}/temp.ZmatMOPAC ${MoleculeName}

# run the AMSOL7.1 calculations:

echo " running AMSOL7.1: SM5.42R (in water solvent) "
timeout $amsoltlimit $SHELL -c "$AMSOLEXE < temp.in-wat > temp.o-wat"
if ( $status ) then
   echo "AMSOL water calculation failed or stalled. Aborting"
   exit -1
endif

echo " running AMSOL7.1: SM5.42R (in hexane solvent) "
timeout $amsoltlimit $SHELL -c "$AMSOLEXE < temp.in-hex > temp.o-hex"
if ( $status ) then
   echo "AMSOL hexadecane calculation failed or stalled. Aborting"
   exit -1
endif

if $?VERBOSE then
echo " extract data from AMSOL7.1 water and hexadecane output files:"
echo " starting process_amsol_mol2.py :"
echo $mol2file
endif
python $DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/process_amsol_mol2.py ${current_path}/temp.o-wat ${current_path}/temp.o-hex ${TEMP_FILE} ${current_path}/output
#python $DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/process_amsol_mol2.py3.py ${current_path}/temp.o-wat ${current_path}/temp.o-hex ${TEMP_FILE} ${current_path}/output

if $?VERBOSE then
echo "process_amsol_mol2.py has finished."
endif

# Ensure (temporary) but expected files are correct
mv -v ${TEMP_FILE} temp-working.mol2
cp -v temp-working.mol2 temp.mol2

#cp temp_AtomTypesFixed.mol2 temp.mol2
