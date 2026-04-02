#!/bin/sh

#THIS IS AN EXAMPLE STUB!!!!!!!!
# EXAMPLE USAGE
# mkdir /tmp/3302344
# cd /tmp/330234
# wget "http://zinc11.docking.org/fget.pl?l=0&z=21995818&f=m"; mv "fget.pl?l=0&z=21995818&f=m" 3302344.mol2 
# $DOCK6BASE/template_pipeline/hdb_lig_gen/generate/build_src/hdb_lig_gen_mol2.sh 3302344.mol2 
# Modified by Trent Balius 


MOL=`basename $1 .gz`
MOL=`basename $MOL .mol2`

$OBABELBASE/bin/obabel -h -imol2 -osmi ${MOL}.mol2 -O ${MOL}_OB.smi
head -1 ${MOL}_OB.smi | awk '{if(NF==1){print "None 0 "$1}; if(NF==2){print $2" 0 "$1}; if(NF>=3){print $2" "$3" "$1}}' > name.txt

#cat name.txt
#exit

#python $DOCK6BASE/template_pipeline/hdb_lig_gen/generate/mol2tosmi_rdkit.py ${MOL}.mol2 ${MOL}.smi 
#echo $1:r $1:r > name.txt
#$DOCK6BASE/template_pipeline/hdb_lig_gen/generate/prepare.py $*
#$DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/calc_solvation_formal_charge_from_mol2.csh $1
$DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/calc_solvation_formal_charge_from_mol2.csh $1

#$SOFT/openbabel/current/bin/obabel -imol2 -osmi output.mol2 -O output.smi

#csh $DOCK6BASE/template_pipeline/hdb_lig_gen/generate/dock6_confgen_db2.csh output.mol2 output.solv name.txt 
csh $DOCK6BASE/template_pipeline/hdb_lig_gen/generate/dock6_confgen_db2.csh output.mol2 output.solv name.txt 

