

## Writen by Trent Balius in Feb 2023 at FNLCR
## mod on Nov 30, 2023

set pwd = `pwd`
set file = $1

source /home/baliuste/zzz.programs/jchem/env.csh

  set CXCALCEXE = `which cxcalc`  
  set MOLCONVERTEXE = `which molconvert` 

$CXCALCEXE stereoisomercount -m 65 $file > stereoisomercount.txt 

#set workdir = $pwd/stereoisomer/

#mkdir $workdir
#cd $workdir

if ! -e $file then
   set file = $pwd/$file
   if ! -e $file then
      echo "error $1 nor $file exist"
      exit
   endif
endif

grep -v "Stereoisomer count" stereoisomercount.txt > temp.txt

paste $file temp.txt | awk '{if($4<=64) {print $1, $2}}' > stereoisomernum_filtered.smi 
paste $file temp.txt | awk '{if($4>64) {print $1, $2}}' > stereoisomernum_filtered2.smi 
rm temp.txt

#cat $file | $CXCALCEXE stereoisomers -m 32 -v true -T | ${MOLCONVERTEXE} smiles -g -T name > ${file:r}_stereoisomers.smi
#cat $file | $CXCALCEXE stereoisomers -f smiles -m 2 -T | $CXCALCEXE stereoisomers -f smiles -T -v true > ${file:r}_stereoisomers.smi
cat stereoisomernum_filtered.smi |  $CXCALCEXE stereoisomers -m 64 -v true -T | ${MOLCONVERTEXE} smiles -g -T name > ${file:r}_stereoisomers.smi 
cat stereoisomernum_filtered2.smi >> ${file:r}_stereoisomers.smi # so that it will at lessed try 1.  

cat ${file:r}_stereoisomers.smi | awk 'BEGIN{count=0}{print $1, $2"_"count;count=count+1}' > temp.smi
mv temp.smi ${file:r}_stereoisomers.smi
#cat $file | $CXCALCEXE stereoisomers -v true -T  > stereoisomers.sdf
#cat stereoisomers.sdf | ${MOLCONVERTEXE} smiles -g  -T name > stereoisomers.smi


