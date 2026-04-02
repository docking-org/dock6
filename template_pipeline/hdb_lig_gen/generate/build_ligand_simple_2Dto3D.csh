#!/bin/csh
# Written by Trent Balius at FNLCR on Feb 6, 2020.  
# This script is to smiplify the Database generation or at less make if easer for me to modify it.   
#

set DOCK6BASE = ${DOCK_HOME}

  if (-e db_build_working) then
      echo "Remove this directory: db_build_working"
      exit
  endif

  mkdir db_build_working
  cd db_build_working
  cp ../$1 .

  # source env
  #source ~/.cshrc.python3
  source ~baliuste/zzz.programs/openbabel/env.csh
  #source /nfs/soft/openbabel/current/env.csh
  source ~baliuste/zzz.programs/jchem/env.csh
  #source /nfs/soft/jchem/current/env.csh
  source /home/baliuste/zzz.programs/corina/env.csh
  #source /nfs/soft/corina/current/env.csh

  #echo "I AM HERE"
  ## Step 1. process smiles with chemAxon to protonate and tautomerize.  
  echo "filename = $1"
#  set smilist = `cat $1 | awk '{print "\""$1"\""}'`   
#  echo "$smilist"
  #set PH = 7.2
  set PH = 7.2
  set TAUTOMER_LIMIT = 30
  #set PROTOMER_LIMIT = 30
  set PROTOMER_LIMIT = 10
  set TAUT_PROT_CUTOFF = 1
  set START = 1 

  set CXCALCEXE = `which cxcalc`  
  set MOLCONVERTEXE = `which molconvert` 

  sed 's/\s\+/ /g' "${1}" | \
        ${CXCALCEXE} -g dominanttautomerdistribution -H "${PH}" -C false -t tautomer-dist | \
        ${MOLCONVERTEXE} sdf -g -c "tautomer-dist>=${TAUTOMER_LIMIT}" | \
        ${CXCALCEXE} -g microspeciesdistribution -H $PH -t protomer-dist | \
        ${MOLCONVERTEXE} smiles -g -c "protomer-dist>=${PROTOMER_LIMIT}" -T name:tautomer-dist:protomer-dist | \
        awk -v "cutoff=${TAUT_PROT_CUTOFF}" -v "start=${START}" '{ if (NR == 1 && start < 2) { print $0, "score" } else { score = ($3 * $4)/100 ; if (score >= cutoff) { print $0, score } } }'  > prot-taut.info
        

 #awk 'BEGIN{count=0}{if(count==0){count=1}else{print $0}}' prot-taut.info | sort -r -n -k2 | awk '{print $1 " " $2}' > prot-taut2.info
  awk 'BEGIN{count=0}{if(count==0){count=1}else{print $0}}' prot-taut.info | awk '{print $1 " " $2}' | sort -u | sort -n -k2 > prot-taut2.info
  
  # step 2. convert smiles to mol2 files.  
  

  split -a 8 -l 1 prot-taut2.info  prot-taut_split_

  #mv prot-taut_split_aaaaaaaa header # move the header
  set count = 0
  #foreach file (`ls prot-taut_split_???????[bcdghijklmnopqrstuvwxyz]`)

  set mountdir = `pwd`
  touch dirlist
  foreach file (`ls prot-taut_split_????????`)
     cd $mountdir
     #set name = `awk -F'\t' '{print $2}' $file`
     set name = `awk  '{print $2}' $file`
     set newname =  "${name}_$count"
     set workdir = ${mountdir}/${name}
     if !(-e $workdir) then 
        mkdir $workdir
        echo $name >> dirlist # remember all of the dir that we make
     endif
     
     echo "$newname"
     if (-e "$workdir/$newname.smi") then
        "Error. "
        exit
     endif
     mv "$file" "$workdir/$newname.smi"
    #awk '{printf "name.txt 0 %s %s | NO_LONG_NAME\n",$2,$1}' "$workdir/$newname.smi" > "$workdir/$newname.name.txt"
     awk '{printf "name.txt 0 %16s %s | NO_LONG_NAME\n",substr($2, 1, 16),$1}' "$workdir/$newname.smi" > "$workdir/$newname.name.txt"
     cd $workdir
     #/home/baliuste/zzz.programs/corina/corina -i t=smiles -o t=mol2 -d rc,flapn,de=6,mc=1,wh $newname.smi $newname.mol2
     corina -i t=smiles -o t=mol2 -d rc,flapn,de=6,mc=1,wh $newname.smi $newname.mol2
      
     @ count = $count + 1
  end
  #awk '{print $1 " " $2}' prot-taut.info > prot-taut.smi
  
  #/home/baliuste/zzz.programs/corina/corina -i t=smiles -o t=mol2 -d rc,flapn,de=6,mc=1,wh prot-taut.smi prot-taut.mol2

  # Step 3. Run Amsol. 
  
  foreach  dir (`cat ${mountdir}/dirlist`)
     cd ${mountdir}/${dir}/
     foreach mol2 (`ls *.mol2`)
        ls -l $mol2
        ${DOCK6BASE}/template_pipeline/hdb_lig_gen/amsol/calc_solvation.csh $mol2 
        mv output.mol2 ${mol2:r}_output.mol2
        mv output.solv ${mol2:r}_output.solv
        cat ${mol2:r}_output.mol2 > ${mol2:r}_output_solv.mol2
        echo "@<TRIPOS>SOLVATION" >> ${mol2:r}_output_solv.mol2
        cat ${mol2:r}_output.solv >> ${mol2:r}_output_solv.mol2
     end
  end
  
