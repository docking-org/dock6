#!/bin/csh
# Written by Trent Balius at FNLCR on Feb 6, 2020.  
# This script is to smiplify the Database generation or at less make if easer for me to modify it.   
# This only generated a mol2 file with one conformation for each molecules in smiles.  
#

set DOCK6BASE = ${DOCK_HOME}

  if (-e db_build_working) then
      echo "Remove this directory: db_build_working"
      exit
  endif

  mkdir db_build_working
  cd db_build_working
  cp ../$1 .

  set ANGLE = $2

  source /home/baliuste/zzz.programs/jchem/env.csh
  source /home/baliuste/zzz.programs/corina/env.csh
  #source /nfs/soft/jchem/current/env.csh
  #source /nfs/soft/corina/current/env.csh
  set TEBSCRIPTS = "/home/baliuste/zzz.github/teb_scripts_programs"
  

  ## Step 1. process smiles with chemAxon to protonate and tautomerize.  
  echo "filename = $1"
  echo "Angle = $2"
#  set smilist = `cat $1 | awk '{print "\""$1"\""}'`   
#  echo "$smilist"
  set PH = 7.2
  set TAUTOMER_LIMIT = 30
  set PROTOMER_LIMIT = 30
  set TAUT_PROT_CUTOFF = 1
  set START = 1 

  set CXCALCEXE = `which cxcalc`  
  set MOLCONVERTEXE = `which molconvert` 


  sed 's/\s\+/ /g' "${1}" | \
        ${CXCALCEXE} -g dominanttautomerdistribution -H "${PH}" -C false -t tautomer-dist | \
        ${MOLCONVERTEXE} sdf -g -c "tautomer-dist>=${TAUTOMER_LIMIT}" | \
        ${CXCALCEXE} -g microspeciesdistribution -H $PH -t protomer-dist | \
        ${MOLCONVERTEXE} smiles -g -c "protomer-dist>=${PROTOMER_LIMIT}" -T name:tautomer-dist:protomer-dist | \
        awk -v "cutoff=${TAUT_PROT_CUTOFF}" -v "start=${START}" '{ if (NR == 1 && start < 2) { print $0, "score" } else { score = ($3 * $4)/100 ; if (score >= cutoff) { print $0, score } } }' > prot-taut.info

  awk 'BEGIN{count=0}{if(count==0){count=1}else{print $0}}' prot-taut.info | awk '{print $1 " " $2}' | sort -u | sort -n -k2 > prot-taut2.info
  
  # step 2. convert smiles to mol2 files.  
  

  #split -a 8 -l 1 prot-taut.info  prot-taut_split_
  split -a 8 -l 1 prot-taut2.info  prot-taut_split_

  #mv prot-taut_split_aaaaaaaa header # move the header
  set count = 0
  #foreach file (`ls prot-taut_split_???????[bcdghijklmnopqrstuvwxyz]`)

  set mountdir = `pwd`
  touch dirlist
  foreach file (`ls prot-taut_split_????????`)
     cd $mountdir
     #set name = `awk -F'\t' '{print $2}' $file`
     set name = `awk '{print $2}' $file`
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
     cd $workdir
     #/home/baliuste/zzz.programs/corina/corina -i t=smiles -o t=mol2 -d rc,flapn,de=6,mc=1,wh $newname.smi $newname.mol2
     corina -i t=smiles -o t=mol2 -d rc,flapn,de=6,mc=1,wh $newname.smi $newname.mol2
      
     @ count = $count + 1
  end
  #awk '{print $1 " " $2}' prot-taut.info > prot-taut.smi
  
  #/home/baliuste/zzz.programs/corina/corina -i t=smiles -o t=mol2 -d rc,flapn,de=6,mc=1,wh prot-taut.smi prot-taut.mol2

  # Step 3. Run Amsol. 
  source /home/baliuste/zzz.programs/openbabel/env.csh
  #source /home/baliuste/.cshrc.python3
  
  foreach  dir (`cat ${mountdir}/dirlist`)
     cd ${mountdir}/${dir}/
     foreach mol2 (`ls *.mol2`)
        #~/zzz.github/DOCK/template_pipeline/hdb_lig_gen/amsol/calc_solvation.csh $mol2 
        $DOCK6BASE/template_pipeline/hdb_lig_gen/amsol/calc_solvation.csh $mol2 
        mv output.mol2 ${mol2:r}_output.mol2
        mv output.solv ${mol2:r}_output.solv
     end
  end

  # Step 4.  convert Si to Du.  Ajust bond angle. 
  #source ~/.bashrc.python3

  #copy files gotten form searching zinc from batch.


  #~/zzz.programs/openbabel/install/bin/obabel -d -ismi BFCARM.smi_sel_92 -osmi -O BFCARM_sel_92_clean.smi

  #python ~/zzz.github/teb_scripts_programs/zzz.scripts/simple_reaction_file.py "[SiH2:5][SiH3].[C:1]=[C:2]-[C:3]=[O:4]>>[SiH3][SiH2:5]-[C:1]-[C:2]-[C:3]=[O:4]" "[SiH2][SiH3]" BFCARM_sel_92_clean.smi temp.smi

  #source ~/.bashrc.python2
  #echo "I AM HERE"
  #source ~/.cshrc.python2
  #echo "I AM HERE"

  foreach  dir (`cat ${mountdir}/dirlist`)
     cd ${mountdir}/${dir}/
     foreach mol2 (`ls *_output.mol2`) 
         echo $mol2
         #echo $mol2:h
         echo $mol2:r
         #echo $mol2:t
         python ${TEBSCRIPTS}/zzz.scripts/mol_covalent_Si_to_Du.py $mol2 ${mol2:r}_Du $ANGLE > mol_covalent_Si_to_Du.log
         #python ~/zzz.github/teb_scripts_programs/zzz.scripts/mol_covalent_Si_to_Du_solv.py $mol2 ${mol2:r}.solv ${mol2:r}_Du $ANGLE > mol_covalent_Si_to_Du.log
     end
  end



  # Step 5. 
