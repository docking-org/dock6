#!/bin/csh -f

# this script is written by Trent E Balius in 2017/03/17
# Modified in 2022. 
# this script will check that all dock directory docked without problems


if ( $#argv != 1 ) then
    echo
    echo "usage: get_not_finished.csh /path/to/docking "
    echo "or"
    echo "usage: get_not_finished.csh ./ "
    echo
    exit 1
endif



set pwd     = `pwd`

set workdir = $1

echo ${workdir}

cd ${workdir}


if !(-e dirlist) then
   echo "Uhoh. there is no dirlist file in $workdir"
   exit
endif

if (-e ${pwd}/dirlist_new) then
   echo "Uhoh. there is already dirlist_new in $workdir"
   
   exit 1
endif

touch $pwd/dirlist_new


foreach dir (`cat dirlist`) 
   if ! (-e $dir/dock.out) then 
      echo "$dir does not have an dock.out" 
      echo $dir >> ${pwd}/dirlist_new
   else if (`tail -1 $dir/dock.out | grep -c "elapsed"` == 0) then
      echo "$dir no elapsed time line in dock.out"
      echo $dir >> ${pwd}/dirlist_new
#   else if (`gzip -t $dir/test.mol2.gz` == 1) then
#      echo "$dir/poses.mol2 is corrupted"
#      echo $dir >> ${pwd}/dirlist_new
   endif
end


if (`cat ${pwd}/dirlist_new | wc -l` == 0) then
  echo "All ok"
  rm ${pwd}/dirlist_new
  exit 0
else
  echo "Wrote `cat ${pwd}/dirlist_new | wc -l` failed directories to ./dirlist_new"
  echo 'to re-submit just the failed ones:'
  echo ''
  echo '   cp dirlist dirlist_original'
  echo '   mv dirlist_new dirlist'
  echo '   ${DOCK_HOME}/template_pipeline/VS_submit_scripts//submit_slurm_array.csh'
  exit 1
endif
