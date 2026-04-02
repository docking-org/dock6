#!/bin/tcsh
#SBATCH -t 24:00:00
#SBATCH --job-name=array_job_test
#SBATCH --partition=norm-oel8
#  #SBATCH --mem-per-cpu=20G
#SBATCH --output=/dev/null

# used mem-per-cpu=40G before fixing the memory leek
if ( $#argv != 0 ) then
    echo
    echo "Run a dock job with SLURM queue."
    echo
    echo "usage: rundock.csh" 
    echo
    exit 1
endif


set uname = `whoami`
if ( $?TMPDIR ) then
    set TMPDIR = /tmp/$uname/
endif

#set pwd = `pwd`

# go to the dir were script was run, this is where the dir file was located 
cd ${SLURM_SUBMIT_DIR} 

set dock = "$DOCK6"

#echo "Copying DOCK to directory: $PWD"
#$DOCKBASE/docking/DOCK/bin/get_dock_files.csh

# switch to subdirectory
#echo "Starting dock in directory: $PWD" >> stderr
set dirarray=`cat dirlist`

if ! $?SLURM_ARRAY_TASK_ID then
    set SLURM_ARRAY_TASK_ID=1
endif

set pth=${SLURM_SUBMIT_DIR}/$dirarray[${SLURM_ARRAY_TASK_ID}]

if ( $?SLURM_JOBID ) then 
    echo "JOB ID: $SLURM_JOBID" >> stderr
endif
if ( $?SLURM_ARRAY_TASK_ID) then
    echo "TASK ID: $SLURM_ARRAY_TASK_ID" >> stderr
endif

set jobdir = ${TMPDIR}/${SLURM_JOBID}/
echo "JOB dir: $jobdir" >> stderr

mkdir -p $jobdir/run 
# now run dock
cd $jobdir/run
#cd $pth
cp $pth/* .
ln -s $pth/../dock6files $jobdir/


touch ligand.db2.gz 
echo "touch ligand.db2.gz" >> stderr 
foreach file ( `cat split_database_index`)
   echo "tar -xzv -O -f $file | gzip -f >> ligand.db2.gz" >> stderr
   tar -xzv -O -f $file | gzip -f >> ligand.db2.gz 
end

echo "HOST: "`hostname` >> stderr

# get the real path with this pwd madness (i.e. dock68 instead of dockenv) 
pushd "$dock:h" > /dev/null
set real_dir = `pwd`
popd > /dev/null

echo "DOCK: $real_dir/$dock:t" >> stderr

#$dock -i dock.in -o dock.out -v >>& stderr
$dock -i dock.in -o dock.out >>& stderr
rm ligand.db2.gz
mv * $pth/
rm ${TMPDIR}/${SLURM_JOBID}/ -r
exit $status
