#!/bin/csh -f

if ( $#argv != 1 ) then
    echo
    echo "Submit jobs to SLURM as array jobs using specified DOCK version."
    echo
    echo "usage: subdock.csh path/to/dock_executable"
    echo
    exit 1
endif

set dock = "$1"

setenv DOCK6 $dock

if ( ! -e dirlist ) then
    echo "Error: Cannot find dirlist, the list of subdirectories!"
    echo "Exiting!"
    exit 1
endif

set dirnum=`cat dirlist | wc -l`
#qsub -t 1-$dirnum  $DOCKBASE/docking/submit/rundock_pbs_array.csh "$dock"
#sbatch --array=1-$dirnum  $DOCKBASE/docking/submit/rundock_slurm_array.csh 
#set scriptpath = /mnt/projects/RAS-CompChem/static/home/work/RAS/sos_pocket/6GJ8/dock6_large_scale_zinc22/submit_scripts_oel8
#set scriptpath = $DOCK6BASE/template_pipeline/VS_submit_scripts
set scriptpath = ${DOCK_HOME}/template_pipeline/VS_submit_scripts
sbatch --array=1-$dirnum  ${scriptpath}/rundock_slurm_array.csh 
