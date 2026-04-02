#!/bin/csh -f

if ( $#argv != 0 ) then
    echo
    echo "Submit jobs to SLURM as array jobs."
    echo
    echo "usage: submit.csh"
    echo
    exit 1
endif
#set scriptpath = /mnt/projects/RAS-CompChem/static/home/work/RAS/sos_pocket/6GJ8/dock6_large_scale_zinc22/submit_scripts_oel8
#set scriptpath = $DOCK6BASE/template_pipeline/VS_submit_scripts
set scriptpath = ${DOCK_HOME}/template_pipeline/VS_submit_scripts
#$DOCKBASE/docking/submit/subdock_slurm_array.csh $DOCKBASE/docking/DOCK/bin/dock.csh
#${scriptpath}/subdock_slurm_array.csh $DOCK6PATH/bin_oel8/dock6
${scriptpath}/subdock_slurm_array.csh ${DOCK_HOME}/bin/dock6
exit $status
