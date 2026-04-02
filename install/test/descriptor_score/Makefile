include ../default_targets.h

test: Des_SGE.dockout Des_SGE+All.dockout \
	Des_DCE+FPS+All.dockout Des_FPS.dockout Des_FPS+DCE.dockout \
	Des_FMS.dockout Des_FMS+SGE.dockout Des_FMS_txt.dockout \
	Des_Tan+SGE.dockout Des_Hun+SGE.dockout Des_LES.dockout \
	Des_VOS.dockout Des_VOS+SGE.dockout \
	MGtests

MGtests: resid Des_MG.dockout Des_MG+All.dockout

#Explicitly check for bmp and nrg files in lieu of inference rules;
#assume if one is found then all are there.
resid:
	@(if [ ! -f "../fp_multigrid_generation/resid_005.bmp" ] ; then \
	    cd ../fp_multigrid_generation && make test ;\
	fi ;\
	)

clean: clean_dock

