include ../default_targets.h

test: gbsa_grid_generation nchemgrid_GB_grid_generation nchemgrid_SA_grid_generation docking

docking: zou_gbsaprimary.dockout verbosezou_gbsaprimary \
	      orientflexiblezou_gbsaprimary.dockout \
	      orientrigidcontactprimaryzou_gbsasecondarymulti.dockout

verbosezou_gbsaprimary:
	@echo
	@echo "Processing zou_gbsa $(DOCK) test $@ in directory `pwd|sed -e's@^.*/@@'`"
	$(DOCK_BIN)/$(DOCK) -i $@.dockin -o $@.dockout -v
	../dockdif -t 4 $@.dockout.save $@.dockout

gbsa_grid_generation:
	cd ../gbsa_grid_generation && make test

nchemgrid_GB_grid_generation:
	cd ../nchemgrid_GB_grid_generation && make test

nchemgrid_SA_grid_generation:
	cd ../nchemgrid_SA_grid_generation && make test

clean: clean_dock
	/bin/rm -f zou_gbsaprimary*mol2 verbosezou_gbsaprimary*mol2 \
	           orientflexiblezou_gbsaprimary*mol2 \
	           orientflexiblecontactprimaryzou_gbsasecondarymulti*mol2
