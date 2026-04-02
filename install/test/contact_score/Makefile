include ../default_targets.h

test: grid_generation contact1.dockout contact2.dockout

grid_generation:
	cd ../grid_generation && make test

clean: clean_dock
	/bin/rm -f contact1*mol2 contact2*mol2
