include ../default_targets.h

test: grid_generation bump1.dockout bump2.dockout

grid_generation:
	cd ../grid_generation && make test

clean: clean_dock
	/bin/rm -f bump1*mol2 bump2*mol2
