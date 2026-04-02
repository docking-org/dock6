include ../default_targets.h

test: grid_generation minimize1.dockout minimize1_restraint.dockout \
	minimize2.dockout minimize3.dockout minimize4.dockout minimize5.dockout \
	minimize6.dockout minimize6_restraint.dockout \
	minimize1_premin.dockout minimize5_premin.dockout \
	minimize5_restraint.dockout minimize5_restraint_premin.dockout

grid_generation:
	cd ../grid_generation && make test

clean: clean_dock
	/bin/rm -f minimize1*mol2 minimize2*mol2 minimize3*mol2 minimize4*mol2 minimize5*mol2 minimize6*mol2
