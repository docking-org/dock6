include ../default_targets.h

test: hawkins_gbsa1.dockout sasa1.dockout \
      gbsa_grid_generation hawkins_gbsa2.dockout sasa2.dockout 

gbsa_grid_generation:
	cd ../gbsa_grid_generation && make test

clean: clean_dock
	/bin/rm -f hawkins_gbsa*mol2 
