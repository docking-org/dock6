include ../default_targets.h

test: 
	make -e DOCK_VERBOSE='-v' db2_grid1.dockout db2_grid2.dockout db2_grid3_gist.dockout db2_grid4_gist_full.dockout

clean: clean_dock
	/bin/rm -f *.mol2
