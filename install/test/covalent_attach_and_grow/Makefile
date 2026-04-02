include ../default_targets.h

test: 
	make -e DOCK_VERBOSE='-v' covalent.dockout covalent_trans_rot.dockout

clean: clean_dock
	/bin/rm -f covalent_output_scored.mol2 covalent.dockout debug.mol2 
