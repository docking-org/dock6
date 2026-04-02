include ../default_targets.h

#test: fresh_solvgrid_generation chemgrid1.dockout chemgrid2.dockout chemgrid3.dockout chemgrid4.dockout
test:
	@# run this test only if the DelPhi executable exists
	@(if [ -x "$(DELPHI_PATH)" ] ; then \
		make fresh_solvgrid_generation ;\
		make chemgrid1.dockout ;\
		make chemgrid2.dockout ;\
		make chemgrid4.dockout ;\
	else \
		echo "" ;\
		echo "The DelPhi executable was not found." ;\
		echo "    Define DELPHI_PATH to the full path of DelPhi." ;\
		echo "    Skipping this test." ;\
		echo "" ;\
	fi ;\
	)


fresh_solvgrid_generation: 
	cd ../solvgrid_generation && make clean && make test

clean: clean_dock
	/bin/rm -f ligand_scored.mol2
