include ../default_targets.h

# run this test only if a grid producing executable exists
test:
	@(if [ -x "$(DELPHI_PATH)" ] ; then \
	        echo "This DelPhi executable was found:" ;\
	        echo "    $(DELPHI_PATH)" ;\
	        make fresh_grid_generation ;\
	        make chemgrid1.dockout ;\
	        make chemgrid1_all.dockout ;\
	        make chemgrid2.dockout ;\
	        make chemgrid3.dockout ;\
	        make chemgrid2_sev.dockout ;\
	        make chemgrid3_sev.dockout ;\
	        make chemgrid4_descriptor.dockout ;\
	elif [ -x "$(QNIFFT_PATH)" ] ; then \
	        echo "This Qnifft executable was found:" ;\
	        echo "    $(QNIFFT_PATH)" ;\
	        make fresh_grid_generation ;\
	        make chemgrid1.dockout ;\
	        make chemgrid1_all.dockout ;\
	        make chemgrid2.dockout ;\
	        make chemgrid3.dockout ;\
	        make chemgrid2_sev.dockout ;\
	        make chemgrid3_sev.dockout ;\
	        make chemgrid4_descriptor.dockout ;\
	else \
	        echo "" ;\
	        echo "An executable to produce grids was not found." ;\
	        echo "    DELPHI_PATH = $(DELPHI_PATH)" ;\
	        echo "    QNIFFT_PATH = $(QNIFFT_PATH)" ;\
	        echo "    Define DELPHI_PATH or QNIFFT_PATH to its full path." ;\
	        echo "    Skipping this test." ;\
	        echo "" ;\
	fi ;\
	)

fresh_grid_generation: 
	cd ../solvgrid_generation && make clean && make test
	cd ../chemgrid_generation && make clean && make test

clean: clean_dock
	/bin/rm -f ligand_scored.mol2
