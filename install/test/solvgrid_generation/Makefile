include ../default_targets.h

# run this test only if a grid producing executable exists
test:
	@(if [ -x "$(DELPHI_PATH)" ] ; then \
	        echo "This DelPhi executable was found:" ;\
	        echo "    $(DELPHI_PATH)" ;\
		make solvgrid_generation ;\
		make sevsolvgrid_generation ;\
		cd hydrogen && make sevsolvgrid_generation_h;\
	elif [ -x "$(QNIFFT_PATH)" ] ; then \
	        echo "This Qnifft executable was found:" ;\
	        echo "    $(QNIFFT_PATH)" ;\
		make solvgrid_generation ;\
		make sevsolvgrid_generation ;\
		cd hydrogen && make sevsolvgrid_generation_h;\
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


# clean: clean_solvgrid_generation clean_chemgrid_generation clean_delphi_grid_generation clean_grid_convert
clean: clean_solvgrid_generation 
