include ../default_targets.h

test:
	@# run this test only if the DelPhi executable exists
	@(if [ -x "$(DELPHI_PATH)" ] ; then \
		make solvgrid_generation ;\
		make delphi_grid_generation ;\
		make grid_convrds ;\
	else \
		echo "" ;\
		echo "The DelPhi executable was not found." ;\
		echo "    Define DELPHI_PATH to the full path of DelPhi." ;\
		echo "    Skipping this test." ;\
		echo "" ;\
	fi ;\
	)

clean: clean_solvgrid_generation clean_delphi_grid_generation clean_grid_convrds
