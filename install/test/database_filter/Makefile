include ../default_targets.h

test: database.dockout

rdkittest:
	@# run this test only if the dock rdkit executable exists
	@(if [ -x "$(DOCK_BIN)/$(DOCK).rdkit" ] ; then \
	    make rdkittest_no_dock6_rdkit_is_an_error ;\
	else \
	    echo "" ;\
	    echo "The dock rdkit executable does not exist." ;\
	    echo "    Skipping this test." ;\
	    echo "" ;\
	fi ;\
	)

rdkittest_no_dock6_rdkit_is_an_error: database.dockrdkitout

clean: clean_dock

