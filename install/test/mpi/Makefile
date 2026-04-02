include ../default_targets.h

test:
	@# run this test only if the parallel dock executable exists
	@(if [ -x "$(DOCK_BIN)/$(DOCK).mpi" ] ; then \
	    make grid_generation ;\
	    make mpi.dockmpiout ;\
	    make residues.dockmpiout ;\
	else \
	    echo "" ;\
	    echo "The parallel dock executable does not exist." ;\
	    echo "    Skipping this test." ;\
	    echo "" ;\
	fi ;\
	)

rdkittest:
	@# run this test only if the parallel dock rdkit executable exists
	@(if [ -x "$(DOCK_BIN)/$(DOCK).rdkit_mpi" ] ; then \
	    make grid_generation ;\
	    make mpi.dockrdkitmpiout ;\
	    make residues.dockrdkitmpiout ;\
	else \
	    echo "" ;\
	    echo "The parallel dock rdkit executable does not exist." ;\
	    echo "    Skipping this test." ;\
	    echo "" ;\
	fi ;\
	)

grid_generation:
	cd ../grid_generation && make test

clean: clean_dock
	/bin/rm -f mpi.dockmpiout mpi.dockmpiout.[1-9]* \
		mpi.dockrdkitmpiout mpi.dockrdkitmpiout.[1-9]* \
		residues.dockmpiout residues.dockmpiout.[1-9]* \
		residues.dockrdkitmpiout residues.dockrdkitmpiout.[1-9]*

