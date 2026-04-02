# DOCK test Makefile.

include ./default_targets.h


# Run all tests silently cleaning first because dependencies are incomplete.
test: silentclean timed_test

# Run all tests without cleaning.
dirtytest: timed_test

# Run all DOCK tests
dock: sampling traditional gbsa scoresofscores filtering \
	parallel chemgrid denovo amberdock

# Run DOCK database filtering tests
filtering:
	cd database_filter && make -e test

# Run DOCK sampling related tests
sampling:
	cd anchor_and_grow && make -e test
	cd bump_filter && make -e test
	cd internal_energy && make -e test
	cd conf_gen_int_score && make -e test
	cd ligand_io && make -e test
	cd ligand_orienting && make -e test
	cd minimization && make -e test

# Run traditional DOCK scores tests
traditional:
	cd contact_score && make -e test
	cd continuous_score && make -e test
	cd grid_score && make -e test
#cd primary_secondary_score && make -e test

# Run DOCK gbsa tests
gbsa:
	cd hawkins_gbsa_score && make -e test
	cd sasa_score && make -e test
	cd zou_gbsa_score && make -e test

# Run DOCK descriptor and footprint scores tests
scoresofscores:
	cd footprint_score && make -e test
	cd fp_multigrid_generation && make -e test
	cd fps_multigrid_score && make -e test
# descriptor_score also depends on fp_multigrid_generation
	cd descriptor_score && make -e test
	cd pharmacophore_score && make -e test

# Run Amber DOCK tests
amberdock:
	cd amber_score_181l && make -e test
	cd amber_score_1a28 && make -e test
	cd amber_score_1lgu && make -e test
	cd amber_score_multiple_ligands && make -e test


# Run chemgrid DOCK tests
chemgrid:
	cd chemgrid_score && make -e test

#Run DOCK denovo test
denovo:
	cd de_novo && make -e test
	cd genetic && make -e test
	cd fraglib_gen && make -e test

# Run DOCK parallel tests
parallel:
	@echo "MPI_HOME is defined to ${MPI_HOME}"
	cd mpi && make -e test

mpi:	parallel

# Run tests on the utilities
utils:
	cd sphere_generation && make test
	cd grid_generation && make test
	cd gbsa_grid_generation && make test
	cd nchemgrid_GB_grid_generation && make test
	cd nchemgrid_SA_grid_generation && make test
	cd solvgrid_generation && make test
	cd chemgrid_generation && make test

timed_test: 
	@echo "Starting validation of "
	@echo ${DOCK_VERSION}
	@echo "at `date`."
	@echo ""
	make utils; make dock
	@echo ""
	@echo "Validation of "
	@echo ${DOCK_VERSION}
	@echo "is complete at `date`."
	@echo ""

# Run all DOCK tests with vanilla valgrind
valgrinddock:
	make -e DOCK_BIN="valgrind --track-origins=yes ../../../bin" dock


# Clean all tests by invoking make clean in all subdirectories that have
# a non-zero-sized Makefile.
clean:
	@( \
	find . -name "Makefile" ! -size 0c -print | grep -v '\./Makefile' | \
	while read path ;\
	do \
	    directory=`dirname $$path` ;\
	    cd $$directory; make clean; cd .. ;\
	done ;\
	)

silentclean:
	@make clean 1> /dev/null 2>&1

