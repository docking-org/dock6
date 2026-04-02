# DOCK default macros, rules and targets for test Makefiles.
# This file is for general purpose targets.
# A target for one specific test should be put in the Makefile for that test.

DOCK=dock6
DOCK_BIN=../../../bin
DOCK_VERSION=`cat ../../src/dock/version.h | grep DOCK_VERSION | cut -d'"' -f2`
RELATIVE_DOCK_HOME=../../..
SHELL = /bin/sh


# This is the default target.
all: test


# Create Amber_Dock input files for a combined receptor ligand complex.
amberize: 
	$(DOCK_BIN)/prepare_amber.pl $(LIG) $(REC)

# there is no dock target;
# instead the .dockin.dockout inference rule is used; see below.

grid:
	# Compute scoring grids
	$(DOCK_BIN)/grid -i grid.in -o grid.out
	../dockdif grid.out.save grid.out

nchemgrid_GB:
	# Calculate GB grid for Zou GB/SA score
	$(DOCK_BIN)/nchemgrid_GB > /dev/null
	../dockdif -w -t 4 OUTCHEM.save OUTCHEM
	../dockdif OUTPARM.save OUTPARM

nchemgrid_SA:
	# Calculate SA grid for Zou GB/SA score
	$(DOCK_BIN)/nchemgrid_SA
	../dockdif -w -t 1 OUTCHEM.save OUTCHEM
	../dockdif OUTPARM.save OUTPARM

rec.prmtop: rec.pdb
	make amberize

showbox:
	# Construct box to enclose spheres
	$(DOCK_BIN)/showbox < box.in > /dev/null
	../dockdif box.pdb.save box.pdb

showsphere:
	# Convert selected spheres into pdb format for viewing
	$(DOCK_BIN)/showsphere < select.in > /dev/null
	../dockdif selected_cluster.pdb.save selected_cluster.pdb

sphgen:
	# Construct sphere set in active site
	$(DOCK_BIN)/sphgen
	../dockdif OUTSPH.save OUTSPH

sphere_selector:
	# Select spheres within 3.0 Ang of ligand 
	$(DOCK_BIN)/sphere_selector struct.sph lig.mol2 3.0
	../dockdif -t 3 selected_spheres.sph.save selected_spheres.sph

chemgrid_generation:
	# generate dock3.5.54 chemgrids for CHEMGRID score
	$(DOCK_BIN)/chemgrid
	../dockdif -w OUTCHEM.save OUTCHEM

delphi_grid_generation:
	# generate DelPhi grid for CHEMGRID score
	$(DOCK_BIN)/make_phimap $(DELPHI_PATH) 1> delphi.out
	../dockdif -w -t 1 delphi.out.save delphi.out

solvgrid_generation:
	# generate dock3.5.54 solvgrids for CHEMGRID score
	../../../bin/solvmap
	../dockdif solvmap.save solvmap
	../dockdif -w -t 1 OUTSOLV.save OUTSOLV
	../../../bin/solvgrid
	../dockdif -w -t 1 OUTRDSL.save OUTRDSL

grid_convert:
	# convert dock3.5.54 chemgrids for CHEMGRID score
	../../../bin/grid-convert -i gconv.in 1> gconv.out 2>&1
	../dockdif gconv.out.save gconv.out

grid_convrds:
	# convert dock3.5.54 chemgrids and solvation grid for CHEMGRID score
	../../../bin/grid-convrds -i gconv.in 1> gconv.out 2>&1
	../dockdif gconv.out.save gconv.out

# Error checking and reporting.
parallel_dock_exists:
	@# abort if the parallel dock executable does not exist
	@(if [ ! -x "$(DOCK_BIN)/$(DOCK).mpi" ] ; then \
	    echo "" ;\
	    echo "The parallel dock executable does not exist." ;\
	    echo "    Skipping this test." ;\
	    echo "" ;\
	    exit 3 ;\
	fi ;\
	)


check:
	@# Find and emit the dif files.
	@(find . -name "*.dif" ! -size 0c -print | \
	while read dif ;\
	do \
	    echo $$dif ;\
	    cat $$dif ;\
	    echo ;\
	done ;\
	)

check_ls:
	@# Find and list the dif files.
	@(find . -name "*.dif" ! -size 0c -print | \
	while read dif ;\
	do \
	    echo $$dif ;\
	done ;\
	)


# Clean targets follow the naming rule clean_program_family_name.
clean_amberize:
	/bin/rm -f amberize_complex.*.out amberize_ligand.*.out amberize_receptor.out \
		*.log *.inpcrd *.prmtop *.frcmod *.amber.pdb lig.*.mopac.out \
		lig.*.mol2 *_scored.mol2 lig.gaff.mol2 *.amber_score.mol2 \
		t4lys_mutant_multiple.1.* t4lys_mutant_multiple.2.* \
		tleap.in tleap.out *.dif

clean_dock:
	/bin/rm -f *.dif *.dockout *.dockrdkitout *rejected.mol2 *_conformers.mol2 \
		*_orients.mol2 *_primary_conformers.mol2 *_primary_ranked.mol2 \
		*_primary_scored.mol2 *_ranked.mol2 *_scored.mol2 \
		*_secondary_conformers.mol2 *_secondary_ranked.mol2 \
		*_secondary_scored.mol2 *_footprint_scored.txt *_hbond_scored.txt\
		*_linker.mol2 *_scaffold.mol2 *_rigid.mol2 *_sidechain.mol2\
		*.denovo_build.mol2 *.prune_dump_layer_*mol2 *.root_layer_*mol2

clean_grid:
	/bin/rm -f grid.bmp grid.cnt grid.nrg grid.out grid.out.dif

clean_nchemgrid_GB:
	/bin/rm -f OUTCHEM OUTPARM PDBPARM NEG_INVA PDBCAV inva screen.para \
		zou_grid.* OUTCHEM.dif OUTPARM.dif

clean_nchemgrid_SA:
	/bin/rm -f OUTCHEM OUTPARM PDBPARM \
		zou_grid.* OUTCHEM.dif OUTPARM.dif

clean_showbox:
	/bin/rm -f box.pdb box.pdb.dif

clean_showsphere:
	/bin/rm -f selected_cluster.pdb selected_cluster.pdb.dif

clean_sphere_selector:
	/bin/rm -f selected_spheres.sph selected_spheres.sph.dif

clean_sphgen:
	/bin/rm -f temp1.ms temp2.sph temp3.atc OUTSPH struct.sph OUTSPH.dif

clean_chemgrid_generation:
	/bin/rm -f OUTCHEM OUTPARM PDBPARM chem.* OUTCHEM.dif

clean_delphi_grid_generation:
	/bin/rm -f rec+sph.phi delphi.out delphi.out.dif

clean_solvgrid_generation:
	/bin/rm -f OUTSOLV solvmap distmap.box OUTRDSL OUTPARM PDBPARM chem.bmp \
		chem.dsl chem.esp chem.vdw solvmap.dif OUTSOLV.dif OUTRDSL.dif

clean_grid_convert:
	/bin/rm -f gconv.out gconv.out.dif \
		chem52.bmp chem52.chm chem52.cmg chem52.dsl chem52.phi

clean_grid_convrds:
	/bin/rm -f gconv.out gconv.out.dif \
		chem52.bmp chem52.chm chem52.cmg chem52.dsl chem52.phi

clean_de_novo:
	/bin/rm -f dn_rec.resid* temp.mol2 *torenv.dat

clean_fraglib:
	/bin/rm	-f	*torenv.dat

.SUFFIXES:
.SUFFIXES:  .dockin .dockout .dif .dockmpiout .dockrdkitin .dockrdkitout \
	    .dockrdkitmpiin .dockrdkitmpiout


# Inference rule to create a DOCK output file of extension .dockout from
# a DOCK input file of extension .dockin
.dockin.dockout: 
	@echo
	@echo "Processing $(DOCK) test $*"
	$(DOCK_BIN)/$(DOCK) $(DOCK_VERBOSE) -i $< -o $@
	../dockdif -t 8 $@.save $@
	@# Compare other $(DOCK) output files if they exist
	@( \
        if [ -f LIG1_anchor1_branch1.mol2 ] ; then \
           ../dockdif -t 8 LIG1_anchor1_branch1.mol2.save LIG1_anchor1_branch1.mol2 ;\
        fi ;\
        if [ -f LIG1_anchor1_branch2.mol2 ] ; then \
           ../dockdif -t 8 LIG1_anchor1_branch2.mol2.save LIG1_anchor1_branch2.mol2 ;\
        fi ;\
        if [ -f LIG1_anchor1_branch3.mol2 ] ; then \
           ../dockdif -t 8 LIG1_anchor1_branch3.mol2.save LIG1_anchor1_branch3.mol2 ;\
        fi ;\
        if [ -f LIG1_anchor1_branch4.mol2 ] ; then \
           ../dockdif -t 8 LIG1_anchor1_branch4.mol2.save LIG1_anchor1_branch4.mol2 ;\
        fi ;\
        if [ -f LIG1_anchor1_branch5.mol2 ] ; then \
           ../dockdif -t 8 LIG1_anchor1_branch5.mol2.save LIG1_anchor1_branch5.mol2 ;\
        fi ;\
	if [ -f $*_conformers.mol2 ] ; then \
	   ../dockdif -t 8 $*_conformers.mol2.save $*_conformers.mol2 ;\
	fi ;\
	if [ -f $*_orients.mol2 ] ; then \
	   ../dockdif -t 8 $*_orients.mol2.save $*_orients.mol2 ;\
	fi ;\
	if [ -f $*_ranked.mol2 ] ; then \
	   ../dockdif -t 8 $*_ranked.mol2.save $*_ranked.mol2 ;\
	fi ;\
	if [ -f $*_scored.mol2 ] ; then \
	   ../dockdif -t 8 $*_scored.mol2.save $*_scored.mol2 ;\
	fi ;\
        if [ -f $*_footprint_scored.txt ] ; then \
           ../dockdif -t 8 $*_footprint_scored.txt.save $*_footprint_scored.txt ;\
        fi ;\
        if [ -f $*_hbond_scored.txt ] ; then \
           ../dockdif -t 8 $*_hbond_scored.txt.save $*_hbond_scored.txt ;\
        fi ;\
        if [ -f $*_primary_conformers.mol2 ] ; then \
	   ../dockdif -t 8 $*_primary_conformers.mol2.save $*_primary_conformers.mol2 ;\
	fi ;\
        if [ -f $*_secondary_conformers.mol2 ] ; then \
	   ../dockdif -t 8 $*_secondary_conformers.mol2.save $*_secondary_conformers.mol2 ;\
	fi ;\
        if [ -f $*_primary_scored.mol2 ] ; then \
	   ../dockdif -t 8 $*_primary_scored.mol2.save $*_primary_scored.mol2 ;\
	fi ;\
        if [ -f $*_secondary_scored.mol2 ] ; then \
	   ../dockdif -t 8 $*_secondary_scored.mol2.save $*_secondary_scored.mol2 ;\
	fi ;\
        if [ -f $*_primary_ranked.mol2 ] ; then \
	   ../dockdif -t 8 $*_primary_ranked.mol2.save $*_primary_ranked.mol2 ;\
	fi ;\
        if [ -f $*_secondary_ranked.mol2 ] ; then \
	   ../dockdif -t 8 $*_secondary_ranked.mol2.save $*_secondary_ranked.mol2 ;\
	fi ;\
        if [ -f $*_branch*.mol2 ] ; then \
           ../dockdif -t 8 $*_branch1.mol2.save $*_branch1.mol2 ;\
        fi ;\
        if [ -f $*.denovo_build.mol2 ] ; then \
           ../dockdif -t 8 $*.denovo_build.mol2.save $*.denovo_build.mol2 ;\
        fi ;\
        if [ -f $*_linker.mol2 ] ; then \
           ../dockdif -t 8 $*_linker.mol2.save $*_linker.mol2 ;\
        fi ;\
        if [ -f $*_rigid.mol2 ] ; then \
           ../dockdif -t 8 $*_rigid.mol2.save $*_rigid.mol2 ;\
        fi ;\
        if [ -f $*_scaffold.mol2 ] ; then \
           ../dockdif -t 8 $*_scaffold.mol2.save $*_scaffold.mol2 ;\
        fi ;\
        if [ -f $*_sidechain.mol2 ] ; then \
           ../dockdif -t 8 $*_sidechain.mol2.save $*_sidechain.mol2 ;\
        fi ;\
        if [ -f $*_torenv.dat ] ; then \
           ../dockdif -w -t 1 $*_torenv.dat.save $*_torenv.dat ;\
        fi ;\
        if [ -f $*.restart0001.mol2.save ] ; then \
           ../dockdif -w -t 3 $*.restart0001.mol2.save $*.restart1.mol2 ;\
        fi ;\
	)

# Inference rule to create a DOCK output file of extension .dockmpiout from
# a DOCK input file of extension .dockin
.dockin.dockmpiout: 
	@echo
	@echo "Processing parallel $(DOCK) test $*"
	@# First assign the number of MPI processes, DOCK_PROCESSES, if needed;
	@# then search in order MPI_HOME, MPICH_HOME, and the path for mpirun.
	@( \
	if [ -z "$(DOCK_PROCESSES)" ] ; then \
	    echo "Environment variable DOCK_PROCESSES is not defined." ;\
	    DOCK_PROCESSES=2 ;\
	fi ;\
	echo "Using $$DOCK_PROCESSES MPI processes for DOCK." ;\
	if [ -n "$(MPI_HOME)" ] ; then \
	    echo "Environment variable MPI_HOME is defined." ;\
	    if [ -x "$(MPI_HOME)/bin/mpirun" ] ; then \
	        echo "$(MPI_HOME)/bin/mpirun -np $$DOCK_PROCESSES "\
                     "$(DOCK_BIN)/$(DOCK).mpi -i $< -o $@ " ;\
	        $(MPI_HOME)/bin/mpirun -np $$DOCK_PROCESSES $(DOCK_BIN)/$(DOCK).mpi -i $< -o $@  ;\
	    else \
	        echo "Error!  $(MPI_HOME)/bin/mpirun is not executable !" ;\
	        exit 1 ;\
	    fi ;\
	elif [ -n "$(MPICH_HOME)" ]; then \
	    echo "Environment variable MPICH_HOME is defined." ;\
	    if [ -x "$(MPICH_HOME)/bin/mpirun" ] ; then \
	        echo "$(MPICH_HOME)/bin/mpirun -np $$DOCK_PROCESSES "\
                     "$(DOCK_BIN)/$(DOCK).mpi -i $< -o $@ " ;\
	        $(MPICH_HOME)/bin/mpirun -np $$DOCK_PROCESSES $(DOCK_BIN)/$(DOCK).mpi -i $< -o $@  ;\
	    else \
	        echo "Error!  $(MPICH_HOME)/bin/mpirun is not executable !" ;\
	        exit 1 ;\
	    fi ;\
	else \
	    echo "Environment variables MPI_HOME and MPICH_HOME are not defined." ;\
	    echo "    Assuming mpirun is in the PATH." ;\
	    echo "mpirun -np $$DOCK_PROCESSES $(DOCK_BIN)/$(DOCK).mpi -i $< -o $@ " ;\
	    mpirun -np $$DOCK_PROCESSES $(DOCK_BIN)/$(DOCK).mpi -i $< -o $@  ;\
	fi ;\
	)
	../dockdif -t 8 $@.save $@
	@# Compare other $(DOCK) output files if they exist
	@( \
	if [ -f $*_conformers.mol2 ] ; then \
	   ../dockdif -t 8 $*_conformers.mol2.save $*_conformers.mol2 ;\
	fi ;\
	if [ -f $*_orients.mol2 ] ; then \
	   ../dockdif -t 8 $*_orients.mol2.save $*_orients.mol2 ;\
	fi ;\
	if [ -f $*_ranked.mol2 ] ; then \
	   ../dockdif -t 8 $*_ranked.mol2.save $*_ranked.mol2 ;\
	fi ;\
	if [ -f $*_scored.mol2 ] ; then \
	   ../dockdif -t 8 $*_scored.mol2.save $*_scored.mol2 ;\
	fi ;\
	if [ -f $*_primary_conformers.mol2 ] ; then \
	   ../dockdif -t 8 $*_primary_conformers.mol2.save $*_primary_conformers.mol2 ;\
	fi ;\
	if [ -f $*_secondary_conformers.mol2 ] ; then \
	   ../dockdif -t 8 $*_secondary_conformers.mol2.save $*_secondary_conformers.mol2 ;\
	fi ;\
	if [ -f $*_primary_scored.mol2 ] ; then \
	   ../dockdif -t 8 $*_primary_scored.mol2.save $*_primary_scored.mol2 ;\
	fi ;\
	if [ -f $*_secondary_scored.mol2 ] ; then \
	   ../dockdif -t 8 $*_secondary_scored.mol2.save $*_secondary_scored.mol2 ;\
	fi ;\
	if [ -f $*_primary_ranked.mol2 ] ; then \
	   ../dockdif -t 8 $*_primary_ranked.mol2.save $*_primary_ranked.mol2 ;\
	fi ;\
	if [ -f $*_secondary_ranked.mol2 ] ; then \
	   ../dockdif -t 8 $*_secondary_ranked.mol2.save $*_secondary_ranked.mol2 ;\
	fi ;\
	)

# Inference rule to create a DOCK output file of extension .dockrdkitout from
# a DOCK input file of extension .dockrdkitin
.dockrdkitin.dockrdkitout: 
	@echo
	@echo "Processing $(DOCK) with RDKit test $*"
	$(DOCK_BIN)/$(DOCK).rdkit $(DOCK_VERBOSE) -i $< -o $@
	../dockdif -t 8 $@.save $@
	@# Compare other $(DOCK) output files if they exist
	@( \
	if [ -f $*_rd.denovo_build.mol2 ] ; then \
           ../dockdif -t 8 $*_rd.denovo_build.mol2.save $*_rd.denovo_build.mol2 ;\
	fi ;\
        if [ -f $*_rd.restart0001.mol2 ] ; then \
           ../dockdif -w -t 3 $*_rd.restart0001.mol2.save $*_rd.restart0001.mol2 ;\
        fi ;\
	if [ -f $*_rd_db_scored.mol2 ] ; then \
	   ../dockdif -t 8 $*_rd_db_scored.mol2.save $*_rd_db_scored.mol2 ;\
	fi ;\
	)

# Inference rule to create a DOCK output file of extension .dockrdkitmpiout from
# a DOCK input file of extension .dockrdkitmpiin
.dockrdkitmpiin.dockrdkitmpiout:
	@echo
	@echo "Processing parallel $(DOCK) with RDKit test $*"
	@# First assign the number of MPI processes, DOCK_PROCESSES, if needed;
	@# then search in order MPI_HOME, MPICH_HOME, and the path for mpirun.
	@( \
	if [ -z "$(DOCK_PROCESSES)" ] ; then \
	    echo "Environment variable DOCK_PROCESSES is not defined." ;\
	    DOCK_PROCESSES=2 ;\
	fi ;\
	echo "Using $$DOCK_PROCESSES MPI processes for DOCK." ;\
	if [ -n "$(MPI_HOME)" ] ; then \
	    echo "Environment variable MPI_HOME is defined." ;\
	    if [ -x "$(MPI_HOME)/bin/mpirun" ] ; then \
	        echo "$(MPI_HOME)/bin/mpirun -np $$DOCK_PROCESSES "\
                     "$(DOCK_BIN)/$(DOCK).rdkit_mpi -i $< -o $@ " ;\
	        $(MPI_HOME)/bin/mpirun -np $$DOCK_PROCESSES $(DOCK_BIN)/$(DOCK).rdkit_mpi -i $< -o $@  ;\
	    else \
	        echo "Error!  $(MPI_HOME)/bin/mpirun is not executable !" ;\
	        exit 1 ;\
	    fi ;\
	elif [ -n "$(MPICH_HOME)" ]; then \
	    echo "Environment variable MPICH_HOME is defined." ;\
	    if [ -x "$(MPICH_HOME)/bin/mpirun" ] ; then \
	        echo "$(MPICH_HOME)/bin/mpirun -np $$DOCK_PROCESSES "\
                     "$(DOCK_BIN)/$(DOCK).rdkit_mpi -i $< -o $@ " ;\
	        $(MPICH_HOME)/bin/mpirun -np $$DOCK_PROCESSES $(DOCK_BIN)/$(DOCK).rdkit_mpi -i $< -o $@  ;\
	    else \
	        echo "Error!  $(MPICH_HOME)/bin/mpirun is not executable !" ;\
	        exit 1 ;\
	    fi ;\
	else \
	    echo "Environment variables MPI_HOME and MPICH_HOME are not defined." ;\
	    echo "    Assuming mpirun is in the PATH." ;\
	    echo "mpirun -np $$DOCK_PROCESSES $(DOCK_BIN)/$(DOCK).rdkit_mpi -i $< -o $@ " ;\
	    mpirun -np $$DOCK_PROCESSES $(DOCK_BIN)/$(DOCK).rdkit_mpi -i $< -o $@  ;\
	fi ;\
	)
	../dockdif -t 8 $@.save $@
	@# Compare other $(DOCK) output files if they exist
	@( \
	if [ -f $*_rd_ranked.mol2 ] ; then \
	   ../dockdif -t 8 $*_rd_ranked.mol2.save $*_rd_ranked.mol2 ;\
	fi ;\
	)


# This special target declares these targets as phony, ie, not based on
# files with that name; it avoids problems when such files do exist.
.PHONY:  all check check_ls clean clean_dock clean_grid \
         clean_nchemgrid_GB clean_nchemgrid_SA clean_showbox clean_showsphere \
         clean_sphere_selector clean_sphgen \
         dock grid nchemgrid_GB nchemgrid_SA showbox showsphere \
         sphere_selector sphgen \
         move parallel realclean test uninstall

