#  Makefile for DOCK

include ../../install/rules.h
include ../../install/config.h

BIN= ../../bin
PROGS= dock6$(DOCK_SUFFIX)

OBJS= amber_typer.o base_grid.o base_mpi.o base_score.o conf_gen_ag.o conf_gen_cg.o conf_gen_dn.o conf_gen_dn_ga.o conf_gen_ga.o conf_gen_hdb.o\
      dock.o dockmol.o filter.o fingerprint.o fraggraph.o fragment.o gasteiger.o grid.o hungarian.o library_file.o \
      master_conf.o master_score.o orient.o sasa.o score.o score_amber.o score_chemgrid.o \
      score_descriptor.o score_footprint.o score_ph4.o score_volume.o score_gist.o score_hbond.o score_internal.o score_multigrid.o score_sasa.o \
      score_solvent.o simplex.o sphere.o trace.o utils.o

ifeq ($(RDKIT),yes)
OBJS += rdtyper.o
endif

all:  $(PROGS)

install:  all
	mv $(PROGS) $(BIN)

dock:  dock6$(DOCK_SUFFIX)

clean:
	/bin/rm -f $(OBJS)
	cd nab && make clean
	-cd gzstream && make cleanall

realclean:  clean
	/bin/rm -f $(PROGS)

# removes all possible dock executables
uninstall:
	cd $(BIN) && /bin/rm -f dock6 dock6.mpi $(PROGS)

dock6$(DOCK_SUFFIX):  configured $(OBJS) nab_objects gzstream_lig
	$(LOAD) $(CFLAGS) -o $@ $(OBJS) nab/*.o gzstream/gzstream.o -lz \
		$(LIBS) $(DOCKBUILDFLAGS)

nab_objects:
	cd nab && make all

# wget https://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz
# tar -xzvf gzstream.tgz
gzstream_lig:
	cd gzstream && make

CONFIG_COMMAND=cd ../../install; ./configure
CONFIG_FILE=../../install/config.h
configured:
	@# abort if the configuration file does not exist
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
	    echo "Error: $(CONFIG_COMMAND) must be executed before make !" ;\
	    exit 2 ;\
	fi ;\
	)

# dependencies generated on hamlet via
# makedepend -Y *.cpp *.C



amber_typer.o: amber_typer.h utils.h dockmol.h fingerprint.h
base_grid.o: base_grid.h utils.h trace.h trace.inl
base_mpi.o: base_mpi.h dockmol.h utils.h trace.h trace.inl
base_score.o: amber_typer.h utils.h base_score.h dockmol.h trace.h trace.inl
conf_gen_ag.o: conf_gen_ag.h dockmol.h utils.h fingerprint.h master_score.h
conf_gen_ag.o: score.h base_score.h score_multigrid.h score_footprint.h
conf_gen_ag.o: score_ph4.h score_volume.h score_hbond.h score_internal.h
conf_gen_ag.o: score_descriptor.h hungarian.h score_gist.h score_amber.h
conf_gen_ag.o: nab/nab.h defreal.h nab/nabtypes.h nab/defreal.h
conf_gen_ag.o: score_chemgrid.h base_grid.h score_solvent.h grid.h
conf_gen_ag.o: score_sasa.h sasa.h simplex.h trace.h trace.inl
conf_gen_cg.o: sphere.h utils.h conf_gen_cg.h dockmol.h conf_gen_ag.h
conf_gen_cg.o: fingerprint.h master_score.h score.h base_score.h
conf_gen_cg.o: score_multigrid.h score_footprint.h score_ph4.h score_volume.h
conf_gen_cg.o: score_hbond.h score_internal.h score_descriptor.h hungarian.h
conf_gen_cg.o: score_gist.h score_amber.h nab/nab.h defreal.h nab/nabtypes.h
conf_gen_cg.o: nab/defreal.h score_chemgrid.h base_grid.h score_solvent.h
conf_gen_cg.o: grid.h score_sasa.h sasa.h simplex.h trace.h trace.inl
conf_gen_dn.o: conf_gen_dn.h amber_typer.h utils.h dockmol.h fragment.h
conf_gen_dn.o: master_score.h score.h base_score.h score_multigrid.h
conf_gen_dn.o: score_footprint.h score_ph4.h score_volume.h score_hbond.h
conf_gen_dn.o: score_internal.h score_descriptor.h fingerprint.h hungarian.h
conf_gen_dn.o: score_gist.h score_amber.h nab/nab.h defreal.h nab/nabtypes.h
conf_gen_dn.o: nab/defreal.h score_chemgrid.h base_grid.h score_solvent.h
conf_gen_dn.o: grid.h score_sasa.h sasa.h fraggraph.h orient.h sphere.h
conf_gen_dn.o: simplex.h conf_gen_ga.h gasteiger.h trace.h trace.inl
conf_gen_dn_ga.o: conf_gen_dn_ga.h amber_typer.h utils.h dockmol.h master_score.h
conf_gen_dn_ga.o: score.h base_score.h score_multigrid.h score_footprint.h score_ph4.h score_descriptor.h
conf_gen_dn_ga.o: score_amber.h nab/nab.h defreal.h nab/nabtypes.h nab/defreal.h
conf_gen_dn_ga.o: score_chemgrid.h base_grid.h score_solvent.h grid.h
conf_gen_dn_ga.o: score_sasa.h sasa.h simplex.h fingerprint.h fraggraph.h fragment.h hungarian.h
conf_gen_dn_ga.o: gasteiger.h orient.h
conf_gen_ga.o: conf_gen_ga.h amber_typer.h utils.h dockmol.h fragment.h
conf_gen_ga.o: master_score.h score.h base_score.h score_multigrid.h
conf_gen_ga.o: score_footprint.h score_ph4.h score_volume.h score_hbond.h
conf_gen_ga.o: score_internal.h score_descriptor.h fingerprint.h hungarian.h
conf_gen_ga.o: score_gist.h score_amber.h nab/nab.h defreal.h nab/nabtypes.h
conf_gen_ga.o: nab/defreal.h score_chemgrid.h base_grid.h score_solvent.h
conf_gen_ga.o: grid.h score_sasa.h sasa.h orient.h sphere.h simplex.h
conf_gen_ga.o: conf_gen_ag.h conf_gen_dn.h fraggraph.h gasteiger.h trace.h
conf_gen_ga.o: trace.inl
conf_gen_hdb.o: conf_gen_hdb.h dockmol.h utils.h master_score.h score.h
conf_gen_hdb.o: base_score.h score_multigrid.h score_footprint.h score_ph4.h
conf_gen_hdb.o: score_volume.h score_hbond.h score_internal.h
conf_gen_hdb.o: score_descriptor.h fingerprint.h hungarian.h score_gist.h
conf_gen_hdb.o: score_amber.h nab/nab.h defreal.h nab/nabtypes.h
conf_gen_hdb.o: nab/defreal.h score_chemgrid.h base_grid.h score_solvent.h
conf_gen_hdb.o: grid.h score_sasa.h sasa.h simplex.h trace.h trace.inl
dock.o: amber_typer.h utils.h master_conf.h dockmol.h conf_gen_cg.h
dock.o: conf_gen_ag.h conf_gen_dn.h fragment.h master_score.h score.h
dock.o: base_score.h score_multigrid.h score_footprint.h score_ph4.h
dock.o: score_volume.h score_hbond.h score_internal.h score_descriptor.h
dock.o: fingerprint.h hungarian.h score_gist.h score_amber.h nab/nab.h
dock.o: defreal.h nab/nabtypes.h nab/defreal.h score_chemgrid.h base_grid.h
dock.o: score_solvent.h grid.h score_sasa.h sasa.h fraggraph.h orient.h
dock.o: sphere.h simplex.h conf_gen_ga.h conf_gen_hdb.h library_file.h
dock.o: base_mpi.h trace.h trace.inl version.h filter.h
dockmol.o: dockmol.h utils.h
filter.o: filter.h dockmol.h utils.h
fingerprint.o: fingerprint.h dockmol.h utils.h
fraggraph.o: dockmol.h utils.h fraggraph.h amber_typer.h master_score.h
fraggraph.o: score.h base_score.h score_multigrid.h score_footprint.h
fraggraph.o: score_ph4.h score_volume.h score_hbond.h score_internal.h
fraggraph.o: score_descriptor.h fingerprint.h hungarian.h score_gist.h
fraggraph.o: score_amber.h nab/nab.h defreal.h nab/nabtypes.h nab/defreal.h
fraggraph.o: score_chemgrid.h base_grid.h score_solvent.h grid.h score_sasa.h
fraggraph.o: sasa.h
fragment.o: dockmol.h utils.h fragment.h amber_typer.h master_score.h score.h
fragment.o: base_score.h score_multigrid.h score_footprint.h score_ph4.h
fragment.o: score_volume.h score_hbond.h score_internal.h score_descriptor.h
fragment.o: fingerprint.h hungarian.h score_gist.h score_amber.h nab/nab.h
fragment.o: defreal.h nab/nabtypes.h nab/defreal.h score_chemgrid.h
fragment.o: base_grid.h score_solvent.h grid.h score_sasa.h sasa.h
gasteiger.o: gasteiger.h dockmol.h utils.h
grid.o: grid.h base_grid.h utils.h trace.h trace.inl
hungarian.o: hungarian.h dockmol.h utils.h
library_file.o: library_file.h base_mpi.h dockmol.h utils.h amber_typer.h
library_file.o: master_score.h score.h base_score.h score_multigrid.h
library_file.o: score_footprint.h score_ph4.h score_volume.h score_hbond.h
library_file.o: score_internal.h score_descriptor.h fingerprint.h hungarian.h
library_file.o: score_gist.h score_amber.h nab/nab.h defreal.h nab/nabtypes.h
library_file.o: nab/defreal.h score_chemgrid.h base_grid.h score_solvent.h
library_file.o: grid.h score_sasa.h sasa.h simplex.h trace.h trace.inl
master_conf.o: master_conf.h dockmol.h utils.h conf_gen_cg.h conf_gen_ag.h
master_conf.o: conf_gen_dn.h amber_typer.h fragment.h master_score.h score.h
master_conf.o: base_score.h score_multigrid.h score_footprint.h score_ph4.h
master_conf.o: score_volume.h score_hbond.h score_internal.h
master_conf.o: score_descriptor.h fingerprint.h hungarian.h score_gist.h
master_conf.o: score_amber.h nab/nab.h defreal.h nab/nabtypes.h nab/defreal.h
master_conf.o: score_chemgrid.h base_grid.h score_solvent.h grid.h
master_conf.o: score_sasa.h sasa.h fraggraph.h orient.h sphere.h simplex.h
master_conf.o: conf_gen_ga.h conf_gen_hdb.h trace.h trace.inl
master_score.o: master_score.h score.h base_score.h utils.h dockmol.h
master_score.o: score_multigrid.h score_footprint.h score_ph4.h
master_score.o: score_volume.h score_hbond.h score_internal.h
master_score.o: score_descriptor.h fingerprint.h hungarian.h score_gist.h
master_score.o: score_amber.h nab/nab.h defreal.h nab/nabtypes.h
master_score.o: nab/defreal.h score_chemgrid.h base_grid.h score_solvent.h
master_score.o: grid.h score_sasa.h sasa.h
orient.o: orient.h dockmol.h utils.h sphere.h
sasa.o: sasa.h score_solvent.h base_score.h utils.h grid.h base_grid.h
sasa.o: score.h dockmol.h
score.o: amber_typer.h utils.h base_score.h dockmol.h grid.h base_grid.h
score.o: score.h
score_amber.o: utils.h dockmol.h amber_typer.h base_score.h score_amber.h
score_amber.o: nab/nab.h defreal.h nab/nabtypes.h nab/defreal.h sphere.h
score_amber.o: nabcode.h nab.h
score_chemgrid.o: dockmol.h utils.h score_chemgrid.h base_grid.h base_score.h
score_chemgrid.o: tnt/tnt.h tnt/tnt_version.h tnt/tnt_math_utils.h
score_chemgrid.o: tnt/tnt_array1d.h tnt/tnt_i_refvec.h tnt/tnt_array2d.h
score_chemgrid.o: tnt/tnt_array3d.h tnt/tnt_array1d_utils.h
score_chemgrid.o: tnt/tnt_array2d_utils.h tnt/tnt_array3d_utils.h
score_chemgrid.o: tnt/tnt_fortran_array1d.h tnt/tnt_fortran_array2d.h
score_chemgrid.o: tnt/tnt_fortran_array3d.h tnt/tnt_fortran_array1d_utils.h
score_chemgrid.o: tnt/tnt_fortran_array2d_utils.h
score_chemgrid.o: tnt/tnt_fortran_array3d_utils.h tnt/tnt_stopwatch.h
score_chemgrid.o: tnt/tnt_subscript.h tnt/tnt_vec.h tnt/tnt_cmat.h
score_descriptor.o: amber_typer.h utils.h base_score.h dockmol.h grid.h
score_descriptor.o: base_grid.h score_descriptor.h score.h score_multigrid.h
score_descriptor.o: score_footprint.h score_ph4.h fingerprint.h hungarian.h
score_descriptor.o: score_volume.h score_gist.h
score_footprint.o: amber_typer.h utils.h base_score.h dockmol.h grid.h
score_footprint.o: base_grid.h score_footprint.h
score_gist.o: amber_typer.h utils.h base_score.h dockmol.h grid.h base_grid.h
score_gist.o: score_gist.h
score_hbond.o: amber_typer.h utils.h base_score.h dockmol.h grid.h
score_hbond.o: base_grid.h score_hbond.h
score_internal.o: amber_typer.h utils.h base_score.h dockmol.h grid.h
score_internal.o: base_grid.h score_internal.h
score_multigrid.o: amber_typer.h utils.h base_score.h dockmol.h grid.h
score_multigrid.o: base_grid.h score.h score_multigrid.h trace.h trace.inl
score_ph4.o: score_ph4.h base_score.h utils.h dockmol.h amber_typer.h
score_sasa.o: amber_typer.h utils.h sasa.h score_solvent.h base_score.h
score_sasa.o: grid.h base_grid.h score.h dockmol.h score_sasa.h
score_solvent.o: amber_typer.h utils.h sasa.h score_solvent.h base_score.h
score_solvent.o: grid.h base_grid.h score.h dockmol.h trace.h trace.inl
score_volume.o: amber_typer.h utils.h base_score.h dockmol.h grid.h
score_volume.o: base_grid.h score_volume.h
simplex.o: master_score.h score.h base_score.h utils.h dockmol.h
simplex.o: score_multigrid.h score_footprint.h score_ph4.h score_volume.h
simplex.o: score_hbond.h score_internal.h score_descriptor.h fingerprint.h
simplex.o: hungarian.h score_gist.h score_amber.h nab/nab.h defreal.h
simplex.o: nab/nabtypes.h nab/defreal.h score_chemgrid.h base_grid.h
simplex.o: score_solvent.h grid.h score_sasa.h sasa.h simplex.h
sphere.o: sphere.h utils.h
trace.o: trace.h trace.inl
utils.o: utils.h
xlogp.o: xlogp.h dockmol.h utils.h
rdtyper.o: rdtyper.h dockmol.h

sffdock6$(DOCK_SUFFIX):  configured $(OBJS) nab_objects
	$(LOAD) $(CFLAGS) -o $@ $(OBJS) \
chirvol.o database.o mme_init.o molio.o molutil.o prm.o regexp.o reslib.o select_atoms.o \
nab/libsff.a \
nab/stubs.o \
nab/wallclock.o \
$(LIBS) $(DOCKBUILDFLAGS)

junk:
	/amber/amber12/AmberTools/src/blas/./dgemm.o \
/amber/amber12/AmberTools/src/blas/./lsame.o \
/amber/amber12/AmberTools/src/blas/./xerbla.o \
/12/lib/libpbsa.a \
	$(LOAD) $(CFLAGS) -o $@ $(OBJS) nab/*.o $(LIBS) $(DOCKBUILDFLAGS)

# DO NOT DELETE
