# Makefile for DOCK grid subdirectory

include ../../install/rules.h
include ../../install/config.h

BIN = ../../bin
PROGS =	grid

OBJS2=	grid.o score.o score_grid.o gridpoint.o define.o \
	io.o io_grid.o io_gridf.o mol_prep.o io_receptor.o \
	io_mol2.o io_pdb.o io_ptr.o io_sph.o \
	label.o label_node.o label_chem.o label_vdw.o label_flex.o \
	mol.o search.o transform.o transformf.o \
	parm.o parm_grid.o rotrans.o \
	utility.o vector.o

all:	$(PROGS)

install: all
	mv $(PROGS) $(BIN)

clean:
	/bin/rm -f $(OBJS2)

realclean: clean
	/bin/rm -f $(PROGS)

uninstall:
	cd $(BIN) && /bin/rm -f $(PROGS)

grid: $(OBJS2)
	$(FC) $(LINK_WITHOUT_FORTRAN_MAIN) $(FFLAGS) $(OBJS2) $(LIBS) -o $@ \
	    $(DOCKBUILDFLAGS)

define.o :	 define.h
grid.o :	 define.h utility.h mol.h global.h score.h \
		 label.h io.h io_receptor.h io_grid.h grid.h \
		 score_grid.h parm_grid.h
io.o :	         define.h utility.h mol.h global.h io.h \
		 io_mol2.h io_pdb.h io_ptr.h io_sph.h transform.h \
		 rotrans.h
io_grid.o :	 define.h utility.h mol.h global.h search.h \
		 label.h score.h io_grid.h io_pdb.h
io_ligand.o :	 define.h utility.h mol.h global.h dock.h \
		 search.h label.h score.h score_dock.h io.h \
		 mol_prep.h io_ligand.h transform.h rotrans.h
io_mol2.o :	 define.h utility.h mol.h global.h io_mol2.h
io_pdb.o :	 define.h utility.h mol.h global.h io_pdb.h
io_ptr.o :	 define.h utility.h mol.h global.h io.h \
		 io_ptr.h vector.h transform.h rotrans.h
io_receptor.o :	 define.h utility.h mol.h global.h search.h \
		 label.h score.h io.h mol_prep.h
io_sph.o :	 define.h utility.h mol.h global.h io_sph.h
label.o :	 define.h utility.h mol.h global.h label.h label_node.h \
		 label_chem.h label_vdw.h label_flex.h
label_chem.o :	 define.h utility.h mol.h global.h label_node.h \
		 label_chem.h
label_flex.o :	 define.h utility.h mol.h global.h search.h \
		 label_node.h label_vdw.h label_flex.h transform.h
label_node.o :	 define.h utility.h mol.h global.h label_node.h
label_vdw.o :	 define.h utility.h mol.h global.h label_node.h \
		 label_vdw.h
screen.o :	 define.h utility.h mol.h global.h label.h 
flex.o :	 define.h utility.h mol.h global.h \
		 label.h transform.h score.h flex.h
match.o :	 define.h utility.h mol.h global.h dock.h \
		 search.h label.h match.h io.h mol_prep.h
mol.o :		 define.h utility.h global.h search.h mol.h vector.h
mol_prep.o :	 define.h utility.h mol.h global.h search.h \
		 label.h mol_prep.h
orient.o :	 define.h utility.h mol.h global.h vector.h \
		 search.h transform.h dock.h label.h score.h \
		 score_dock.h match.h orient.h
parm.o :	 define.h utility.h global.h parm.h
parm_dock.o :	 define.h utility.h mol.h global.h dock.h \
		 search.h label.h io.h score.h score_dock.h \
		 match.h orient.h parm.h parm_dock.h
parm_grid.o :	 define.h utility.h mol.h global.h parm.h \
		 label.h score.h score_grid.h grid.h io.h
rank.o :	 define.h utility.h mol.h global.h dock.h \
		 search.h label.h score.h score_dock.h io.h \
		 io_ligand.h transform.h rank.h
rotrans.o :	 define.h vector.h rotrans.h
score.o :	 define.h utility.h mol.h global.h label.h \
		 score.h io_grid.h io_receptor.h vector.h
score_dock.o :	 define.h utility.h mol.h global.h dock.h \
		 search.h label.h score.h score_dock.h simplex.h \
		 io_ligand.h io_grid.h transform.h
score_grid.o :	 define.h mol.h global.h search.h label.h \
		 score.h grid.h score_grid.h
search.o :	 define.h utility.h global.h search.h
simplex.o :	 define.h utility.h global.h
transform.o :	 define.h utility.h mol.h global.h vector.h \
		 transform.h
utility.o :	 define.h
vector.o :	 define.h vector.h

