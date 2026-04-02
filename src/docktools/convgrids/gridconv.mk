include ../../../install/rules.h
include ../../../install/config.h

PROGS =	grid-convert
BIN =	../../../bin

OBJS=	grid.o \
	io_grid.o io_gridf.o \
	label.o label_node.o label_chem.o label_vdw.o \
	parm.o parm_grid.o \
	utility.o 

all:	$(PROGS)

install: clean all
	mv $(PROGS) $(BIN)

clean:
	/bin/rm -f $(OBJS)

realclean: clean
	/bin/rm -f $(PROGS)

uninstall:
	cd $(BIN); /bin/rm -f $(PROGS)

grid-convert: $(OBJS)
	$(FC) $(LINK_WITHOUT_FORTRAN_MAIN) $(FFLAGS) $(OBJS) $(LIBS) -o $@ \
	    $(DOCKBUILDFLAGS)

grid.o :	 define.h utility.h mol.h global.h score.h \
		 label.h io.h io_receptor.h io_grid.h grid.h \
		 score_grid.h parm_grid.h parm.h
