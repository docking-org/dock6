# Makefile for DOCK chemgrid directory

include ../../../install/rules.h
include ../../../install/config.h

BIN = ../../../bin
PROGS =  nchemgrid_GB

OBJS1 = nchemgrid.o dist.o nddist.o ndconst.o ngrdout.o nparmrec.o \
       npolariz.o second.o

all: $(PROGS)

install: all
	mv $(PROGS) $(BIN)

clean:
	/bin/rm -f $(OBJS1)

realclean: clean
	/bin/rm -f $(PROGS)

uninstall:
	cd $(BIN) && /bin/rm -f $(PROGS)

nchemgrid_GB: $(OBJS1)
	$(FC) $(FFLAGS) $(OBJS1) -o $@

nchemgrid.o: nchemgrid.h
nddist.o: nchemgrid.h
ndconst.o:  nchemgrid.h
ngrdout.o:  nchemgrid.h
nparmrec.o: nparmrec.h
nparmrec.sybH.o: nparmrec.h
npolariz.o: nchemgrid.h
npolariz.debug.o: nchemgrid.h
