# Makefile for DOCK chemgrid directory

include ../../../install/rules.h
include ../../../install/config.h

BIN = ../../../bin
PROGS =  nchemgrid_SA

OBJS = nchemgrid.o dist.o nddist.o ndconst.o ngrdout.o nparmrec.o \
       npolariz.o nsas.o nsphgrid.o second.o

all: $(PROGS)

install: all
	mv $(PROGS) $(BIN)

clean:
	/bin/rm -f $(OBJS)

realclean: clean
	/bin/rm -f $(PROGS)

uninstall:
	cd $(BIN) && /bin/rm -f $(PROGS)

nchemgrid_SA: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@

nchemgrid.o: nchemgrid.h
nddist.o: nchemgrid.h
ndconst.o:  nchemgrid.h
ngrdout.o:  nchemgrid.h
nparmrec.o: nparmrec.h
npolariz.o: nchemgrid.h
nsas.o: nchemgrid.h
nsphgrid.o: nchemgrid.h
