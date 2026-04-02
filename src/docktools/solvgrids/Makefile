# Makefile for DOCK chemgrid directory

include ../../../install/rules.h
include ../../../install/config.h

BIN = ../../../bin
PROGS =  solvgrid solvmap

OBJS1 = chemgrid.o dist.o ddist.o dconst.o grdout.o parmrec.o 

all: $(PROGS)

install: all
	mv $(PROGS) $(BIN)

clean:
	/bin/rm -f $(OBJS1) $(OBJS2)

realclean: clean
	/bin/rm -f $(PROGS)

uninstall:
	cd $(BIN) && /bin/rm -f $(PROGS)

solvgrid: $(OBJS1)
	$(FC) $(FFLAGS) $(OBJS1) -o $@

chemgrid.o: chemgrid.h
ddist.o: chemgrid.h
dconst.o:  chemgrid.h
grdout.o:  chemgrid.h
parmrec.o: parmrec.h
