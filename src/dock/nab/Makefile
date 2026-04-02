# Makefile for Dock's NAB
# This has extra targets to handle LEX related problems:
# is_LEX_available and copy.lex.mm_options.c.

include ../../../install/rules.h
include ../../../install/config.h

OBJS=molio.o prm.o sff.o traceback.o select_atoms.o \
memutil.o reslib.o jacobi.o molutil.o \
chirvol.o regexp.o binpos.o time.o rand2.o conjgrad.o database.o \
tleap.o 

all:  $(OBJS)
	$(AR) nab.a $(OBJS)
	$(RANLIB) nab.a

clean:
	/bin/rm -f *.o nab.a lex.mm_options.c

is_LEX_available:
	@# abort if no LEX executable is in the PATH.
	@( $(LEX) < mm_options.l 1> /dev/null 2>&1 || \
	    { echo "Error: $(LEX) was not found !" ;\
	      echo "  Verify $(CONFIG_COMMAND) and examine the LEX definition" ;\
	      echo "  in the configuration file, $(CONFIG_FILE)." ;\
	      echo "  If a LEX cannot be found then try this:" ;\
	      echo "  cd ../src/dock/nab" ;\
	      echo "  make clean" ;\
	      echo "  make copy.lex.mm_options.c" ;\
	      echo "  make -e DOCKBUILDFLAGS='-Dflex' all" ;\
	      echo "  cd ../../../install" ;\
	      echo "  make install" ;\
	      exit 2 ;\
	    } ;\
	    /bin/rm -f lex.yy.c ;\
	)

lex.mm_options.c: mm_options.l
	make is_LEX_available
	$(LEX) -t mm_options.l | sed 's/yy/mmo/g' > lex.mm_options.c

# See header comments and target is_LEX_available.
copy.lex.mm_options.c:
	cp -p flex.lex.mm_options.c lex.mm_options.c
	touch lex.mm_options.c

molio.o:  molio.c
	$(CC) -c $(OCFLAGS) $(DOCKBUILDFLAGS) molio.c

sff.o:  sff.c
	$(CC) -c $(OCFLAGS) -DDOCK -D$(LEX) $(DOCKBUILDFLAGS) sff.c

conjgrad.o:  conjgrad.c
	$(CC) -c $(OCFLAGS) $(DOCKBUILDFLAGS) conjgrad.c

rand2.o:  rand2.c
	$(CC) -c $(OCFLAGS) $(DOCKBUILDFLAGS) rand2.c

chirvol.o:	\
	nab.h
conjgrad.o:  \
	memutil.h
jacobi.o:	\
	memutil.h
memutil.o:	\
	nab.h
molio.o:	\
	nab.h	\
	errormsg.h
molutil.o:	\
	nab.h	\
	errormsg.h	\
	memutil.h
prm.o: \
	nab.h
rand2.o: \
	defreal.h
select_atoms.o:	\
	nab.h
sff.o: \
	nab.h \
	memutil.h \
	lex.mm_options.c \
	rattle.c \
	eff.c \
	debug.h
traceback.o: \
	errormsg.h
