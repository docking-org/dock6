include ../../install/rules.h
include ../../install/config.h

PROGRAM = resp 
BINDIR = ../../bin

OBJS= amopen.o lapack.o mexit.o resp.o  
SRCS= amopen.F lapack.F mexit.F resp.F

all:		$(PROGRAM)

$(PROGRAM):	$(OBJS)
		$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

install:	$(PROGRAM)
		mv $(PROGRAM) $(BINDIR)

clean:
	/bin/rm -f $(OBJS)

realclean: clean
	/bin/rm -f $(PROGRAM)

uninstall:
	cd $(BINDIR) && /bin/rm -f $(PROGRAM)

