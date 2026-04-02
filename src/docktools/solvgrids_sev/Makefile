# Makefile for DOCK chemgrid directory

include ../../../install/rules.h
include ../../../install/config.h

BIN = ../../../bin
#PROGS =  solvgrid solvmap
PROGS = sevsolv

OBJS1 = sas.o cube.o memalloc.o solvmap.o

#CFLAGS = -O3 -DLINUX 
#CFLAGS = -m32 -g -DLINUX 
#CFLAGS = -g -DLINUX 
CFLAGS = -O3 -g -DLINUX 
#CFLAGS = -m32 -DLINUX 
#CFLAGS =  -DLINUX 
#CFLAGS =  -O3 -DLINUX 
#FFLAGS2 = "$(FFLAGS) -fcray-pointer"
#FFLAGS2 = -m32 -DLINUX -fcray-pointer -fdefault-integer-8 -fdefault-real-8
#FFLAGS2 = -DLINUX -fcray-pointer 
FFLAGS2 = -O3 -DLINUX -fcray-pointer -fdefault-integer-8 -fdefault-real-8
#FFLAGS2 = -DLINUX -fcray-pointer -fdefault-integer-8 -fdefault-real-8
#FFLAGS2 = -O3 -DLINUX -fcray-pointer -fdefault-integer-8 -fdefault-real-8
#FFLAGS3 = -m32 -Bstatic
#FFLAGS2 =  -DLINUX -fcray-pointer
#FFLAGS3 =  -Bstatic
FFLAGS3 = 
FFLAGS=  $(FFLAGS2)

#-r--r----- 1 baliuste baliuste   357 Mar 29 14:54 acc.h
#-r--r----- 1 baliuste baliuste 12592 Mar 29 14:54 cube.f
#-r--r----- 1 baliuste baliuste   370 Mar 29 14:54 cubic.h
#-r--r----- 1 baliuste baliuste   500 Mar 29 14:54 dist.h
#-r--r----- 1 baliuste baliuste  1290 Mar 29 14:54 memalloc.c
#-r--r----- 1 baliuste baliuste   221 Mar 29 14:54 pointer.h

#-r--r----- 1 baliuste baliuste 31733 Mar 29 14:54 sas.f
#-r--r----- 1 baliuste baliuste 21961 Mar 29 14:54 solvmap.f

# cube.f:   include 'cubic.h'
# memalloc.c:#include <stdlib.h>
# memalloc.c:#include <string.h>
# sas.f:         include 'acc.h'
# sas.f:         include 'cubic.h'
# sas.f:          include 'pointer.h'
# sas.f:         include 'pointer.h'
# solvmap.f:       include 'acc.h'
# solvmap.f:       include 'dist.h'
# solvmap.f:        include 'pointer.h'
# solvmap.f:       include 'pointer.h'
#

# gfortran -O3 -DLINUX -fcray-pointer -c  sas.f
# gfortran -O3 -DLINUX -fcray-pointer -c  cube.f
# /usr/bin/gcc -O3 -DLINUX -c  memalloc.c
# gfortran -O3 -DLINUX -fcray-pointer -c  solvmap.f
# gfortran -O3 -DLINUX -fcray-pointer -Bstatic -o sevsolv sas.o cube.o memalloc.o solvmap.o



all: $(PROGS)

install: all
	mv $(PROGS) $(BIN)

clean:
	/bin/rm -f $(OBJS1) $(OBJS2)

realclean: clean
	/bin/rm -f $(PROGS)

uninstall:
	cd $(BIN) && /bin/rm -f $(PROGS)

sevsolv: $(OBJS1)
	$(FC) $(FFLAGS) $(FFLAGS3) $(OBJS1) -o $@


memalloc.o: 
sas.o: pointer.h 
cube.o: cubic.h pointer.h memalloc.o 
solvmap.o: acc.h dist.h pointer.h  



