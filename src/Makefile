# DOCK Makefile for source directories.
# This is for building executables from source files.
# It merely invokes with the appropriate target, for example install
# or clean, in all the specified source subdirectories.

include ../install/rules.h

all: install

install: dock6 utils

dock6:
	$(MAKE) -C dock install

utils:
	$(MAKE) -C accessories install
	$(MAKE) -C gbsa_grids  install
	$(MAKE) -C grid  install
	$(MAKE) -C tleap  install
	$(MAKE) -C antechamber  install
	$(MAKE) -C mopac6/src  install
	$(MAKE) -C resp  install
	$(MAKE) -C docktools  install

clean:
	$(MAKE) -C dock  clean
	$(MAKE) -C accessories  clean
	$(MAKE) -C gbsa_grids  clean
	$(MAKE) -C grid  clean
	$(MAKE) -C tleap  clean
	$(MAKE) -C antechamber  clean
	$(MAKE) -C mopac6/src  clean
	$(MAKE) -C resp  clean
	$(MAKE) -C docktools  clean

realclean: 
	$(MAKE) -C dock  realclean
	$(MAKE) -C accessories  realclean
	$(MAKE) -C gbsa_grids  realclean
	$(MAKE) -C grid  realclean
	$(MAKE) -C tleap  clean
	$(MAKE) -C antechamber  realclean
	$(MAKE) -C mopac6/src  realclean
	$(MAKE) -C resp  realclean
	$(MAKE) -C docktools  realclean

uninstall:
	$(MAKE) -C dock  uninstall
	$(MAKE) -C accessories  uninstall
	$(MAKE) -C gbsa_grids  uninstall
	$(MAKE) -C grid  uninstall
	$(MAKE) -C tleap  uninstall
	$(MAKE) -C antechamber  uninstall
	$(MAKE) -C mopac6/src  uninstall
	$(MAKE) -C resp  uninstall
	$(MAKE) -C docktools  uninstall

