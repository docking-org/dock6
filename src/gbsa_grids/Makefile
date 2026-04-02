# DOCK gbsa_grids Makefile.
# This is for building the gbsa_grids executables from source files.
# It merely invokes make with the appropriate target, for example install
# or clean, in all the specified source subdirectories.

all: install

install:
	$(MAKE) -C schemgrid_GB  install
	$(MAKE) -C schemgrid_SA  install

clean:
	$(MAKE) -C schemgrid_GB  clean
	$(MAKE) -C schemgrid_SA  clean

realclean:
	$(MAKE) -C schemgrid_GB  realclean
	$(MAKE) -C schemgrid_SA  realclean

uninstall:
	$(MAKE) -C schemgrid_GB  uninstall
	$(MAKE) -C schemgrid_SA  uninstall

