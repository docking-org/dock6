# DOCK docktools Makefile.
# This is for building the docktools executables from source files.
# Converts chemgrids as well as receptor desolvation grids. ifdefs in 
# the code to choose one over the other and hence two makefiles (a bit
# sloppy but works for now)

all: install

install:
	make -f gridconv.mk install
	make -f gridrdsl.mk install

clean:
	make -f gridconv.mk clean
	make -f gridrdsl.mk clean

realclean:
	make -f gridconv.mk realclean
	make -f gridrdsl.mk realclean

uninstall:
	make -f gridconv.mk uninstall
	make -f gridrdsl.mk uninstall
