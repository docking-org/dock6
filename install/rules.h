# DOCK default macros and rules for source Makefiles.
# This should not contain targets so that default targets will work.
# This should be included in Makefiles before the configuration
# include file to ensure that the defaults can be overridden.
# 
# Note that DOCKBUILDFLAGS provides a hook into the build process for
# installers; for example, to build debug versions of the DOCK programs:
# make -e DOCKBUILDFLAGS="-DDEBUG -O0 -g"
# or, for example, to build profiling versions of the DOCK programs:
# make clean; make -e DOCKBUILDFLAGS="-pg"
# which for dock6 would then be analysed via
# dock6 -i bla; gprof -C dock6 |less


CONFIG_COMMAND=./configure
CONFIG_FILE=config.h
DOCK_VERSION=`cat ../src/dock/version.h | grep DOCK_VERSION | cut -d'"' -f2`
# This is used to put whitespace into a macro.
EMPTY=
SHELL=/bin/sh


# This target ensures that the right target (all) is the default target.
# Originally the idea was to avoid problems due to including this file
# in the wrong place.  But all known makes scan the whole makefile.
# However, AIX's make does not understand .PHONY; consequently aix make
# would treat all the targets in .PHONY as the default target !
default: all


# Append to the list of known suffixes
.SUFFIXES:  .cpp

# C++ compilation rule
.cpp.o:
	$(CXX) -c $(CXXFLAGS) -o $@ $(DOCKBUILDFLAGS) $<

# C compilation rule
.c.o:
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $(DOCKBUILDFLAGS) $<

# Fortran compilation rule
.f.o:
	$(FC) -c $(FFLAGS) -o $@ $(DOCKBUILDFLAGS) $<

# Fortran compilation rule with preprocessing
.F.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -o $@ $(DOCKBUILDFLAGS) $<


# This special target declares these targets as phony, ie, not based on
# files with that name; it avoids problems when such files do exist.
.PHONY:  all archive_export clean cleandock configured distclean dock \
         dockclean install rdkitclean realclean superclean test uninstall utils

