# Makefile for DOCK accessories subdirectory

include ../../install/rules.h
include ../../install/config.h

BIN = ../../bin

OBJS = \
	showbox.o \
	showsphere.o \
	sphgen.o

PROGS = \
	showbox \
	showsphere \
	sphgen \
	sphere_selector

SCRIPTS = \
	amberize_complex amberize_ligand amberize_receptor \
	prepare_amber.pl prepare_rna_amber.pl \
	mol2.py mol2.pyc multigrid_fp_gen.py mol2bild.py combine_torenv.py 

#############################################################################

all: $(PROGS)

install: all
	mv $(PROGS) $(BIN)
	sed "s?DOCKHOME?$(DOCKHOME)?g" prepare_amber.pl > $(BIN)/prepare_amber.pl
	sed "s?DOCKHOME?$(DOCKHOME)?g" prepare_rna_amber.pl > $(BIN)/prepare_rna_amber.pl
	sed "s?DOCKHOME?$(DOCKHOME)?g" amberize_ligand > $(BIN)/amberize_ligand
	sed "s?DOCKHOME?$(DOCKHOME)?g" amberize_complex > $(BIN)/amberize_complex
	sed "s?DOCKHOME?$(DOCKHOME)?g" amberize_receptor > $(BIN)/amberize_receptor
	cd $(BIN) && chmod +x amberize* prepare_amber.pl prepare_rna_amber.pl
	cp -p mol2.py $(BIN)/mol2.py
	cp -p combine_torenv.py $(BIN)/combine_torenv.py
	cp -p mol2bild.py $(BIN)/mol2bild.py
	sed "s?DOCKHOME?$(DOCKHOME)?g" multigrid_fp_gen.py > $(BIN)/multigrid_fp_gen.py
	cd $(BIN) && chmod +x multigrid_fp_gen.py

clean:
	/bin/rm -f $(OBJS)

realclean: clean
	/bin/rm -f $(PROGS) 

uninstall:
	cd $(BIN) && /bin/rm -f $(PROGS) $(SCRIPTS) __pycache__/mol2.cpython-36.pyc
	-cd $(BIN) && rmdir __pycache__

#############################################################################

#_____________________________________________________________________________
# Programs

showbox: showbox.o
	$(FC) $(FFLAGS) -o $@ $(DOCKBUILDFLAGS) $?

showsphere: showsphere.o
	$(FC) $(FFLAGS) -o $@ $(DOCKBUILDFLAGS) $?

sphgen: sphgen.o
	$(FC) $(FFLAGS) -o $@ $(DOCKBUILDFLAGS) $?

sphere_selector: sphere_selector.cpp sphere_selector.h
	$(CXX) $(CXXFLAGS) sphere_selector.cpp -lm -o $@ $(DOCKBUILDFLAGS)

