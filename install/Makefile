# DOCK top-level installation Makefile.
# This is mostly for building and testing executables from source files.
# It merely invokes make with the appropriate target, for example install
# or clean, in all the specified directories.

include ./rules.h


install: configured
	@echo "Starting installation of "
	@echo ${DOCK_VERSION}
	@echo "at `date`."
	@echo ""
	$(MAKE) -C ../src  install
	@echo ""
	@echo "Installation of "
	@echo ${DOCK_VERSION}
	@echo "is complete at `date`."
	@echo ""

all: install

dock: configured
	@echo "Starting installation of "
	@echo ${DOCK_VERSION}
	@echo "program dock at `date`."
	@echo ""
	$(MAKE) -C ../src/dock  install
	@echo ""
	@echo "Installation of "
	@echo ${DOCK_VERSION}
	@echo "program dock is complete at `date`."
	@echo ""

parallel: configuredparallel
	make cleandock
	make dock

rdkit:
	$(MAKE) -C ../src/rdkit install

utils: configured
	@echo "Starting installation of "
	@echo ${DOCK_VERSION}
	@echo "utilities at `date`."
	@echo ""
	$(MAKE) -C ../src utils
	@echo ""
	@echo "Installation of "
	@echo ${DOCK_VERSION}
	@echo "utilities is complete at `date`."
	@echo ""

test: 
	$(MAKE) -C  test  test
	$(MAKE) -C  test check

testparallel: 
	$(MAKE) -C test  silentclean
	$(MAKE) -C test  parallel
	$(MAKE) -C test  check

check: 
	$(MAKE) -C test  check

# This target removes all object files.
clean: 
	$(MAKE) -C ../src  clean

# This target removes dock object files.
cleandock: 
	$(MAKE) -C ../src/dock  clean

dockclean: cleandock

rdkitclean:
	$(MAKE) -C ../src/rdkit rdkitclean

# This target removes installed executables from bin.
uninstall:
	$(MAKE) -C ../src  uninstall

# This target cleans and removes executables from the source directories.
realclean:
	$(MAKE) -C ../src  realclean

# This target realcleans, removes the configuration file and cleans the tests.
distclean: realclean
	$(MAKE) -C ./test  clean
	/bin/rm -f $(CONFIG_FILE)

# This target cleans and removes all executables and the configuration file.
superclean: uninstall distclean
	$(MAKE) -C ../src/rdkit rdkitclean

configured:
	@# abort if the configuration file does not exist
	@( if [ ! -f $(CONFIG_FILE) ] ; then \
	    echo "Error: $(CONFIG_COMMAND) must be executed before make !" ;\
	    exit 2 ;\
	fi ;\
	)

configuredparallel: configured
	@# abort if the configuration file does not support MPI
	@( if ! grep MPI $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: the DOCK configuration file does not support MPI !" ;\
	    echo "       Rerun $(CONFIG_COMMAND) specifying parallel !" ;\
	    exit 2 ;\
	fi ;\
	)

help:
	@echo "configure is the first step in the installation.  Execute this: "
	@echo "./configure -help"


#
# The targets below are for DOCK developers.
#

# This target exports and archives the DOCK directory tree.
# It is a template for a bourne script taking the tag name and destination
# as arguments.
DISTRIBUTION_TAG=dock6-1-release
DISTRIBUTION_TAG=dock6-2-release
DISTRIBUTION_TAG=dock6-3-release
DISTRIBUTION_TAG=dock6-4-release
DISTRIBUTION_TAG=dock6-5-release
DISTRIBUTION_TAG=dock6-6-release
DISTRIBUTION_TAG=dock6-7-release
DISTRIBUTION_TAG=dock6-8-release
DISTRIBUTION_TAG=dock6-9-beta
DISTRIBUTION_TAG=dock6-9-rc1
DISTRIBUTION_TAG=dock6-9-release
DISTRIBUTION_TAG=dock6-10-beta
DISTRIBUTION_TAG=dock6-10-release
DISTRIBUTION_TAG=dock6-11-beta
DISTRIBUTION_TAG=dock6-11-release
#DOCK_HOME=..
DOCK_MODULE=dock6
DOCK_ARCHIVE=$(DOCK_MODULE).$(DISTRIBUTION_TAG)
archive_export:
	# git archive isn't flexible enough, so clone and remove
	cd /tmp && git clone git@github.com:docking-org/dock6.git --depth 1 $(DOCK_MODULE)
	# plan to use lightweight tags for prereleases, eg; git tag dock6-10-rc1
	# and branches after releases, eg: git checkout -b dock6-10-release
	cd /tmp/$(DOCK_MODULE) && git tag $(DISTRIBUTION_TAG)
	cd /tmp/$(DOCK_MODULE) && rm -f .gitmodules .gitmodules  # for example
	cd /tmp/$(DOCK_MODULE) && mkdir bin  # because git fumbled empty directories
	cd /tmp/$(DOCK_MODULE) && mkdir tutorials/amber_scoring_demo/t4lys  # ditto
	cd /tmp/$(DOCK_MODULE) && mkdir tutorials/ligand_sampling_demo/1_struct # ditto
	cd /tmp/$(DOCK_MODULE) && rm -f README.md   # non-public
	cd /tmp/$(DOCK_MODULE) && rm -rf developer  # non-public
	cd /tmp/$(DOCK_MODULE) && rm -rf tutorials/bks_fxn_demo  # broken
	cd /tmp/$(DOCK_MODULE) && rm -rf install/intel.ip.mkl  # redundant?
	cd /tmp/$(DOCK_MODULE)/install && make populatedocktutorials
	cd /tmp/$(DOCK_MODULE)/install && make webarchiveexport
	cd /tmp/$(DOCK_MODULE) && rm -rf docs  # docs is only released on the web
	cd /tmp && tar cf $(DOCK_ARCHIVE).tar $(DOCK_MODULE)
	cd /tmp && tar --delete "$(DOCK_MODULE)/\.git/" -f $(DOCK_ARCHIVE).tar
	cd /tmp && gzip -9 $(DOCK_ARCHIVE).tar
	mv /tmp/$(DOCK_ARCHIVE).tar.gz /tmp/$(DOCK_ARCHIVE).tgz
	#rm -rf /tmp/$(DOCK_MODULE)
	-chmod 700 /tmp/$(DOCK_MODULE)
	-cd /tmp && scp $(DOCK_ARCHIVE).tgz $(DESTINATION)
	rm -i /tmp/$(DOCK_ARCHIVE).tgz
	-chmod 600 /tmp/$(DOCK_ARCHIVE).tgz
	@echo "Exporting of $(DOCK_ARCHIVE) is completed."

# intended to be made after target updatedocktutorials
webarchiveexport: generating_grid.tar.gz lig_sampling.tar.gz \
		site_files.tar.gz struct_files.tar.gz
	cd ../docs/ && gtar zvcf /tmp/webarchiveexport.$(DISTRIBUTION_TAG).tgz .
	chmod 600 /tmp/webarchiveexport.$(DISTRIBUTION_TAG).tgz
	#-cd ../docs/ && scp webarchiveexport.$(DISTRIBUTION_TAG).tgz /tmp
	#cd /var/www/html/kuntz/DOCK_6 && extract

# update dock input and output files in docs/tutorials/ by
# running script_demo in tutorials/
updatedocktutorials: populatedocktutorials
	#If things don't look correct after making this once then make it again.
	#only amber does not need populatedocktutorials for updating;  12-2012 ?
	#this target explicitly updates only dock program inputs and outputs.
	#for amber only readmes and script*
	#exist in cvs in /tutorials; everything else is in cvs in /docs/tutorials.
	#run dock first to update dock input file.
	#error ignore is for repeated making.
	-cd ../docs/tutorials/amber_score && \
	  ../../../bin/dock6 -i dock.in -o dock.out
	#run tutorial script to update the rest.
	cd ../docs/tutorials/amber_score && script_clean && script_demo
	#for ligand_sampling_demo only readmes and script*
	#exist in cvs in /tutorials; everything else is in cvs in /docs/tutorials.
	cd ../tutorials/ligand_sampling_demo/4_dock && \
	  ../../../bin/dock6 -i rigid.in -o dock.out
	cp -p ../tutorials/ligand_sampling_demo/4_dock/dock.out \
	  ../docs/tutorials/ligand_sampling_dock/rigid.out
	cd ../tutorials/ligand_sampling_demo/4_dock && \
	  ../../../bin/dock6 -i anchor_and_grow.in -o dock.out
	cp -p ../tutorials/ligand_sampling_demo/4_dock/dock.out \
	  ../docs/tutorials/ligand_sampling_dock/anchor_and_grow.out
	cd ../tutorials/ligand_sampling_demo && \
	  script_clean && script_demo 1> demo.out 2>&1
	#demo.out contains two dock output files.
	cd ../tutorials/ligand_sampling_demo/4_dock && \
	  cp -p anchor_and_grow.in rigid.in *.mol2 \
	  ../../../docs/tutorials/ligand_sampling_dock 
	#for this tutorial everything in cvs exists in /tutorials.
	cd ../tutorials/solvent_scoring_demo && \
	  script_clean && script_demo 1> demo.out 2>&1
	#for this tutorial everything in cvs exists in /tutorials.
	#for mpi, only the input file is distributed, but run everything to verify.
	#error ignore is for repeated making; in particular, this can be ignored:
	##Could not open ../3_grid/1VRT.rec.nrg for reading. Program will terminate.
	-cd ../tutorials/mpi_demo/4_dock && ../../../bin/dock6 -i mpi.in
pdatedocktutorials: 
	cd ../tutorials/mpi_demo && \
	  script_clean && script_demo 1> demo.out 2>&1

generating_grid.tar.gz:
	cd ../docs/tutorials/grid_generation && \
	gtar zvcf generating_grid.tar.gz box.in grid.bmp grid.in grid.nrg \
	  grid.out rec_box.pdb

lig_sampling.tar.gz:
	cd ../docs/tutorials/ligand_sampling_dock && \
	gtar zvcf lig_sampling.tar.gz anchor_and_grow.in anchor_and_grow.out \
	  anchor_and_grow_scored.mol2 rigid.in rigid.out rigid_scored.mol2

site_files.tar.gz:
	cd ../docs/tutorials/sphere_generation && \
	gtar zvcf site_files.tar.gz INSPH OUTSPH rec.ms rec.sph \
	  selected_spheres.in selected_spheres.pdb selected_spheres.sph \
	  sphgen_cluster.in sphgen_cluster.pdb

struct_files.tar.gz:
	cd ../docs/tutorials/struct_prep && \
	gtar zvcf struct_files.tar.gz 1ABE.pdb rec_noH.pdb rec_charged.mol2 \
	  lig_charged.mol2

# Copy input files from docs/tutorials to tutorials.
# In the future consider whether to eliminate directory tutorials from cvs
# and distribute a complete mirror of the docs/tutorials
# which is also the web site tutorials.
# The effort needed is mainly to untangle directory names and scripts,
# but the consideration should include simple addition of new tutorials.
populatedocktutorials:
	[ -d ../tutorials/amber_scoring_demo/t4lys ] || \
	  mkdir ../tutorials/amber_scoring_demo/t4lys
	cd ../docs/tutorials/amber_score && \
	  cp -p 1lgu.pdb dock.in lig.mol2 script_clean script_demo amber_score.htm \
	  README Introductory_Tutorial_for_Amber_score_DOCK6.* \
	  ../../../tutorials/amber_scoring_demo/t4lys
	[ -d ../tutorials/ligand_sampling_demo/1_struct ] || \
	  mkdir ../tutorials/ligand_sampling_demo/1_struct
	cd ../docs/tutorials/struct_prep && \
	  cp -p 1ABE.pdb lig_charged.mol2 rec_charged.mol2 rec_noH.pdb 0README \
	  ../../../tutorials/ligand_sampling_demo/1_struct
	cd ../docs/tutorials/sphere_generation && \
	  cp -p INSPH rec.ms selected_spheres.in sphgen_cluster.in \
	  ../../../tutorials/ligand_sampling_demo/2_site
	cd ../docs/tutorials/grid_generation && \
	  cp -p box.in grid.in rec_box.pdb \
	  ../../../tutorials/ligand_sampling_demo/3_grid
	cd ../docs/tutorials/ligand_sampling_dock && \
	  cp -p anchor_and_grow.in rigid.in \
	  ../../../tutorials/ligand_sampling_demo/4_dock

