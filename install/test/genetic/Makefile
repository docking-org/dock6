include ../default_targets.h

test:   combine_torsions tournament.GA.dockout roulette.GA.dockout \
	elitism.GA_max.dockout  

combine_torsions:
	#Processing dock6 test combine_torsions
	 $(DOCK_BIN)/combine_torenv.py ./fraglib_torenv_large.dat ./fraglib_torenv_small.dat
	../dockdif unique_full_sorted_fraglib.dat.save unique_full_sorted_fraglib.dat

clean: clean_dock
	/bin/rm -f *restart*mol2 *crossover*mol2 full_fraglib.dat \
		full_sorted_fraglib.dat unique_full_sorted_fraglib.dat
