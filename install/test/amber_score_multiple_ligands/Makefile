include ../default_targets.h

REC = 1lgu.pdb
LIG = t4lys_mutant_multiple.mol2


test: amberize  move_ligand.dockout 

benz:
	$(DOCK_BIN)/prepare_amber.pl benz_binders.mol2 $(REC)
	make benz_move_nothing.dockout 

benzclean:
	/bin/rm -f  benz_binders.*.gaff.mol2 benz_binders.*.mol2 \
	            benz_binders.*.mopac.out benz_binders.*.log \
	            benz_binders.*.frcmod benz_binders.*.amber.pdb \
	            benz_binders.*.inpcrd benz_binders.*.prmtop \
	            amberize_rec.out amberize_ligand.*.out amberize_complex.*.out \
	            1lgu.benz_binders.*.log 1lgu.benz_binders.*.amber.pdb \
	            1lgu.benz_binders.*.inpcrd 1lgu.benz_binders.*.prmtop 

clean: clean_amberize clean_dock
