include ../default_targets.h

# prepare_amber.pl arguments:
REC = 1lgu.pdb
LIG = lig.mol2

test: work_in_progress_amberize move_nothing.dockout move_ligand.dockout \
      move_by_atom_expression.dockout move_by_distance.dockout \
      move_everything.dockout verbose_move_by_distance

verbose_move_by_distance:
	@# assume amberize has been made
	@echo
	@echo "Processing test $@"
	$(DOCK_BIN)/$(DOCK) -i move_by_distance.dockin -v 1> $@.dockout
	../dockdif -t 8 $@.dockout.save $@.dockout
	../dockdif -t 2 1lgu.lig.1.final_pose.amber.pdb.save 1lgu.lig.1.final_pose.amber.pdb
	../dockdif -t 2 1lgu.final_pose.amber.pdb.save 1lgu.final_pose.amber.pdb
	../dockdif -t 2 lig.1.final_pose.amber.pdb.save lig.1.final_pose.amber.pdb
	../dockdif -t 4 1lgu.lig.1.final_pose.amber.restart.save 1lgu.lig.1.final_pose.amber.restart
	../dockdif -t 2 1lgu.final_pose.amber.restart.save 1lgu.final_pose.amber.restart
	../dockdif -t 4 lig.1.final_pose.amber.restart.save lig.1.final_pose.amber.restart

work_in_progress_amberize: amberize
	@# receptor prep
	../dockdif 1lgu.amber.pdb.save 1lgu.amber.pdb
	@# ligand prep
	../dockdif -t 3 lig.1.gaff.mol2.save lig.1.gaff.mol2
	@# complex prep
	@# prepare_amber
	../dockdif lig.amber_score.mol2.save lig.amber_score.mol2

clean: clean_amberize clean_dock
	/bin/rm -f 1lgu.lig.1.final_pose.amber.pdb \
		1lgu.final_pose.amber.pdb \
		lig.1.final_pose.amber.pdb \
		1lgu.lig.1.final_pose.amber.restart \
		1lgu.final_pose.amber.restart \
		lig.1.final_pose.amber.restart

