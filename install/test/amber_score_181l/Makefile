include ../default_targets.h

REC = 181l.pdb
LIG = lig.mol2

test: amberize move_ligand.dockout 

clean: clean_amberize clean_dock
	/bin/rm -f 181l.lig.1.final_pose.amber.restart \
		lig.1.final_pose.amber.restart

