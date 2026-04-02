include ../default_targets.h

test: fps_grids 

clean: clean_fps_grids clean_dock

fps_grids: important_residues.dockout
	# Compute grids for important residues and the one grid for the remainder residues.
	# getrange_mkgrides
	$(DOCK_BIN)/multigrid_fp_gen.py rec.mol2 resid grid.in \
	    `grep -A 1 range_union important_residues.dockout | sed -e'/range_union/d' -e'/--/d' -e's/ //g' | tr ',' '\n' | sort -nu`
	../dockdif resid_005.out.save           resid_005.out
	../dockdif resid_013.out.save           resid_013.out
	../dockdif resid_020.out.save           resid_020.out
	../dockdif resid_022.out.save           resid_022.out
	../dockdif resid_048.out.save           resid_048.out
	../dockdif resid_049.out.save           resid_049.out
	../dockdif resid_050.out.save           resid_050.out
	../dockdif resid_051.out.save           resid_051.out
	../dockdif resid_052.out.save           resid_052.out
	../dockdif resid_053.out.save           resid_053.out
	../dockdif resid_054.out.save           resid_054.out
	../dockdif resid_055.out.save           resid_055.out
	../dockdif resid_071.out.save           resid_071.out
	../dockdif resid_078.out.save           resid_078.out
	../dockdif resid_remaining.out.save     resid_remaining.out

clean_fps_grids:
	/bin/rm -f *.bmp *.cnt *.nrg *.out *.out.dif resid*.mol2 temp.mol2
