include ../default_targets.h

test: 
	make -e DOCK_VERBOSE='-v' blurry_gist_trilinear.dockout  blurry_gist_trilinear_flex.dockout  descriptor_score_bgist_flex.dockout  descriptor_score_grid_bgist_flex.dockout  full_gist_displacement.dockout full_gist_displacement_united.dockout full_gist_displacement_blurry.dockout 

clean: clean_dock
	/bin/rm -f LIG1_anchor1_branch*.mol2
