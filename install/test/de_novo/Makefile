include ../default_targets.h

test:   MkGrids dn_generate_fragments.dockout dn_de_novo_graph.dockout \
	dn_de_novo_ex.dockout dn_de_novo_rand.dockout dn_refinement.dockout

MkGrids: 
	$(DOCK_BIN)/multigrid_fp_gen.py dn_rec.mol2 dn_rec.resid \
	dn_generate_multigrid.gridin `cat dn_primary_residues.dat | sed -e 's/\n/ /g'`

clean: clean_dock clean_de_novo

