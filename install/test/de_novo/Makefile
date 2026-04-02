include ../default_targets.h

test:   MkGrids dn_generate_fragments.dockout dn_de_novo_graph.dockout \
		dn_generate_fragments_biasing.dockout \
	dn_de_novo_ex.dockout dn_de_novo_rand.dockout dn_refinement.dockout \
	dn_de_novo_graph_no_pruning.dockout dn_de_novo_graph_pruning_no_sampling.dockout \
	dn_de_novo_graph_biasing.dockout dn_de_novo_rand_biasing.dockout

rdkittest:
	@# run this test only if the dock rdkit executable exists
	@(if [ -x "$(DOCK_BIN)/$(DOCK).rdkit" ] ; then \
	    make rdkittest_no_dock6_rdkit_is_an_error ;\
	else \
	    echo "" ;\
	    echo "The dock rdkit executable does not exist." ;\
	    echo "    Skipping this test." ;\
	    echo "" ;\
	fi ;\
	)

rdkittest_no_dock6_rdkit_is_an_error:   MkGrids \
	dn_generate_fragments.dockrdkitout \
	dn_de_novo_graph.dockrdkitout \
	dn_de_novo_ex.dockrdkitout \
	dn_de_novo_rand.dockrdkitout

MkGrids:
	$(DOCK_BIN)/multigrid_fp_gen.py dn_rec.mol2 dn_rec.resid \
	dn_generate_multigrid.gridin `cat dn_primary_residues.dat | sed -e 's/\n/ /g'`

clean: clean_dock clean_de_novo

