# DOCK_GA regression case matrix

A suggested golden-master suite whose cases, together, exercise every surface the
`conf_gen_ga` refactor touched. Use it to build the first golden references (see
[`hpc_golden_verification.md`](hpc_golden_verification.md) for how to capture/diff, and
the [regression harness README](../src/dock/tests/regression/README.md) for the case
directory layout).

Each case is a copy of a shared **base** `dock.in` with a few keys overridden. Keep them
small and fast (few generations, small ensemble) so the whole suite reruns after every
refactor step.

---

## Coverage matrix

| Case dir | Refactored surface exercised | Key overrides vs base |
|---|---|---|
| `sel_elitism` | selection registry → **elitism** handler (`ga_selection`) | *(base; `ga_selection_method elitism`)* |
| `sel_tournament` | selection registry → **tournament** handler | `ga_selection_method tournament` |
| `sel_roulette` | selection registry → **roulette** handler | `ga_selection_method roulette` |
| `mut_deletion` | mutation type-pool → **deletion** arm (`ga_mutation`) | only `ga_mutate_deletion yes` |
| `mut_addition` | mutation type-pool → **addition** arm | only `ga_mutate_addition yes` |
| `mut_substitution` | mutation type-pool → **substitution** (delete+add) arm | only `ga_mutate_substitution yes` |
| `mut_replacement` | mutation type-pool → **replacement** arm + segment-type pool | only `ga_mutate_replacement yes` |
| `xover_random` | `breeding_rand` path | *(base; `ga_xover_sampling_method_rand yes`)* |
| `xover_exhaustive` | `breeding_exhaustive` + `check_exhaustive` path | `ga_xover_sampling_method_rand no`, `ga_check_overlap yes` |
| `filter_mw_hard` | `ga_filters` hard MW + descriptors | `ga_mol_wt_cutoff_type hard`, tight MW bounds |
| `filter_mw_soft` | `ga_filters` **soft** MW (the RNG-injected path) | `ga_mol_wt_cutoff_type soft`, tight bounds, set std-dev |
| `filter_tight` | `ga_filters` rot-bond / HA / HD / charge predicates all firing | tight `ga_constraint_*` limits |

`ga_descriptors` (MW/rot/HA/HD/charge + covalent radius) and `ga_naming` (titles) are
exercised by **every** case, so no dedicated case is needed for them — but do eyeball a
few molecule titles in the golden output to confirm the `_gNNNN_iNNN` naming (incl. the
gen-vs-loc padding asymmetry) is what you expect.

Twelve cases is thorough; if you want a minimal first pass, the four that cover the most
distinct code are: `sel_tournament`, `sel_roulette`, `mut_replacement`, `filter_mw_soft`
(elitism + random-xover + hard-MW + the other mutations are all in the base already).

---

## Base `dock.in` — the GA block

This is only the **GA-specific** block. A working `dock.in` also needs the usual DOCK
sections your production runs already use: the receptor grid / footprint / scoring setup
(`grid_score_*`, `*_score_secondary`, footprint/pharmacophore files as applicable),
`orient_ligand` + sphere/box files, the AMBER/vdw parameter files, and the `simplex_*`
minimization block. Reuse those from an existing run. The two non-negotiables for a
reproducible case are in **bold** below.

```text
##### molecule / de-novo GA #####
conf_search_type            de_novo_ga        # (or however your build enters DOCK_GA)
ga_molecule_file            inputs/start.mol2 # small starting ligand (multi-segment — see caveats)
ga_utilities                no
ga_charge_parent_gasteiger  no

# fragment libraries (needed for addition / substitution / replacement)
ga_fraglib_scaffold_file    inputs/scaffold.mol2
ga_fraglib_linker_file      inputs/linker.mol2
ga_fraglib_sidechain_file   inputs/sidechain.mol2

# small + fast so the suite reruns quickly
ga_max_generations          3
ga_ensemble_size            20
ga_max_num_gen_with_no_crossover  25

##### crossover #####
ga_xover_on                 yes
ga_xover_sampling_method_rand  yes            # -> breeding_rand ; "no" + check_overlap -> exhaustive
ga_xover_max                50
ga_bond_tolerance           0.5
ga_angle_cutoff             0.14
ga_check_overlap            no

##### mutation (rates over ENABLED types must sum to 100) #####
ga_mutations                yes
ga_mutate_addition          yes
ga_mutate_deletion          yes
ga_mutate_substitution      yes
ga_mutate_replacement       yes
ga_add_rate                 25
ga_del_rate                 25
ga_sub_rate                 25
ga_rep_rate                 25
ga_mutate_parents           no
ga_omut_rate                0.7
ga_max_mut_cycles           5
ga_mut_sampling_method      rand
ga_num_random_picks         15
ga_max_root_size            5

##### pruning / filters #####
ga_energy_cutoff            100
ga_heur_unmatched_num       1
ga_heur_matched_rmsd        2.0
ga_mol_wt_cutoff_type       hard
ga_constraint_upper_mol_wt  500.0
ga_constraint_lower_mol_wt  0
ga_constraint_mol_wt_std_dev  35.0
ga_constraint_rot_bon       10
ga_constraint_H_accept      10
ga_constraint_H_don         5
ga_constraint_formal_charge 2

##### selection #####
ga_selection_method         elitism           # elitism | tournament | roulette
ga_elitism_combined         yes
ga_elitism_option           max               # percent | number | max
ga_tournament_p_vs_c        yes
ga_roulette_separate        yes

##### output / naming #####
ga_name_identifier          ga
ga_output_prefix            ga_output

##### REQUIRED for reproducibility #####
minimize_ligand             yes               # ** RNG is only seeded when the minimizer is on **
simplex_random_seed         12345             # ** fixed seed **
```

---

## Per-case overrides

Copy the base into each `cases/<name>/dock.in` and apply only these lines.

**Selection** — change one key:
```text
# sel_tournament
ga_selection_method   tournament
# sel_roulette
ga_selection_method   roulette
# sel_elitism = base (no change)
```

**Mutation** — enable exactly one type so the type-pool is that type only. When only one
type is enabled its rate defaults to 100, but set it explicitly to be safe:
```text
# mut_deletion
ga_mutate_addition no
ga_mutate_deletion yes
ga_mutate_substitution no
ga_mutate_replacement no
ga_del_rate 100

# mut_addition   -> ga_mutate_addition yes, others no, ga_add_rate 100
# mut_substitution -> ga_mutate_substitution yes, others no, ga_sub_rate 100
# mut_replacement  -> ga_mutate_replacement yes, others no, ga_rep_rate 100
```

**Crossover**:
```text
# xover_exhaustive
ga_xover_sampling_method_rand no
ga_check_overlap              yes
# xover_random = base (no change)
```

**Filters** — make the constraints bite so the predicates actually reject molecules
(pick numbers relative to your ligands' real MW/rotatable-bond counts):
```text
# filter_mw_hard  (tighten the window so some children fall outside it)
ga_mol_wt_cutoff_type       hard
ga_constraint_upper_mol_wt  350.0
ga_constraint_lower_mol_wt  150.0

# filter_mw_soft  (same window, soft probabilistic accept -> exercises the RNG-injected path)
ga_mol_wt_cutoff_type       soft
ga_constraint_upper_mol_wt  350.0
ga_constraint_lower_mol_wt  150.0
ga_constraint_mol_wt_std_dev  20.0

# filter_tight  (force rot/HA/HD/charge rejects)
ga_constraint_rot_bon       3
ga_constraint_H_accept      3
ga_constraint_H_don         1
ga_constraint_formal_charge 0
```

---

## `manifest.txt` per case

List the scientific output files to compare (basenames the run writes into its working
dir). Exact names depend on `ga_output_prefix` and how many generations run; a typical
set:
```text
# restart / ensemble snapshots per generation
restart0000.mol2
restart0001.mol2
restart0002.mol2
restart0003.mol2
# whatever else your run emits (pruned / filtered / rejected / scores)
# ga_output... (add the real basenames after one baseline run)
```
Easiest path: run the baseline once, `ls` the output dir, and list the files you care
about. Or omit `manifest.txt` entirely to compare **every** produced file (then use
`--ignore-regex 'seconds'` to skip wall-clock log lines).

---

## Caveats that make or break a case

- **Ligand must have enough segments.** Deletion requires `num_segments > 2`; substitution
  and replacement need multiple attachment-point segments (sidechain/linker/scaffold). A
  tiny rigid ligand will skip those operators and the case won't actually exercise them —
  pick a starting molecule with several rotatable-bond-separated pieces.
- **Fragment libraries must be present and non-trivial** for addition/substitution/
  replacement to produce children; `ga_fraglib_*` files must contain compatible fragments.
- **Reproducibility gate.** Every case needs `minimize_ligand yes` + a fixed
  `simplex_random_seed`. Validate each case by running the **baseline twice** and diffing
  (`compare.py`) before capturing its golden — if two baseline runs differ, the case is not
  seed-deterministic and any "regression" it reports later is noise.
- **Rates sum to 100** over the *enabled* mutation types, else the run `exit(0)`s at parse.
- **soft-MW draws RNG**, so `filter_mw_soft` only reproduces if the seed is fixed (it is,
  via the gate above) — this case is the end-to-end check that the injected-`rand()`
  refactor of `mw_cutoff` preserved the draw sequence.
- Keep generations/ensemble small but **> 0 real work**: at least a couple of generations
  so selection and mutation actually run and the golden captures their effect.
```
