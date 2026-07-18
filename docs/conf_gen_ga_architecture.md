# DOCK_GA Architecture Map — `src/dock/conf_gen_ga.{cpp,h}`

**Phase 0 deliverable** for the `conf_gen_ga` refactor (see `conf_gen_ga_refactor_spec.md`).
This is a **map, not a refactor** — no source was modified to produce it.

- Target file: [`src/dock/conf_gen_ga.cpp`](../src/dock/conf_gen_ga.cpp) — **14,750 lines**, 601 KB
- Header: [`src/dock/conf_gen_ga.h`](../src/dock/conf_gen_ga.h) — **593 lines**
- Version characterized: **v6.13.1** (tag `v6.13.1`)
- Everything lives in **one class, `GA_Recomb`** (~90 methods) plus two small helper classes (`Tor_Env`, `SEGMENTS`).

> **How to read this doc.** §1 is the symbol inventory. §2 clusters those symbols into the concerns the refactor will split along. §3 is the state/coupling map, including the full RNG-consumption enumeration (the determinism linchpin). §4 is build wiring. §5 is the extension-point audit that drives Phase 3. §6 is the risk register. §7 records how to build/run (tests come in Phase 1). §8 lists decisions the user must confirm before Phase 2.

---

## 0. Executive summary

`conf_gen_ga.cpp` implements **DOCK_GA**, the molecular-evolution genetic algorithm: it evolves candidate ligands in the binding site via crossover ("breeding"), mutation, fitness pruning, and selection, reusing DOCK6's typing, minimization, and scoring machinery. The code is a **single ~15k-line god-class** with a de-facto "blackboard" architecture: ~90 methods communicate through ~40 shared mutable member vectors (`parents`, `children`, `scored_generation`, `tmp_*`, `mutants`, …) and ~120 `ga_*` configuration members rather than through parameters.

Findings that shape the whole refactor:

1. **Determinism rests on one seed and ~38 ordered `rand()` sites.** RNG is C stdlib `rand()`, seeded **once** at [`conf_gen_ga.cpp:967`](../src/dock/conf_gen_ga.cpp) (`srand(simplex.random_seed)`). Every later draw consumes from that single global stream, so *any* change to the count or order of draws anywhere shifts all downstream results. This is the #1 constraint on every Phase 2 step.
2. **No MPI in this file.** Contrary to the spec's general worry, `conf_gen_ga.{cpp,h}` contain **zero** MPI code (the two `mpi` grep hits are the word "compilation" in comments). MPI lives in `dock.cpp`/`base_mpi.*`. This file is purely serial; MPI safety here means "don't break the MPI *build*," not "preserve rank logic."
3. **The extension points the spec targets are all "edit-a-scattered-set-of-sites" today** — selection is an if/else ladder over 5 booleans spread across 7 sites in 2 files; mutation is a 645-line dispatcher; filters are a 5-check block copy-pasted (and already drifted) across 3 functions; crossover is a boolean ladder feeding triplicated geometry code. Every Goal-#1 seam is genuinely painful today (details in §5).
4. **These files carry no license header.** Siblings (e.g. `dock.cpp`) carry a UCSF copyright-grant block; the repo `LICENSE` is BSD-3-Clause. New sibling files should match the *sibling* header convention — see decision **G** in §8.
5. **Significant dead / disabled code** coexists with live code behind `#if 0` and forced-false flags (RDKit descriptor-drive, niching, extinction, SUS/Metropolis selection, `get_stddev_changes`). Do **not** delete without confirming across build configs (§6).

---

## 1. Symbol inventory

All symbols are methods of `class GA_Recomb` unless noted. Line numbers are definition start lines in `conf_gen_ga.cpp`. "Size" is approximate. **Bold = god-function / high-complexity hotspot.**

### 1.1 Largest / most complex functions (refactor hotspots)

| Function | Lines | ~Size | Why it's a hotspot |
|---|---:|---:|---|
| **`max_breeding`** | 947 | ~836 | The entire generational control loop; seeds RNG; drives every phase |
| **`breeding_rand`** | 2701 | ~800 | Random-pair crossover pipeline; 5-deep nested loops; stack VLA |
| **`mutation_selection`** | 5386 | ~645 | Central mutation dispatcher (type pick → segment classify → operator → tag/filter) |
| **`breeding_exhaustive`** | 2057 | ~636 | Full pairwise crossover pipeline; duplicate of the `breeding_rand` geometry |
| **`replace_scaffold`** | 7774 | ~520 | Fragment-replacement attachment engine |
| **`selection_elite`** | 10053 | ~396 | Elitism with percent/number/max × combined × extinction × niching |
| **`add_H`** | 8693 | ~390 | Deletion/substitution H-capping + mol reconstruction |
| **`deletion_mutation`** | 6039 | ~360 | Deletion operator with deep aps/size branching |
| **`additions`** | 6406 | ~350 | Addition operator; drives de novo growth |
| **`rand_replacement`** | 6951 | ~330 | Fragment-library draw + 6-stage compatibility guarding |
| **`replace_combine_mols`** | 4471 | ~322 | Replacement geometry (duplicate of recomb/attach) |
| **`iso_replace`** | 14383 | ~280 | Isoswap scaffold replacement |
| **`check_exhaustive`** | 1796 | ~253 | Overlap detection (5-deep nested loop) |
| **`selection_tournament`** | 10456 | ~252 | Tournament selection; 5 `rand()` sites |
| **`roulette`** | 11022 | ~173 | Roulette + SUS wheel (SUS branch flagged "THIS HAS TO BE REDONE") |
| **`niche_sharing` / `niche_crowding`** | 11507 / 12210 | ~172 / ~133 | 8-objective niching pipelines |

### 1.2 Full symbol table by definition order

| Line | Symbol | Responsibility |
|---:|---|---|
| 43 | `GA_Recomb()` / `~GA_Recomb()` | ctor/dtor |
| 97 | `Tor_Env()` / `~Tor_Env()` | torsion-environment helper class ctor/dtor |
| 110 | `input_parameters` | Parse ~90 `ga_*` parameters from `Parameter_Reader` |
| 663 | `input_parameters_selection` | Parse per-selection-method sub-parameters |
| 727 | `initialize` | Load parent + fragment libraries, set DNM/iso flags |
| 813 | `initialize_fraglib` | Compute per-fragment size/ring/aps-distance matrices |
| 886 | `initialize_internal_energy_parms` | Copy internal-energy params into members |
| 907 | `read_library` | Read a mol2 file into a `vector<DOCKMol>` |
| **947** | **`max_breeding`** | **Main GA control loop (see §2.1)** |
| 1785 | `print_molecule_counts` | Print per-generation counters |
| 1796 | `check_exhaustive` | Identify parent pairs with overlapping rotatable bonds |
| 2057 | `breeding_exhaustive` | Full pairwise crossover + minimize + prune + mutate |
| 2701 | `breeding_rand` | Random-pair crossover + minimize + prune + mutate |
| 3515 | `similarity_compare` | Hungarian-RMSD test of two candidate halves |
| 3625 | `activate_half_mol` | Activate one connected half of a molecule for recombination |
| 3738 | `calc_cov_radius` | Sybyl atom type → covalent radius (if/else ladder) |
| 3786 | `recomb_mols` | Translate/rotate/attach two halves → child (legacy) |
| 3953 | `recomb_mols_xover` | DNM-aware version of `recomb_mols` (duplicate geometry) |
| 4166 | `rotate` | Rotate a mol to align bond vectors (manual 3×3 matmul) |
| 4255 | `attach` | Build combined child DOCKMol from two active halves |
| 4471 | `replace_combine_mols` | Recomb+attach for fragment replacement |
| 4801 | `switch_active_halves` | Reset flags, activate opposite halves → second child |
| 4839 | `mol_to_frag` | Wrap a DOCKMol into a `Fragment` |
| 4865 | `prepare_mol_torenv` | Activate origin/target of last bond for torenv check |
| 4893 | `prepare_mut_torenv` | Same for a caller-supplied bond |
| 4922 | `prepare_torenv_indices` | Collect rotatable-bond pairs into `torenv_recheck_indices` |
| 4953 | `master_mut_exhaustive` | Mutate every molecule once until quota/cycle cap |
| 5038 | `master_mut_counter` | Tally & print add/del/sub/repl success counts |
| 5154 | `master_mut_rand` | Randomly pick molecules to mutate per rate |
| 5267 | `segment_id` | Amber-type + AG conformer search → segment topology |
| **5386** | **`mutation_selection`** | **Central mutation dispatcher (see §2.3)** |
| 6039 | `deletion_mutation` | Deletion operator (inactivate smaller side, cap with H) |
| 6406 | `additions` | Addition operator (H→Du, de novo growth) |
| 6765 | `replacement` | Replacement operator wrapper (linker vs scaffold lib) |
| 6836 | `prepare_replacement_segment` | Build segment size/ring/aps metadata for replacement |
| 6951 | `rand_replacement` | Draw library fragments + 6-stage compatibility guard |
| 7294 | `compare_bond_types` | ref.aps × frag.aps bond-type match matrix |
| 7332 | `allowed_bond_combos` | Classify all-same / all-diff / mixed aps bond types |
| 7500 | `orient_frag_to_ref` | Sphere-match orient fragment onto ref segment |
| 7721 | `overlapping_aps` | Du–Du RMSD overlap test between aps |
| 7774 | `replace_scaffold` | Attach oriented fragment to ref (replacement engine) |
| 8302 | `calc_norm_ref` | Ring-normal of ref when its 2 aps are colinear |
| 8331 | `compare_norm` | Compare ref vs frag ring-plane normals |
| 8400 | `calc_norm` | Averaged normalized ring-normal from coords |
| 8500 | `exit_vector` | Cosine between two bond vectors |
| 8542 | `compare_rmsd_bt` | Pick frag aps closest to ref aps respecting bond type |
| 8598 | `rand_H_to_Du` | Randomly convert one hydrogen to a dummy atom |
| 8656 | `H_to_Du` | Set an atom to type Du, charge 0 |
| 8693 | `add_H` | Cap a cut site with H, rebuild + minimize |
| 9087 | `cleanup_mutation_selection` | Clear mutation scratch state (load-bearing) |
| 9109 | `prepare_internal_energy` | Configure scorer for internal-energy calc |
| 9145 | `minimize_children` | Type/charge/minimize/score a child |
| 9216 | `uniqueness_prune_mut` | Dedup children vs parents (records pruned parents) |
| 9336 | `uniqueness_prune` | Dedup children vs parents (no pruned-parent tracking) |
| 9459 | `fitness_pruning` | Prune ensemble by score/IE/similarity cutoffs |
| 9566 | `hard_filter` | Batch filter a mol vector (MW/rot/HA/HD/charge) |
| 9671 | `hard_filter_mol_HA_HD` | Single-mol filter incl. HA/HD |
| 9749 | `hard_filter_mol` | Single-mol filter (HA/HD commented out — drifted copy) |
| 9832 | `calc_descriptors` | Dispatcher for the four descriptor calculators |
| 9850 | `calc_mol_wt` | Sum atomic weights → `mol.mol_wt` |
| 9911 | `calc_rot_bonds` | Count rotatable bonds → `mol.rot_bonds` |
| 9936 | `num_HA_HD` | Count H-bond acceptors/donors |
| 9969 | `calc_formal_charge` | Sum partial charges → `mol.formal_charge` |
| **9998** | **`selection_method`** | **Top-level selection dispatcher (see §2.4)** |
| 10053 | `selection_elite` | Elitism selection |
| 10456 | `selection_tournament` | Tournament selection |
| 10715 | `tournament_trad` | Traditional pair comparison (by score) |
| 10736 | `tournament_niche` | Niche pair comparison (by rank/crowding) |
| 10769 | `selection_roulette` | Roulette-wheel selection |
| 10913 | `selection_sus` | Stochastic Universal Sampling (disabled at input) |
| 11022 | `roulette` | Cumulative-fitness wheel core (roulette + SUS) |
| 11206 | `selection_metropolis` | Metropolis MC selection (disabled at input) |
| 11330 | `mc_metropolis` | Metropolis accept (older, buggy variant) |
| 11423 | `true_mc_metropolis` | Corrected Metropolis accept |
| 11507 | `niche_sharing` | Fitness sharing over up to 8 objectives |
| 11686 | `distance_matrix` | Pairwise \|Δscore\| matrix per objective |
| 11817 | `niche_radius` | σ from mean nearest-neighbor distance |
| 11865 | `sharing_function` | Shared fitness = score / niche count |
| 12068 | `score_scaling` | Min-max rescale each objective |
| 12210 | `niche_crowding` | NSGA-II rank + crowding distance |
| 12350 | `increment_rank` | Add positional index to rank |
| 12367 | `crowding_dist` | NSGA-II crowding distance per objective |
| 12527 | `deactivate_mol` / 12550 `activate_mol` / 12573 `activate_vector` | Toggle atom/bond active flags |
| 12748 | `renumber_atom_numbering` | Rename atoms per element (H1, C1, …) |
| 12878 | `is_active_vec` / 12909 `is_active` | Debug dumps of active flags |
| 12935 | `print_torenv` | Debug dump of torsion-environment table |
| 12965 | `time_seconds` | Wall-clock seconds |
| 12985 | `score_parents` | Re-type/charge/rescore each parent |
| 13045 | `unique_parents` | Dedup parents via Hungarian RMSD |
| 13116 | `naming_function` | Build unique molecule title strings |
| 13218 | `calc_pairwise_distance` | Debug warn on clashing atom pairs |
| 13249 | `pairwise_tanimoto_calc` | Pairwise Tanimoto → CSV, optional prune |
| 13348 | `pairwise_hms_calc` | Pairwise Hungarian score → txt matrix |
| 13396 | `set_DNM_bools` / 13415 `set_DNM_bools_start` | Set "do-not-mutate" flags from tags |
| 13444 | `frag_erase` | Erase an element from a `vector<Fragment>` |
| 13461 | `print_molecules` | Write mol2 + descriptor header block |
| **13624** | **`mw_cutoff`** | **Soft/hard molecular-weight gate (6.13 feature)** |
| 13678 | `read_frag_library` | Read a mol2 fragment library |
| 14087 | `prepare_isolibrary` | Read/register isosteric fragment library |
| 14212 | `iso_combine_fragments` | Geometrically attach two iso fragments |
| 14287 | `iso_addition_sidechain` | Iso-based sidechain expansion |
| 14383 | `iso_replace` | Iso-based scaffold replacement |
| 14667 | `get_stddev_changes` | Ramp std-dev for RDKit drive — **entirely `#if 0` (dead)** |

---

## 2. Concern clustering

The ~90 methods group into eleven concerns. The right-hand column is the candidate Phase 2 translation unit (names indicative, to be confirmed in decision **B**).

| Concern | Key symbols | Candidate unit |
|---|---|---|
| **GA control loop / generations** | `max_breeding`, `print_molecule_counts`, `time_seconds` | stays in `conf_gen_ga.*` (orchestrator) |
| **Selection strategies** | `selection_method`, `selection_elite`, `selection_tournament`, `tournament_trad/niche`, `selection_roulette`, `selection_sus`, `roulette`, `selection_metropolis`, `mc_metropolis`, `true_mc_metropolis` | `ga_selection.*` |
| **Niching (sharing/crowding)** | `niche_sharing`, `distance_matrix`, `niche_radius`, `sharing_function`, `score_scaling`, `niche_crowding`, `increment_rank`, `crowding_dist` | `ga_niching.*` (or under selection) |
| **Crossover / breeding** | `check_exhaustive`, `breeding_exhaustive`, `breeding_rand`, `similarity_compare`, `activate_half_mol`, `switch_active_halves` | `ga_crossover.*` |
| **Molecule construction / geometry** | `recomb_mols`, `recomb_mols_xover`, `rotate`, `attach`, `replace_combine_mols`, `calc_cov_radius`, `calc_norm*`, `exit_vector`, `mol_to_frag` | `ga_geometry.*` |
| **Mutation operators** | `mutation_selection`, `master_mut_*`, `deletion_mutation`, `additions`, `replacement`, `rand_replacement`, `replace_scaffold`, `add_H`, `rand_H_to_Du`, `H_to_Du`, `segment_id`, `cleanup_mutation_selection`, `prepare_*_torenv` | `ga_mutation.*` |
| **Fragment compatibility guards** | `compare_bond_types`, `allowed_bond_combos`, `compare_rmsd_bt`, `compare_norm`, `overlapping_aps`, `orient_frag_to_ref`, `prepare_replacement_segment` | `ga_fragment_match.*` |
| **Isoswap** | `prepare_isolibrary`, `iso_combine_fragments`, `iso_addition_sidechain`, `iso_replace` | `ga_isoswap.*` (currently dead — `ga_num_iso_picks=0`) |
| **Fitness / scoring interface** | `minimize_children`, `prepare_internal_energy`, `score_parents`, `fitness_pruning`, `uniqueness_prune*` | `ga_fitness.*` |
| **Filters & cutoffs / descriptors** | `hard_filter`, `hard_filter_mol*`, `calc_descriptors`, `calc_mol_wt`, `calc_rot_bonds`, `num_HA_HD`, `calc_formal_charge`, `mw_cutoff` | `ga_filters.*` |
| **Population / ensemble bookkeeping & I/O** | `activate*`/`deactivate_mol`, `renumber_atom_numbering`, `naming_function`, `print_molecules`, `read_library`, `read_frag_library`, `pairwise_*_calc`, `set_DNM_bools*`, `frag_erase`, `unique_parents`, debug dumps | `ga_population.*` / `ga_io.*` |
| **Parameter parsing** | `input_parameters`, `input_parameters_selection` | `ga_params.*` (or stays with orchestrator) |

### 2.1 The control loop (`max_breeding`, 947–1782)

Seeds RNG once (`srand`, 967), builds/scores/sorts the parent generation (Gen 0 setup, ~1069–1280), then `for (i = 0; i < max_generations; i++)`:

1. **Breeding dispatch** (~1299–1402): optional `check_exhaustive` → `breeding_rand` **or** `breeding_exhaustive`. These internally run crossover → torenv filter → `minimize_children` → `fitness_pruning` → `hard_filter` → mutations → divergent-molecule prune.
2. **RDKit filtering** (~1404–1522): **`#if 0` — disabled.**
3. **Offspring stats** (~1548–1574): average-score reductions.
4. **Selection** (~1624–1665): optional niching rescore, extinction handling (disabled), then `selection_method`.
5. **Output** (~1684–1745): write `restart####`, `pruned`, `filtered_`, `rejected` mol2 files.
6. **Timing / mutation stats** (~1752–1764).

### 2.2 Crossover dispatch

Three boolean parameters (`ga_check_only`, `ga_check_overlap`, `ga_xover_sampling_method_rand`) select among `check_exhaustive` (identify overlapping pairs only, no children, no RNG), `breeding_exhaustive` (enumerate all pairs, no RNG), and `breeding_rand` (sample pairs via `rand()`). The three overlap-detection loops and the recomb/attach geometry are **triplicated / duplicated** across these functions.

### 2.3 Mutation dispatch (`mutation_selection`)

Four operators, keyed by macros in the header (`DELETION_TYPE 0`, `ADDITION_TYPE 1`, `SUBSTITUTION_TYPE 2`, `REPLACEMENT_TYPE 3`):

- **Deletion** → `deletion_mutation` (6039)
- **Addition** → `additions` (6406) → de novo `DN_GA_Build`
- **Substitution** → *no dedicated function*; it is `deletion_mutation(…true…)` **then** `additions(…true…)` (a `subst` bool woven through both)
- **Replacement** → `replacement` (6765) → `rand_replacement` (6951) → `replace_scaffold` (7774)

Dispatch is a weighted-vector-then-random-pick (type chosen at 5434) followed by an `if / else-if` chain on `mutation_type` at ~5819/5857/5887/5937. Segments are classified by attachment-point count into SIDECHAIN(1)/LINKER(2)/SCAFFOLD(≥3)/RIGID(0), which constrains which operator can run on which segment. **Isoswap** (`iso_*`) is an orthogonal expansion invoked inside addition/replacement when `ga_num_iso_picks > 0` (currently 0 → dead).

### 2.4 Selection dispatch (`selection_method`, 9998)

An `if / else-if` chain over **five separate booleans** (`ga_selection_method_{elitism,tournament,roulette,sus,metropolis}`, 10020–10039), with **no default/error arm**. Only elitism/tournament/roulette are reachable — the input validator (`ga_selection_method` allowed-values string, ~554) advertises only those three, so `sus`/`metropolis` dispatch arms are dead. **Niching** is a cross-cutting modifier gated by `ga_niching` (hard-disabled, forced false) that every strategy re-implements the same way.

---

## 3. State & coupling map

### 3.1 The "blackboard"

`GA_Recomb` is a de-facto global blackboard. Methods rarely pass molecule collections as parameters; instead they read and mutate shared members. The heavy shared mutable containers (all `std::vector<DOCKMol>` unless noted, declared in `conf_gen_ga.h` ~233–275):

`parents`, `xover_parents`, `tmp_parents`, `saved_parents`, `mutated_parents`, `children`, `tmp_children`, `pruned_children`, `scored_generation`, `mutants`, `divergent_children`, `pruned_parents`, `filtered_mols`, `tmp_parent1/2`, `new_tmp_parent1/2`, plus fragment libraries `tmp_scaffolds/linkers/sidechains` and `orig_segments`.

Segment/mutation scratch state that operators read and must remember to clear via `cleanup_mutation_selection`: `num_segments`, `num_segs_removed`, `no_seg_bonds`, `bond_btwn_seg`, `terminal_seg`, `seg_exclude_indices`, `dnm_encountered`, `dnm_enabled`, and the `success_add/del/sub/replace`, `total_muts`, `total_mut_attempts` counters.

Run-state counters: `current_generation`, `generated_molecule_counter`, `valid_torenv_molecule_counter`, `unpruned_molecule_counter`, `new_temp`.

**Only true `static` state** is three `static const` formatting constants (`DELIMITER`, `FLOAT_WIDTH`, `STRING_WIDTH`, header ~106–108). Everything else is per-instance — but because a single `GA_Recomb` instance runs the whole GA, the members behave like globals. This coupling is the central obstacle to unit testing (Phase 1B) and to splitting files (Phase 2 step 2).

### 3.2 Configuration state

~120 `ga_*` members hold parsed configuration (files, crossover/mutation/selection/filter parameters). They are set once in `input_parameters` / `input_parameters_selection` and read everywhere. Several are **hardcoded rather than parsed** (forced values), which creates large dead branches: `ga_num_iso_picks=0`, `ga_niching=false`, `ga_selection_extinction=false`, `ga_use_dn_roulette=false`, `ga_check_only=false`, `ga_use_torenv_table=true`.

### 3.3 RNG consumption — full enumeration (**determinism linchpin**)

**Generator:** C stdlib `rand()` / `srand()` (`<cstdlib>`). **Not** a C++ `<random>` engine.
**Seeding:** exactly one call, `srand(simplex.random_seed)` at [`conf_gen_ga.cpp:967`](../src/dock/conf_gen_ga.cpp), inside `max_breeding`, before the generation loop. The seed reaches the file **only** through `simplex.random_seed` (the `Simplex_Minimizer` passed into `max_breeding`). A comment at ~6982 ("rand has been seeded iff the minimizer is turned on") implies seeding effectively depends on the minimizer being active — **verify** this invariant before touching RNG plumbing.
**Historical instability:** per-generation reseed lines (`//srand((simplex.random_seed + i)*10)` at 1064, 1776) are commented out — evidence that reseeding was tried and abandoned.

**35 bare `rand()` sites + `random_shuffle` usage.** By concern:

| Site(s) | Function | Purpose |
|---|---|---|
| 2789, 2790 | `breeding_rand` | Pick two parent indices per breeding iteration |
| 5209 | `master_mut_rand` | Pick a molecule to mutate |
| 5434 | `mutation_selection` | Pick mutation **type** from weighted vector |
| 5601 | `mutation_selection` | Pick a segment **type** to mutate |
| 5642, 5659, 5687, 5702, 5723, 5738 | `mutation_selection` | Pick a segment index (variable, data-dependent draw count in DNM reselect loops) |
| 6058 | `deletion_mutation` | Pick a bond (guarded by `if (bonds.size() >> 1)` — a bit-shift, not a comparison; likely a latent bug) |
| 6983 | `rand_replacement` | Pick a library fragment (draw count depends on rejects before a hit) |
| 8623 | `rand_H_to_Du` | Pick a hydrogen to convert to dummy |
| 10521, 10522, 10549, 10550, 10573, 10574, 10619, 10620, 10650, 10651 | `selection_tournament` | Pick tournament pair indices |
| 11068, 11069, 11105 | `roulette` | Roulette/SUS wheel draws (`rand()/RAND_MAX`) |
| 11400 | `mc_metropolis` | Metropolis accept draw (dead path) |
| 11476 | `true_mc_metropolis` | Metropolis accept draw |
| 13637, 13649 | `mw_cutoff` | Soft-MW probabilistic accept |
| 14335, 14337, 14418, 14420 | `iso_addition_sidechain`, `iso_replace` | Iso tail picks (dead — `ga_num_iso_picks=0`) |

**Draw-count is data-dependent in several places** (DNM reselect loops in `mutation_selection`; reject-before-hit loops in `rand_replacement`, `rand_H_to_Du`, iso picks). This means the RNG stream position after a mutation depends on molecule content, not just call order — so even a refactor that "preserves call order" can still shift the stream if it changes how many candidates are examined. Phase 2 step 5 (centralizing RNG behind a thin `rand()` wrapper) must be done alone and proven bit-identical.

### 3.4 Cross-file (public) surface

Anything declared in `conf_gen_ga.h` that other TUs use must stay stable. `conf_gen_ga.h` is included by `dock.cpp`, `conf_gen_dn.cpp`, `master_conf.cpp` (per the Makefile dependency lines). The externally-visible entry is the `GA_Recomb` class and its `input_parameters` / `initialize` / `max_breeding` entry methods. Internal helpers (the ~85 other methods) are free to move. `dock.h` is included transitively; nothing in `conf_gen_ga.*` appears to *add* to `dock.h`.

---

## 4. Build wiring

- **Compiler/flags** come from `install/config.h` + `install/rules.h`, generated by `install/configure` from one of the platform profiles in `install/` (e.g. `gnu`, `homebrew`, `intel`, `gnu.rdkit`, `intel.intelmpi.parallel`, …). `src/dock/Makefile` does `include ../../install/rules.h` and `../../install/config.h`.
- **Object list:** `conf_gen_ga.o` is listed in `OBJS` in [`src/dock/Makefile`](../src/dock/Makefile) (line ~10) and built by the implicit `.cpp.o` rule using `CXXFLAGS`.
- **Explicit header deps:** lines ~97–105 list `conf_gen_ga.o`'s prerequisite headers (amber_typer, dockmol, fragment, master_score, the score_* family, hungarian, fingerprint, conf_gen_ag, conf_gen_dn, fraggraph, gasteiger, trace, nab, etc.).
- **Link:** all `OBJS` link into `dock6$(DOCK_SUFFIX)`; `make install` moves it to `../../bin`. Targets present: `all`, `install`, `dock`, `clean`. **There is no `test` target** — Phase 1 adds an optional one.

### 4.1 Conditional-compilation macros affecting this file

| Macro | Effect on `conf_gen_ga.*` | Set by |
|---|---|---|
| `BUILD_DOCK_WITH_RDKIT` | Enables RDKit includes (`rdtyper.h`), RDKit descriptor-drive members/params, and several `#ifdef` blocks (18, 166, 9646, 9737, 9815, 13530, header 35, 527–568). **Most of the descriptor-drive body is *additionally* `#if 0`, so even with RDKit on, much stays off.** | `install/gnu.rdkit`, `install/homebrew.rdkit` (`-DBUILD_DOCK_WITH_RDKIT`) |
| `BUILD_DOCK_WITH_MPI` | **No effect on this file** — no `#ifdef BUILD_DOCK_WITH_MPI` anywhere in `conf_gen_ga.*`. Affects `dock.cpp`/`base_mpi.*`. | `install/intel.intelmpi.parallel`, `install/sgi.parallel`, `install/bluegene*` (`-DBUILD_DOCK_WITH_MPI`) |
| `#if 0` blocks (in-file) | Disable RDKit drive (~167–234, 1404–1522), `get_stddev_changes` (~14665), and other experimental paths regardless of RDKit. | — |

**Build configurations to keep green:** serial (`gnu`/`homebrew`), MPI (`*.parallel`), and ±RDKit (`*.rdkit`). Because this file has no MPI code, MPI "safety" reduces to *not breaking compilation* under `-DBUILD_DOCK_WITH_MPI`; there is no rank logic here to preserve.

---

## 5. Extension-point audit (drives Phase 3)

For each Goal-#1 seam: how you add one **today**, and why it hurts.

### 5.1 Selection strategy — **most painful; #1 target**

**To add one today, edit ~7 sites across 2 files:** (1) new `bool ga_selection_method_<new>` + secondary flag in the header; (2) initialize both to false in the ctor and in `input_parameters`; (3) extend the `query_param` allowed-values string (~554), else the validator rejects it; (4) add a `.compare()` arm to the primary decode cascade (~563–576) and the secondary cascade (~618–631); (5) add sub-parameter parsing in `input_parameters_selection`; (6) add an `else if` arm to `selection_method` (~10020–10039); (7) implement `selection_<new>` and hand-replicate the ~100–250 lines of extinction/niching/separate-vs-combined scaffolding every existing strategy copies.

**Why it hurts:** method identity is 5 parallel booleans (invariant "exactly one true" is unenforced; dispatcher has no default), duplicated into a parallel "secondary/extinction" set; each strategy re-implements the same boilerplate; the niching-hook pattern is copy-pasted into every strategy.

### 5.2 Mutation operator — **645-line monolith**

**To add one today:** add a `#define NEW_TYPE 4` + tag macro + `int new_co` weight + `success_new` vector + method decl in the header; add an enable/weight block, a "trying type" print, segment-eligibility rules, a new `else if (mutation_type == NEW_TYPE)` dispatch arm, and tagging/`success_new` bookkeeping inside `mutation_selection`; add tally loops in **both** the parents and offspring branches of `master_mut_counter` (duplicated); add input parsing.

**Why it hurts:** no operator interface/registry — the operator is grafted into a 645-line function plus 3–4 others in lockstep; operators communicate through mutable scratch globals with manual cleanup; substitution is already an implicit delete+add, blurring operator boundaries; many magic constants are inlined.

### 5.3 Crossover / breeding

**To add one today:** add `ga_*` flag(s) + parsing; add an `else if` arm to the boolean ladder in `max_breeding`; write a new `breeding_*` that almost certainly copy-pastes the triplicated overlap-detection loop and the duplicated `recomb_mols`/`recomb_mols_xover`/`replace_combine_mols` geometry.

**Why it hurts:** boolean-ladder dispatch (no strategy/enum); triplicated overlap loops; duplicated geometry; reliance on shared scratch members instead of passed buffers.

### 5.4 Fitness / scoring option

The scoring boundary is `minimize_children` / `prepare_internal_energy` / `score_parents` calling `Master_Score`. Toggles like internal-energy and (dead) RDKit descriptor-drive are threaded as `if (use_internal_energy)` / `#if 0` branches rather than composed. **To add a scoring adjustment today:** thread a new flag through `minimize_children`, `fitness_pruning` (cutoffs), and possibly the selection fitness basis. **Why it hurts:** no single scoring-adjustment seam; ligand-efficiency-style adjustments would branch in several places.

### 5.5 Filter / cutoff

**To add one today:** add a field to `DOCKMol`; write a `calc_*` function; call it in `calc_descriptors`; add a constraint member + `query_param` read; add the inline `if (field > constraint)` check + counter + `cout` in **all three** `hard_filter*` copies — which have **already drifted** (`hard_filter_mol` has HA/HD commented out).

**Why it hurts:** the 5-check block is triplicated and inconsistent; constraints are ungrouped scalars; no filter registry.

**Worked example already in tree — the soft MW cutoff (6.13):** `mw_cutoff` (13624) implements a Metropolis-style soft accept: for a mol over the upper (or under the lower) bound it computes `Z = excessMW / ga_constraint_mol_wt_std_dev`, `acceptRate = exp(-Z²)`, and rejects only if `acceptRate < rand()/100`. **Clarification for the spec:** `mw_cutoff` uses a *fixed* std-dev (`ga_constraint_mol_wt_std_dev`, default 35.0). `get_stddev_changes` (14667) — which the spec pairs with it — is a **separate, `#if 0` dead** helper that only ramps std-dev for the RDKit *drive* descriptors, not for MW. The "soft cutoff with ramping std-dev" concept is not wired into the live build.

---

## 6. Risk register

Ordered by blast radius. These are the regions where a "harmless" refactor most easily changes scientific output.

1. **Single global RNG stream (highest).** One `srand` at 967; ~38 ordered draws; several with **data-dependent draw counts** (§3.3). Any reordering, added/removed draw, or changed candidate-examination count anywhere shifts every later result. *Mitigation:* golden-master regression before/after every step; do RNG centralization (Phase 2 step 5) in isolation with a bit-identical proof.
2. **`std::sort` instability / tie-breaking.** Selection and niching sort by `current_score` / `rank` / `fitness` with comparators (`mol_sort`, `compare_rank`, `fitness_sort`) that are not guaranteed strict-weak or unique-keyed. Equal-key molecules get an implementation-defined order, which then decides who breeds. `fitness_pruning` even sorts *inside* a loop. *Mitigation:* treat comparators as behavior; if they must change, prove identical ordering on goldens; consider `stable_sort` only as an explicit, approved, separately-reviewed change (it *can* alter output).
3. **Float accumulation order.** Naive sequential single-precision sums: average-score reductions (control loop), roulette `fitness_sum`/wheel, `niche_radius`, `sharing_function`, `calc_norm`, distance matrices, geometry (`sin_theta = sqrt(1-cos²)`, manual 3×3 matmul in `rotate`). Near-threshold comparisons (`ga_bond_tolerance`, `ga_angle_cutoff`, hardcoded 0.5/1.0/0.7071 cutoffs) can flip on any reordering. *Mitigation:* preserve summation order; never parallelize these without approval + quantified tolerance.
4. **Erase-while-iterating / order-dependent containers.** `fitness_pruning`, tournament winner-erase, `roulette` erase, `master_mut_*` erase-then-push_back all reorder collections mid-loop; downstream sampling depends on the resulting order. *Mitigation:* preserve exact traversal/erase semantics on extraction.
5. **Stack VLAs.** `int ppairs[parents.size()][parents.size()]` in `breeding_rand` (2766) and similar — large-N stack risk; also a portability wart. *Note but do not "fix" during a pure refactor* unless approved.
6. **Dead/disabled code that is not obviously dead.** RDKit descriptor-drive (`#if 0` **and** `#ifdef RDKIT`), niching, extinction, SUS/Metropolis selection, `get_stddev_changes`, iso-swap (`ga_num_iso_picks=0`). Some is reachable only in build configs not exercised here. *Mitigation:* per non-goal in the spec, do **not** delete without confirming unreachable in **all** configs; prefer leaving in place during Phase 2.
7. **Latent bugs to preserve-or-fix-deliberately (not silently).** `if (bonds.size() >> 1)` bit-shift guard (6058); `selection_metropolis` combined branch passing `tmp_parents` instead of `temp` (~11287); `mc_metropolis` off-by-one start (`initial_num-1`); `score_scaling` objective-6 max using `score_nrg` instead of `score_hun` (~12173); `if` vs `else if` for objectives 6/7 in `distance_matrix`/`sharing_function`/`crowding_dist`; `prepare_mol_torenv`/`prepare_mut_torenv` set bond flags `true` despite "false" comments. These change output if "fixed"; each is a separate, flagged, approved change per the working agreement — never bundled into a refactor.
8. **No MPI here, but keep the MPI build green.** Adding sibling `.cpp` files or headers must compile under `-DBUILD_DOCK_WITH_MPI`.

---

## 7. Building, testing & running

### 7.1 Normal DOCK build (unchanged by this work)

```sh
cd install
./configure gnu           # or: homebrew, gnu.rdkit, intel.intelmpi.parallel, ...
cd ../src/dock
make                      # builds dock6; make install moves it to ../../bin
```

- Serial: `./configure gnu` (Linux) / `./configure homebrew` (macOS; wants Homebrew `gcc-11`/`gfortran-11`).
- RDKit: `./configure gnu.rdkit` (needs `$RDBASE`, `$BOOST`).
- MPI: `./configure intel.intelmpi.parallel` (or another `*.parallel` profile).

The test infrastructure below is **purely additive** — `make` and `make install` do not build it, require it, or gain any new dependency or step.

### 7.2 Test infrastructure (Phase 1 — landed)

Framework: **doctest v2.4.11**, vendored as a single header at `src/dock/tests/doctest.h` (MIT), nothing to install (decision **F**). Full details in [`src/dock/tests/README.md`](../src/dock/tests/README.md).

Three layers of safety net:

```sh
# 1. Unit suite (doctest). Integrated (uses your configured compiler):
cd src/dock && make test
#    Standalone (no ./configure needed — works even where DOCK can't fully build):
cd src/dock/tests && make run

# 2. Compile-check gate — fast syntax-only check of conf_gen_ga.cpp after each refactor step:
cd src/dock && make check-compile
#    or directly:  src/dock/tests/compile_check.sh   (add --rdkit if $RDBASE/$BOOST are set)

# 3. End-to-end golden-master regression (run on your cluster — see below):
src/dock/tests/regression/check_case.sh <case_dir> /path/to/dock6
```

- `make test` from `src/dock` needs a configured tree (`install/config.h`) like every target there; the **standalone** `cd tests && make run` path needs no configure and is how tests run in the refactor dev environment.
- The unit suite currently holds framework + determinism **smoke tests**; real characterization/seam tests are added *as Phase 2 extracts low-dependency units* (spec §5.3). Adding-a-test instructions are in the tests README.

### 7.3 End-to-end regression is the user's pre-merge responsibility

A full `dock6` build needs a **Fortran toolchain** (nab/grid/score modules), absent in the refactor dev environment, and **no de-novo/GA case ships in `tutorials/`** (they are amber/bks/ligand-sampling/solvent/mpi demos). So end-to-end golden capture and diffing must run **on your cluster**, against real receptor grids + fragment libraries + seeded `dock.in` files. The harness (`run_case.sh` / `capture_golden.sh` / `check_case.sh` / `compare.py`) is built and self-tested here with a stub binary — it catches score drift and count/ordering changes exactly — you supply the cases and the built binary. Workflow and coverage checklist: [`src/dock/tests/regression/README.md`](../src/dock/tests/regression/README.md).

In the dev environment, the standing safety net after every Phase 2 step is therefore: **`make check-compile` + `make test`**, with the end-to-end goldens run on the cluster before merge (decision **D**). A full HPC how-to (build baseline + refactored, capture goldens, diff, SLURM example) is in [`hpc_golden_verification.md`](hpc_golden_verification.md).

---

## 8. Decisions (accepted)

The spec's Appendix (A–F) plus one new decision (G) surfaced by discovery. **The user accepted the conservative defaults for all of A–G**, with one refinement to D: end-to-end golden verification runs on the user's cluster at the final build, and the standing per-step gate in the dev environment is **compile-check + unit tests** (full-speed local verification). Defaults, for the record:

- **A. Reproducibility target.** *Default: exact same-seed bit-for-bit.* Confirm the lab requires exact reproducibility (recommended given the single global RNG stream), vs. tolerating small documented float differences across builds.
- **B. File splitting.** *Default: allowed*, into sibling `.cpp/.h` wired into `src/dock/Makefile` like existing siblings, with `conf_gen_ga.*` remaining the orchestrator. Confirm, or restrict Phase 2 to in-file extraction only. (Candidate units in §2.)
- **C. C++ standard.** *Default: whatever the profile sets* (the RDKit profiles use `-std=c++11`; serial profiles set none explicitly). Add nothing requiring a newer standard.
- **D. MPI test access.** *Note:* this file has no MPI code, so serial characterization covers its logic fully; MPI verification is only "does the `-DBUILD_DOCK_WITH_MPI` build still compile/link." Confirm the environment can at least compile an MPI config, or accept that as a pre-merge responsibility.
- **E. Priority tilt.** *Default: extensibility > readability > testability > performance.* Discovery supports leading with the **selection** and **mutation** seams (most painful).
- **F. Test framework & wiring.** *Default: vendored single-header **doctest** under `src/dock/tests/`, optional `make test`, normal build untouched.* Confirm framework choice and that local `make test` (no CI) is acceptable for now.
- **G. License header for new files (new).** `conf_gen_ga.{cpp,h}` currently carry **no** header; siblings (e.g. `dock.cpp`) carry a UCSF copyright-grant block and the repo `LICENSE` is BSD-3-Clause (2026, Regents of UC). Confirm new sibling files should copy the **sibling UCSF header block** (recommended), and whether to also add that header to `conf_gen_ga.*` (a separate, flagged change).

---

## 9. Changelog (structural change → verification)

Each Phase 2/3 change is behavior-preserving and listed with how it was verified. Dev-env
verification is `make check-compile` (conf_gen_ga.cpp: 0 errors) + `make test` (unit suite green);
bit-for-bit end-to-end verification is the user's cluster responsibility (§7.3).

| # | Phase | Change | New/edited files | Verification |
|---|---|---|---|---|
| 1 | 2 (extract) | Extracted the four pure descriptor computations (`calc_mol_wt`, `calc_rot_bonds`, `num_HA_HD`, `calc_formal_charge`) from `GA_Recomb` into free functions in `ga_descriptors.{h,cpp}`; the methods now delegate. Float+double accumulation order/type preserved exactly (incl. the legacy `DN_GA_Build` warning text). | +`ga_descriptors.h`, +`ga_descriptors.cpp`, +`tests/test_ga_descriptors.cpp`; edited `conf_gen_ga.cpp` (4 bodies + include), `Makefile` (OBJS + deps), `tests/Makefile` | compile-check PASS (0 errors, 13 pre-existing warnings, no new); unit suite 10 cases / 49 assertions PASS incl. a bit-for-bit MW-accumulation test; `ga_descriptors.cpp` compiles to `.o` standalone |
| 2 | 2 (extract) | Extracted the MW cutoff (`mw_cutoff`: hard + 6.13 soft path) and the pure threshold predicates (rot bonds, HA, HD, formal-charge range) into `ga_filters.{h,cpp}`. `mw_cutoff` now delegates with `rand()` **injected** as a callback, so the soft path's data-dependent draws occur at the exact same points/count (RNG order preserved); the vector `hard_filter` routes its threshold checks through the predicates. Exact float/double arithmetic and comparison-conversion types preserved. The two drifted single-mol `hard_filter_mol*` twins keep their inline one-liners for the Phase 3 filter registry (their MW check already flows through the delegated `mw_cutoff`). | +`ga_filters.h`, +`ga_filters.cpp`, +`tests/test_ga_filters.cpp`; edited `conf_gen_ga.cpp` (`mw_cutoff` body, `hard_filter` checks, include), `Makefile`, `tests/Makefile` | compile-check PASS (0 errors, 13 warnings, no new); unit suite 17 cases / 78 assertions PASS incl. deterministic soft-cutoff decisions via a stub RNG and draw-count assertions; `ga_filters.cpp` compiles to `.o` standalone |
| 3 | 3 (seam) | **Selection seam.** Replaced the 5-way if/else-if ladder in `selection_method` with a table-driven registry (`{label, enabling-flag member-ptr, handler member-ptr}` rows) + the pure `ga_selection::first_enabled` decision. Adding a strategy no longer touches the dispatcher. Behavior-identical: same rows in the same order, first enabled wins, same labels/handlers, no-op if none set. | +`ga_selection.h`, +`ga_selection.cpp`, +`tests/test_ga_selection.cpp`; edited `conf_gen_ga.cpp` (dispatcher + include), `Makefile`, `tests/Makefile` | compile-check PASS (0 errors, 13 warnings, no new); unit suite 20 cases / 87 assertions PASS incl. first-match-wins / none-set contract; `ga_selection.cpp` compiles to `.o` standalone |
| 4 | 2 (extract) | Extracted STEP 1 of `mutation_selection` — the probability-weighted mutation-type pool — into `ga_mutation::build_weighted_type_pool` (a `{enabled, coefficient, type_code}` table). Push order (deletion, addition, substitution, replacement) and per-type multiplicity preserved exactly; the `rand()%size` pick and `exit(1)`-on-empty stay in the method, so RNG order and side effects are unchanged. The non-uniform operator **dispatch was deliberately left as explicit if/else** (a uniform registry would be forced abstraction — see §10.2). | +`ga_mutation.h`, +`ga_mutation.cpp`, +`tests/test_ga_mutation.cpp`; edited `conf_gen_ga.cpp` (STEP 1 + include), `Makefile`, `tests/Makefile` | compile-check PASS (0 errors, 13 warnings, no new); unit suite 25 cases / 92 assertions PASS; `ga_mutation.cpp` compiles to `.o` standalone |
| 5 | 2 (extract/reuse) | Extracted `calc_cov_radius` into the pure `ga_descriptors::covalent_radius` lookup (values identical to the header `COV_RADII_*`; unknown → 0.71 fallback, warning kept verbatim in the wrapper). Also reused `ga_mutation::build_weighted_type_pool` (coefficient 1 = "offer once if available") for the segment-type candidate pool in `mutation_selection`, replacing the three conditional `push_back`s; order (sidechain, linker, scaffold) and the subsequent `rand()%size` pick unchanged. | edited `ga_descriptors.{h,cpp}`, `tests/test_ga_descriptors.cpp`, `conf_gen_ga.cpp` (2 sites) | compile-check PASS (0 errors, 13 warnings, no new); unit suite 27 cases / 116 assertions PASS incl. covalent-radius recognized/fallback cases; `ga_descriptors.cpp` compiles to `.o` standalone |
| 6 | 3 (seam) | **Completed the selection seam:** consolidated the 5-arm token→bool parse cascade in `input_parameters` into a `{token, enabling-flag}` table + loop, paired with the dispatch registry (change #3). Adding a selection strategy now touches a row in each of two tables instead of an if/else-if arm in both the parser and the dispatcher. Behavior identical: exactly one flag true for a valid token, same `"You chose...poorly."`+`exit(0)` on an unknown token. | edited `conf_gen_ga.cpp` (parse cascade) | compile-check PASS (0 errors, 13 warnings, no new); unit suite 27 cases / 116 assertions PASS |
| 7 | 2 (extract) | Extracted the molecule-title formatter from `naming_function` into `ga_naming::build_molecule_title`. Preserves the exact zero-padding, **including the deliberate gen-vs-loc asymmetry** (loc not padded to width 4 for loc ≥ 100) — pinned by tests so it can't be silently "fixed". The method keeps the (dead-but-preserved) parent-name `list` computation and the `mol.energy != PARENT_TAG` assignment guard. | +`ga_naming.h`, +`ga_naming.cpp`, +`tests/test_ga_naming.cpp`; edited `conf_gen_ga.cpp` (2 blocks + include), `Makefile`, `tests/Makefile` | compile-check PASS (0 errors, 13 warnings, no new); unit suite 30 cases / 129 assertions PASS; `ga_naming.cpp` compiles to `.o` standalone |

**Safe-target batch complete.** Remaining un-refactored logic in `conf_gen_ga.cpp` is dominated by float-sensitive geometry (`exit_vector`, `calc_norm`, `compare_norm`, `rotate`, niching/roulette reductions) and stateful, DOCKMol-coupled routines. Those should be verified against cluster golden references (below) before/after any change, since local unit tests cannot catch floating-point bit drift.

*Phase 0 map complete and accepted; Phase 1 test net complete; Phase 2/3 underway.*

---

## 10. Extension guides ("how to add one")

As each seam lands, its worked "how to add one" recipe is recorded here.

### 10.1 Add a new selection strategy

The dispatch is table-driven (changelog #3), so you no longer edit `selection_method`.

1. **Declare the enabling flag** in `conf_gen_ga.h` next to the existing ones:
   `bool ga_selection_method_<name>;` (initialize it `false` in the constructor).
2. **Parse it** in `GA_Recomb::input_parameters`: extend the `ga_selection_method`
   allowed-values string (currently `"elitism | tournament | roulette"`) and add one row
   `{ "<name>", &GA_Recomb::ga_selection_method_<name> }` to the `selection_options[]`
   table (no new `.compare()` arm). Add any strategy-specific sub-parameters to
   `input_parameters_selection`.
3. **Write the handler** as a `GA_Recomb` method with the strategy signature
   `void selection_<name>( std::vector<DOCKMol> &, Master_Score &, AMBER_TYPER & );`
   (declare it in the header, define it in `conf_gen_ga.cpp`). Reuse the shared
   scaffolding the other strategies use.
4. **Register one row** in the `strategies[]` table inside `selection_method`:
   `{ "<Label>", &GA_Recomb::ga_selection_method_<name>, &GA_Recomb::selection_<name> },`.
   Rows are tried in order; the first enabled flag wins.
5. **Add unit tests**: `tests/test_ga_selection.cpp` for any pure decision logic, plus a
   test of your handler's pure pieces (extract them as free functions where practical, and
   test RNG-dependent parts with an injected stub generator, as `ga_filters` does).

So adding a strategy is: one flag + one parse-table row + one dispatch-registry row + the
handler. Neither the parse loop nor the dispatch loop (nor `ga_selection::first_enabled`)
changes.

*(Possible further tidy: the parse table and the dispatch registry are two tables listing
the same strategies; a single shared registry — e.g. a private static `GA_Recomb`
accessor holding token, label, flag-ptr, and handler-ptr — would make it one. Deferred as
lower-value; the two-table form is already a strict improvement over the old 7 edit sites.)*

### 10.2 Add a new mutation operator

Unlike selection, the mutation operators are **not** behind a uniform registry — and
deliberately so. Deletion / addition / substitution / replacement have different
signatures, tagging (`-d` / `-a` / `-s` / `-r`), DNM handling, and post-processing
(substitution is literally deletion-then-addition), so a common interface would be
forced abstraction that obscures rather than helps. Adding an operator today:

1. **Define its type code + tag** in `conf_gen_ga.h` (next to `DELETION_TYPE` … and
   `DELETION_TAG` …), and add an enable flag `ga_mutate_<name>` + coefficient `<name>_co`
   (parsed in `input_parameters`, and included in the rate-sum-to-100 validation).
2. **Register its weight** in the `type_weights[]` table at STEP 1 of `mutation_selection`
   (`{ ga_mutate_<name>, <name>_co, <NAME>_TYPE }`) so it can be drawn.
3. **Add its dispatch arm** in STEP 4 of `mutation_selection` (the explicit if/else on
   `mutation_type`), with the operator's own tagging/filter/success-bookkeeping, and add
   its tally to `master_mut_counter` (both the parents and offspring branches).
4. **Write the operator** as a `GA_Recomb` method, and **unit-test** its pure pieces by
   extracting them as free functions (test any RNG-dependent part with an injected stub,
   as `ga_filters` does).

The probability-weighting math itself (`ga_mutation::build_weighted_type_pool`) does not
change — only the table row.
