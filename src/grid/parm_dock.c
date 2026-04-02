/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/
#include "define.h"
#include "utility.h"
#include "mol.h"
#include "global.h"
#include "dock.h"
#include "search.h"
#include "label.h"
#include "io.h"
#include "score.h"
#include "score_dock.h"
#include "match.h"
#include "orient.h"
#include "parm.h"
#include "parm_dock.h"

void get_parameters
(
  DOCK *dock,
  ORIENT *orient,
  SCORE *score,
  LABEL *label
)
{
  int i, j;

  PARM parm = {0};

  STRING40 parameter_name;
  STRING100 parameter_value;

  enum FILE_FORMAT format;

/*
* Read parameters into buffer
* 8/96 te
*/
  read_parameters (&parm);

/*
* Begin reading in individual parameters
* 8/96 te
*/
  fprintf (global.outfile, 
    "\n________________General_Parameters________________\n");

  get_parameter
  (
    (void *) &label->flex.flag,
    &parm, Boolean, "flexible_ligand",
    "no",
    TRUE
  );

  get_parameter
  (
    (void *) &orient->flag,
    &parm, Boolean, "orient_ligand",
    "no",
    TRUE
  );

  dock->multiple_orients =
    dock->multiple_orients || orient->flag;

  get_parameter
  (
    (void *) &score->flag,
    &parm, Boolean, "score_ligand",
    "no",
    TRUE
  );

  get_parameter
  (
    (void *) &score->minimize.flag,
    &parm, Boolean, "minimize_ligand",
    "no",
    score->flag
  );

  get_parameter
  (
    (void *) &dock->multiple_ligands,
    &parm, Boolean, "multiple_ligands",
    "no",
    TRUE
  );

  get_parameter
  (
    (void *) &label->chemical.screen.flag,
    &parm, Boolean, "chemical_screen",
    "no",
    !orient->flag && !score->flag && dock->multiple_ligands
  );

  label->chemical.flag = label->chemical.screen.flag;

  label->vdw.flag =
    orient->flag || label->flex.flag || score->flag ||
    dock->multiple_ligands || label->chemical.screen.flag;

  get_parameter
  (
    (void *) &dock->parallel.flag,
    &parm, Boolean, "parallel_jobs",
    "no",
    dock->multiple_ligands
  );

  get_parameter
  (
    (void *) &i,
    &parm, Integer, "random_seed",
    "0",
    label->flex.flag || score->minimize.flag || orient->flag
  );

  srand (i);


  if (label->flex.flag || (global.output_volume == 'v'))
    fprintf (global.outfile,
      "\n____________Flexible_Ligand_Parameters____________\n");

  get_parameter
  (
    (void *) &label->flex.anchor_flag,
    &parm, Boolean, "anchor_search",
    "no",
    label->flex.flag && score->flag
  );

  get_parameter
  (
    (void *) &label->flex.multiple_anchors,
    &parm, Boolean, "multiple_anchors",
    "no",
    label->flex.anchor_flag
  );

  get_parameter
  (
    (void *) &label->flex.anchor_size,
    &parm, Integer, "anchor_size",
    label->flex.multiple_anchors ? "10" : "0",
    label->flex.multiple_anchors
  );

  get_parameter
  (
    (void *) &label->flex.periph_flag,
    &parm, Boolean, "peripheral_search",
    label->flex.anchor_flag ? "yes" : "no",
    FALSE    /* Control over this parameter has been disabled 3/98 te */
  );

  get_parameter
  (
    (void *) &label->flex.write_flag,
    &parm, Boolean, "write_partial_structures",
    "no",
    !dock->multiple_ligands && label->flex.periph_flag
  );

  get_parameter
  (
    (void *) &label->flex.drive_flag,
    &parm, Boolean, "torsion_drive",
    "no",
    label->flex.flag && (!label->flex.anchor_flag || label->flex.periph_flag)
  );

  dock->multiple_conforms = label->flex.drive_flag;

  get_parameter
  (
    (void *) &label->flex.clash_overlap,
    &parm, Real, "clash_overlap",
    label->flex.drive_flag ? "0.5" : "0",
    label->flex.drive_flag
  );

  get_parameter
  (
    (void *) &label->flex.max_conforms,
    &parm, Integer,
    label->flex.anchor_flag
      ? "configurations_per_cycle|peripheral_seeds"
      : "conformation_cutoff_factor|maximum_conformations",
    label->flex.drive_flag ? (label->flex.anchor_flag ? "25" : "5") : "1",
    label->flex.drive_flag
  );

  if ((label->flex.max_conforms <= 0) || (label->flex.max_conforms >= INT_MAX))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "Value for %s not acceptable.\n",
      label->flex.anchor_flag ?
        "configurations_per_cycle" : "conformation_cutoff_factor"));

  get_parameter
  (
    (void *) &label->flex.minimize_flag,
    &parm, Boolean, "torsion_minimize",
    "no",
    label->flex.flag && score->minimize.flag &&
      (!label->flex.anchor_flag || label->flex.periph_flag)
  );

  if (label->flex.anchor_flag &&
    !label->flex.drive_flag &&
    !label->flex.minimize_flag)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: torsion_drive or torsion_minimize must be \n"
      "  selected with anchor_search.\n"));

  get_parameter
  (
    (void *) &label->flex.reminimize_layers,
    &parm, Integer, "reminimize_layer_number",
    label->flex.minimize_flag && label->flex.periph_flag ? "2" : "0",
    label->flex.minimize_flag && label->flex.periph_flag
  );

  get_parameter
  (
    (void *) &label->flex.minimize_anchor_flag,
    &parm, Boolean, "minimize_anchor",
    score->minimize.flag ? "yes" : "no",
    score->minimize.flag && label->flex.periph_flag
  );

  get_parameter
  (
    (void *) &label->flex.reminimize_anchor_flag,
    &parm, Boolean, "reminimize_anchor",
    score->minimize.flag && label->flex.periph_flag ? "yes" : "no",
    score->minimize.flag && label->flex.periph_flag
  );

  get_parameter
  (
    (void *) &label->flex.reminimize_ligand_flag,
    &parm, Boolean, "reminimize_ligand",
    score->minimize.flag && label->flex.periph_flag ? "yes" : "no",
    score->minimize.flag && label->flex.periph_flag
  );

  get_parameter
  (
    (void *) &label->flex.max_torsions,
    &parm, Integer, "flexible_bond_maximum",
    label->flex.flag && dock->multiple_ligands ? "10" : "<infinity>",
    label->flex.flag && dock->multiple_ligands
  );


  if (orient->flag || (global.output_volume == 'v'))
    fprintf (global.outfile,
      "\n_____________Orient_Ligand_Parameters_____________\n");

  get_parameter
  (
    (void *) &orient->match.flag,
    &parm, Boolean, "match_receptor_sites",
    label->chemical.screen.pharmaco_flag &&
      label->chemical.screen.fold_flag ? "yes" : "no",
    orient->flag
  );

  get_parameter
  (
    (void *) &orient->random_flag,
    &parm, Boolean, "random_search",
    "no",
    orient->flag
  );

  if (orient->flag &&
    !orient->match.flag &&
    !orient->random_flag)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: No orient_ligand options selected\n"));

  get_parameter
  (
    (void *) &orient->match.centers_flag,
    &parm, Boolean, "ligand_centers|match_ligand_centers",
    "no",
    orient->flag && !dock->multiple_ligands
  );

  get_parameter
  (
    (void *) &orient->match.auto_flag,
    &parm, Boolean, "automated_matching|uniform_sampling",
    orient->match.flag && !dock->multiple_ligands ? "yes" : "no",
    orient->match.flag
  );

  if ((label->flex.periph_flag == TRUE) && (orient->flag == TRUE))
    fprintf (global.outfile,
      "ATTENTION get_parameters: "
      "For this run, maximum_orientations value\n"
      "  used to allocate memory for array.  So don't set it too large.\n");

  get_parameter
  (
    (void *) &orient->max,
    &parm, Integer, "maximum_orientations|total_orientations",
    orient->flag ? (orient->match.auto_flag ? "500" : "5000") : "1",
    orient->flag
  );

  if ((orient->max <= 0) || (orient->max >= INT_MAX))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "Value for maximum_orientations not acceptable.\n"));

  get_parameter
  (
    (void *) &dock->write_orients,
    &parm, Boolean,
    label->flex.drive_flag ?
      (orient->flag ? "write_configurations|write_orientations"
      : "write_conformations|write_orientations") :
      "write_orientations",
    "no",
    orient->flag || label->flex.drive_flag
  );

  get_parameter
  (
    (void *) &dock->rank_orients,
    &parm, Boolean,
    label->flex.drive_flag ?
      (orient->flag ? "rank_configurations|rank_orientations"
      : "rank_conformations|rank_orientations") :
      "rank_orientations",
    score->flag && dock->write_orients ? "yes" : "no",
    score->flag && dock->write_orients &&
      orient->flag && !label->flex.drive_flag
  );

  get_parameter
  (
    (void *) &dock->rank_orient_total,
    &parm, Integer,
    label->flex.drive_flag ?
      (orient->flag ? "write_configuration_total|rank_orientation_total"
      : "write_conformation_total|rank_orientation_total")
      : "rank_orientation_total",
    dock->rank_orients ? "100" : "1",
    dock->rank_orients
  );

  if ((dock->rank_orient_total <= 0) || (dock->rank_orient_total >= INT_MAX))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: Value for %s not acceptable.\n",
      label->flex.drive_flag ?
        (orient->flag ? "write_configuration_total" :
          "write_conformation_total") :
        "rank_orientation_total"));


  if (label->flex.periph_flag == TRUE)
  {
    dock->rank_anchors = TRUE;
    dock->rank_anchor_total = orient->flag ? orient->max : 1;
  }

  else
  {
    dock->rank_anchors = dock->rank_orients;
    dock->rank_anchor_total = dock->rank_orient_total;
  }


  if ((orient->match.flag && !orient->match.auto_flag) ||
    (global.output_volume == 'v'))
    fprintf (global.outfile,
      "\n_________________Match_Parameters_________________\n");

  get_parameter
  (
    (void *) &orient->match.clique_size_min,
    &parm, Integer, "nodes_minimum",
    orient->match.flag ? "3" : "0",
    orient->match.flag && !orient->match.auto_flag
  );

  get_parameter
  (
    (void *) &orient->match.clique_size_max,
    &parm, Integer, "nodes_maximum",
    orient->match.flag ? "10" : "0",
    orient->match.flag && !orient->match.auto_flag
  );

  get_parameter
  (
    (void *) &orient->match.distance_tolerance,
    &parm, Real, "distance_tolerance",
    orient->match.flag ? "0.25" : "0",
    orient->match.flag && !orient->match.auto_flag
  );

  get_parameter
  (
    (void *) &orient->match.distance_minimum,
    &parm, Real, "distance_minimum",
    orient->flag && orient->match.flag ? "2.0" : "0",
    orient->flag && orient->match.flag && !orient->match.auto_flag
  );

  get_parameter
  (
    (void *) &orient->match.degeneracy_flag,
    &parm, Boolean, "check_degeneracy",
    "no",
    orient->flag && orient->match.flag && !orient->match.auto_flag
  );

  get_parameter
  (
    (void *) &orient->match.reflect_flag,
    &parm, Boolean, "reflect_ligand",
    "no",
    orient->flag && orient->match.flag && !orient->match.auto_flag &&
      (orient->match.clique_size_max > 3)
  );

  get_parameter
  (
    (void *) &orient->match.critical_flag,
    &parm, Boolean, "critical_points",
    "no",
    orient->flag && orient->match.flag && !orient->match.auto_flag
  );

  get_parameter
  (
    (void *) &orient->match.multiple_flag,
    &parm, Boolean, "multiple_points",
    "no",
    orient->match.critical_flag
  );

  get_parameter
  (
    (void *) &orient->match.chemical_flag,
    &parm, Boolean, "chemical_match",
    label->chemical.screen.pharmaco_flag ? "yes" : "no",
    orient->flag && orient->match.flag && !orient->match.auto_flag
  );

  label->chemical.flag =
    label->chemical.flag || orient->match.chemical_flag;



  if (score->flag || (global.output_volume == 'v'))
    fprintf (global.outfile,
      "\n________________Scoring_Parameters________________\n");

  get_parameter
  (
    (void *) &score->intra_flag,
    &parm, Boolean, "intramolecular_score",
    label->flex.flag && score->flag ? "yes" : "no",
    label->flex.flag && score->flag
  );

  get_parameter
  (
    (void *) &score->inter_flag,
    &parm, Boolean, "intermolecular_score",
    score->flag ? "yes" : "no",
    score->flag
  );

  if (score->flag && !score->intra_flag && !score->inter_flag)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: No scoring options selected\n"));

  get_parameter
  (
    (void *) &score->grid.flag,
    &parm, Boolean, "gridded_score",
    score->inter_flag ? "yes" : "no",
    score->inter_flag
  );

  get_parameter
  (
    (void *) &score->grid.version,
    &parm, Real, "grid_version",
    "4.0",
    score->grid.flag
  );

  if ((score->grid.version < 3.0) || (score->grid.version > 4.0))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: grid_version selection not supported\n"));

  get_parameter
  (
    (void *) &score->grid.size,
    &parm, Integer, "grid_points",
    (score->grid.version < 4) ? "1000000" : "0",
    score->grid.version < 4
  );

  get_parameter
  (
    (void *) &score->grid.spacing,
    &parm, Real, "receptor_atom_grid_spacing",
    score->inter_flag && !score->grid.flag ? "2.0" : "0.0",
    FALSE
  );

  if (score->inter_flag && !score->grid.flag && (score->grid.spacing <= 0.0))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: Innappropriate grid_spacing value.\n"));

  get_parameter
  (
    (void *) &score->bump.flag,
    &parm, Boolean, "bump_filter",
    "no",
    orient->flag && score->grid.flag
  );

  get_parameter
  (
    (void *) &score->bump.maximum,
    &parm, Integer, "bump_maximum",
    "0",
    score->bump.flag
  );

  score->type[NONE].flag = !score->flag;
  strcpy (score->type[NONE].name, "none");
  strcpy (score->type[NONE].abbrev, "out");

  get_parameter
  (
    (void *) &score->contact.flag,
    &parm, Boolean, "contact_score",
    "no",
    score->flag &&
      (!score->grid.flag || (score->grid.version >= 4.0))
  );

  score->type[CONTACT].flag = score->contact.flag;
  strcpy (score->type[CONTACT].name, "contact");
  strcpy (score->type[CONTACT].abbrev, "cnt");

  get_parameter
  (
    (void *) &score->contact.distance,
    &parm, Real, "contact_cutoff_distance",
    score->contact.flag && (!score->grid.flag || score->intra_flag) ?
      "4.5" : "0.0",
    score->contact.flag && (!score->grid.flag || score->intra_flag)
  );

  get_parameter
  (
    (void *) &score->contact.clash_overlap,
    &parm, Real, "contact_clash_overlap",
    score->contact.flag && (score->intra_flag || !score->grid.flag)
      ? "0.75" : "0.0",
    score->contact.flag && (score->intra_flag || !score->grid.flag)
  );

  get_parameter
  (
    (void *) &score->contact.clash_penalty,
    &parm, Real, "contact_clash_penalty",
    score->contact.flag ? "50" : "0",
    score->contact.flag
  );

  get_parameter
  (
    (void *) &score->chemical.flag,
    &parm, Boolean, "chemical_score",
    "no",
    score->flag &&
      (!score->grid.flag || (score->grid.version >= 4.0))
  );

  label->chemical.flag =
    label->chemical.flag || score->chemical.flag;

  score->type[CHEMICAL].flag = score->chemical.flag;
  strcpy (score->type[CHEMICAL].name, "chemical");
  strcpy (score->type[CHEMICAL].abbrev, "chm");

  get_parameter
  (
    (void *) &score->energy.flag,
    &parm, Boolean, "energy_score",
    "no",
    score->flag
  );

  score->type[ENERGY].flag = score->energy.flag;
  strcpy (score->type[ENERGY].name, "energy");
  strcpy (score->type[ENERGY].abbrev, "nrg");

  get_parameter
  (
    (void *) &score->energy.distance,
    &parm, Real, "energy_cutoff_distance",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag) ? "10.0" : "0.0",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag)
  );

  get_parameter
  (
    (void *) &score->energy.distance_dielectric,
    &parm, Boolean, "distance_dielectric",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag) ? "yes" : "no",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag)
  );

  get_parameter
  (
    (void *) &score->energy.dielectric_factor,
    &parm, Real, "dielectric_factor",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag) ? "4.0" : "0.0",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag)
  );

  get_parameter
  (
    (void *) &score->energy.attractive_exponent,
    &parm, Integer, "attractive_exponent",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag) ? "6" : "0",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag)
  );

  get_parameter
  (
    (void *) &score->energy.repulsive_exponent,
    &parm, Integer, "repulsive_exponent",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag) ? "12" : "0",
    (!score->grid.flag || score->intra_flag) &&
      (score->energy.flag || score->chemical.flag)
  );

  get_parameter
  (
    (void *) &score->energy.atom_model,
    &parm, Character, "atom_model",
    score->energy.flag || score->chemical.flag ? "UNITED all" : "0",
    score->energy.flag || score->chemical.flag
  );

  get_parameter
  (
    (void *) &score->energy.scale_vdw,
    &parm, Real, "vdw_scale",
    score->energy.flag || score->chemical.flag ? "1" : "0",
    score->energy.flag || score->chemical.flag
  );

  get_parameter
  (
    (void *) &score->energy.scale_electro,
    &parm, Real, "electrostatic_scale",
    score->energy.flag || score->chemical.flag ? "1" : "0",
    score->energy.flag || score->chemical.flag
  );

  get_parameter
  (
    (void *) &score->energy.decomp_flag,
    &parm, Boolean, "output_atom_scores",
    "no",
    FALSE
  );

  get_parameter
  (
    (void *) &score->rmsd.flag,
    &parm, Boolean, "rmsd_score",
    "no",
    score->inter_flag && !score->intra_flag &&
      !score->grid.flag && !orient->flag
  );

  score->type[RMSD].flag = score->rmsd.flag;
  strcpy (score->type[RMSD].name, "rmsd");
  strcpy (score->type[RMSD].abbrev, "rmsd");

  if
  (
    score->inter_flag &&
    !score->contact.flag &&
    !score->chemical.flag &&
    !score->energy.flag &&
    !score->rmsd.flag
  )
    exit (fprintf (global.outfile,
      "ERROR get_parameters: no intermolecular scoring options selected.\n"));

  for (i = 1; i < SCORE_TOTAL; i++)
  {
    sprintf (parameter_name, "%s_maximum", score->type[i].name);

    get_parameter
    (
      (void *) &score->type[i].maximum,
      &parm, Real, parameter_name,
      "0",
      score->type[i].flag &&
        ((dock->write_orients && !dock->rank_orients) ||
        (!dock->multiple_orients &&
          dock->multiple_ligands && !dock->rank_ligands))
    );

    sprintf (parameter_name, "%s_size_penalty", score->type[i].name);

    get_parameter
    (
      (void *) &score->type[i].size_penalty,
      &parm, Real, parameter_name,
      "0",
      score->type[i].flag && dock->rank_ligands && (i != RMSD)
    );
  }

  get_parameter
  (
    (void *) &score->rmsd_override,
    &parm, Real, "rmsd_override",
    "0.0",
    score->flag && !dock->rank_ligands &&
      dock->write_orients && !dock->rank_orients
  );


  if (score->minimize.flag || (global.output_volume == 'v'))
    fprintf (global.outfile, 
      "\n______________Minimization_Parameters_____________\n");

  for (i = 1, j = 0; i < SCORE_TOTAL; i++)
  {
    sprintf (parameter_name, "%s_minimize", score->type[i].name);

    get_parameter
    (
      (void *) &score->type[i].minimize,
      &parm, Boolean, parameter_name,
      "no",
      score->type[i].flag && score->minimize.flag
    );

    if (score->type[i].minimize == TRUE) j++;
  }

  if
  (
    (score->minimize.flag == TRUE) &&
    ((j == 0) ||
    ((score->inter_flag == FALSE) &&
      (label->flex.minimize_flag == FALSE)))
  )
    exit (fprintf (global.outfile,
      "ERROR get_parameters: no minimization selected\n"));

  get_parameter
  (
    (void *) &score->minimize.translation,
    &parm, Real, "initial_translation",
    score->inter_flag && score->minimize.flag ? "1.0" : "0.0",
    score->inter_flag && score->minimize.flag
  );

  get_parameter
  (
    (void *) &score->minimize.rotation,
    &parm, Real, "initial_rotation",
    score->inter_flag && score->minimize.flag ? "0.1" : "0.0",
    score->inter_flag && score->minimize.flag
  );

  get_parameter
  (
    (void *) &score->minimize.torsion,
    &parm, Real, "initial_torsion",
    label->flex.minimize_flag ? "10.0" : "0.0",
    label->flex.minimize_flag
  );

  get_parameter
  (
    (void *) &score->minimize.iteration,
    &parm, Integer, "maximum_iterations",
    score->minimize.flag ? "100" : "0",
    score->minimize.flag
  );

  if ((score->minimize.flag == TRUE) && (score->minimize.iteration < 1))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: maximum_iterations < 1\n"));

  for (i = 1; i < SCORE_TOTAL; i++)
  {
    sprintf (parameter_name, "%s_convergence", score->type[i].name);

    get_parameter
    (
      (void *) &score->type[i].convergence,
      &parm, Real, parameter_name,
      score->type[i].minimize && (score->minimize.iteration > 1) ? "0.1" : "0",
      score->type[i].minimize && (score->minimize.iteration > 1)
    );
  }

  get_parameter
  (
    (void *) &score->minimize.cycle,
    &parm, Integer, "maximum_cycles",
    score->minimize.flag ? "1" : "0",
    score->minimize.flag
  );

  if ((score->minimize.flag == TRUE) && (score->minimize.cycle < 1))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: maximum_cycles < 1\n"));

  get_parameter
  (
    (void *) &score->minimize.cycle_converge,
    &parm, Real, "cycle_convergence",
    (score->minimize.cycle > 1) ? "1.0" : "0",
    score->minimize.cycle > 1
  );

  if ((score->minimize.cycle > 1) && (score->minimize.cycle_converge <= 0))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: cycle_convergence <= 0\n"));

  for (i = 1; i < SCORE_TOTAL; i++)
  {
    sprintf (parameter_name, "%s_termination", score->type[i].name);

    get_parameter
    (
      (void *) &score->type[i].termination,
      &parm, Real, parameter_name,
      score->type[i].minimize && (score->minimize.cycle > 1) ? "1.0" : "0",
      score->type[i].minimize && (score->minimize.cycle > 1)
    );
  }


  if (label->chemical.screen.flag || (global.output_volume == 'v'))
    fprintf (global.outfile,
      "\n____________Chemical_Screen_Parameters____________\n");

  get_parameter
  (
    (void *) &label->chemical.screen.construct_flag,
    &parm, Boolean, "construct_screen",
    "no",
    label->chemical.screen.flag
  );

  get_parameter
  (
    (void *) &label->chemical.screen.process_flag,
    &parm, Boolean, "screen_ligands",
    label->chemical.screen.flag ?
      (label->chemical.screen.construct_flag ? "no" : "yes") : "no",
    label->chemical.screen.flag && !label->chemical.screen.construct_flag
  );

  if (label->chemical.screen.flag &&
    !label->chemical.screen.process_flag &&
    !label->chemical.screen.construct_flag)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: No chemical screen options selected!\n"));

  get_parameter
  (
    (void *) &label->chemical.screen.pharmaco_flag,
    &parm, Boolean, "pharmacophore_screen",
    "no",
    label->chemical.screen.process_flag
  );

  if (label->chemical.screen.pharmaco_flag)
    orient->match.chemical_flag = TRUE;

  get_parameter
  (
    (void *) &label->chemical.screen.similar_flag,
    &parm, Boolean, "similarity_screen",
    label->chemical.screen.process_flag ?
      (label->chemical.screen.pharmaco_flag ? "no" : "yes") : "no",
    label->chemical.screen.process_flag &&
      !label->chemical.screen.pharmaco_flag
  );

  get_parameter
  (
    (void *) &label->chemical.screen.fold_flag,
    &parm, Boolean, "fold_keys",
    "yes",
    FALSE
  );

  if (label->chemical.screen.process_flag &&
    !label->chemical.screen.pharmaco_flag &&
    !label->chemical.screen.similar_flag)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: If screen_ligands selected, then either\n"
      "  pharmacophore_screen or similarity_screen must also be selected!\n"));

  get_parameter
  (
    (void *) &label->chemical.screen.dissimilar_maximum,
    &parm, Real, "dissimilarity_maximum",
    "0.25",
    label->chemical.screen.similar_flag
  );

  get_parameter
  (
    (void *) &label->chemical.screen.distance_minimum,
    &parm, Real, "distance_begin",
    label->chemical.screen.flag ? "2" : "0",
    label->chemical.screen.flag
  );

  get_parameter
  (
    (void *) &label->chemical.screen.distance_maximum,
    &parm, Real, "distance_end",
    label->chemical.screen.flag ? "17" : "0",
    label->chemical.screen.flag
  );

  get_parameter
  (
    (void *) &label->chemical.screen.distance_interval,
    &parm, Real, "distance_interval",
    label->chemical.screen.flag ? "0.5" : "0",
    label->chemical.screen.flag
  );

  label->chemical.screen.interval_total = (int)
    ((label->chemical.screen.distance_maximum -
    label->chemical.screen.distance_minimum) /
    label->chemical.screen.distance_interval) + 2;

  if (label->chemical.screen.interval_total > MASK_LENGTH)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: Chemical screen intervals > %d!\n",
      MASK_LENGTH - 2));


  if (dock->parallel.flag || (global.output_volume == 'v'))
    fprintf (global.outfile,
      "\n______________Parallel_Job_Parameters_____________\n");

  get_parameter
  (
    (void *) &dock->parallel.server,
    &parm, Boolean, "parallel_server",
    "no",
    dock->parallel.flag
  );

/*
* Make sure only a single type of scoring has been requested
* 10/96 te
*/
  if (dock->parallel.server)
  {
    for (i = j = 0; i < SCORE_TOTAL; i++)
      if (score->type[i].flag) j++;

    if (j != 1)
      exit (fprintf (global.outfile,
        "ERROR get_parameters: "
        "Parallel job server must perform single type of scoring.\n"));
  }

  get_parameter
  (
    (void *) &dock->parallel.server_name,
    &parm, String, "server_name",
    dock->parallel.server ? global.job_name : "server",
    dock->parallel.flag
  );

  get_parameter
  (
    (void *) &dock->parallel.client_total,
    &parm, Integer, "client_total",
    dock->parallel.server ? "5" : (dock->parallel.flag ? "1" : "0"),
    dock->parallel.server
  );

  if (dock->parallel.flag)
  {
    if (dock->parallel.client_total < 1)
      exit (fprintf (global.outfile,
        "ERROR get_parameters: Inappropriate value for client_total.\n"));

    ecalloc
    (
      (void **) &dock->parallel.client_name,
      dock->parallel.client_total,
      sizeof (STRING20),
      "parallel client names",
      global.outfile
    );

    for (i = 0; i < dock->parallel.client_total; i++)                 
    {
      if (dock->parallel.server)
      {
        sprintf (parameter_name, "client_name_%d", i + 1);
        sprintf (parameter_value, "client%d", i + 1);
      }

      else
      {
        strcpy (parameter_name, "client_name");
        strcpy (parameter_value, global.job_name);
      }

      get_parameter
      (
        (void *) &dock->parallel.client_name[i],
        &parm, String, parameter_name,
        parameter_value,
        TRUE
      );
    }
  }


  if (dock->multiple_ligands || (global.output_volume == 'v'))
    fprintf (global.outfile,
      "\n____________Multiple_Ligand_Parameters____________\n");
                                                                      
  get_parameter
  (
    (void *) &dock->max_ligands,
    &parm, Integer, "ligands_maximum",
    dock->multiple_ligands ? "1000" : "1",
    dock->multiple_ligands
  );

  get_parameter
  (
    (void *) &dock->initial_skip,
    &parm, Integer, "initial_skip",
    "0",
    dock->multiple_ligands
  );

  get_parameter
  (
    (void *) &dock->interval_skip,
    &parm, Integer, "interval_skip",
    "0",
    dock->multiple_ligands
  );

  get_parameter
  (
    (void *) &dock->min_heavies,
    &parm, Integer, "heavy_atoms_minimum",
    "0",
    dock->multiple_ligands
  );

  get_parameter
  (
    (void *) &dock->max_heavies,
    &parm, Integer, "heavy_atoms_maximum",
    orient->flag && dock->multiple_ligands ? "100" : "<infinity>",
    dock->multiple_ligands
  );

  get_parameter
  (
    (void *) &dock->rank_ligands,
    &parm, Boolean,
    "rank_ligands",
    "no",
    dock->multiple_ligands &&
      (score->flag || label->chemical.screen.similar_flag) &&
      (!dock->write_orients || dock->rank_orients)
  );

  get_parameter
  (
    (void *) &dock->rank_ligand_total,
    &parm, Integer,
    "rank_ligand_total",
    dock->rank_ligands ? "100" : "1",
    dock->rank_ligands
  );

  get_parameter
  (
    (void *) &dock->restart_interval,
    &parm, Integer, "restart_interval",
    dock->rank_ligands ? "100" : "0",
    dock->rank_ligands
  );


  fprintf (global.outfile,
    "\n____________________File_Input____________________\n");

  get_parameter
  (
    (void *) &dock->ligand_file_name,
    &parm, String, "ligand_atom_file",
    dock->multiple_ligands ?
      (label->chemical.screen.process_flag ? "database.ptr" : "database.mol2")
      : "ligand.mol2",
    !dock->parallel.flag || dock->parallel.server
  );

  format = check_file_extension (dock->ligand_file_name, TRUE);

  if (format == Unknown)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "Unrecognized file extension of ligand_atom_file\n"));

  else if (label->chemical.screen.process_flag && (format != Ptr))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "PTR file required to chemical screen ligands.\n"));

  get_parameter
  (
    (void *) &orient->match.ligand_file_name,
    &parm, String, "ligand_center_file",
    "ligand_center.sph",
    orient->match.centers_flag
  );

  if (check_file_extension (orient->match.ligand_file_name, TRUE) == Unknown)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "Unrecognized file extension of ligand_center_file\n"));

  get_parameter
  (
    (void *) &orient->match.receptor_file_name,
    &parm, String, "receptor_site_file",
    label->chemical.screen.process_flag ? "target.mol2" : "receptor_site.sph",
    orient->flag || label->chemical.screen.process_flag
  );

  if (check_file_extension (orient->match.receptor_file_name, TRUE) == Unknown)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "Unrecognized file extension of receptor_site_file\n"));

  get_parameter
  (
    (void *) &score->grid.file_prefix,
    &parm, String, "score_grid_prefix",
    "score_grid",
    score->inter_flag && score->grid.flag
  );

  get_parameter
  (
    (void *) &score->grid.receptor_file_name,
    &parm, String, "receptor_atom_file",
    "receptor.mol2",
    score->inter_flag && !score->grid.flag
  );

  if (!check_file_extension (score->grid.receptor_file_name, TRUE))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "Unrecognized file extension of receptor_atom_file\n"));

  get_parameter
  (
    (void *) &label->vdw.file_name,
    &parm, String, "vdw_definition_file",
    PARAMETER_PATH "vdw.defn",
    label->vdw.flag
  );

  get_parameter
  (
    (void *) &label->chemical.file_name,
    &parm, String, "chemical_definition_file",
    PARAMETER_PATH "chem.defn",
    label->chemical.flag
  );

  get_parameter
  (
    (void *) &label->chemical.match_file_name,
    &parm, String, "chemical_match_file",
    PARAMETER_PATH "chem_match.tbl",
    orient->match.chemical_flag
  );

  get_parameter
  (
    (void *) &label->chemical.score_file_name,
    &parm, String, "chemical_score_file",
    PARAMETER_PATH "chem_score.tbl",
    score->chemical.flag
  );

  get_parameter
  (
    (void *) &label->chemical.screen_file_name,
    &parm, String, "chemical_screen_file",
    PARAMETER_PATH "chem_screen.tbl",
    label->chemical.screen.similar_flag
  );

  get_parameter
  (
    (void *) &label->flex.file_name,
    &parm, String, "flex_definition_file",
    PARAMETER_PATH "flex.defn",
    label->flex.flag
  );

  get_parameter
  (
    (void *) &label->flex.search_file_name,
    &parm, String, "flex_drive_file",
    PARAMETER_PATH "flex_drive.tbl",
    label->flex.drive_flag
  );

  strcat (strcpy (parameter_value, global.job_name), ".quit");

  get_parameter
  (
    (void *) &dock->quit_file_name,
    &parm, String, "quit_file",
    parameter_value,
    dock->rank_ligands || (dock->multiple_ligands && dock->parallel.flag)
  );

  strcat (strcpy (parameter_value, global.job_name), ".dump");

  get_parameter
  (
    (void *) &dock->dump_file_name,
    &parm, String, "dump_file",
    parameter_value,
    dock->rank_ligands
  );


  fprintf (global.outfile,
    "\n____________________File_Output___________________\n");

  for (i = 0; i < SCORE_TOTAL; i++)
  {
    sprintf (parameter_name, "ligand_%s_file", i ? score->type[i].name : "out");
    sprintf (parameter_value, "%s_%s.%s", global.job_name,
      score->type[i].abbrev,
      label->chemical.screen.flag || dock->parallel.flag ? "ptr" : "mol2");

    get_parameter
    (
      (void *) &score->type[i].file_name,
      &parm, String, parameter_name,
      parameter_value,
      score->type[i].flag && !dock->parallel.server
    );

    if (score->type[i].flag)
    {
      if (dock->parallel.flag && !dock->parallel.server)
      {
        if (check_file_extension (score->type[i].file_name, TRUE) != Ptr)
          exit (fprintf (global.outfile,
            "ERROR get_parameters: "
            "Parallel client must write in PTR format\n"));
      }

      else if (label->chemical.screen.construct_flag)
      {
        if (check_file_extension (score->type[i].file_name, TRUE) != Ptr)
          exit (fprintf (global.outfile,
            "ERROR get_parameters: "
            "Constructed screen must be written in PTR format\n"));
      }

      else
      {
        if (check_file_extension (score->type[i].file_name, TRUE) == Unknown)
          exit (fprintf (global.outfile,
            "ERROR get_parameters: "
            "Unrecognized file extension of ligand output file\n"));
      }
    }
  }

  strcat (strcpy (parameter_value, global.job_name), ".info");

  get_parameter
  (
    (void *) &dock->info_file_name,
    &parm, String, "info_file",
    parameter_value,
    score->flag && dock->multiple_ligands
  );

  strcat (strcpy (parameter_value, global.job_name), ".rst");

  get_parameter
  (
    (void *) &dock->restart_file_name,
    &parm, String, "restart_file",
    parameter_value,
    dock->rank_ligands
  );

/*
* Flush output and clean up memory
* 11/96 te
*/
  fprintf (global.outfile, "\n\n");
  fflush (global.outfile);

  if (global.infile != stdin)
  {
    fclose (global.infile);
    efree ((void **) &parm.line);
  }
}

/* ////////////////////////////////////////////////////////////// */
int process_commands
(
  DOCK *dock,
  int argc,
  char *argv[]
)
{
  int i = 1, j;
  int molecule_io = FALSE;
  FILE_NAME infile_name = "", outfile_name = "";

/*
* Extract name of global.executable from command line,
* ignoring any path description if present
* 6/95 te
*/
  strcpy (global.executable,
    strrchr (argv[0], '/') ? strrchr (argv[0], '/') + 1 : argv[0]);
  strcpy (global.job_name, "dock");

/*
* Step through each command line argument
* 6/95 te
*/
  while (i < argc)
  {

/*
*   Look for command flag prefix (-)
*   6/95 te
*/
    if (argv[i][0] == '-')
    {

/*
*     For "-i" flag read in next field as the input file name
*     6/95 te
*/
      if (argv[i][1] == 'i')
      {
        if ((i + 1 < argc) && (argv[i + 1][0] != '-'))
        {
          strncpy (infile_name, argv[i + 1], 80);
          i++;

/*
*         Set the job name to the root name of the the "*.in" file
*         10/95 te
*/
          if (strrchr (infile_name, '.') &&
            !strcmp (strrchr (infile_name, '.'), ".in"))
          {
            for (j = 0; j < strlen (infile_name); j++)
            {
              if (!strcmp (&infile_name[j], ".in"))
              {
                global.job_name[j] = 0;
                break;
              }

              else
                global.job_name[j] = infile_name[j];
            }
          }

          else
            strcpy (global.job_name, infile_name);
        }

        else
          sprintf (infile_name, "%s.in", global.job_name);
      }

/*
*     For "-o" flag read in next field as the output file name.
*     If the next field does not exist or is another flag, then
*     use default output file name.
*     6/95 te
*/
      else if (argv[i][1] == 'o')
      {
        if ((i + 1 < argc) && (argv[i + 1][0] != '-'))
        {
          strncpy (outfile_name, argv[i + 1], 80);
          i++;
        }

        else
          sprintf (outfile_name, "%s.out", global.job_name);

        if (!strcmp (infile_name, ""))
          break;
      }

/*
*     For "-p" flag, set performance flag
*     1/97 te
*/
      else if (argv[i][1] == 'p')
        dock->performance_flag = TRUE;

/*
*     For "-r" flag, turn on flag to perform restart run
*     6/95 te
*/
      else if (argv[i][1] == 'r')
      {
        if (strcmp (infile_name, ""))
          dock->restart = TRUE;

        else
          break;
      }

/*
*     For "-s" flag, allow all input and output through standard streams
*     6/95 te
*/
      else if (argv[i][1] == 's')
      {
        if ((strcmp (infile_name, "")) || (strcmp (outfile_name, "")))
          break;
      }

/*
*     For "-t" flag, set output volume to terse
*     6/95 te
*/
      else if (argv[i][1] == 't')
        global.output_volume = 't';

/*
*     For "-v" flag, set output volume to verbose
*     6/95 te
*/
      else if (argv[i][1] == 'v')
        global.output_volume = 'v';

      else
        break;

      i++;
    }

    else
      break;
  }

/*
* Check if no arguments were given, to conform to old style input/output
* 10/95 te
*/
  if (argc == 1)
  {
/*
*   If INDOCK file exists, then assume user wants output directed to OUTDOCK
*   10/95 te
*/
    if (global.infile = fopen ("INDOCK", "r"))
    {
      fclose (global.infile);

      if (global.infile = fopen ("OUTDOCK", "r"))
      {
        fclose (global.infile);
        fprintf (stderr,
          "ERROR get_parameters: \n"
          "       Presence of INDOCK file implies OUTDOCK to receive output,\n"
          "       but OUTDOCK file already exists.  For reverse compatibility\n"
          "       reasons, this process will stop.\n");
        exit (EXIT_FAILURE);
      }

      strcpy (infile_name, "INDOCK");
      strcpy (outfile_name, "OUTDOCK");
    }

/*
*   If it doesn't, then trigger an error condition
*   10/95 te
*/
    else
      i = 0;
  }

/*
* If any error was made in the command line, then quit
* 6/95 te
*/
  if (i < argc)
  {
    fprintf (stderr,
      "  \n"
      "[33mUsage: %s [-i [input_file]] [-o [output_file]] ...\n"
      "  [-restart] [-standard_i/o] [-terse] [-verbose][0m\n"
      "  \n"
      "  -i: read from %s.in or input_file, standard_in otherwise\n"
      "  -o: write to %s.out or output_file (-i required), \n"
      "      standard_out otherwise\n"
      "  -p: monitor performance\n"
      "  -r: restart run (-i required)\n"
      "  -s: read from and write to standard streams (-i and/or -o illegal)\n"
      "  -t: terse program output\n"
      "  -v: verbose program output\n" 
      "  \n"
      "  \n"
      "To create or update an input file interactively just type:\n"
      "  %s -i your_name.in\n"
      "  \n",
      global.executable, global.executable, global.executable,
      global.executable);
    exit (EXIT_FAILURE);
  }

/*
* Open the input file
* 6/95 te
*/
  if (strcmp (infile_name, ""))
  {
    if (strcmp (outfile_name, ""))
    {
/*
*     Check to see if user has supplied molecule files as input/output
*     6/95 te
*/
      if ((check_file_extension (infile_name, FALSE) != Unknown) &&
        (check_file_extension (outfile_name, FALSE) != Unknown))
      {
        molecule_io = TRUE;

/*
*       Open up scratch file and write necessary parameters to it
*       6/95 te
*/
        global.infile = tmpfile ();

        fprintf (global.infile, "orient_ligand             no\n");
        fprintf (global.infile, "flexible_ligand           no\n");
        fprintf (global.infile, "score_ligand              no\n");
        fprintf (global.infile, "multiple_ligands          yes\n");
        fprintf (global.infile, "chemical_screen           no\n");
        fprintf (global.infile, "parallel_jobs             no\n");
        fprintf (global.infile, "ligands_maximum           <infinity>\n");
        fprintf (global.infile, "initial_skip              0\n");
        fprintf (global.infile, "interval_skip             0\n");
        fprintf (global.infile, "atoms_maximum             <infinity>\n");
        fprintf (global.infile, "heavy_atoms_minimum       0\n");
        fprintf (global.infile, "heavy_atoms_maximum       <infinity>\n");
        fprintf (global.infile, "ligand_atom_file          %s\n", infile_name);
        fprintf (global.infile, "vdw_definition_file       %s\n",
          PARAMETER_PATH "vdw.defn");
        fprintf (global.infile, "ligand_out_file          %s\n", outfile_name);

        rewind (global.infile);
      }

      else
        global.infile = efopen (infile_name, "r", stdout);
    }

    else
      global.infile = efopen (infile_name, "a+", stdout);
  }

  else
    global.infile = stdin;

/*
* Open the output file
* 6/95 te
*/
  if ((strcmp (outfile_name, "")) && (!molecule_io))
  {
    if (dock->restart)
      global.outfile = efopen (outfile_name, "a", stdout);

    else
      global.outfile = efopen (outfile_name, "w", stdout);
  }

  else
    global.outfile = stdout;

  return TRUE;
}

