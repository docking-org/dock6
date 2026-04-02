/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
12/96
*/
#include "define.h"
#include "utility.h"
#include "mol.h"
#include "global.h"
#include "dock.h"
#include "search.h"
#include "label.h"
#include "score.h"
#include "score_dock.h"
#include "flex.h"
#include "simplex.h"
#include "io_ligand.h"
#include "io_grid.h"
#include "transform.h"
#include "rank.h"
#include "vector.h"

typedef struct simplex_struct
{
  DOCK		*dock;
  LABEL		*label;
  SCORE		*score;
  MOLECULE	*mol_ref;
  MOLECULE	*mol_ori;
  MOLECULE	*mol_score;
  MOLECULE	mol_min;
  MOLECULE	mol_best;
  LIST		*list;
  int		rigid_flag;
  int		layer_inner;
  int		layer_outer;
  int		cycle;

} SIMPLEX;


/* ///////////////////////////////////////////////////////////////

Subroutine to calculate and optimize the score for a ligand orientation.
The best scoring orientation is updated.

Return values:
	TRUE	molecule orientation was scored
	FALSE	molecule orientation not scored

1/97 te

/////////////////////////////////////////////////////////////// */

int get_anchor_score
(
  DOCK		*dock,
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*mol_ref,
  MOLECULE	*mol_init,
  MOLECULE	*mol_ori,
  MOLECULE	*mol_score,
  LIST		*best_anchors,
  int		orient_id,
  int		write_flag
)
{
  int type;			/* Score type iterator */

/*
* Read/construct scoring grids
* 11/96 te
*/
  if (!orient_id)
  {
    if (score->inter_flag && !score->grid.init_flag)
    {
      if (score->grid.flag)
        read_grids
        (
          &score->grid,
          &score->bump,
          &score->contact,
          &score->chemical,
          &score->energy,
          &label->chemical
        );

      else
        make_receptor_grid
        (
          &score->grid,
          &score->energy,
          label
        );
    }
  }

  else
    reset_score (&mol_ori->score);

/*
* Evaluate whether anchor orientation bumps with receptor
* 12/96 te
*/
  if
  (
    (score->inter_flag == TRUE) &&
    (score->bump.flag == TRUE) &&
    (check_bump (&score->grid, &score->bump, label, mol_ori) >
      score->bump.maximum)
  )
    return FALSE;

  if (mol_ori->transform.flag)
    mol_ori->transform.rmsd =
      calc_rmsd (mol_ori, mol_init);

  copy_transform (mol_score, mol_ori);
  copy_score (&mol_score->score, &mol_ori->score);
  copy_coords (mol_score, mol_ori);
  copy_torsions (mol_score, mol_ori);
  copy_segments (mol_score, mol_ori);
  copy_layers (mol_score, mol_ori);

/*
* Loop through all scoring types
* 11/96 te
*/
  for (type = 0; type < SCORE_TOTAL; type++)
  {
    if (!score->type[type].flag)
      continue;

    mol_score->score.type = type;

/*
*   Either minimize this orientation or just score it
*   11/96 te
*/
    if (score->type[type].minimize && label->flex.minimize_anchor_flag)
    {
      minimize_ligand
        (dock, label, score, mol_ref, mol_ori, mol_score, TRUE, 0, 0);
      mol_score->transform.rmsd =
        calc_rmsd (mol_score, mol_init);
    }

    else
      calc_score (label, score, mol_score, 0);

/*
*   Update the list of best anchors
*   11/96 te
*/
    update_list (best_anchors, mol_score);

/*
*   Write out the orientation if flagged, and if:
*
*   1. Ligand has moved AND rmsd is within override, OR
*   2. No scoring is performed, OR
*   3. The score is below maximum cutoff
*
*   3/96 te
*/
    if
    (
      write_flag &&
      (
        ((mol_score->transform.flag) &&
          (mol_score->transform.rmsd <=
            score->rmsd_override)) ||
        (mol_score->score.total <= score->type[type].maximum)
      )
    )
    {
      write_ligand
      (
        dock,
        score,
        mol_score,
        score->type[type].file_name,
        score->type[type].file
      );
    }

/*
*   Reset coordinates previous to minimization
*   6/96 te
*/
    if (score->type[type].minimize)
    {
      copy_transform (mol_score, mol_ori);
      copy_torsions (mol_score, mol_ori);
      copy_coords (mol_score, mol_ori);
    }
  }

  return TRUE;
}


/* /////////////////////////////////////////////////////////////// */

void get_peripheral_score
(
  DOCK		*dock,
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*mol_ref,
  MOLECULE	*mol_init,
  MOLECULE	*mol_conf,
  MOLECULE	*mol_ori,
  MOLECULE	*mol_score,
  LIST		*best_anchors,
  LIST		*best_orients,
  int		write_flag,
  int		anchor
)
{
  int type;			/* Score type iterator */
  int seed;			/* Seed conformation iterator */
  int layer;			/* Segment layer iterator */
  int segment_id;		/* Segment iterator */
  int segment;			/* Segment id */
  static int allocate_flag = TRUE; /* Flag for whether lists allocated */
  static LIST seed_conform;	/* List of seed conformations */
  static LIST best_conform;	/* List of best conformations */

/*
* Initialize conformation storage lists
* 12/96 te
*/
  if (allocate_flag == TRUE)
  {
    allocate_lists (score, &seed_conform, label->flex.max_conforms, FALSE);
    allocate_flag = FALSE;
  }

  else
    reset_lists (score, &seed_conform);

  copy_lists (score, &best_conform, best_anchors);

/*
* Loop through all scoring types
* 12/96 te
*/
  for (type = 0; type < SCORE_TOTAL; type++)
  {
    if (!score->type[type].flag)
      continue;

/*
*   Loop through each shell of flexible segments
*   3/97 te
*/
    for
    (
      layer = 0;
      (layer < mol_conf->total.layers) &&
        (best_conform.total[type] > 0);
      layer++
    )
    {
/*
*     Loop through all segments in current layer
*     Except first layer, just do first segment and only do pruning of it
*     2/98 te
*/
      for
      (
        segment_id = 0;
        ((layer == 0) && (segment_id < 1)) ||
        ((layer != 0) && (segment_id < mol_conf->layer[layer].segment_total))
          && (best_conform.total[type] > 0);
        segment_id++
      )
      {
        if (layer > 0)
        {
/*
*         Allocate sufficient space for best conformations
*         1/98 te
*/
          segment = mol_conf->layer[layer].segment[segment_id];

          best_conform.total[type] =
            seed_conform.total[type] * mol_conf->segment[segment].conform_total;
          reallocate_list (&best_conform, type);
          reset_list (&best_conform, type);

/*
*         Loop through all seed conformations
*         3/97 te
*/
          for (seed = 0; seed < seed_conform.total[type]; seed++)
          {
              copy_member (FALSE, mol_ori, seed_conform.member[type][seed]);
              mol_ori->layer[layer].active_flag = TRUE;
              mol_ori->segment[segment].active_flag = TRUE;
              mol_ori->segment[segment].min_flag = TRUE;
              mol_ori->segment[segment].conform_count = 0;

/*
*             Drive through conformations of this segment
*             3/97 te
*/
              while
              (
                get_segment_conformation
                (
                  label,
                  score,
                  mol_ori,
                  layer,
                  segment_id
                ) != EOF
              )
              {
                copy_molecule (mol_score, mol_ori);

/*
*               Either minimize this conformation or just score it
*               If minimize, allow rigid+outer layer OR outer two layers
*               2/97 te
*/
                if (score->type[type].minimize)
                  minimize_ligand
                    (dock, label, score, mol_ref, mol_ori, mol_score,
                    label->flex.reminimize_anchor_flag,
                    MAX (0, layer - label->flex.reminimize_layers),
                    layer);

                else
                  calc_score
                    (label, score, mol_score, layer);

                if (global.output_volume == 'v')
                  mol_score->transform.rmsd =
                    calc_rmsd (mol_score, mol_init);

                mol_score->segment[segment].min_flag = FALSE;
                update_list (&best_conform, mol_score);

              } /* End of conformer loop */
          } /* End of seed loop */
        } /* End of if (layer > 0) */

        if (global.output_volume == 'v')
        {
          fprintf (global.outfile,
            "\nAnchor-first search: anchor %d, layer %d, segment %d, "
            "configurations %d\n", anchor + 1,
            layer + 1, segment_id + 1, best_conform.total[type]);

          print_list (score, &best_conform, type, global.outfile);
          fflush (global.outfile);
        }

/*
*       Write out the conformationally expanded ensemble
*       2/98 te
*/
        if (label->flex.write_flag == TRUE)
          write_periph_structures
          (
            dock,
            score,
            &best_conform,
            mol_ref,
            mol_ori,
            type,
            anchor,
            layer,
            segment_id,
            FALSE
          );

/*
*       Identify the most different best conforms as seeds for next cycle
*       2/98 te
*/
        shrink_list
        (
          &best_conform,
          &seed_conform,
          type,
          label->flex.max_conforms
        );

        if (global.output_volume == 'v')
        {
          fprintf (global.outfile,
            "\nAnchor-first search: anchor %d, layer %d, segment %d, "
            "prunings %d\n", anchor + 1,
            layer + 1, segment_id + 1, seed_conform.total[type]);

          print_list (score, &seed_conform, type, global.outfile);
          fflush (global.outfile);
        }

/*
*       Write out the pruned ensemble
*       2/98 te
*/
        if (label->flex.write_flag == TRUE)
          write_periph_structures
          (
            dock,
            score,
            &seed_conform,
            mol_ref,
            mol_ori,
            type,
            anchor,
            layer,
            segment_id,
            TRUE
          );

      } /* End of segment count loop */
    } /* End of layer loop */

/*
*   Process the top scoring configurations
*   2/97 te
*/
    if (seed_conform.total[type] > 0)
    {
/*
*     Re-minimize the top scoring configurations
*     2/97 te
*/
      if
      (
        score->type[type].minimize &&
        label->flex.reminimize_ligand_flag &&
        (mol_ori->total.layers > 2)
      )
      {
        for (seed = 0; seed < seed_conform.total[type]; seed++)
        {
          copy_molecule (mol_score, seed_conform.member[type][seed]);

          for
          (
            segment_id = 0;
            segment_id < mol_score->total.segments;
            segment_id++
          )
            mol_score->segment[segment_id].min_flag = TRUE;

          minimize_ligand
            (dock, label, score, mol_ref, seed_conform.member[type][seed],
            mol_score, TRUE, 0, mol_score->total.layers - 1);
          mol_score->transform.rmsd = calc_rmsd (mol_score, mol_init);

          copy_molecule (seed_conform.member[type][seed], mol_score);
        }

        sort_list (&seed_conform, type);

        if (global.output_volume == 'v')
        {
          fprintf (global.outfile,
            "Anchor-first search: anchor %d, reminimize_ligand results\n",
            anchor + 1);
          print_list (score, &seed_conform, type, global.outfile);
        }
      }

      merge_lists (score, best_orients, &seed_conform);

/*
*     Write out the best configurations to file, if requested
*     2/97 te
*/
      if (write_flag)
      {
        for (seed = 0; seed < seed_conform.total[type]; seed++)
        {
          if
          (
            ((seed_conform.member[type][seed]->transform.flag) &&
              (seed_conform.member[type][seed]->transform.rmsd <=
                score->rmsd_override)) ||
            (seed_conform.member[type][seed]->score.total <=
              score->type[type].maximum)
          )
          {
            write_ligand
            (
              dock,
              score,
              seed_conform.member[type][seed],
              score->type[type].file_name,
              score->type[type].file
            );
          }
        } /* End of seed loop */
      } /* End of write IF */
    } /* End of process IF */
  } /* End of score type loop */
}


/* ////////////////////////////////////////////////////////////////// */

void write_periph_structures
(
  DOCK		*dock,
  SCORE		*score,
  LIST		*list,
  MOLECULE	*mol_ref,
  MOLECULE	*mol_ori,
  int		type,
  int		anchor,
  int		layer,
  int		segment,
  int		shrink_flag
)
{
  FILE_NAME original_name;	/* File name of output file */
  FILE *original_file;		/* File pointer of output file */
  int number_written;		/* Number of structures written to file */

/*
* Store, then edit name of output file
* 2/98 te
*/
  strcpy (original_name, score->type[type].file_name);
  original_file = score->type[type].file;
  number_written = score->type[type].number_written;
  score->type[type].number_written = 0;

  sprintf (strrchr (score->type[type].file_name, '.') , "-%d-%02d-%d-%d%s",
    anchor + 1, layer + 1, segment + 1, shrink_flag + 1,
    strrchr (original_name, '.'));
  score->type[type].file =
    efopen (score->type[type].file_name, "w", global.outfile);

  write_topscorers
  (
    dock,
    score,
    list,
    mol_ref,
    mol_ori
  );

  score->type[type].number_written = number_written;
  efclose (&score->type[type].file);
  strcpy (score->type[type].file_name, original_name);
  score->type[type].file = original_file;
}



/* ////////////////////////////////////////////////////////////////// */

void minimize_ligand
(
  DOCK		*dock,
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*mol_ref,
  MOLECULE	*mol_ori,
  MOLECULE	*mol_score,
  int		rigid_flag,
  int		layer_inner,
  int		layer_outer
)
{
  int i, j;			/* Iteration variable */
  int torsion;			/* Torsion id */
  int segment;			/* Segment id */
  int layer;			/* Layer id */
  int neighbor;			/* Segment neighbor id */
  int flex;			/* Flex label id */
  SIMPLEX simplex = {0};	/* Simplex data */
  int iteration_total;		/* Number of simplex iterations */
  int vertex_total = 0;		/* Total simplex vertices */
  float *vertex = NULL;		/* Simplex vertices */
  float delta_cycle;		/* Change in score during a cycle */
  float delta_total;		/* Total change in score */
  float distance;		/* Simplex vertex distance */

  if ((score->intra_flag) && (label->flex.minimize_flag == TRUE))
  {
    if
    (
      (layer_inner < 0) ||
      (layer_inner > layer_outer) ||
      (layer_outer >= mol_score->total.layers)
    )
      exit (fprintf (global.outfile,
        "ERROR minimize_ligand: molecule layer selection inappropriate\n"));

/*
*   Identify inner neighbors of minimizing segments
*   4/97 te
*/
    for (layer = layer_outer; layer >= layer_inner; layer--)
    {
      for (i = 0; i < mol_score->layer[layer].segment_total; i++)
      {
        segment = mol_score->layer[layer].segment[i];

        if (mol_score->segment[segment].min_flag)
        {
          if
          (
            ((torsion = mol_score->segment[segment].torsion_id) != NEITHER) &&
            ((flex = mol_score->torsion[torsion].flex_id) != NEITHER) &&
            (label->flex.member[flex].minimize_flag == TRUE)
          )
          {
            vertex_total++;
            mol_score->transform.tors_flag = TRUE;
          }

          if (layer > layer_inner)
          {
            for (j = 0; j < mol_score->segment[segment].neighbor_total; j++)
            {
              neighbor = mol_score->segment[segment].neighbor[j].id;

              if
              (
                (mol_score->segment[neighbor].layer_id < layer) &&
                (mol_score->segment[neighbor].layer_id >= layer_inner)
              )
                mol_score->segment[neighbor].min_flag = TRUE;
            }
          }
        }
      }
    }

/*
*   If no inner bonds can be minimized, then perform rigid-body minimization
*   4/97 te
*/
    if (vertex_total == 0)
      rigid_flag = TRUE;

/*
*   Identify outer neighbors of minimizing segments (not already identified)
*   4/97 te
*/
    for (layer = layer_inner + 1; layer <= layer_outer; layer++)
    {
      for (i = 0; i < mol_score->layer[layer].segment_total; i++)
      {
        segment = mol_score->layer[layer].segment[i];

        if (mol_score->segment[segment].min_flag != TRUE)
        {
          for (j = 0; j < mol_score->segment[segment].neighbor_total; j++)
          {
            neighbor = mol_score->segment[segment].neighbor[j].id;

            if
            (
              (mol_score->segment[neighbor].layer_id < layer) &&
              (mol_score->segment[neighbor].layer_id >= layer_inner) &&
              (mol_score->segment[neighbor].min_flag == TRUE)
            )
            {
              if
              (
                ((torsion = mol_score->segment[segment].torsion_id)
                  != NEITHER) &&
                ((flex = mol_score->torsion[torsion].flex_id) != NEITHER) &&
                (label->flex.member[flex].minimize_flag == TRUE)
              )
              {
                vertex_total++;
                mol_score->transform.tors_flag = TRUE;
              }

              mol_score->segment[segment].min_flag = TRUE;
              break;
            }
          }
        }
      }
    }
  }

  if ((score->inter_flag) && (rigid_flag == TRUE))
  {
    simplex.rigid_flag = TRUE;
    vertex_total += 6;
    mol_score->transform.trans_flag =
      mol_score->transform.rot_flag = TRUE;
  }

  if (vertex_total > 0)
  {
    ecalloc
    (
      (void **) &vertex,
      vertex_total,
      sizeof (float),
      "simplex vertex array",
      global.outfile
    );

    simplex.dock = dock;
    simplex.score = score;
    simplex.label = label;
    simplex.mol_ref = mol_ref;
    simplex.mol_ori = mol_ori;
    simplex.mol_score = mol_score;
    simplex.layer_inner = layer_inner;
    simplex.layer_outer = layer_outer;

    copy_molecule (&simplex.mol_min, mol_score);
    copy_molecule (&simplex.mol_best, mol_score);
    simplex.mol_best.score.total = INITIAL_SCORE;
    delta_total = 0;

    for
    (
      simplex.cycle = 0, distance = FLT_MAX;
      (simplex.cycle < score->minimize.cycle) &&
        (distance > score->minimize.cycle_converge) &&
        ((simplex.cycle == 0) ||
          (simplex.mol_best.score.total <
          score->type[mol_score->score.type].termination));
      simplex.cycle++
    )
    {
      simplex_optimize
      (
        (void *) &simplex,
        vertex,
        vertex_total,
        score->type[mol_score->score.type].convergence,
        &iteration_total,
        score->minimize.iteration,
        simplex_score,
        &delta_cycle
      );

      copy_molecule (mol_score, &simplex.mol_best);

/*
*     Compute distance that vertex has moved (and update initial vertex)
*     1/97 te
*/
      for (i = distance = 0; i < vertex_total; i++)
        distance += SQR (vertex[i]);

      distance = sqrt (distance) / (float) (simplex.cycle + 1);
      memset (vertex, 0, vertex_total * sizeof (float));

/*
*     Update minimizer statistics
*     1/97 te
*/
      delta_total += delta_cycle;

      score->minimize.iteration_total += iteration_total;
      score->minimize.iteration_max =
        MAX (score->minimize.iteration_max, iteration_total);
      score->minimize.iteration_min =
        MIN (score->minimize.iteration_min, iteration_total);
    }

    score->minimize.call_sub_total++;

    score->minimize.vertex_total += vertex_total;
    score->minimize.vertex_max =
      MAX (score->minimize.vertex_max, vertex_total);
    score->minimize.vertex_min =
      MIN (score->minimize.vertex_min, vertex_total);

    score->minimize.cycle_total += simplex.cycle;
    score->minimize.cycle_max =
      MAX (score->minimize.cycle_max, simplex.cycle);
    score->minimize.cycle_min =
      MIN (score->minimize.cycle_min, simplex.cycle);

    score->minimize.delta_total += delta_total;
    score->minimize.delta_max =
      MAX (score->minimize.delta_max, delta_total);
    score->minimize.delta_min =
      MIN (score->minimize.delta_min, delta_total);

    free_molecule (&simplex.mol_min);
    free_molecule (&simplex.mol_best);
    efree ((void **) &vertex);
  }

  else
    calc_score
    (
      label,
      score,
      mol_score,
      layer_outer
    );

/*
* Turn off all minimize flags
* 3/97 te
*/
  if (label->flex.minimize_flag == TRUE)
    for (layer = layer_inner; layer <= layer_outer; layer++)
      for (i = 0; i < mol_score->layer[layer].segment_total; i++)
      {
        segment = mol_score->layer[layer].segment[i];
        mol_score->segment[segment].min_flag = FALSE;
      }
}

/* ////////////////////////////////////////////////////////////////// */

float simplex_score (void *simplex_input, float *vertex)
{
  int i;
  int flex;			/* Flex label id */
  int torsion;			/* Torsion id */
  int segment;			/* Segment id */
  int layer;			/* Layer id */
  int intra_current_flag;	/* Flag for whether score up-to-date */
  int vertex_count = 0;
  SIMPLEX *simplex;
  float score;

  simplex = (SIMPLEX *) simplex_input;

/*
* Extract transformation variables from simplex array if zeroth layer
* 12/96 te
*/
  if ((simplex->score->inter_flag) && (simplex->rigid_flag == TRUE))
  {
    for (i = 0; i < 3; i++)
    {
      simplex->mol_min.transform.translate[i] =
        simplex->mol_score->transform.translate[i] +
        vertex[i] * simplex->score->minimize.translation /
        (float) (simplex->cycle + 1);

      simplex->mol_min.transform.rotate[i] =
        simplex->mol_score->transform.rotate[i] +
        vertex[i + 3] * simplex->score->minimize.rotation /
        (float) (simplex->cycle + 1);
    }

    transform_molecule (&simplex->mol_min, simplex->mol_ref);

    vertex_count = 6;

/*
*   Flag all segment intermolecular scores as out-of-date
*   3/97 te
*/
    for (segment = 0; segment < simplex->mol_min.total.segments; segment++)
      simplex->mol_min.segment[segment].score.inter.current_flag = FALSE;
  }

  else
    vertex_count = 0;

/*
* Extract torsion variables from simplex array
* 3/97 te
*/
  if ((simplex->score->intra_flag) &&
    (simplex->label->flex.minimize_flag == TRUE))
  {
    intra_current_flag = TRUE;

    for (layer = simplex->layer_inner; layer <= simplex->layer_outer; layer++)
    {
      for (i = 0; i < simplex->mol_min.layer[layer].segment_total; i++)
      {
        segment = simplex->mol_min.layer[layer].segment[i];

        if ((simplex->mol_min.segment[segment].min_flag == TRUE) &&
          ((torsion = simplex->mol_min.segment[segment].torsion_id) != NEITHER))
        {
          if
          (
            ((flex = simplex->mol_min.torsion[torsion].flex_id) != NEITHER) &&
            (simplex->label->flex.member[flex].minimize_flag == TRUE)
          )
          {
            simplex->mol_min.torsion[torsion].target_angle =
              simplex->mol_score->torsion[torsion].target_angle +
              vertex[vertex_count++] *
                simplex->score->minimize.torsion * PI / 180.0 /
                (float) (simplex->cycle + 1) /
                (float) (simplex->layer_outer - layer + 1);

            torsion_transform (&simplex->mol_min, torsion);
          }

          intra_current_flag = FALSE;
          simplex->mol_min.segment[segment].score.inter.current_flag = FALSE;
        }
      }

      if (intra_current_flag == FALSE)
        for (i = 0; i < simplex->mol_min.layer[layer].segment_total; i++)
        {
          segment = simplex->mol_min.layer[layer].segment[i];
          simplex->mol_min.segment[segment].score.intra.current_flag = FALSE;
        }
    }
  }

/*
* Evaluate the score of this ligand position
* 12/96 te
*/
  score =
    calc_score
    (
      simplex->label,
      simplex->score,
      &simplex->mol_min,
      simplex->layer_outer
    );

  if (score < simplex->mol_best.score.total)
    copy_molecule (&simplex->mol_best, &simplex->mol_min);

  return score;
}


/* /////////////////////////////////////////////////////////////////// */

float calc_score
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  int		layer_outer
) 
{
  int li, lj;			/* Layer id */
  int sli, slj;			/* Segment id in layer list */
  int si, sj;			/* Segment id */
  int asi, asj;			/* Atom id in segment list */
  int ai, aj;			/* Atom id */

/*
* Update near flag array which speeds up intramolecular score calculation
* 3/97 te
*/
  if (score->intra_flag)
    initialize_near (&score->near, molecule);

/*
* Initialize score data
* 1/97 te
*/
  reset_score_parts (&molecule->score.intra);
  reset_score_parts (&molecule->score.inter);
  molecule->score.total = 0;

/*
* Loop through all layers of interest
* 1/97 te
*/
  for (li = layer_outer; li >= 0; li--)
  {
/*
*   Loop through layer segments
*   3/97 te
*/
    for (sli = molecule->layer[li].segment_total - 1; sli >= 0; sli--)
    {
      si = molecule->layer[li].segment[sli];
      if (!molecule->segment[si].active_flag) continue;

/*
*     Calculate intermolecular score
*     3/97 te
*/
      if (score->inter_flag)
      {
        if (!molecule->segment[si].score.inter.current_flag ||
          (molecule->segment[si].score.type != molecule->score.type))
        {
          reset_score_parts (&molecule->segment[si].score.inter);
          molecule->segment[si].score.inter.current_flag = TRUE;

/*
*         Loop through segment atoms
*         3/97 te
*/
          for (asi = molecule->segment[si].atom_total - 1; asi >= 0; asi--)
          {
            ai = molecule->segment[si].atom[asi];

            calc_inter_score
            (
              label,
              score,
              molecule,
              ai,
              &molecule->segment[si].score.inter
            );
          }
        }

        sum_score
        (
          molecule->score.type,
          &molecule->score.inter,
          &molecule->segment[si].score.inter
        );
      }

/*
*     Calculate intramolecular score
*     3/97 te
*/
      if (score->intra_flag)
      {
        if (!molecule->segment[si].score.intra.current_flag ||
          (molecule->segment[si].score.type != molecule->score.type))
        {
          reset_score_parts (&molecule->segment[si].score.intra);
          molecule->segment[si].score.intra.current_flag = TRUE;

/*
*         Loop from current layer to all inner layers
*         3/97 te
*/
          for (lj = li; lj >= 0; lj--)
          {
/*
*           Loop through layer segments
*             For current layer, only consider previous segments
*             For inner layers, consider all segments
*           3/97 te
*/
            for
            (
              slj =
                (li == lj
                  ? sli - 1
                  : molecule->layer[lj].segment_total - 1);
              slj >= 0;
              slj--
            )
            {
              sj = molecule->layer[lj].segment[slj];

/*
*             Loop through segment atoms
*             3/97 te
*/
              for (asi = molecule->segment[si].atom_total - 1; asi >= 0; asi--)
              {
                ai = molecule->segment[si].atom[asi];

                for
                  (asj = molecule->segment[sj].atom_total - 1; asj >= 0; asj--)
                {
                  aj = molecule->segment[sj].atom[asj];

                  if (!score->near.flag[ai][aj])
                    calc_intra_score
                    (
                      label,
                      score,
                      molecule,
                      ai,
                      aj,
                      &molecule->segment[si].score.intra
                    );

                } /* End of asj loop */
              } /* End of asi loop */
            } /* End of slj loop */
          } /* End of lj loop */
        } /* End of current_flag IF block */

        sum_score
        (
          molecule->score.type,
          &molecule->score.intra,
          &molecule->segment[si].score.intra
        );

      } /* End of intra-scoring IF block */

      molecule->segment[si].score.type = molecule->score.type;
      molecule->segment[si].score.total =
        molecule->segment[si].score.intra.total +
        molecule->segment[si].score.inter.total;

    } /* End of sli loop */
  } /* End of li loop */

  molecule->score.total =
    molecule->score.intra.total +
    molecule->score.inter.total;

  return molecule->score.total;
}


/* ////////////////////////////////////////////////////////////////// */

void calc_intra_score
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  int		atomi,
  int		atomj,
  SCORE_PART	*intra
)
{
  if (molecule->score.type == CONTACT)
    calc_intra_contact
    (
      &score->contact,
      label,
      molecule,
      atomi,
      atomj,
      intra
    );

  else if (molecule->score.type == ENERGY)
    calc_intra_energy
    (
      &score->energy,
      label,
      molecule,
      atomi,
      atomj,
      intra
    );

  else if (molecule->score.type == CHEMICAL)
    calc_intra_energy
    (
      &score->energy,
      label,
      molecule,
      atomi,
      atomj,
      intra
    );

  else
    exit (fprintf (global.outfile,
      "ERROR update_score_pair: illegal score type requested\n"));
}


/* ////////////////////////////////////////////////////////////////// */

void calc_inter_score
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
)
{
  if (molecule->score.type == CONTACT)
    calc_inter_contact
    (
      &score->grid,
      &score->bump,
      &score->contact,
      label,
      molecule,
      atom,
      inter
    );

  else if (molecule->score.type == CHEMICAL)
    calc_inter_chemical
    (
      &score->grid,
      &score->energy,
      &score->chemical,
      label,
      molecule,
      atom,
      inter
    );

  else if (molecule->score.type == ENERGY)
    calc_inter_energy
    (
      &score->grid,
      &score->energy,
      label,
      molecule,
      atom,
      inter
    );

  else if (molecule->score.type == RMSD)
    calc_inter_rmsd
    (
      &score->grid,
      molecule,
      atom,
      inter
    );
}


/* ////////////////////////////////////////////////////////////////// */

void sum_score
(
  int		score_type,
  SCORE_PART	*sum,
  SCORE_PART	*increment
)
{
  if (score_type == CONTACT)
    sum_contact (sum, increment);

  else if (score_type == ENERGY)
    sum_energy (sum, increment);

  else if (score_type == CHEMICAL)
    sum_chemical (sum, increment);

  else if (score_type == RMSD)
    sum_rmsd (sum, increment);

  else
    exit (fprintf (global.outfile, "ERROR sum_score: unknown score type\n"));
}

/* ////////////////////////////////////////////////////////////////// */

void output_score_info
(
  DOCK		*dock,
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  LIST		*list,
  int		ligand_id,
  float		time
)
{
  int i;
  STRING80 line;

  if (global.output_volume == 't')
  {
    if (ligand_id == 0)
    {
      fprintf (global.outfile, "%s", "Result_header");
      fprintf (global.outfile, " %s", "name");

      if (label->flex.flag)
      {
        if (label->flex.multiple_anchors)
          fprintf (global.outfile, " %s", "anchors");

        if (label->flex.drive_flag && !label->flex.anchor_flag)
          fprintf (global.outfile, " %s", "conforms");
      }

      if (dock->multiple_orients)
      {
        fprintf (global.outfile, " %s", "matches");

        if (score->bump.flag)
          fprintf (global.outfile, " %s", "orients");
      }

      for (i = 0; i < SCORE_TOTAL; i++)
        if (score->type[i].flag)
        {
          fprintf (global.outfile, " %s", score->type[i].name);

          if (score->intra_flag && score->intra_flag)
            fprintf (global.outfile, " %s %s", "intra", "inter");

          if (dock->multiple_orients || score->type[i].minimize)
            fprintf (global.outfile, " %s", "rmsd");
        }

      fprintf (global.outfile, " %s\n", "time");
    }

    fprintf (global.outfile, "%s", "Results");
    fprintf (global.outfile, " %s", molecule->info.name);

    if (label->flex.flag)
    {
      if (label->flex.multiple_anchors)
        fprintf (global.outfile, " %d", dock->anchor_total);

      if (label->flex.drive_flag && !label->flex.anchor_flag)
        fprintf (global.outfile, " %d", dock->conform_total);
    }

    if (dock->multiple_orients)
    {
      fprintf (global.outfile, " %d", dock->orient_total);

      if (score->bump.flag)
        fprintf (global.outfile, " %d", dock->score_total);
    }

    for (i = 1; i < SCORE_TOTAL; i++)
      if (score->type[i].flag)
      {
        fprintf (global.outfile, " %.2f",
          list->member[i][0]->score.total);

        if (score->intra_flag && score->intra_flag)
          fprintf (global.outfile, " %.2f %.2f",
            list->member[i][0]->score.intra.total,
            list->member[i][0]->score.inter.total);

        if (dock->multiple_orients || score->type[i].minimize)
          fprintf (global.outfile, " %.2f", list->member[i][0]->transform.rmsd);
      }

    fprintf (global.outfile, " %.2f\n", time);
  }

  else
  {
    fprintf (global.outfile,
      "\n_______________________Docking_Results_______________________\n");

    fprintf (global.outfile, "%-12s : %s\n", "Name", molecule->info.name);
    fprintf (global.outfile, "%-12s : %s\n", "Description",
      molecule->info.comment);

    if (label->flex.flag)
    {
      if (label->flex.multiple_anchors)
        fprintf (global.outfile, "%-50s : %10d\n", "Total anchors",
          dock->anchor_total);

      if (label->flex.drive_flag && !label->flex.anchor_flag)
        fprintf (global.outfile, "%-50s : %10d\n", "Total conformations",
          dock->conform_total);
    }

    if (dock->multiple_orients)
    {
      fprintf (global.outfile, "%-50s : %10d\n", "Orientations tried",
        dock->orient_total);

      if (score->bump.flag)
        fprintf (global.outfile, "%-50s : %10d\n", "Orientations scored",
          dock->score_total);
    }

    fprintf (global.outfile, "\n");

    for (i = 1; i < SCORE_TOTAL; i++)
      if (score->type[i].flag)
      {
        if (score->intra_flag || score->inter_flag)
        {
          sprintf (line, "Best %s%s score",
            score->intra_flag
              ? (score->inter_flag ? "" : "intramolecular ")
              : (score->inter_flag ? "intermolecular " : ""),
            score->type[i].name);

          fprintf (global.outfile, "%-50s : %10.2f\n", line,
            list->member[i][0]->score.total);

          if (score->intra_flag && score->inter_flag)
          {
            sprintf (line, "  Intramolecular %s score", score->type[i].name);
            fprintf (global.outfile, "%-50s : %10.2f\n", line,
              list->member[i][0]->score.intra.total);

            sprintf (line, "  Intermolecular %s score", score->type[i].name);
            fprintf (global.outfile, "%-50s : %10.2f\n", line,
              list->member[i][0]->score.inter.total);
          }
        }

        if (dock->multiple_orients|| score->type[i].minimize)
        {
          sprintf (line, "RMSD of best %s scorer (A)",
            score->type[i].name);
          fprintf (global.outfile, "%-50s : %10.2f\n", line,
            list->member[i][0]->transform.rmsd);
        }

        if
        (
          !dock->rank_ligands &&
          dock->write_orients &&
          !dock->rank_orients
        )
          fprintf (global.outfile, "%-50s : %10d\n", "Orientations written",
            score->type[i].number_written);

        fprintf (global.outfile, "\n");
      }

    fprintf (global.outfile,
      "%-50s : %10.2f\n", "Elapsed cpu time (sec)", time);

    fprintf (global.outfile, "\n\n");
  }

  fflush (global.outfile);
}


/* //////////////////////////////////////////////////////////////////

Top score list processing routines

////////////////////////////////////////////////////////////////// */

void allocate_lists
(
  SCORE		*score,
  LIST		*list,
  int		length,
  int		lite_flag
)
{
  int i;

  list->lite_flag = lite_flag;

  for (i = 0; i < SCORE_TOTAL; i++)
    if (score->type[i].flag == TRUE)
    {
      list->total[i] = length;
      allocate_list (list, i);
      list->total[i] = 0;
    }
}

/* //////////////////////////////////////////////////////////////////// */

void allocate_list
(
  LIST		*list,
  int		type
)
{
  int i;

  ecalloc
  (
    (void **) &list->member[type],
    list->total[type],
    sizeof (MOLECULE *),
    "molecule list",
    global.outfile
  );

  for (i = 0; i < list->total[type]; i++)
  {
    ecalloc
    (
      (void **) &list->member[type][i],
      1,
      sizeof (MOLECULE),
      "molecule list",
      global.outfile
    );

    reset_molecule (list->member[type][i]);
  }

  list->max[type] = list->total[type];
}

/* //////////////////////////////////////////////////////////////////// */

void reset_lists (SCORE *score, LIST *list)
{
  int i;

  for (i = 0; i < SCORE_TOTAL; i++)
    if (score->type[i].flag == TRUE)
      reset_list (list, i);
}

/* //////////////////////////////////////////////////////////////////// */

void reset_list (LIST *list, int type)
{
  int i;

  for (i = 0; i < list->max[type]; i++)
    reset_molecule (list->member[type][i]);

  list->total[type] = 0;
}

/* //////////////////////////////////////////////////////////////////// */

void free_lists (SCORE *score, LIST *list)
{
  int i;

  for (i = 0; i < SCORE_TOTAL; i++)
    if (score->type[i].flag == TRUE)
      free_list (list, i);
}


/* //////////////////////////////////////////////////////////////////// */

void free_list (LIST *list, int type)
{
  int i;

  for (i = 0; i < list->max[type]; i++)
  {
    free_molecule (list->member[type][i]);
    efree ((void **) &list->member[type][i]);
  }

  efree ((void **) &list->member[type]);
  list->max[type] = 0;
}


/* //////////////////////////////////////////////////////////////////// */

void reallocate_lists (SCORE *score, LIST *list)
{
  int i;

  for (i = 0; i < SCORE_TOTAL; i++)
    if (score->type[i].flag == TRUE)
      reallocate_list (list, i);
}


/* //////////////////////////////////////////////////////////////////// */

void reallocate_list (LIST *list, int type)
{
  if (list->max[type] < list->total[type])
  {
    free_list (list, type);
    list->max[type] = list->total[type];
    allocate_list (list, type);
  }
}

/* //////////////////////////////////////////////////////////////////// */

void copy_lists (SCORE *score, LIST *copy, LIST *original)
{
  int i;

  for (i = 0; i < SCORE_TOTAL; i++)
    if (score->type[i].flag == TRUE)
      copy_list (copy, original, i);
}

/* //////////////////////////////////////////////////////////////////// */

void copy_list (LIST *copy, LIST *original, int type)
{
  int i;

  copy->total[type] = original->total[type];
  reallocate_list (copy, type);

  for (i = 0; i < original->total[type]; i++)
    copy_member
      (original->lite_flag, copy->member[type][i], original->member[type][i]);

  copy->total[type] = original->total[type];
}

/* //////////////////////////////////////////////////////////////////// */

void inter_lists (SCORE *score, LIST *list, MOLECULE *reference)
{
  int i, j;

  for (i = 0; i < SCORE_TOTAL; i++)
    if (score->type[i].flag == TRUE)
    {
      for (j = 0; j < list->total[i]; j++)
        list->member[i][j]->score.total = 
          list->member[i][j]->score.inter.total +
          list->member[i][j]->transform.heavy_total *
            score->type[i].size_penalty;

      sort_list (list, i);

      if (i != RMSD)
      {
        for (j = 0; j < list->total[i]; j++)
          if (list->member[i][j]->score.total > 0)
            list->total[i] = j;

        if (list->total[i] == 0)
        {
          reference->score.type = i;
          reference->score.total = 0;
          update_list (list, reference);
        }
      }
    }
}

/* //////////////////////////////////////////////////////////////////// */

void merge_lists (SCORE *score, LIST *list, LIST *add_list)
{
  int i, j;

  for (i = 0; i < SCORE_TOTAL; i++)
    if (score->type[i].flag == TRUE)
      for (j = 0; (j < add_list->total[i]) && (j < list->max[i]); j++)
        update_list (list, add_list->member[i][j]);
}

/* //////////////////////////////////////////////////////////////////// */

void update_list
(
  LIST *list,
  MOLECULE *molecule
)
{
  int type;

  type = molecule->score.type;

  if ((type < 0) || (type >= SCORE_TOTAL))
  {
    fprintf (global.outfile, "ERROR update_list: Accessing illegal list.\n");
    exit (EXIT_FAILURE);
  }

  if (list->total[type] < list->max[type])
    copy_member
    (
      list->lite_flag,
      list->member[type][list->total[type]++],
      molecule
    );

  else if (type == 0)
    copy_member
    (
      list->lite_flag,
      list->member[type][list->max[type] - 1],
      molecule
    );

  else if (molecule->score.total <
    list->member[type][list->max[type] - 1]->score.total)
    copy_member
    (
      list->lite_flag,
      list->member[type][list->max[type] - 1],
      molecule
    );

  else
    return;

  sort_list (list, type);
}

/* //////////////////////////////////////////////////////////////////// */

int compare_score (MOLECULE *reference, MOLECULE *candidate)
{
  if
  (
    (reference->info.name != NULL) &&
    (candidate->info.name != NULL) &&
    !strcmp (reference->info.name, candidate->info.name)
  )
  {
    if (candidate->score.total < reference->score.total)
      return TRUE;
  }

  else
  {
    if (candidate->score.inter.total < reference->score.inter.total)
      return TRUE;
  }

  return FALSE;
}

/* //////////////////////////////////////////////////////////////////// */

void sort_list (LIST *list, int type)
{
  if (list->total[type] > 1)
    qsort
    (
      list->member[type], 
      list->total[type], 
      sizeof (MOLECULE *),
      (int (*)()) compare_member
    );
}

/* //////////////////////////////////////////////////////////////////// */

void copy_member (int lite_flag, MOLECULE *old, MOLECULE *new)
{
  if (lite_flag == TRUE)
  {
    copy_info (old, new);
    copy_transform (old, new);
    copy_score (&old->score, &new->score);
    copy_torsions (old, new);
    copy_keys (old, new);
  }

  else
    copy_molecule (old, new);
}

/* //////////////////////////////////////////////////////////////////// */

int compare_member (MOLECULE **mol1, MOLECULE **mol2)
{
  return ((*mol1)->score.total == (*mol2)->score.total) ? 0 :
    (((*mol1)->score.total > (*mol2)->score.total) ? 1 : -1);
}
 
/* //////////////////////////////////////////////////////////////////// */

int print_lists (SCORE *score, LIST *list, FILE *file)
{
  int i;

  if (file == NULL)
    return FALSE;

  for (i = 0; i < SCORE_TOTAL; i++)
    if (score->type[i].flag == TRUE)
      print_list (score, list, i, file);

  return TRUE;
}

 
/* //////////////////////////////////////////////////////////////////// */

int print_list (SCORE *score, LIST *list, int type, FILE *file)
{
  int i;

  if (file == NULL)
    return FALSE;

  if (score->type[type].flag == TRUE)
  {
    fprintf (file, "%-4s %8s %8s\n", "rank", score->type[type].name, "rmsd");

    for (i = 0; i < list->total[type]; i++)
      fprintf (file, "%4d %8.2f %8.2f\n",
        i + 1, list->member[type][i]->score.total,
        list->member[type][i]->transform.rmsd);
  }

  else
    return FALSE;

  return TRUE;
}

/* //////////////////////////////////////////////////////////////////// */

int save_lists (SCORE *score, LIST *list, FILE *file)
{
  int i, j;

  for (i = 0; i < SCORE_TOTAL; i++)
  {
    if (score->type[i].flag == TRUE)
    {
      efwrite (&list->total[i], sizeof (int), 1, file);

      for (j = 0; j < list->total[i]; j++)
        save_molecule (list->member[i][j], file);
    }
  }

  return TRUE;
}

/* //////////////////////////////////////////////////////////////////// */

int load_lists (SCORE *score, LIST *list, FILE *file)
{
  int i, j;

  for (i = 0; i < SCORE_TOTAL; i++)
  {
    if (score->type[i].flag == TRUE)
    {
      efread (&list->total[i], sizeof (int), 1, file);

      if (list->total[i] > list->max[i])
        return FALSE;

      for (j = 0; j < list->total[i]; j++)
        load_molecule (list->member[i][j], file);
    }
  }

  return TRUE;
}


/* //////////////////////////////////////////////////////////////////// */

void shrink_list
(
  LIST		*full,
  LIST		*shrunk,
  int		type,
  float		rank_rmsd_cutoff
)
{
  int i, j;
  int remainder;
  static int *discard = NULL;
  static int size = 0;
  int discard_member;
  float rmsd;

  if (full->total[type] <= 0)
  {
    fprintf (global.outfile,
      "WARNING shrink_list: No members provided\n");
    return;
  }

/*
* Eliminate members with extremely repulsive scores
* 1/98 te
*/
  for (i = 0; i < full->total[type]; i++)
    if (full->member[type][i]->score.total >= INITIAL_SCORE)
      full->total[type] = i;
  
/*
* Allocate arrays
* 3/97 te
*/
  if (size < full->total[type])
  {
    efree ((void **) &discard);

    size = full->total[type];

    ecalloc
    (
      (void **) &discard,
      size,
      sizeof (int),
      "shrink discard array",
      global.outfile
    );
  }

  else
    memset (discard, 0, size * sizeof (int));

/*
* Discard members outside rank/rmsd cutoff of better scoring members
* 1/98 te
*/
  for (i = 0; i < full->total[type]; i++)
    if (!discard[i])
      for (j = i + 1; j < full->total[type]; j++)
        if (!discard[j])
        {
          rmsd = calc_layer_rmsd (full->member[type][i], full->member[type][j]);
          rmsd = MAX (rmsd, 0.001);

          if ((float) j / rmsd > rank_rmsd_cutoff)
            discard[j] = TRUE;
        }

/*
* Transfer remainder over to shrunk list
* 3/97 te
*/
  if (shrunk->total[type] < full->total[type])
  {
    shrunk->total[type] = full->total[type];
    reallocate_list (shrunk, type);
    shrunk->total[type] = 0;
  }

  else
    reset_list (shrunk, type);

  for (i = 0; i < full->total[type]; i++)
    if (!discard[i])
      update_list (shrunk, full->member[type][i]);
}


/* /////////////////////////////////////////////////////////////////// */

void initialize_near
(
  NEAR		*near,
  MOLECULE	*molecule
) 
{
  int i, j, k;
  int ni, nni;

/*
* Check if near data is up-to-date
* 3/97 te
*/
  if
  (
    (near->total == molecule->total.atoms) &&
    (near->name != NULL) &&
    (molecule->info.name != NULL) &&
    !strcmp (near->name, molecule->info.name)
  )
    return;

/*
* Either reallocate space for the flags or reset them
* 3/97 te
*/
  if (near->total < molecule->total.atoms)
  {
    free_near (near);
    near->total = molecule->total.atoms;
    allocate_near (near);
  }

  else
    for (i = 0; i < near->total; i++)
      memset (near->flag[i], 0, near->total * sizeof (int));

  vstrcpy (&near->name, molecule->info.name);

  for (i = 0; i < molecule->total.atoms; i++)
    for (j = 0; j < molecule->atom[i].neighbor_total; j++)
    {
      ni = molecule->atom[i].neighbor[j].id;
      near->flag[i][ni] = near->flag[ni][i] = TRUE;

      for (k = 0; k < molecule->atom[ni].neighbor_total; k++)
      {
        nni = molecule->atom[ni].neighbor[k].id;
        near->flag[i][nni] = near->flag[nni][i] = TRUE;
      }
    }
}


/* /////////////////////////////////////////////////////////////////// */

void allocate_near (NEAR *near)
{
  int i;

  ecalloc
  (
    (void **) &near->flag,
    near->total,
    sizeof (int *),
    "near flag array",
    global.outfile
  );

  for (i = 0; i < near->total; i++)
    ecalloc
    (
      (void **) &near->flag[i],
      near->total,
      sizeof (int *),
      "near flag array",
      global.outfile
    );
}


/* /////////////////////////////////////////////////////////////////// */

void free_near (NEAR *near)
{
  int i;

  for (i = 0; i < near->total; i++)
    efree ((void **) &near->flag[i]);

  efree ((void **) &near->flag);
  efree ((void **) &near->name);
}


/* /////////////////////////////////////////////////////////////////// */

void free_scores (LABEL *label, SCORE *score)
{
  free_near (&score->near);

  if (score->grid.flag)
    free_grids
    (
      &score->grid,
      &score->bump,
      &score->contact,
      &score->chemical,
      &score->energy,
      &label->chemical
    );

  if (score->energy.flag || score->chemical.flag)
    free_vdw_energy (&score->energy);

  if (score->chemical.flag)
    free_table (&label->chemical, &label->chemical.score_table);
}
