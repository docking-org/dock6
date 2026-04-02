/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/
#include "define.h"
#include "mol.h"
#include "global.h"
#include "search.h"
#include "label.h"
#include "score.h"
#include "grid.h"
#include "score_grid.h"

void make_grids
(
  SCORE_GRID *grid,
  SCORE_BUMP *bump,
  SCORE_CONTACT *contact,
  SCORE_CHEMICAL *chemical,
  SCORE_ENERGY *energy,
  LABEL *label,
  MOLECULE *receptor
)
{
  int i, j, k, l, atomi, index;
  int ilo, ihi, jlo, jhi, klo, khi;
  int contact_minimum;
  int chem_id;
  int inside_grid;
  float dist, dist_sq;
  float dist_inv;
  float dist_min;
  float dist_sq_min;
  float dist_power;
  int report_increment;

  XYZ origin_coord, gridpt_coord;
  int grid_cutoff, grid_coord[3];
  float square_distance (XYZ, XYZ);

  dist_min = DISTANCE_MIN;
  dist_sq_min = SQR (DISTANCE_MIN);

  if (!energy->vdw_init_flag)
    initialize_vdw_energy (energy, &label->vdw);

/*
  if (energy->flag)
    energy->dielectric_factor = 322.0 / energy->dielectric_factor;
*/

/*
* Determine the longest interaction distance to consider.
* Take the square of each radii in preparation for grid calculation.
*/
  grid->distance = 0.0;

  if (bump->flag)
  {
    for (i = 0, bump->distance = 0.0; i < label->vdw.total; i++)
      bump->distance = MAX
        (bump->distance,
        2.0 * bump->clash_overlap * label->vdw.member[i].radius);

    grid->distance = MAX (grid->distance, bump->distance);
  }

  if (contact->flag)
    grid->distance = MAX (grid->distance, contact->distance);

  if (energy->flag)
    grid->distance = MAX (grid->distance, energy->distance);

  grid_cutoff = (int) (grid->distance / grid->spacing + 1.0);
  // trent e. balius:  for multigrid docking 
  // we may wish to have grids for single residues
  // 
  if (receptor->total.atoms >= 10)
     report_increment = receptor->total.atoms / 10;
  else{
     fprintf (global.outfile,"Warning: receptor has less than 10 atoms\n");
     report_increment = 1;
  }
  for (atomi = 0; atomi < receptor->total.atoms; atomi++)
  {
    inside_grid = TRUE;
    for (i = 0; i < 3; i++)
    {
      origin_coord[i] = receptor->coord[atomi][i] - grid->origin[i];
      grid_coord[i] = NINT (origin_coord[i] / grid->spacing);

      if (grid_coord[i] >= (grid->span[i] + grid_cutoff))
        inside_grid = FALSE;
      if (grid_coord[i] < (0 - grid_cutoff))
        inside_grid = FALSE;
    }

    if (inside_grid)
    {
      ilo = MAX (0, (grid_coord[0] - grid_cutoff));
      ihi = MIN (grid->span[0], (grid_coord[0] + grid_cutoff + 1));
      jlo = MAX (0, (grid_coord[1] - grid_cutoff));
      jhi = MIN (grid->span[1], (grid_coord[1] + grid_cutoff + 1));
      klo = MAX (0, (grid_coord[2] - grid_cutoff));
      khi = MIN (grid->span[2], (grid_coord[2] + grid_cutoff + 1));

      chem_id = receptor->atom[atomi].chem_id;

      for (i = ilo; i < ihi; i++)
      {
        gridpt_coord[0] =
          ((float) (i)) *
          grid->spacing + grid->origin[0];

        for (j = jlo; j < jhi; j++)
        {
          gridpt_coord[1] =
            ((float) (j)) *
            grid->spacing + grid->origin[1];

          for (k = klo; k < khi; k++)
          {
            gridpt_coord[2] =
              ((float) (k)) *
              grid->spacing + grid->origin[2];

            index =
              grid->span[0] * grid->span[1] * k +
              grid->span[0] * j + i;

            dist_sq = square_distance (receptor->coord[atomi], gridpt_coord);

            if (dist_sq > dist_sq_min)
              dist = sqrt (dist_sq);

            else
            {
              dist = dist_min;
              dist_sq = dist_sq_min;
            }

            if (dist <= grid->distance)
            {
/*
*             Increment bump grids
*/
              if
              (
                bump->flag &&
                (receptor->atom[atomi].heavy_flag == TRUE) &&
                (dist <= bump->distance)
              )
                for (l = 0; l < label->vdw.total; l++)
                  if
                  (
                    (label->vdw.member[l].heavy_flag == TRUE) &&
                    (label->vdw.member[l].bump_id < bump->grid[index]) &&
                    (dist < bump->clash_overlap *
                      (label->vdw.member[receptor->atom[atomi].vdw_id].radius +
                      label->vdw.member[l].radius))
                  )
                    bump->grid[index] = label->vdw.member[l].bump_id;

/*
*             Increment contact grids
*/
              if
              (
                contact->flag &&
                (dist <= contact->distance) &&
                (receptor->atom[atomi].heavy_flag == TRUE)
              )
                contact->grid[index] -= 1;

/*
*             Increment energy and chemical grids
*/
              if
              (
                energy->flag &&
                (label->vdw.member[receptor->atom[atomi].vdw_id].well_depth !=
                  0.0) &&
                (dist <= energy->distance)
              )
              {
                dist_inv = 1.0 / dist;

                POWER (dist_inv, energy->repulsive_exponent, dist_power);
                energy->avdw[index] +=
                  energy->vdwA[receptor->atom[atomi].vdw_id] *
                  dist_power;

                POWER (dist_inv, energy->attractive_exponent, dist_power);
                energy->bvdw[index] +=
                  energy->vdwB[receptor->atom[atomi].vdw_id] *
                  dist_power;

                if (chemical->flag)
                  chemical->grid[chem_id][index] +=
                    energy->vdwB[receptor->atom[atomi].vdw_id] *
                    dist_power;

                if (energy->distance_dielectric)
                  energy->es[index] += energy->dielectric_factor *
                    receptor->atom[atomi].charge * SQR (dist_inv);

                else
                  energy->es[index] += energy->dielectric_factor *
                    receptor->atom[atomi].charge * dist_inv;
              }
            }
          }
        }
      }
    }

    if (!(atomi % report_increment))
    {
      fprintf (global.outfile, "%-40s: %8d\n",
        "Percent of protein atoms processed",
        atomi / report_increment * 10);
      fflush (global.outfile);
    }
  }
}


