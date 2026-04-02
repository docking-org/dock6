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
#include "screen.h"
#include "match.h"
#include "io.h"
#include "mol_prep.h"
#include "rotrans.h"


/* ////////////////////////////////////////////////////////////// 

Routine to manage matching procedure
3/96 te

////////////////////////////////////////////////////////////// */

int get_match
(
  DOCK		*dock,
  MATCH		*match,
  LABEL		*label,
  MOLECULE      *mol_conf,
  int           conformation_id,
  int           orientation_id
)
{
  if (orientation_id == 0)
  {
    match->total = 0;
    match->max = 0;
    match->cycle = 0;
    match->chiral_flag = TRUE;

    if (!match->init_flag)
      init_match_site (match, label);

    if (!match->centers_flag)
      get_ligand_centers (match, mol_conf);

    if (conformation_id == 0)
      init_match_ligand (match, label);
  }

  if (match->auto_flag)
  {
    while (match->total >= match->max)
    {
      match->cycle++;
      initialize_adjacency (match, label);

      if (compute_matches (match) == EOF)
        return EOF;

      match->total = 0;
    }

    copy_clique (&match->clique, &match->clique_sort[match->total]);
    match->total++;
    return TRUE;
  }

  else
  {
    if (orientation_id == 0)
    {
      match->cycle++;
      initialize_adjacency (match, label);
    }

    if (match_distance (match) != EOF)
    {
      match->total++;
      return TRUE;
    }
   
    else
    {
      if (!dock->multiple_ligands || !dock->multiple_conforms)
        output_match_info (match);

      return EOF;
    }
  }
}


/* ////////////////////////////////////////////////////////////// 

Routine to free matching arrays
5/97 te

////////////////////////////////////////////////////////////// */

void free_matches
(
  LABEL		*label,
  MATCH		*match
)
{
  free_clique_list (match);
  free_clique_link (match);
  free_match (match);
  free_match_site (match, label);
  free_ligand_centers (match);
}


/* ////////////////////////////////////////////////////////////// 

Routine to allocate space for matching arrays
3/96 te

////////////////////////////////////////////////////////////// */

void init_match_site
(
  MATCH		*match,
  LABEL		*label
)
{
/*
* Construct match table
* 11/96 te
*/
  if (match->chemical_flag && !label->chemical.match_table)
    get_table
    (
      &label->chemical,
      label->chemical.match_file_name,
      &label->chemical.match_table
    );

/*
* Read in site points for matching
* 8/96 te
*/
  get_site (match, label);
  get_centers (match, label);

/*
* Update the critical cluster identifiers of the site points
* 7/95 te
*/
  id_site_points (match);

/*
* Reorder site points so that critical clusters come first
* 7/95 te
*/
  if (match->critical_flag)
    order_site_points (match);

/*
* Allocate memory for the critical cluster filters
* 7/95 te
*/
  ecalloc
  (
    (void **) &match->critical_filter,
    match->cluster_total + 1,
    sizeof (char *),
    "matching filter",
    global.outfile
  );

/*
* Construct distance matrix for receptor site points
* 6/95 te
*/
  match->bin_length =
    calculate_distances
    (
      &match->receptor_site,
      &match->receptor_distances,
      &match->receptor_distance_size
    );

/*
* Make keys for pharmacophore search
* 1/97 te
*/
  if (label->chemical.screen.pharmaco_flag)
  {
    if (label->chemical.screen.fold_flag)
    {
      update_folded_keys (&label->chemical, &match->receptor_site, 0);
      return;
    }

    else
      update_unfolded_keys (&label->chemical, &match->receptor_site, 0);
  }

  ecalloc
  (
    (void **) &match->clique_adjacency_total,
    match->clique_size_max + 1,
    sizeof (int),
    "clique adjacency list total",
    global.outfile
  );

  ecalloc
  (
    (void **) &match->clique_adjacency_list,
    match->clique_size_max + 1,
    sizeof (EDGE *),
    "clique adjacency list",
    global.outfile
  );

  ecalloc
  (
    (void **) &match->adjacent,
    match->clique_size_max + 1,
    sizeof (int),
    "adjacent iterator list",
    global.outfile
  );

  match->clique.edge_max = match->clique_size_max;
  allocate_clique (&match->clique);

  ecalloc
  (
    (void **) &match->degenerate,
    match->clique_size_max,
    sizeof (int),
    "degenerate clique list",
    global.outfile
  );

  emalloc
  (
    (void **) &match->histogram,
    match->clique_size_max * sizeof (int),
    "clique histogram",
    global.outfile
  );

/*
* Construct matrix needed to calculate chirality
* 10/95 te
*/
  calculate_vectors
  (
    &match->receptor_site,
    &match->receptor_vectors,
    &match->receptor_vector_size
  );

/*
* Construct receptor distance bins
* 1/97 te
*/
  make_distance_bins (match);

  allocate_clique_atoms (match);

  match->init_flag = TRUE;
}


/* ////////////////////////////////////////////////////////////// 

Routine to free space for matching arrays
4/97 te

////////////////////////////////////////////////////////////// */

void free_match_site
(
  MATCH		*match,
  LABEL		*label
)
{
/*
* Free match table
* 5/97 te
*/
  if (match->chemical_flag && label->chemical.match_table)
    free_table
    (
      &label->chemical,
      &label->chemical.match_table
    );

/*
* Free site points for matching
* 5/97 te
*/
  free_molecule (&match->receptor_site);

  free_distances
  (
    &match->receptor_distances,
    &match->receptor_distance_size
  );

  efree ((void **) &match->clique_adjacency_total);
  efree ((void **) &match->clique_adjacency_list);
  efree ((void **) &match->adjacent);
  free_clique (&match->clique);
  efree ((void **) &match->degenerate);
  efree ((void **) &match->histogram);

  if (match->critical_flag)
    free_site_points (match);

  efree ((void **) &match->critical_filter);

  free_vectors
  (
    &match->receptor_vectors,
    &match->receptor_vector_size
  );

  free_distance_bins (match);
  free_molecule (&match->ligand_clique);
  free_molecule (&match->receptor_clique);
}


/* ////////////////////////////////////////////////////////////// 

Routine to allocate space for matching arrays
3/96 te

////////////////////////////////////////////////////////////// */

void allocate_match (MATCH *match)
{
  int i;

  ecalloc
  (
    (void **) &match->adjacency_total,
    match->node_max,
    sizeof (int),
    "adjacency list total",
    global.outfile
  );

  ecalloc
  (
    (void **) &match->adjacency_max,
    match->node_max,
    sizeof (int),
    "adjacency list max",
    global.outfile
  );

  ecalloc
  (
    (void **) &match->adjacency_list,
    match->node_max,
    sizeof (EDGE *),
    "adjacency matrix",
    global.outfile
  );

  for (i = 0; i < match->clique_size_max + 1; i++)
    ecalloc
    (
      (void **) &match->clique_adjacency_list[i],
      match->node_max,
      sizeof (EDGE),
      "clique adjacency matrix",
      global.outfile
    );

  for (i = 0; i < match->node_max; i++)
  {
    match->clique_adjacency_list[0][i].node = i;
    match->clique_adjacency_list[0][i].residual = 0;
  }

/*
* Allocate memory for and initialize the chemical matching filter
* 7/95 te
*/
  emalloc
  (
    (void **) &match->chemical_filter,
    match->node_max * sizeof (char),
    "chemical label filter",
    global.outfile
  );

  memset (match->chemical_filter, 1, match->node_max);

/*
* Allocate memory for the critical cluster filters
* 7/95 te
*/
  for (i = 0; i < match->cluster_total + 1; i++)
    ecalloc
    (
      (void **) &match->critical_filter[i],
      match->node_max,
      sizeof (char),
      "matching filter",
      global.outfile
    );
}


/* ////////////////////////////////////////////////////////////// 

Routine to reset match variables.
11/96 te

////////////////////////////////////////////////////////////// */

void reset_match (MATCH *match)
{
  int i;

  memset (match->adjacency_total, 0, match->node_max * sizeof (int));

  memset
    (match->clique_adjacency_total, 0,
    (match->clique_size_max + 1) * sizeof (int));

  for (i = 1; i < match->clique_size_max + 1; i++)
    memset
      (match->clique_adjacency_list[i], 0,
      match->node_max * sizeof (EDGE));

  reset_clique (&match->clique);
  memset (match->adjacent, 0, (match->clique_size_max + 1) * sizeof (int));
  memset (match->degenerate, 0, match->clique_size_max * sizeof (int));
  memset (match->histogram, 0, match->clique_size_max * sizeof (int));
}

/* ////////////////////////////////////////////////////////////// 

Routine to free match variables.
11/96 te

////////////////////////////////////////////////////////////// */

void free_match (MATCH *match)
{
  int i;

  efree ((void **) &match->adjacency_total);
  efree ((void **) &match->adjacency_max);

  for (i = 0; i < match->node_max; i++)
    efree ((void **) &match->adjacency_list[i]);

  efree ((void **) &match->adjacency_list);

  if (match->clique_adjacency_list)
    for (i = 0; i < match->clique_size_max + 1; i++)
      efree ((void **) &match->clique_adjacency_list[i]);

  efree ((void **) &match->chemical_filter);

  if (match->critical_filter)
    for (i = 0; i < match->cluster_total + 1; i++)
      efree ((void **) &match->critical_filter[i]);

  match->node_max = 0;
}


/* ////////////////////////////////////////////////////////////// 

Routine to reallocate match variables.
3/97 te

////////////////////////////////////////////////////////////// */

void reallocate_match (MATCH *match)
{
  if (match->node_max < match->node_total)
  {
    free_match (match);
    match->node_max = match->node_total;
    allocate_match (match);
  }
}


/* ///////////////////////////////////////////////////////////// */

void allocate_clique_atoms (MATCH *match)
{
  match->ligand_clique.max.atoms = match->clique_size_max;
  allocate_molecule (&match->ligand_clique);

  match->receptor_clique.max.atoms = match->clique_size_max;
  allocate_molecule (&match->receptor_clique);
}


/* ////////////////////////////////////////////////////////////// */

void allocate_clique (CLIQUE *clique)
{
  if (clique->edge_max > 0)
    ecalloc
    (
      (void **) &clique->edge,
      clique->edge_max,
      sizeof (EDGE),
      "clique edges",
      global.outfile
    );
}
    
/* ////////////////////////////////////////////////////////////// */

void reallocate_clique (CLIQUE *clique)
{
  if (clique->edge_max < clique->edge_total)
  {
    free_clique (clique);
    clique->edge_max = clique->edge_total;
    allocate_clique (clique);
  }
}

/* ////////////////////////////////////////////////////////////// */

void free_clique (CLIQUE *clique)
{
  if (clique->edge_max > 0)
    efree ((void **) &clique->edge);

  clique->edge_max = 0;
}

/* ////////////////////////////////////////////////////////////// */

void reset_clique (CLIQUE *clique)
{
  if (clique->edge_max > 0)
    memset (clique->edge, 0, clique->edge_max * sizeof (EDGE));

  clique->edge_total = 0;
  clique->residual = 0;
  clique->reflect_flag = FALSE;
}

/* ////////////////////////////////////////////////////////////// */

void copy_clique (CLIQUE *copy, CLIQUE *original)
{
  int i;

  copy->edge_total = original->edge_total;
  reallocate_clique (copy);

  for (i = 0; i < original->edge_total; i++)
    copy->edge[i] = original->edge[i];

  copy->residual = original->residual;
  copy->reflect_flag = original->reflect_flag;
}

/* ////////////////////////////////////////////////////////////// 

Routine to set up arrays used during matching.
7/95 te

////////////////////////////////////////////////////////////// */

int init_match_ligand
(
  MATCH		*match,
  LABEL		*label
)
{
  match->node_total =
    match->ligand_center.total.atoms * match->receptor_site.total.atoms;

  reallocate_match (match);

/*
* If chemical matching is performed, construct new chemical filter
* 7/95 te
*/
  if (match->chemical_flag)
    make_chemical_filter (&label->chemical, match);

/*
* If critical cluster matching is performed, construct new filter
* to guide matching
* 6/95 te
*/
  if (match->critical_flag)
    make_critical_filter (match);

  return TRUE;
}


/* ////////////////////////////////////////////////////////////// 

Routine to set up arrays used during matching.
7/95 te

////////////////////////////////////////////////////////////// */

int initialize_adjacency
(
  MATCH		*match,
  LABEL		*label
)
{
  reset_match (match);

  if (label->chemical.screen.process_flag)
    compute_chemical_adjacency (match, label);

  else
    compute_adjacency (match);

  match->clique_adjacency_total[0] = match->node_total;

  return TRUE;
}


/* ////////////////////////////////////////////////////////////// 

Routine to get site points
11/96 te

////////////////////////////////////////////////////////////// */

void get_site
(
  MATCH		*match,
  LABEL		*label
)
{
  FILE *file;

/*
* Read in receptor site points
* 3/96 te
*/
  file = efopen (match->receptor_file_name, "r", global.outfile);

  if (read_molecule
    (&match->receptor_site, NULL, match->receptor_file_name, file, TRUE)
    != TRUE)
    exit (fprintf (global.outfile,
      "ERROR get_site: Unable to read receptor site points from %s\n",
      match->receptor_file_name));

  if (match->clique_size_min > match->receptor_site.total.atoms)
    exit (fprintf (global.outfile,
      "ERROR get_site: clique_size_min > number of site points.\n"));

  if (prepare_molecule
    (
      &match->receptor_site,
      match->receptor_file_name,
      file,
      label,
      'a',
      FALSE,
      match->chemical_flag,
      label->chemical.screen.flag
    ) != TRUE)
    exit (fprintf (global.outfile,
      "ERROR get_site: Unable to prepare receptor site points.\n"));

  fclose (file);
}


/* ////////////////////////////////////////////////////////////// 

Routine to get (or just allocate) ligand centers
3/96 te

////////////////////////////////////////////////////////////// */

void get_centers
(
  MATCH		*match,
  LABEL		*label
)
{
  FILE *file;

/*
* If requested, read in ligand centers
* 3/96 te
*/
  if (match->centers_flag)
  {
    file = efopen (match->ligand_file_name, "r", global.outfile);

    if (read_molecule
      (
        &match->ligand_center,
        NULL,
        match->ligand_file_name,
        file,
        TRUE
      ) != TRUE)
      exit (fprintf (global.outfile,
        "ERROR get_centers: Unable to read ligand site points from %s\n",
        match->ligand_file_name));

    if (prepare_molecule
      (
        &match->ligand_center,
        match->ligand_file_name,
        file,
        label,
        'a',
        FALSE,
        match->chemical_flag,
        FALSE
      ) != TRUE)
      exit (fprintf (global.outfile,
        "ERROR get_centers: Unable to prepare ligand site points.\n"));

    center_of_mass
    (
      match->ligand_center.coord,
      match->ligand_center.total.atoms,
      match->ligand_center.transform.com
    );

    fclose (file);

    calculate_distances
    (
      &match->ligand_center,
      &match->ligand_distances,
      &match->ligand_distance_size
    );

    calculate_vectors
    (
      &match->ligand_center,
      &match->ligand_vectors,
      &match->ligand_vector_size
    );
  }
}


/* ////////////////////////////////////////////////////////////// 

Routine to fix the critical cluster identifiers of the site points
so that the zeroth cluster (unclustered site points) becomes the
last cluster.
7/95 te

////////////////////////////////////////////////////////////// */

void id_site_points (MATCH *match)
{
  int i;

  if (match->critical_flag)
  {
    match->cluster_total = match->receptor_site.total.substs;

/*
*   Make sure the largest cluster is not greater than the largest
*   allowed clique size
*   7/95 te
*/
    if (match->cluster_total > match->clique_size_max)
      exit (fprintf (global.outfile,
        "ERROR id_site_points: "
        "Critical cluster identifiers cannot exceed nodes_maximum,\n"
        "      or if a zeroth cluster is present, nodes_maximum - 1.\n"));
  }

/*
* If critical cluster filtering is not to be performed, then
* set all cluster identifiers to zero.
* 7/95 te
*/
  else
  {
    for (i = 0; i < match->receptor_site.total.atoms; i++)
      match->receptor_site.atom[i].subst_id = 0;

    match->cluster_total = 0;
  }
}


/* ////////////////////////////////////////////////////////////// 

Routine to reshuffle order of the site points to get the members of
critical clusters first
7/95 te

////////////////////////////////////////////////////////////// */

void order_site_points (MATCH *match)
{
  int i, j;
  MOLECULE temporary = {0};
  int compare_site_points ();

/*
* Allocate space for the cluster lists
* 11/96 te
*/
  ecalloc
  (
    (void **) &match->cluster,
    match->cluster_total,
    sizeof (int *),
    "cluster list",
    global.outfile
  );

  ecalloc
  (
    (void **) &match->cluster_size,
    match->cluster_total,
    sizeof (int),
    "cluster list",
    global.outfile
  );

/*
* Determine the size of each cluster
* 11/96 te
*/
  for (i = 0; i < match->receptor_site.total.atoms; i++)
    match->cluster_size[match->receptor_site.atom[i].subst_id]++;
  
  for (i = 0; i < match->cluster_total; i++)
    ecalloc
    (
      (void **) &match->cluster[i],
      match->cluster_size[i],
      sizeof (int),
      "cluster list",
      global.outfile
    );

/*
* Allocate memory for space to reorder site points
* 7/95 te
*/
  temporary.max.atoms = match->receptor_site.total.atoms;
  vstrcpy (&temporary.info.name, "temporary");
  allocate_molecule (&temporary);

/*
* Copy over site points in critical clusters
* 7/95 te
*/
  for (i = temporary.total.atoms = 0; i < match->cluster_total; i++)
    for (j = match->cluster_size[i] = 0;
      j < match->receptor_site.total.atoms; j++)
      if (match->receptor_site.atom[j].subst_id == i)
      {
        copy_atom
          (&temporary.atom[temporary.total.atoms],
          &match->receptor_site.atom[j]);
        copy_coord
          (temporary.coord[temporary.total.atoms],
          match->receptor_site.coord[j]);

        match->cluster[i][match->cluster_size[i]++] = temporary.total.atoms++;
      }

/*
* Check to make the correct number of atoms were copied
* 7/95 te
*/
  if (match->receptor_site.total.atoms != temporary.total.atoms)
    exit (fprintf (global.outfile,
      "ERROR order_site_points: "
      "Critical cluster identifiers must be in the range "
      "[0, nodes_maximum].\n"));

/*
* Copy the reordered site points back into molecular structure
* 7/95 te
*/
  copy_atoms (&match->receptor_site, &temporary);
  free_molecule (&temporary);

/*
* Check to make sure the first cluster is identified as number 1
* 7/95 te
*/
  if (match->receptor_site.atom[0].subst_id != 0)
    exit (fprintf (global.outfile,
      "ERROR order_site_points: "
      "Critical cluster identifiers must begin at 1.\n"));

/*
* Check to make sure the clusters are numbered sequentially with no gaps
* 7/95 te
*/
  for (i = 1; i < match->receptor_site.total.atoms; i++)
    if ((match->receptor_site.atom[i].subst_id !=
      match->receptor_site.atom[i - 1].subst_id) &&
      (match->receptor_site.atom[i].subst_id !=
      match->receptor_site.atom[i - 1].subst_id + 1))
      exit (fprintf (global.outfile,
        "ERROR order_site_points: "
        "Critical cluster identifiers cannot have gaps.\n"));

/*
* Print out critical clusters
* 7/95 te
*/
  fprintf (global.outfile, "____Critical_Clusters____\n");

  fprintf (global.outfile, "%5s %5s %5s\n", "clust", "name", "type");
  for (i = 0; (i < match->receptor_site.total.atoms) &&
    (match->receptor_site.atom[i].subst_id <= match->cluster_total); i++)
    fprintf (global.outfile, "%5d %5s %5s\n",
      match->receptor_site.atom[i].subst_id + 1,
      match->receptor_site.atom[i].name, match->receptor_site.atom[i].type);

  fprintf (global.outfile, "\n\n");
}


/* ////////////////////////////////////////////////////////////// 

Routine to free critical clusters
5/97 te

////////////////////////////////////////////////////////////// */

void free_site_points (MATCH *match)
{
  int i;

/*
* Free space for the cluster lists
* 11/96 te
*/
  for (i = 0; i < match->cluster_total; i++)
    efree ((void **) &match->cluster[i]);

  efree ((void **) &match->cluster);
  efree ((void **) &match->cluster_size);
}



/* ////////////////////////////////////////////////////////////// 

Routine to construct a vector matrix: a N by N matrix in which each
element is a difference vector between the positions of each of N atoms
or site points.  This matrix is constructed to speed up the clique
chirality checking algorithms.
7/95 te

////////////////////////////////////////////////////////////// */

void calculate_vectors (MOLECULE *molecule, XYZ ***matrix, int *size)
{
  int i, j, k;

/*
* Either allocate space for difference vectors, or reset space
* 3/97 te
*/
  if (*size < molecule->total.atoms)
  {
    free_vectors (matrix, size);

    *size = molecule->total.atoms;

    ecalloc
    (
      (void **) matrix,
      *size,
      sizeof (XYZ *),
      "vectors matrix",
      global.outfile
    );

    for (i = 0; i < *size; i++)
      ecalloc
      (
        (void **) &(*matrix)[i],
        *size,
        sizeof (XYZ),
        "vectors matrix",
        global.outfile
      );
  }

  else
    for (i = 0; i < *size; i++)
      memset ((*matrix)[i], 0, *size * sizeof (XYZ));

/*
* Loop through all pairs of atoms/site points
* 3/97 te
*/
  for (i = 0; i < molecule->total.atoms; i++)
    for (j = 0; j < molecule->total.atoms; j++)
      for (k = 0; k < 3; k++)
        (*matrix)[i][j][k] =
          molecule->coord[j][k] -
          molecule->coord[i][k];
}


/* ////////////////////////////////////////////////////////////// */

void free_vectors (XYZ ***matrix, int *size)
{
  int i;

  for (i = 0; i < *size; i++)
    efree ((void **) &(*matrix)[i]);

  efree ((void **) matrix);
  *size = 0;
}


/* ////////////////////////////////////////////////////////////// 

Routine to construct the receptor distance bins used in match_driver.
1/97 te

////////////////////////////////////////////////////////////// */

void make_distance_bins (MATCH *match)
{
  int i, j, k;
  SLINT2 **current = NULL;
  SLINT2 **previous = NULL;

/*
* Allocate space for receptor distance bins
* 1/97 te
*/
  match->bin_width = 0.5;
  match->bin_total = NINT (match->bin_length / match->bin_width) + 1;

  ecalloc
  (
    (void **) &match->bin,
    match->bin_total,
    sizeof (SLINT2 *),
    "receptor distance bins",
    global.outfile
  );

  ecalloc
  (
    (void **) &previous,
    match->bin_total,
    sizeof (SLINT2 *),
    "receptor distance bins",
    global.outfile
  );

  ecalloc
  (
    (void **) &current,
    match->bin_total,
    sizeof (SLINT2 *),
    "receptor distance bins",
    global.outfile
  );

/*
* Divide receptor distances into bins
* 1/97 te
*/
  for (i = 0; i < match->receptor_site.total.atoms; i++)
    for (j = i + 1; j < match->receptor_site.total.atoms; j++)
    {
      k = NINT (match->receptor_distances[i][j] / match->bin_width);

      ecalloc
      (
        (void **) &current[k],
        1,
        sizeof (SLINT2),
        "next bin occupant",
        global.outfile
      );

      current[k]->i = i;
      current[k]->j = j;

      if (match->bin[k] != NULL)
        previous[k]->next = current[k];

      else
        match->bin[k] = current[k];

      previous[k] = current[k];
      current[k] = NULL;
    }

  efree ((void **) &current);
  efree ((void **) &previous);
}


/* ////////////////////////////////////////////////////////////// 

Routine to free the receptor distance bins used in match_driver.
5/97 te

////////////////////////////////////////////////////////////// */

void free_distance_bins (MATCH *match)
{
  int i;
  SLINT2 *prev = NULL;
  SLINT2 *next = NULL;

  for (i = 0; i < match->bin_total; i++)
    for (prev = match->bin[i]; prev != NULL; prev = next)
    {
      next = prev->next;
      efree ((void **) &prev);
    }

  match->bin_total = 0;
  efree ((void **) &match->bin);
}


/* ////////////////////////////////////////////////////////////// 

Routine to extract ligand centers from a ligand
3/96 te

////////////////////////////////////////////////////////////// */

int get_ligand_centers
(
  MATCH		*match,
  MOLECULE	*molecule
)
{
  int atom;
  int segment;
  int layer;

  match->ligand_center.total.atoms = molecule->total.atoms;
  reallocate_atoms (&match->ligand_center);

/*
* Allocate space for the ligand atom key
* 3/97 te
*/
  if (match->ligand_key_total < match->ligand_center.total.atoms)
  {
    efree ((void **) &match->ligand_key);
    match->ligand_key_total = match->ligand_center.total.atoms;

    ecalloc
    (
      (void **) &match->ligand_key,
      match->ligand_key_total,
      sizeof (int),
      "match ligand key",
      global.outfile
    );
  }

/*
* Copy the anchor heavy atoms over
* 3/96 te
*/
  for (atom = match->ligand_center.total.atoms = 0;
    atom < molecule->total.atoms; atom++)
  {
    if
    (
      ((segment = molecule->atom[atom].segment_id) != NEITHER) &&
      (segment < molecule->total.segments) &&
      ((layer = molecule->segment[segment].layer_id) != NEITHER) &&
      (layer < molecule->total.layers) &&
      (layer != 0)
    )
      continue;

    if (molecule->atom[atom].heavy_flag != TRUE)
      continue;

/*
    fprintf (global.outfile, "Lig cent %s, seg %d, layer %d\n",
      molecule->atom[atom].name, segment, layer);

    if (match->ligand_center.total.atoms >= match->ligand_center.max.atoms)
      exit (fprintf (global.outfile,
        "ERROR get_ligand_centers: Too many centers in molecule\n"));
*/

    copy_atom
      (&match->ligand_center.atom[match->ligand_center.total.atoms],
      &molecule->atom[atom]);
    copy_coord
      (match->ligand_center.coord[match->ligand_center.total.atoms],
      molecule->coord[atom]);

    match->ligand_key[match->ligand_center.total.atoms] = atom;
    match->ligand_center.total.atoms++;
  }

  copy_transform (&match->ligand_center, molecule);

/*
* Calculate distance matrix
* 3/97 te
*/
  calculate_distances
  (
    &match->ligand_center,
    &match->ligand_distances,
    &match->ligand_distance_size
  );

/*
* Compute vectors for chirality checking
* 1/97 te
*/
  calculate_vectors
  (
    &match->ligand_center,
    &match->ligand_vectors,
    &match->ligand_vector_size
  );

  return TRUE;
}

/* ////////////////////////////////////////////////////////////// 

Routine to free ligand centers
5/97 te

////////////////////////////////////////////////////////////// */

void free_ligand_centers (MATCH *match)
{
  free_molecule (&match->ligand_center);
  efree ((void **) &match->ligand_key);

  free_distances
  (
    &match->ligand_distances,
    &match->ligand_distance_size
  );

  free_vectors
  (
    &match->ligand_vectors,
    &match->ligand_vector_size
  );
}

/* ////////////////////////////////////////////////////////////// 

Routine to extract ligand keys from a ligand
3/96 te

////////////////////////////////////////////////////////////// */

int get_ligand_keys
(
  MATCH		*match,
  MOLECULE	*molecule
)
{
  int i;

  if (molecule->transform.fold_flag == TRUE)
    exit (fprintf (global.outfile,
      "ERROR get_ligand_keys: ligand keys are folded\n"));

  if (molecule->total.atoms > match->ligand_center.max.atoms)
    exit (fprintf
      (global.outfile, "ERROR get_ligand_keys: too many ligand atoms\n"));

  match->ligand_center.total.atoms = molecule->total.atoms;

  for (i = 0; i < molecule->total.atoms; i++)
    match->ligand_center.atom[i].chem_id = molecule->atom[i].chem_id;

  copy_keys (&match->ligand_center, molecule);

  return TRUE;
}

/* ////////////////////////////////////////////////////////////// 

Routine to make the chemical matching filter.
7/95 te

////////////////////////////////////////////////////////////// */

void make_chemical_filter
(
  LABEL_CHEMICAL *label_chemical,
  MATCH		*match
)
{
  int i, j, node;

  for (i = node = 0; i < match->receptor_site.total.atoms; i++)
    for (j = 0; j < match->ligand_center.total.atoms; j++, node++)
      match->chemical_filter[node] = (char) label_chemical->match_table
        [match->ligand_center.atom[j].chem_id]
        [match->receptor_site.atom[i].chem_id];
}
    

/* ////////////////////////////////////////////////////////////// 

Routine to construct the critical cluster filter used during matching.
7/95

////////////////////////////////////////////////////////////// */

void make_critical_filter (MATCH *match)
{
  int i, j;

/*
* Reset critical filters
* 7/95 te
*/
  for (i = 0; i < match->cluster_total; i++)
    memset (match->critical_filter[i], 0, match->node_max);

  memset (match->critical_filter[match->cluster_total], 1, match->node_max);

/*
* Turn ON selected elements in filters
* 7/95 te
*/
  for (i = 0; i < match->receptor_site.total.atoms; i++)
    for (j = i * match->ligand_center.total.atoms; 
      j < (i + 1) * match->ligand_center.total.atoms; j++)
      match->critical_filter[match->receptor_site.atom[i].subst_id][j] = 1;

/*
* If user wants MULTIPLE selections from each critical cluster,
* then make each screen INCLUDE selection from all preceding clusters
* 7/95 te
*/
  if (match->multiple_flag)
    for (i = 0; i < match->cluster_total - 1; i++)
      for (j = 0; j < match->node_max; j++)
        match->critical_filter[i + 1][j] |=
          match->critical_filter[i][j];
}

/* ////////////////////////////////////////////////////////////// 

Routine to construct the node adjacency matrices used to guide clique
formation during the matching routine.  
7/95 te

////////////////////////////////////////////////////////////// */

void compute_adjacency (MATCH *match)
{
  int li, lj;		/* Ligand atom iterators */
  int ri, rj;		/* Receptor point iterators */
  int inode, jnode;	/* Node values */

  int bin_id;		/* Bin position */
  int bin_tolerance;	/* Tolerance in bin position */
  int bin_center;	/* Center bin to query */
  int bin_initial;	/* Initial bin to query */
  int bin_final;	/* Final bin to query */
  SLINT2 *bin;		/* Bin pointer */

  float residual;	/* Difference in distances */

  int compare_nodes ();

  bin_tolerance = (int)
    (match->cycle * match->distance_tolerance / match->bin_width + .9999);

/*
* Loop through all possible ligand center pairs
* 7/95 te
*/
  for (li = 0; li < match->ligand_center.total.atoms; li++)
    for (lj = li + 1; lj < match->ligand_center.total.atoms; lj++)
    {
/*
*     Check first geometric criteria for adjacency
*     7/95 te
*/
      if (match->ligand_distances[li][lj] < match->distance_minimum)
        continue;

/*
*     Determine which receptor distance bins to check
*     1/97 te
*/
      bin_center =
        NINT (match->ligand_distances[li][lj] / match->bin_width);
      bin_initial = MAX (0, bin_center - bin_tolerance);
      bin_final = MIN (match->bin_total, bin_center + bin_tolerance + 1);

/*
*     Loop through all retrieved site point pairs
*     1/97 te
*/
      for (bin_id = bin_initial; bin_id < bin_final; bin_id++)
        for (bin = match->bin[bin_id]; bin != NULL; bin = bin->next)
        {
/*
*         Extract the site point id
*         1/97 te
*/
          ri = bin->i;
          rj = bin->j;

/*
*         Determine if the two receptor site points can join the two
*         ligand centers to form two adjacent nodes.  Node adjacency is
*         determined by the following geometric criteria:
*
*           1) satisfies *distance_tolerance*
*           2) satisfies *distance_minimum*
*
*         1/97 te
*/
          residual =
            match->ligand_distances[li][lj] -
            match->receptor_distances[ri][rj];

          residual = ABS (residual);

          if ((residual > match->cycle * match->distance_tolerance) ||
            (match->receptor_distances[ri][rj] < match->distance_minimum))
            continue;

/*
*         Compute the first two of four possible nodes and update
*         adjacency matrices.
*         7/95 te
*/
          inode = ri * match->ligand_center.total.atoms + li;
          jnode = rj * match->ligand_center.total.atoms + lj;

          append_adjacency (match, inode, jnode, residual);

          if (match->degeneracy_flag)
            append_adjacency (match, jnode, inode, residual);

/*
*         Compute the second two of four possible nodes and update
*         adjacency matrices.
*         7/95 te
*/
          inode = ri * match->ligand_center.total.atoms + lj;
          jnode = rj * match->ligand_center.total.atoms + li;

          append_adjacency (match, inode, jnode, residual);

          if (match->degeneracy_flag)
            append_adjacency (match, jnode, inode, residual);
        }
    }

/*
* Sort each adjacency list by ascending order
* 7/95 te
*/
  for (li = 0; li < match->node_total; li++)
  {
    if (!match->chemical_flag || match->chemical_filter[li])
      qsort
      (
        match->adjacency_list[li],
        match->adjacency_total[li],
        sizeof (EDGE),
        compare_nodes
      );
  }

/*
  for (li = 0; li < match->node_total; li++)
  {
    fprintf (global.outfile, "node %d:", li);
    for (lj = 0; lj < match->adjacency_total[li]; li++)
      fprintf (global.outfile, " %d", match->adjacency_list[li][lj].node);
    fprintf (global.outfile, "\n");
  }
*/
}


/* ////////////////////////////////////////////////////////////// 

Routine to construct the node adjacency matrices used to guide clique
formation during the matching routine (for chemical screening).  
7/95 te

////////////////////////////////////////////////////////////// */

void compute_chemical_adjacency (MATCH *match, LABEL *label)
{
  int li, lj;		/* Ligand atom iterators */
  int ri, rj;		/* Receptor point iterators */
  int inode, jnode;	/* Node values */
  int lil, ljl;		/* Label id of ligand iterators */
  int ril, rjl;		/* Label id of receptor iterators */

  int compare_nodes ();

  update_uncertainty
  (
    &label->chemical,
    match->distance_tolerance,
    &match->ligand_center
  );

/*
* Loop through all possible ligand center pairs
* 11/96 te
*/
  for (ri = 0; ri < match->receptor_site.total.atoms; ri++)
  {
    ril = match->receptor_site.atom[ri].chem_id;

    for (rj = ri + 1; rj < match->receptor_site.total.atoms; rj++)
    {
      if (match->receptor_distances[ri][rj] < match->distance_minimum)
        continue;

      rjl = match->receptor_site.atom[rj].chem_id;

      for (li = 0; li < match->ligand_center.total.atoms; li++)
      {
        lil = match->ligand_center.atom[li].chem_id;

        if (!label->chemical.match_table[ril][lil])
          continue;

        for (lj = li + 1; lj < match->ligand_center.total.atoms; lj++)
        {
          ljl = match->ligand_center.atom[lj].chem_id;

          if (!label->chemical.match_table[rjl][ljl])
            continue;

          if (!(match->receptor_site.key[ri][rj].distance &
            match->ligand_center.key[li][lj].distance))
            continue;

          inode = ri * match->ligand_center.total.atoms + li;
          jnode = rj * match->ligand_center.total.atoms + lj;

          append_adjacency (match, inode, jnode, 0);
        }
      }
    }
  }

/*
* Sort each adjacency list by ascending order
* 11/96 te
*/
  for (li = 0; li < match->node_total; li++)
  {
    if (match->chemical_filter[li])
      qsort
      (
        match->adjacency_list[li],
        match->adjacency_total[li],
        sizeof (EDGE),
        compare_nodes
      );
  }
}


/* ////////////////////////////////////////////////////////////// 

Routine used by qsort to determine the order of the adjacency lists.
7/95 te

////////////////////////////////////////////////////////////// */

int compare_nodes
(
  EDGE *pta,
  EDGE *ptb
)
{
  return (pta->node == ptb->node) ? 0 : ((pta->node > ptb->node) ? 1 : -1);
}


/* ////////////////////////////////////////////////////////////// 

Routine to add a node to an adjacency list.
11/96 te

////////////////////////////////////////////////////////////// */
void append_adjacency (MATCH *match, int list, int node, float residual)
{
  if (match->adjacency_total[list] + 1 >
    match->adjacency_max[list])
  {
    match->adjacency_max[list] += 10;

    erealloc
    (
      (void **) &match->adjacency_list[list],
      match->adjacency_max[list] * sizeof (EDGE),
      "adjacency list",
      global.outfile
    );
  }

  match->adjacency_list
    [list][match->adjacency_total[list]].node = node;
  match->adjacency_list
    [list][match->adjacency_total[list]].residual = residual;
  match->adjacency_total[list]++;
}

/* ////////////////////////////////////////////////////////////// 

Routine to assemble all matches for a given distance tolerance

Return values:
	TRUE	matches recorded
	EOF	insufficient matches recorded (terminate orienting)

1/97 te

////////////////////////////////////////////////////////////// */

int compute_matches (MATCH *match)
{
  int i;
  int previous_max;		/* Record of the number of previous matches */
  SLCLIQUE *current = NULL;	/* Pointer to current clique */
  SLCLIQUE *previous = NULL;	/* Pointer to previous clique */

  int compare_cliques ();

/*
* Loop through all matches
* 1/97 te
*/
  for
  (
    previous_max = match->max, match->max = 0;
    match_distance (match) != EOF;
    match->max++
  )
  {
    ecalloc
    (
      (void **) &current,
      1,
      sizeof (SLCLIQUE),
      "clique link list",
      global.outfile
    );

    copy_clique (&current->clique, &match->clique);

    if (match->clique_link != NULL)
      previous->next = current;

    else
      match->clique_link = current;

    previous = current;
    current = NULL;
  }

/*
* Copy cliques over to linear array and sort them
* 1/97 te
*/
  if (match->clique_sort_total < match->max)
  {
    free_clique_list (match);
    match->clique_sort_total = match->max;

    ecalloc
    (
      (void **) &match->clique_sort,
      match->clique_sort_total,
      sizeof (CLIQUE),
      "clique sort list",
      global.outfile
    );
  }

  for
  (
    i = 0, current = match->clique_link;
    (i < match->max) && (current != NULL);
    i++, current = current->next
  )
    copy_clique (&match->clique_sort[i], &current->clique);

  if (i != match->max)
    exit (fprintf (global.outfile,
      "ERROR compute_matches: Insufficient cliques copied\n"));

  qsort
  (
    match->clique_sort,
    match->max,
    sizeof (CLIQUE),
    compare_cliques
  );

  free_clique_link (match);

/*
* Check if sufficient matches were made
* 1/97 te
*/
  if
  (
    (match->cycle >= 10) &&
    ((match->max == 0) ||
    (match->max <= previous_max))
  )
    return EOF;

  else
    return TRUE;
}


/* //////////////////////////////////////////////////////////////////// */

void free_clique_list (MATCH *match)
{
  int i;

  for (i = 0; i < match->clique_sort_total; i++)
    free_clique (&match->clique_sort[i]);

  efree ((void **) &match->clique_sort);
  match->clique_sort_total = 0;
}


/* //////////////////////////////////////////////////////////////////// */

void free_clique_link (MATCH *match)
{
  SLCLIQUE *current = NULL;	/* Pointer to current clique */
  SLCLIQUE *previous = NULL;	/* Pointer to previous clique */

/*
* Free up linked list
* 3/97 te
*/
  for
  (
    previous = match->clique_link;
    previous != NULL;
    previous = current
  )
  {
    current = previous->next;
    free_clique (&previous->clique);
    efree ((void **) &previous);
  }

  match->clique_link = NULL;
}

int compare_cliques
(
  CLIQUE *pta,
  CLIQUE *ptb
)
{
  return (pta->residual == ptb->residual) ? 0 :
    ((pta->residual > ptb->residual) ? 1 : -1);
}


/* ////////////////////////////////////////////////////////////// 

Routine to generate ligand orientations by finding subsets of
ligand atoms to superimpose on receptor site points.

1.  At start, no nodes are in clique, and all nodes are available
    in clique adjacency list.

2.  Successive nodes are taken from the initial clique adjacency list
    to seed the clique.

3.  A new adjacency list is formed by the intersection of the previous
    adjacency list and the new node adjacency list.

4.  The clique is grown in a depth-first fashion until no more nodes
    are available in the adjacency list.  If the clique is large enough,
    it is returned as a match.

5.  Backtracking is implemented by successively removing nodes from
    the clique to try alternative nodes from the previous adjacency list.

6.  When backtracking has exhausted all nodes in the clique and all
    nodes in the initial clique adjacency list, then matching terminates.

Return values:
	TRUE:	suitable match found
	EOF:	no more matches can be formed

11/96 te

////////////////////////////////////////////////////////////// */

int match_distance (MATCH *match)
{
  int i, j;			/* Counter variable */
  EDGE edge;			/* Edge under consideration */
  EDGE iedge, jedge;		/* Edge comparisons */
  int valid_clique;		/* Flag for whether an valid_clique found */

  int match_chirality (MATCH *, int);

/*
* If a prior clique existed, then erase the last node and
* decrement the clique size
* 11/96 te
*/
  if (match->clique.edge_total > 0)
  {
    match->clique.edge_total--;
    match->clique.edge[match->clique.edge_total].node = 0;
    match->clique.edge[match->clique.edge_total].residual = 0;
  }

/*
* Perform clique detection until valid clique is found
* 11/96 te
*/
  for (valid_clique = FALSE; !valid_clique;)
  {
/*
*   Can we expand the current clique?  Check if:
*   1. The clique is smaller than the maximum size
*   2. There are unused nodes in the clique adjacency list
*   11/96 te
*/
    if
    (
      (match->clique.edge_total < match->clique_size_max) &&
      (match->adjacent[match->clique.edge_total] <
        match->clique_adjacency_total[match->clique.edge_total])
    )
    {
/*
*     Identify the next adjacent node and increment the adjacency list
*     11/96 te
*/
      edge = match->clique_adjacency_list
        [match->clique.edge_total]
        [match->adjacent[match->clique.edge_total]++];

/*
*     STAGE 1: The node is geometrically feasible, if
*     1. The node passes the critical cluster filter:
*        a. The filter is not requested, OR
*        b. The filter indicates valid node, AND
*     2. The clique has proper chirality:
*        a. Chirality checking is not necessary, OR
*        b. The node is not being added to a three-node clique, OR
*        c. The chirality of the new clique is correct, OR
*        d. The ligand is allowed to be reflected
*     11/96 te
*/
      if
      (
        (
          !match->critical_flag ||
          (((match->clique.edge_total == 0) &&
            match->critical_filter[0][edge.node]) ||
          ((match->clique.edge_total > 0) &&
            match->critical_filter
              [match->receptor_site.atom
                [match->clique.edge[match->clique.edge_total - 1].node /
                  match->ligand_center.total.atoms].subst_id + 1][edge.node]))
        ) &&
        (
          !match->chiral_flag ||
          (match->clique.edge_total != 3) ||
          !(match->clique.reflect_flag = match_chirality (match, edge.node)) ||
          match->reflect_flag
        )
      )
      {
/*
*       STAGE 2: The node is worth adding, if:
*       1. The current node isn't forming a redundant clique:
*          a. We don't care, OR
*          b. It the first in the clique, OR
*          c. The node is numerically above the previously added node, AND
*       2. The current node passes the chemical filter
*          a. The filter is not requested, OR
*          a. The filter indicates a valid node
*       11/96 te
*/
        if
        (
          (
            !match->degeneracy_flag ||
            !match->clique.edge_total ||
            (edge.node > match->clique.edge[match->clique.edge_total - 1].node)
          ) &&
          (
            !match->chemical_flag ||
            match->chemical_filter[edge.node]
          )
        )
        {
/*
*         Record the node in the clique and increment the clique size.
*         Also, document the finding of a clique of this size,
*         and erase the record of the next larger one.
*         11/96 te
*/
          match->clique.edge[match->clique.edge_total++] = edge;

          if (match->degeneracy_flag)
          {
            match->degenerate[match->clique.edge_total - 1] = TRUE;
            match->degenerate[match->clique.edge_total] = FALSE;
          }

/*
*         Construct an adjacency list for the new clique
*         11/96 te
*/
          for
          (
            i = j = match->adjacent[match->clique.edge_total] = 
              match->clique_adjacency_total[match->clique.edge_total] = 0;
            (i < match->clique_adjacency_total[match->clique.edge_total - 1]) &&
              (j < match->adjacency_total[edge.node]);
          )
          {
            iedge = match->clique_adjacency_list
              [match->clique.edge_total - 1][i];
            jedge = match->adjacency_list[edge.node][j];

            if (jedge.node < iedge.node)
              j++;

            else if (iedge.node < jedge.node)
              i++;

            else
            {
              match->clique_adjacency_list
                [match->clique.edge_total]
                [match->clique_adjacency_total[match->clique.edge_total]++] =
                (iedge.residual > jedge.residual) ? iedge : jedge;

              i++, j++;
            }
          }
        }

/*
*       If this node fails STAGE2 check, then for the sake of degeneracy
*       checking, record that this clique *could* have existed
*       11/96 te
*/
        else if (match->degeneracy_flag)
          match->degenerate[match->clique.edge_total] = TRUE;
      }
    }

/*
*   If we cannot expand the current clique, then process it.
*   10/96 te
*/
    else if (match->clique.edge_total > 0)
    {
/*
*     Compute residual of clique and check its validity
*     1/97 te
*/
      for (i = match->clique.residual = 0; i < match->clique.edge_total; i++)
        match->clique.residual =
          MAX (match->clique.residual, match->clique.edge[i].residual);

/*
*     This clique is sufficient for generating an orientation, if:
*     1. The clique is large enough, AND
*     2. The clique passes the degeneracy filter
*        a. The filter is not requested, OR
*        b. The filter indicates the cliqe is not degenerate
*     11/96 te
*/
      if
      (
        (match->clique.edge_total >= match->clique_size_min) &&
        (!match->degeneracy_flag ||
          !match->degenerate[match->clique.edge_total]) &&
        (match->clique.residual >
          ((match->cycle - 1) * match->distance_tolerance))
      )
        valid_clique = TRUE;

/*
*     Otherwise, erase the last node and decrement the clique size
*     11/96 te
*/
      else
      {
        match->clique.edge_total--;
        match->clique.edge[match->clique.edge_total].node = 0;
        match->clique.edge[match->clique.edge_total].residual = 0;
      }
    }

    else
      return EOF;
  }

/*
* Update match statistics
* 11/96 te
*/
  match->histogram[match->clique.edge_total - 1]++;

/*
* Reset ligand reflection state
* 1/97 te
*/
  if (match->clique.edge_total < 4)
    match->clique.reflect_flag = FALSE;

  return TRUE;
}


/* ////////////////////////////////////////////////////////////// 

Routine to check if the ligand and receptor portions of a clique
have the same or opposite chirality.

Return values:
	TRUE:	The ligand must be reflected
	FALSE:	The ligand must not be reflected
7/95 te

////////////////////////////////////////////////////////////// */

int match_chirality (MATCH *match, int node)
{
  int i;
  static int atoms[4], sites[4];

  float calc_chirality (int *, XYZ **);
  
  for (i = 0; i < 3; i++)
  {
    atoms[i] = match->clique.edge[i].node % match->ligand_center.total.atoms;
    sites[i] = match->clique.edge[i].node / match->ligand_center.total.atoms;
  }

  atoms[3] = node % match->ligand_center.total.atoms;
  sites[3] = node / match->ligand_center.total.atoms;

  if (calc_chirality (atoms, match->ligand_vectors) *
    calc_chirality (sites, match->receptor_vectors) < 0.0)
    return TRUE;

  else
    return FALSE;
}

/* ////////////////////////////////////////////////////////////// 

Routine to calculate the whether the fourth member is above or
below the plane defined by the first three members of a group.
7/95 te

////////////////////////////////////////////////////////////// */

float calc_chirality (int atoms[4], XYZ **vector)
{
  static XYZ normal;
  void vcrossv3 (XYZ, XYZ, XYZ);
  float vdotv3 (XYZ, XYZ);

  vcrossv3 (vector[atoms[0]][atoms[1]], vector[atoms[0]][atoms[2]], normal);
  return vdotv3 (vector[atoms[0]][atoms[3]], normal);
}



/* ////////////////////////////////////////////////////////////// 

Routine to output match routine statistics
7/95 te

////////////////////////////////////////////////////////////// */

void output_match_info (MATCH *match)
{
  int i;

/*
* Output histogram of clique formation
* 6/95 te
*/
  if ((global.output_volume != 't') && (match->total > 0))
  {
    fprintf (global.outfile, "\n_______Matching_Histogram______\n");

    for (i = 0; i < match->clique_size_max; i++)
      if (match->histogram[i] > 0)
        fprintf (global.outfile, "%7d node%s: %9d\n",
          i + 1, (i == 0 ? " " : "s"),
          match->histogram[i]);

    fprintf (global.outfile, "\n");
  }
}

/* ///////////////////////////////////////////////////////////// */

void extract_clique (MATCH *match, MOLECULE *molecule)
{
  int i, lig, rec;

  molecule->transform.refl_flag = match->clique.reflect_flag;

  for (i = 0; i < match->clique.edge_total; i++)
  {
    rec = match->clique.edge[i].node / match->ligand_center.total.atoms;

    copy_atom
      (&match->receptor_clique.atom[i],
      &match->receptor_site.atom[rec]);
    copy_coord
      (match->receptor_clique.coord[i],
      match->receptor_site.coord[rec]);
  }

  for (i = 0; i < match->clique.edge_total; i++)
  {
    lig = match->clique.edge[i].node % match->ligand_center.total.atoms;

    copy_atom
      (&match->ligand_clique.atom[i],
      &match->ligand_center.atom[lig]);
    copy_coord
      (match->ligand_clique.coord[i],
      match->ligand_center.coord[lig]);
  }

  match->ligand_clique.total.atoms =
    match->receptor_clique.total.atoms = match->clique.edge_total;
}


/* ////////////////////////////////////////////////////////////// 

Routine to construct a random match
3/97 te

////////////////////////////////////////////////////////////// */

int get_random_match
(
  MATCH		*match,
  LABEL		*label,
  MOLECULE      *mol_conf,
  int           molecule_id,
  int           conformation_id,
  int           orientation_id
)
{
  int i, j;
  int unique;

  if (orientation_id == 0)
  {
    if ((molecule_id == 0) && (conformation_id == 0))
    {
      get_site (match, label);
      get_centers (match, label);

      match->clique.edge_max = match->clique_size_min;
      allocate_clique (&match->clique);
      match->clique.edge_total = match->clique_size_min;

      allocate_clique_atoms (match);
    }
    
    if (!match->centers_flag)
      get_ligand_centers (match, mol_conf);

    if (conformation_id == 0)
      match->node_total =
        match->ligand_center.total.atoms * match->receptor_site.total.atoms;
  }

/*
* Assign nodes to clique
* 3/97 te
*/
  for (i = 0; i < match->clique.edge_total;)
  {
    match->clique.edge[i].node = rand () % match->node_total;

    for (j = 0, unique = TRUE; j < i; j++)
      if ((match->clique.edge[i].node % match->ligand_center.total.atoms ==
        match->clique.edge[j].node % match->ligand_center.total.atoms) ||
        (match->clique.edge[i].node / match->ligand_center.total.atoms ==
        match->clique.edge[j].node / match->ligand_center.total.atoms))
        unique = FALSE;

    if (unique)
      i++;
  }

  return TRUE;
}


/* ////////////////////////////////////////////////////////////// 

Routine to free random matching arrays
5/97 te

////////////////////////////////////////////////////////////// */

void free_random_matches (MATCH *match)
{
  free_clique (&match->clique);
  free_ligand_centers (match);
  free_molecule (&match->receptor_site);
  free_molecule (&match->ligand_clique);
  free_molecule (&match->receptor_clique);
}


/* ////////////////////////////////////////////////////////////// 

Routine to manage chemical screening
3/96 te

////////////////////////////////////////////////////////////// */

int check_screen
(
  MATCH		*match,
  LABEL		*label,
  MOLECULE      *molecule,
  int           molecule_id
)
{
  static MOLECULE temporary = {0};

  match->chiral_flag = FALSE;

  if (label->chemical.screen.pharmaco_flag)
  {
    if (molecule_id == 0)
    {
      init_match_site (match, label);
      update_keys (&label->chemical, &match->receptor_site, 0);
    }

    if (label->chemical.screen.fold_flag)
      return
        check_pharmacophore
        (
          &label->chemical,
          match->distance_tolerance,
          &match->receptor_site,
          molecule
        );

    else
    {
      get_ligand_keys (match, molecule);
      init_match_ligand (match, label);
      initialize_adjacency (match, label);

      if (match_distance (match) != EOF)
        return TRUE;

      else
        return FALSE;
    }
  }

  else if (label->chemical.screen.similar_flag)
  {
    if (molecule_id == 0)
      init_similar_site (match, label);

    molecule->score.total =
      molecule->score.inter.total =
      check_dissimilarity (&label->chemical, &match->receptor_site, molecule);

    if (molecule->score.total <=
      label->chemical.screen.dissimilar_maximum)
      return TRUE;

    else
      return FALSE;
  }

  else if (label->chemical.screen.fold_flag)
  {
    if (!label->chemical.init_flag)
      get_chemical_labels (&label->chemical);

    fold_keys (&label->chemical, &temporary, molecule);
    copy_keys (molecule, &temporary);
    return TRUE;
  }

  else
    return FALSE;
}



/* ////////////////////////////////////////////////////////////// 

Routine to free chemical screening arrays
5/97 te

////////////////////////////////////////////////////////////// */

void free_screen
(
  MATCH		*match,
  LABEL		*label
)
{
  if (label->chemical.screen.pharmaco_flag)
  {
    free_match_site (match, label);
/*
    if (label->chemical.screen.fold_flag)
      return
        check_pharmacophore
        (
          &label->chemical,
          match->distance_tolerance,
          &match->receptor_site,
          molecule
        );

    else
    {
      get_ligand_keys (match, molecule);
      init_match_ligand (match, label);
      initialize_adjacency (match, label);

      if (match_distance (match) != EOF)
        return TRUE;

      else
        return FALSE;
    }
*/
  }


  else if (label->chemical.screen.similar_flag)
  {
    free_molecule (&match->receptor_site);
    free_table (&label->chemical, &label->chemical.screen_table);
  }

/*
  else if (label->chemical.screen.fold_flag)
  {
    if (!label->chemical.init_flag)
      get_chemical_labels (&label->chemical);

    fold_keys (&label->chemical, &temporary, molecule);
    copy_keys (molecule, &temporary);
    return TRUE;
  }
*/
}


/* ////////////////////////////////////////////////////////////// 

Routine to prepare for similarity search
3/96 te

////////////////////////////////////////////////////////////// */

void init_similar_site
(
  MATCH *match,
  LABEL *label
)
{
  FILE *file;

  get_chemical_labels (&label->chemical);

/*
* Read in reference molecule for similarity search
* 11/96 te
*/
  file = efopen (match->receptor_file_name, "r", global.outfile);

  if (read_molecule
    (&match->receptor_site, NULL, match->receptor_file_name, file, FALSE)
    != TRUE)
    exit (fprintf (global.outfile,
      "ERROR init_similar_site: Error reading in reference molecule.\n"));

  fclose (file);
}



