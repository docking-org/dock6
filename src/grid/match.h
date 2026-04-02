/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
3/96
*/

typedef struct distance_list_struct
{
  float distance;
  int i, j;

} DISTANCE_LIST;

typedef struct edge_struct
{
  float residual;
  int node;

} EDGE;

typedef struct clique_struct
{
  EDGE *edge;
  int edge_total;
  int edge_max;
  float residual;
  int reflect_flag;

} CLIQUE;

typedef struct singly_linked_clique_struct
{
  CLIQUE clique;
  struct singly_linked_clique_struct *next;

} SLCLIQUE;

typedef struct match_struct
{
  int flag;			/* Flag for matching */
  int auto_flag;		/* Flag for automated matching */
  int random_flag;		/* Flag for random matching */
  int centers_flag;		/* Flag for ligand centers from file */
  int chiral_flag;		/* Flag to check clique chirality */
  int reflect_flag;		/* Flag to reflect ligand */
  int degeneracy_flag;		/* Flag for checking if clique is subclique */
  int critical_flag;		/* Flag for critical spheres to prune search */
  int multiple_flag;		/* Flag for multiple selections from cluster */
  int chemical_flag;		/* Flag for chemical labels to prune search */

  int clique_size_min;		/* Minimum nodes in match for orient */
  int clique_size_max;		/* Maximum nodes in match for orient */
  float distance_minimum;	/* Minimum internal distance to match */
  float distance_tolerance;	/* Difference in internal distances to match */

  int total;			/* Number of matches */
  int max;			/* Maximum number of matches */

  FILE_NAME ligand_file_name;	/* File containing ligand centers */
  FILE_NAME receptor_file_name;	/* File containing receptor site points */

  float **receptor_distances;	/* Receptor distance matrix */
  int receptor_distance_size;	/* Receptor distance matrix size */

  float **ligand_distances;	/* Ligand distance matrix */
  int ligand_distance_size;	/* Ligand distance matrix size */

  XYZ **receptor_vectors;	/* Receptor vector matrix (for chirality) */
  int receptor_vector_size;	/* Receptor vector matrix size */

  XYZ **ligand_vectors;		/* Ligand vector matrix */
  int ligand_vector_size;	/* Ligand vector matrix size */

  SLINT2 **bin;			/* Bin array of receptor distances */
  float bin_length;		/* Longest distance in bin array */
  float bin_width;		/* Width of each bin */
  int bin_total;		/* Number of bins */

  int node_max;			/* Maximum number of nodes */
  int node_total;		/* Current number of nodes */

  int init_flag;		/* Flag for whether match info initialized */
  CLIQUE clique;		/* Current clique */
  int *adjacent;		/* Node under consideration in adjacency list */
  int *degenerate;		/* Record of whether clique is degenerate */
  SLCLIQUE *clique_link;	/* Linked list of all cliques generated */
  CLIQUE *clique_sort;		/* Truncated, sorted list of cliques */
  int clique_sort_total;	/* Length of list of sorted cliques */
  int cycle;			/* Number of times matching performed */

  EDGE **adjacency_list;	/* Adjacency lists */
  int *adjacency_total;		/* Length of each adjacency list */
  int *adjacency_max;		/* Maximum length of each adjacency list */
  EDGE **clique_adjacency_list;	/* Adjacency lists for cliques */
  int *clique_adjacency_total;	/* Length of clique adjacency lists */

  char *chemical_filter;	/* Array encoding chemical nodes */
  char **critical_filter;	/* Array encoding critical clusters */
  int cluster_total;		/* Number of critical clusters */
  int **cluster;		/* Members of each cluster */
  int *cluster_size;		/* Number of members in each cluster */
  int *histogram;		/* Histogram of clique populations */

  int *ligand_key;		/* List of which atoms used as centers */
  int ligand_key_total;		/* Length of list of atoms as centers */
  MOLECULE ligand_center;	/* Ligand centers for matching */
  MOLECULE receptor_site;	/* Receptor site points for matching */

  MOLECULE ligand_clique;	/* Ligand centers in clique */
  MOLECULE receptor_clique;	/* Receptor site points in clique */

} MATCH;

/*
Routines defined in match.c, that are called by outside functions
*/

int get_match
(
  DOCK          *dock,
  MATCH         *match,
  LABEL         *label,
  MOLECULE      *mol_conf,
  int           conformation_id,
  int           orientation_id
);

int get_random_match
(
  MATCH         *match,
  LABEL         *label,
  MOLECULE      *mol_conf,
  int           molecule_id,
  int           conformation_id,
  int           orientation_id
);

void init_match_site
(
  MATCH         *match,
  LABEL         *label
);

void free_matches
(
  LABEL         *label,
  MATCH         *match
);

void free_random_matches (MATCH *match);

void free_match_site
(
  MATCH         *match,
  LABEL         *label
);

int init_match_ligand
(
  MATCH         *match,
  LABEL         *label
);

int initialize_adjacency
(
  MATCH         *match,
  LABEL         *label
);

void calculate_vectors (MOLECULE *molecule, XYZ ***matrix, int *size);
void free_vectors (XYZ ***matrix, int *size);

void get_site
(
  MATCH         *match,
  LABEL         *label
);

void get_centers
(
  MATCH         *match,
  LABEL         *label
);

void id_site_points (MATCH *match);
void order_site_points (MATCH *match);
void free_site_points (MATCH *match);
void make_distance_bins (MATCH *match);
void free_distance_bins (MATCH *match);
void allocate_clique_atoms (MATCH *match);

void allocate_match (MATCH *match);
void reallocate_match (MATCH *match);
void free_match (MATCH *match);
void reset_match (MATCH *match);

void allocate_clique (CLIQUE *clique);
void reallocate_clique (CLIQUE *clique);
void free_clique (CLIQUE *clique);
void reset_clique (CLIQUE *clique);
void copy_clique (CLIQUE *copy, CLIQUE *original);

int get_ligand_centers
(
  MATCH         *match,
  MOLECULE      *molecule
);

void free_ligand_centers (MATCH *match);

int get_ligand_keys
(
  MATCH         *match,
  MOLECULE      *molecule
);

int screen_match
(
  DOCK          *dock,
  MATCH         *match,
  LABEL         *label,
  MOLECULE      *molecule
);

void make_chemical_filter
(
  LABEL_CHEMICAL *label_chemical,
  MATCH         *match
);

void compute_adjacency (MATCH *match);
void compute_chemical_adjacency (MATCH *match, LABEL *label);
void append_adjacency (MATCH *match, int list, int node, float residual);
void reset_adjacency (MATCH *match);
void free_adjacency (MATCH *match);
void reset_match (MATCH *match);
void free_match (MATCH *match);
void make_critical_filter (MATCH *match);

int compute_matches (MATCH *);
void free_clique_list (MATCH *);
void free_clique_link (MATCH *);
int match_distance (MATCH *);

void output_match_info (MATCH *match);
void extract_clique (MATCH *match, MOLECULE *molecule);

int check_screen
(
  MATCH         *match,
  LABEL         *label,
  MOLECULE      *mol_init,
  int           molecule_id
);

void free_screen
(
  MATCH         *match,
  LABEL         *label
);

void init_similar_site
(
  MATCH         *match,
  LABEL         *label
);

