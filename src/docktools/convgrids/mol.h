/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
5/97
*/

/*
Structures to hold molecule data
5/97 te
*/

typedef struct info_struct
{
  int	input_id;		/* id from input file */
  int	output_id;		/* id for output file */
  char	*name;			/* name string */
  char	*comment;		/* comment string */
  char	*molecule_type;		/* type string (from sybyl) */
  char	*charge_type;		/* charge string (from sybyl) */
  char	*status_bits;		/* status string (from sybyl) */
  char	*source_file;		/* name of source file (ptr format) */
  long	source_position;	/* file position in source (ptr format) */
  long	file_position;		/* position in input file */

  int	allocated;		/* flag for array allocation */
  int	assign_chem;		/* flag for chemical label assignment */
  int	assign_vdw;		/* flag for vdw label assignment */
  int	assign_flex;		/* flag for flex label assignment */

} INFO;

typedef struct transform_struct
{
  int	flag;			/* flag for transformation */
  int	trans_flag;		/* flag for rigid body translation */
  int	rot_flag;		/* flag for rigid body rotation */
  int	refl_flag;		/* flag for chiral reflection */
  int	tors_flag;		/* flag for bond rotation */

  XYZ	translate;		/* translation from reference */
  XYZ	rotate;			/* rotation from reference */
  int	conf_total;		/* Number of conformations */
  int	anchor_segment;		/* Anchor segment */
  int	anchor_atom;		/* Anchor atom */
  int	active_layer;		/* Current active layer */
  int	fold_flag;		/* flag for chemical key folding */

  XYZ	com;			/* center of mass of current position */
  float	rmsd;			/* RMS deviation from reference */
  int	heavy_total;		/* Number of non-hydrogen atoms */

} TRANSFORM;

typedef struct score_part_struct
{
  int	current_flag;		/* Flag for whether score up-to-date */
  float electro;		/* Electrostatic component of score */
  float vdw;			/* Van der Waals component of score */
  float total;			/* Total score */

} SCORE_PART;

typedef struct mol_score_struct
{
  int	type;			/* type of score method */
  int	bumpcount;		/* number of bumps */
  SCORE_PART	intra;		/* intramolecular score */
  SCORE_PART	inter;		/* intermolecular score */
  float	total;			/* total score */

} MOL_SCORE;

typedef struct arraysize_struct
{
  int	atoms;			/* Number of atoms */
  int	bonds;			/* Number of bonds */
  int	substs;			/* Number of substructures */
  int	torsions;		/* Number of torsions */
  int	sets;			/* Number of sets */
  int	segments;		/* Number of segments */
  int	layers;			/* Number of layers */
  int	keys;			/* Number of keys (along one dimension) */

} ARRAYSIZE;

typedef struct link_struct
{
  int	id;			/* Link id */
  int	bond_id;		/* Bond id */
  int	out_flag;		/* Flag for whether outward or inward link */

} LINK;

typedef struct atom_struct
{
  int	number;			/* Atom number */
  int	subst_id;		/* Substructure id */
  int	segment_id;		/* Molecular segment id */
  int	chem_id;		/* Chemical label id */
  int	vdw_id;			/* VDW id */
  int	heavy_flag;		/* Heavy atom flag */
  int	centrality;		/* Proximity of most distant atom */
  int	flag;			/* General purpose flag */
  float	charge;			/* Atomic partial charge */
  char	*name;			/* Atom name */
  char	*type;			/* Atom type */

  LINK	*neighbor;		/* Atom neighbor list */
  int	neighbor_total;		/* Number of neighbors */
  int	neighbor_max;		/* Maximum number of neighbors */

} ATOM;

typedef struct bond_struct
{
  int	id;			/* Bond id */
  int	origin;			/* Atom at bond origin */
  int	target;			/* Atom at bond target */
  int	ring_flag;		/* Molecular ring id */
  int	flex_id;		/* Flexible label id */
  char	*type;			/* Bond type */

} BOND;

typedef struct subst_struct
{
  int	number;			/* Substructure number */
  int	root_atom;		/* First atom in substructure */
  int	dict_type;		/* Dictionary type (sybyl) */
  int	inter_bonds;		/* Number of inter-substructure bonds */
  char	*name;			/* Substructure name */
  char	*type;			/* Substructure type */
  char	*chain;			/* Substructure chain identifier */
  char	*sub_type;		/* Sub-type */
  char	*status;		/* Status string (sybyl) */

} SUBST;

typedef struct torsion_struct
{
  int	flex_id;		/* Flexible label id */
  int	bond_id;		/* Bond id */
  int	segment_id;		/* Segment id */
  int	flag;			/* General purpose flag */
  int	periph_flag;		/* Flag for peripheral torsion (eg. hydroxyl) */
  int	reverse_flag;		/* Flag for reversal of origin/target */
  int	origin;			/* Atom at bond origin */
  int	target;			/* Atom at bond target */
  int	origin_neighbor;	/* Neighbor to atom at bond origin */
  int	target_neighbor;	/* Neighbor to atom at bond origin */
  float	current_angle;		/* Current torsion angle */
  float	target_angle;		/* Target torsion angle */

} TORSION;

typedef struct set_struct
{
  char	*name;			/* Set name */
  char	*type;			/* Set type */
  char	*obj_type;		/* Set object type (ATOMS, BONDS, SUBSTS) */
  char	*sub_type;		/* Set subtype (Sybyl) */
  char	*status;		/* Set statis bits (Sybyl) */
  char	*comment;		/* Set comment */
  int	member_total;		/* Number of members in the set */
  int	*member;		/* Set members */

} SET;

typedef struct key_struct
{
  int	count;			/* Number of distances stored (folded) */
  MASK	distance;		/* Distance keys between labels */

} KEY;

typedef struct segment_struct
{
  int	id;			/* Segment id */
  int	torsion_id;		/* Torsion in this segment */
  int	layer_id;		/* Layer to which segment belongs */
  int	periph_flag;		/* Flag for peripheral segment */
  int	heavy_total;		/* Number of heavy atoms */
  int	active_flag;		/* Flag for whether to score, etc. */
  int	min_flag;		/* Flag for whether to minimize */
  int	conform_total;		/* Number of conformations in segment */
  int	conform_count;		/* Number of current conformation */
  int	conform_seed;		/* Random seed for conformation search */

  MOL_SCORE	score;		/* Score for segment */

  int	*atom;			/* Atoms in this segment */
  int	atom_max;		/* Maximum number of atoms in segment */
  int	atom_total;		/* Number of atoms in segment */

  LINK	*neighbor;		/* Neighboring segments */
  int	neighbor_max;		/* Maximum number of neighboring segments */
  int	neighbor_total;		/* Number of neighboring segments */

} SEGMENT;

typedef struct layer_struct
{
  int	id;			/* Layer id */
  int	active_flag;		/* Flag for whether to score, etc. */
  int	conform_total;		/* Total number of conformations for layer */
  MOL_SCORE	score;		/* Score for layer */

  int	*segment;		/* Segments in this layer */
  int	segment_max;		/* Maximum number of segments in layer */
  int	segment_total;		/* Number of segments in layer */

} LAYER;


typedef struct molecule_struct
{
  INFO		info;		/* General info */
  TRANSFORM	transform;	/* Transformation info */
  MOL_SCORE	score;		/* Molecule score info */
  ARRAYSIZE	total;		/* Current size of arrays */
  ARRAYSIZE	max;		/* Allocated size of arrays */

  ATOM		*atom;		/* Atom info */
  XYZ		*coord;		/* Atom coordinates */
  BOND		*bond;		/* Bond info */
  SUBST		*subst;		/* Sybyl substructure info */
  SET		*set;		/* Sybyl set info */
  TORSION	*torsion;	/* Rotatable torsion info */
  SEGMENT	*segment;	/* Rigid segment info */
  LAYER		*layer;		/* Segment layer info */
  KEY		**key;		/* Chemical keys */

} MOLECULE;

typedef struct linked_molecule
{
  struct linked_molecule *next_head, *next_member;
  MOLECULE *molecule;

} LINKED_MOLECULE;


/*
Routines used to manipulate molecule data structures
12/96 te
*/

void allocate_molecule			(MOLECULE *);
void allocate_info			(MOLECULE *);
void allocate_transform			(MOLECULE *);
void allocate_score			(MOL_SCORE *);
void allocate_atoms			(MOLECULE *);
void allocate_atom_neighbors		(ATOM     *);
void allocate_bonds			(MOLECULE *);
void allocate_torsions			(MOLECULE *);
void allocate_substs			(MOLECULE *);
void allocate_segments			(MOLECULE *);
void allocate_segment_atoms		(SEGMENT  *);
void allocate_segment_neighbors		(SEGMENT  *);
void allocate_layers			(MOLECULE *);
void allocate_layer_torsions		(LAYER    *);
void allocate_layer_segments		(LAYER    *);
void allocate_sets			(MOLECULE *);
void allocate_keys			(MOLECULE *);

void reset_molecule			(MOLECULE *);
void reset_info				(MOLECULE *);
void reset_transform			(MOLECULE *);
void reset_score			(MOL_SCORE *);
void reset_score_parts			(SCORE_PART *);
void reset_atoms			(MOLECULE *);
void reset_atom				(ATOM     *);
void reset_atom_neighbors		(ATOM     *);
void reset_bonds			(MOLECULE *);
void reset_bond				(BOND     *);
void reset_torsions			(MOLECULE *);
void reset_torsion			(TORSION  *);
void reset_substs			(MOLECULE *);
void reset_subst			(SUBST    *);
void reset_segments			(MOLECULE *);
void reset_segment			(SEGMENT  *);
void reset_segment_neighbors		(SEGMENT  *);
void reset_layers			(MOLECULE *);
void reset_layer			(MOLECULE *, int);
void reset_sets				(MOLECULE *);
void reset_set				(SET      *);
void reset_keys				(MOLECULE *);

void free_molecule			(MOLECULE *);
void free_info				(MOLECULE *);
void free_atoms				(MOLECULE *);
void free_atom				(ATOM     *);
void free_atom_neighbors		(ATOM     *);
void free_bonds				(MOLECULE *);
void free_bond				(BOND     *);
void free_torsions			(MOLECULE *);
void free_substs			(MOLECULE *);
void free_subst				(SUBST    *);
void free_segments			(MOLECULE *);
void free_segment_atoms			(SEGMENT  *);
void free_segment_neighbors		(SEGMENT  *);
void free_layers			(MOLECULE *);
void free_layer_segments		(LAYER    *);
void free_sets				(MOLECULE *);
void free_set				(SET      *);
void free_keys				(MOLECULE *);

void reallocate_molecule		(MOLECULE *);
void reallocate_atoms			(MOLECULE *);
void reallocate_atom_neighbors		(ATOM     *);
void reallocate_bonds			(MOLECULE *);
void reallocate_torsions		(MOLECULE *);
void reallocate_substs			(MOLECULE *);
void reallocate_segments		(MOLECULE *);
void reallocate_segment_atoms		(SEGMENT  *);
void reallocate_segment_neighbors	(SEGMENT  *);
void reallocate_layers			(MOLECULE *);
void reallocate_layer_segments		(LAYER    *);
void reallocate_sets			(MOLECULE *);
void reallocate_keys			(MOLECULE *);

void copy_molecule			(MOLECULE *, MOLECULE *);
void copy_info				(MOLECULE *, MOLECULE *);
void copy_transform			(MOLECULE *, MOLECULE *);
void copy_score				(MOL_SCORE *, MOL_SCORE *);
void copy_atoms				(MOLECULE *, MOLECULE *);
void copy_atom				(ATOM     *, ATOM     *);
void copy_atom_neighbors		(ATOM     *, ATOM     *);
void copy_coords			(MOLECULE *, MOLECULE *);
void copy_coord				(XYZ       , XYZ       );
void copy_bonds				(MOLECULE *, MOLECULE *);
void copy_bond				(BOND     *, BOND     *);
void copy_torsions			(MOLECULE *, MOLECULE *);
void copy_torsion			(TORSION  *, TORSION  *);
void copy_substs			(MOLECULE *, MOLECULE *);
void copy_subst				(SUBST    *, SUBST    *);
void copy_segments			(MOLECULE *, MOLECULE *);
void copy_segment			(SEGMENT  *, SEGMENT  *);
void copy_layers			(MOLECULE *, MOLECULE *);
void copy_layer				(LAYER    *, LAYER    *);
void copy_sets				(MOLECULE *, MOLECULE *);
void copy_set				(SET      *, SET      *);
void copy_keys				(MOLECULE *, MOLECULE *);

void save_molecule			(MOLECULE *, FILE *);
void save_info				(MOLECULE *, FILE *);
void save_transform			(MOLECULE *, FILE *);
void save_score				(MOL_SCORE *, FILE *);
void save_atoms				(MOLECULE *, FILE *);
void save_atom				(ATOM     *, FILE *);
void save_atom_neighbors		(ATOM     *, FILE *);
void save_bonds				(MOLECULE *, FILE *);
void save_bond				(BOND     *, FILE *);
void save_torsions			(MOLECULE *, FILE *);
void save_substs			(MOLECULE *, FILE *);
void save_subst				(SUBST    *, FILE *);
void save_segments			(MOLECULE *, FILE *);
void save_segment			(SEGMENT  *, FILE *);
void save_layers			(MOLECULE *, FILE *);
void save_layer				(LAYER    *, FILE *);
void save_sets				(MOLECULE *, FILE *);
void save_set				(SET      *, FILE *);
void save_keys				(MOLECULE *, FILE *);

void load_molecule			(MOLECULE *, FILE *);
void load_info				(MOLECULE *, FILE *);
void load_transform			(MOLECULE *, FILE *);
void load_score				(MOL_SCORE *, FILE *);
void load_atoms				(MOLECULE *, FILE *);
void load_atom				(ATOM     *, FILE *);
void load_atom_neighbors		(ATOM     *, FILE *);
void load_bonds				(MOLECULE *, FILE *);
void load_bond				(BOND     *, FILE *);
void load_torsions			(MOLECULE *, FILE *);
void load_substs			(MOLECULE *, FILE *);
void load_subst				(SUBST    *, FILE *);
void load_segments			(MOLECULE *, FILE *);
void load_segment			(SEGMENT  *, FILE *);
void load_layers			(MOLECULE *, FILE *);
void load_layer				(LAYER    *, FILE *);
void load_sets				(MOLECULE *, FILE *);
void load_set				(SET      *, FILE *);
void load_keys				(MOLECULE *, FILE *);

void save_string			(char **, FILE *);
void load_string			(char **, FILE *);

int get_identifier			(void);

void atom_neighbors			(MOLECULE *);
void revise_atom_neighbors		(MOLECULE *);

void get_torsion_neighbors		(MOLECULE *);
void reverse_torsion			(TORSION  *);
float calculate_distances		(MOLECULE *, float ***, int *);
void free_distances			(float ***, int *);

int get_atom_neighbor
(
  void *atom,
  int atom_id,
  int neighbor_id
);

void flag_atom_neighbor
( 
  void  *atom_in,
  int   atom_id,
  int   neighbor,
  int   flag   
);

void flag_segment_neighbor
( 
  void  *segment_in,
  int   segment_id,
  int   neighbor,
  int   flag   
);

int get_segment_neighbor
(
  void *segment,
  int segment_id,
  int neighbor_id
);

int get_atom_centrality
(
  MOLECULE *molecule,
  int atom_id
);

