/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

/* Structures to store atom labelling definitions */

typedef struct node_struct
{
  STRING5 type;
  int include, next_total;
  int vector_atom, multiplicity;
  float weight;
  struct node_struct *next[6];
} NODE;


/* Routines to manipulate atom label structures */

NODE *create_node (void);
int assign_node   (NODE *node, int include);
void free_node    (NODE *node);
void print_node   (NODE *node, int level);
int check_type    (char *candidate, char *reference);

int check_atom
(
  MOLECULE *molecule,
  int      current_atom,
  NODE     *node
);

int check_bonded_atoms
(
  MOLECULE *molecule,
  int      current_atom,
  int      previous_atom,
  NODE     *node
);

