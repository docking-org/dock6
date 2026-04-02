/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Structures and routines to perform searches of linked elements
(eg. bonded atoms) in either depth-first or breadth-first manner.

Written by Todd Ewing
11/96
*/

typedef struct search_struct
{
  int complete_flag;
  int total_size;
  int max_size;
  int thread_total;
  int radius;
  int **target;
  int **origin;
  int **thread;
  int *total;
  int *count;
  int *neighbor_id;
  int *log;
} SEARCH;


/*
Routines to perform searches
11/96 te
*/

void allocate_search	(SEARCH *search);
void reset_search	(SEARCH *search);
void free_search	(SEARCH *search);
void reallocate_search	(SEARCH *search);

int breadth_search
(
  SEARCH        *search,
  void      	*array,
  int           total,
  int           get_neighbor (void *, int, int),
  void          flag_neighbor (void *, int, int, int),
  int           *seed,
  int           seed_total,
  int           avoid,
  int		iteration
);

int depth_search
(
  SEARCH        *search,
  void      	*array,
  int           total,
  int           get_neighbor (void *, int, int),
  void          flag_neighbor (void *, int, int, int),
  int           *seed,
  int           seed_total,
  int           avoid,
  int		iteration
);

int get_search_radius
(
  SEARCH        *search,
  int           target,
  int           origin
);

int get_search_origin
(
  SEARCH        *search,
  int           target,
  int           radius,
  int           count
);

int get_search_target
(
  SEARCH        *search,
  int           origin,
  int           radius,
  int           count
);

