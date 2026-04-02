/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

void read_grids
(
  SCORE_GRID *,
  SCORE_BUMP *,
  SCORE_CONTACT *,
  SCORE_CHEMICAL *,
  SCORE_ENERGY *,
  LABEL_CHEMICAL *
);

void write_grids
(
  SCORE_GRID *,
  SCORE_BUMP *,
  SCORE_CONTACT *,
  SCORE_CHEMICAL *,
  SCORE_ENERGY *,
  LABEL_CHEMICAL *
);

void free_grids
(
  SCORE_GRID *,
  SCORE_BUMP *,
  SCORE_CONTACT *,
  SCORE_CHEMICAL *,
  SCORE_ENERGY *,
  LABEL_CHEMICAL *
);

int read_box
(
  FILE_NAME file_name,
  XYZ center_of_mass,
  XYZ dimension
);

