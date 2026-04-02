/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

int write_info
(
  DOCK *dock,
  SCORE *score,
  LIST *list,
  int ligand_read_num,
  int ligand_dock_num,
  int ligand_skip_num,
  float time
);

int write_topscorers
(
  DOCK *dock,
  SCORE *score,
  LIST *list,
  MOLECULE *mol_ref,
  MOLECULE *mol_ori
);

int write_restartinfo
(
  DOCK *dock,
  SCORE *score,
  LIST *list,
  int ligand_read_num,
  int ligand_dock_num,
  int ligand_skip_num,
  float clock_elapsed
);

int read_restartinfo
(
  DOCK *dock,
  SCORE *score,
  LIST *list,
  int *ligand_read_num,
  int *ligand_dock_num,
  int *ligand_skip_num,
  float *clock_elapsed
);

