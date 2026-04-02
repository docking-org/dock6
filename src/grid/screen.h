/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

int  update_keys		(LABEL_CHEMICAL *, MOLECULE *, int);
int  update_folded_keys		(LABEL_CHEMICAL *, MOLECULE *, int);
int  update_unfolded_keys	(LABEL_CHEMICAL *, MOLECULE *, int);
MASK get_mask			(CHEMICAL_SCREEN *, float);

void fold_keys
(
  LABEL_CHEMICAL        *label_chemical,
  MOLECULE              *copy,
  MOLECULE              *original
);

int check_pharmacophore
(
  LABEL_CHEMICAL        *label_chemical,
  float                 uncertainty,
  MOLECULE              *target,
  MOLECULE              *candidate
);

void update_equivalency
(
  LABEL_CHEMICAL        *label_chemical,
  MOLECULE              *molecule
);

void update_uncertainty
(
  LABEL_CHEMICAL        *label_chemical,
  float         	uncertainty,
  MOLECULE      	*molecule
);

int mask_keys
(
  MOLECULE      *target,
  MOLECULE      *candidate
);

float check_dissimilarity
(
  LABEL_CHEMICAL        *label_chemical,
  MOLECULE              *target,
  MOLECULE              *candidate
);

float compare_keys
(
  LABEL_CHEMICAL        *label_chemical,
  MOLECULE      	*target,
  MOLECULE      	*candidate
);

int bit_count (MASK mask);
