/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

typedef struct parm_struct
{
  STRING200 *line;
  int total;

} PARM;

enum VARIABLE_TYPE {Boolean, Integer, Real, Character, String};

void read_parameters (PARM *parm);

void get_parameter
(
  void			*variable,
  PARM			*parm,
  enum VARIABLE_TYPE	variable_type,
  char			name[],
  char			suggest[],
  int			query
);

int yes_or_no (char *value);
