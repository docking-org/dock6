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
#include "parm.h"
#include "label.h"
#include "score.h"
#include "score_grid.h"
#include "grid.h"
#include "io.h"

int get_parameters
(
  GRID *grid,
  SCORE_GRID *score_grid,
  SCORE_BUMP *score_bump,
  SCORE_CONTACT *score_contact,
  SCORE_CHEMICAL *score_chemical,
  SCORE_ENERGY *score_energy,
  LABEL *label
)
{
  PARM parm = {0};

  read_parameters (&parm);

  fprintf (global.outfile, 
    "\n________________General_Parameters________________\n");

  get_parameter
  (
    (void *) &score_grid->flag,
    &parm, Boolean, "compute_grids",
    "NO yes",
    TRUE
  );

  if (score_grid->flag)
    label->vdw.flag = TRUE;

  get_parameter
  (
    (void *) &score_grid->spacing,
    &parm, Real, "grid_spacing",
    score_grid->flag ? "0.3" : "0.0",
    score_grid->flag
  );

  get_parameter
  (
    (void *) &grid->output_molecule,
    &parm, Boolean, "output_molecule",
    "NO yes",
    TRUE
  );

  if (!score_grid->flag && !grid->output_molecule)
    exit (fprintf (global.outfile,
      "WARNING get_parameters: "
      "No processing requested.  Execution terminated.\n"));

  if (score_grid->flag)
    fprintf (global.outfile, 
      "\n________________Scoring_Parameters________________\n");

  get_parameter
  (
    (void *) &score_contact->flag,
    &parm, Boolean, "contact_score",
    "NO yes",
    score_grid->flag
  );

  get_parameter
  (
    (void *) &score_contact->distance,
    &parm, Real, "contact_cutoff_distance",
    score_contact->flag ? "4.5" : "0.0",
    score_contact->flag
  );

  get_parameter
  (
    (void *) &score_chemical->flag,
    &parm, Boolean, "chemical_score",
    "NO yes",
    score_grid->flag
  );

  if (score_chemical->flag)
    label->chemical.flag = TRUE;

  get_parameter
  (
    (void *) &score_energy->flag,
    &parm, Boolean, "energy_score",
    "NO yes",
    score_grid->flag
  );

  if (score_chemical->flag && !score_energy->flag)
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "energy scoring is required with chemical scoring.\n"));

  get_parameter
  (
    (void *) &score_energy->distance,
    &parm, Real, "energy_cutoff_distance",
    score_energy->flag ? "10.0" : "0.0",
    score_energy->flag
  );

  get_parameter
  (
    (void *) &score_energy->atom_model,
    &parm, Character, "atom_model",
    score_energy->flag ? "UNITED all" : "0",
    score_energy->flag
  );

  if (score_energy->flag &&
    (score_energy->atom_model != 'a') &&
    (score_energy->atom_model != 'u'))
    return FALSE;

  get_parameter
  (
    (void *) &score_energy->attractive_exponent,
    &parm, Integer, "attractive_exponent",
    score_energy->flag ? "6" : "0",
    score_energy->flag
  );

  get_parameter
  (
    (void *) &score_energy->repulsive_exponent,
    &parm, Integer, "repulsive_exponent",
    score_energy->flag ? "12" : "0",
    score_energy->flag
  );

  if (score_energy->flag &&
    (score_energy->attractive_exponent >= score_energy->repulsive_exponent))
    exit (fprintf (global.outfile,
      "ERROR get_parameters: Require EXP(repulse) > EXP(attract).\n"));

  get_parameter
  (
    (void *) &score_energy->distance_dielectric,
    &parm, Boolean, "distance_dielectric",
    score_energy->flag ? "YES no" : "NO",
    score_energy->flag
  );

  get_parameter
  (
    (void *) &score_energy->dielectric_factor,
    &parm, Real, "dielectric_factor",
    score_energy->flag ? "4.0" : "0.0",
    score_energy->flag
  );

  get_parameter
  (
    (void *) &score_bump->flag,
    &parm, Boolean, "bump_filter",
    "NO yes",
    score_grid->flag
  );

  if
  (
    score_grid->flag &&
    !score_bump->flag &&
    !score_contact->flag &&
    !score_chemical->flag &&
    !score_energy->flag
  )
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "At least ONE type of scoring must be selected.\n"));

  if
  (
    !score_bump->flag &&
    (score_contact->flag ||
      score_chemical->flag)
  )
    exit (fprintf (global.outfile,
      "ERROR get_parameters: "
      "Contact and chemical scoring require a bump grid.\n"));

  get_parameter
  (
    (void *) &score_bump->clash_overlap,
    &parm, Real, "bump_overlap",
    score_bump->flag ? "0.75" : "0.0",
    score_bump->flag
  );


  fprintf (global.outfile, 
    "\n____________________File_Input____________________\n");

  get_parameter
  (
    (void *) &grid->in_file_name,
    &parm, String, "receptor_file",
    "receptor.mol2",
    TRUE
  );

/*  
  if (!check_file_extension (grid->in_file_name, TRUE))
    return FALSE;
*/ 

  get_parameter
  (
    (void *) &grid->box_file_name,
    &parm, String, "box_file",
    "site_box.pdb",
    score_grid->flag
  );

  get_parameter
  (
    (void *) &label->vdw.file_name,
    &parm, String, "vdw_definition_file",
    PARAMETER_PATH "vdw.defn",
    label->vdw.flag
  );

  get_parameter
  (
    (void *) &label->chemical.file_name,
    &parm, String, "chemical_definition_file",
    PARAMETER_PATH "chem.defn",
    label->chemical.flag
  );


  fprintf (global.outfile, 
    "\n____________________File_Output___________________\n");

  get_parameter
  (
    (void *) &score_grid->file_prefix,
    &parm, String, "score_grid_prefix",
    global.job_name,
    score_grid->flag
  );

  get_parameter
  (
    (void *) &grid->out_file_name,
    &parm, String, "receptor_out_file",
    "receptor_out.mol2",
    grid->output_molecule
  );

/*
  if (!check_file_extension (grid->out_file_name, TRUE))
    return FALSE;
*/    

  if (global.infile != stdin)
    fclose (global.infile);

  fprintf (global.outfile, "\n\n");
  return TRUE;
}


/* ////////////////////////////////////////////////////////////// */

int process_commands
(
  int argc,
  char *argv[]
)
{
  int i = 1, j;
  FILE_NAME infile_name = "", outfile_name = "";
  FILE *tmpfile ();

/*
* Extract name of global.executable from command line,
* ignoring any path description if present
* 6/95 te
*/
  strcpy (global.executable,
    strrchr (argv[0], '/') ? strrchr (argv[0], '/') + 1 : argv[0]);
  strcpy (global.job_name, global.executable);

/*
* Step through each command line argument
* 6/95 te
*/
  while (i < argc)
  {

/*
*   Look for command flag prefix (-)
*   6/95 te
*/
    if (argv[i][0] == '-')
    {

/*
*     For "-i" flag read in next field as the input file name
*     6/95 te
*/
      if (argv[i][1] == 'i')
      {
        if ((i + 1 < argc) && (argv[i + 1][0] != '-'))
        {
          strncpy (infile_name, argv[i + 1], 80);
          i++;

/*
*         Set the job name to the root name of the the "*.in" file
*         10/95 te
*/
          if (strrchr (infile_name, '.') &&
            !strcmp (strrchr (infile_name, '.'), ".in"))
          {
            for (j = 0; j < strlen (infile_name); j++)
            {
              if (!strcmp (&infile_name[j], ".in"))
              {
                global.job_name[j] = 0;
                break;
              }

              else
                global.job_name[j] = infile_name[j];
            }
          }

          else
            strcpy (global.job_name, infile_name);
        }

        else
          sprintf (infile_name, "%s.in", global.job_name);
      }

/*
*     For "-o" flag read in next field as the output file name.
*     If the next field does not exist or is another flag, then
*     use default output file name.
*     6/95 te
*/
      else if (argv[i][1] == 'o')
      {
        if ((i + 1 < argc) && (argv[i + 1][0] != '-'))
        {
          strncpy (outfile_name, argv[i + 1], 80);
          i++;
        }

        else
          sprintf (outfile_name, "%s.out", global.job_name);

        if (!strcmp (infile_name, ""))
          break;
      }

/*
*     For "-s" flag, allow all input and output through standard streams
*     6/95 te
*/
      else if (argv[i][1] == 's')
      {
        if ((strcmp (infile_name, "")) || (strcmp (outfile_name, "")))
          break;
      }

/*
*     For "-t" flag, set output volume to terse
*     6/95 te
*/
      else if (argv[i][1] == 't')
        global.output_volume = 't';

/*
*     For "-v" flag, set output volume to verbose
*     6/95 te
*/
      else if (argv[i][1] == 'v')
        global.output_volume = 'v';

      else
        break;

      i++;
    }

    else
      break;
  }

/*
* Check if no arguments were given, to conform to old style input/output
* 10/95 te
*/
  if (argc == 1)
  {
    strcpy (infile_name, "INCHEM");

/*
*   If INCHEM file exists, then assume user wants output directed to OUTCHEM
*   10/95 te
*/
    if (global.infile = fopen (infile_name, "r"))
    {
      fclose (global.infile);
      strcpy (outfile_name, "OUTCHEM");
    }

/*
*   If it doesn't, then trigger an error condition
*   10/95 te
*/
    else
      i = 0;
  }

/*
* If any error was made in the command line, then quit
* 6/95 te
*/
  if (i < argc)
  {
    fprintf (stderr,
      "  \n"
      "[33mUsage: %s [-i [input_file]] [-o [output_file]] ...\n"
      "  [-s[0mtandard_i/o[33m] [-t[0merse[33m] [-v[0merbose[33m][0m\n"
      "  \n"
      "  -i: read from %s.in or input_file, standard_in otherwise\n"
      "  -o: write to %s.out or output_file (-i required), \n"
      "      standard_out otherwise\n"
      "  -s: read from and write to standard streams (-i and/or -o illegal)\n"
      "  -t: terse program output\n"
      "  -v: verbose program output\n" 
      "  \n",
      global.executable, global.executable, global.executable);
    exit (EXIT_FAILURE);
  }

/*
* Open the input file
* 6/95 te
*/
  if (strcmp (infile_name, ""))
    global.infile = efopen (infile_name, "a+", stdout);

  else
    global.infile = stdin;

/*
* Open the output file
* 6/95 te
*/
  if (strcmp (outfile_name, ""))
      global.outfile = efopen (outfile_name, "w", stdout);

  else
    global.outfile = stdout;

  return TRUE;
}

