/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "define.h"
#include "utility.h"
#include "global.h"
#include "parm.h"

/* /////////////////////////////////////////////////////////// */

void read_parameters (PARM *parm)
{
  STRING400 line;

  if (global.infile == stdin)
    return;

/*
* Count the number of lines in the parameter file
* 3/96 te
*/
  fseek (global.infile, 0, SEEK_SET);
  for
  (
    parm->total = 0;
    fgets (line, sizeof (STRING400), global.infile) != NULL;
    parm->total++
  );

  fseek (global.infile, 0, SEEK_SET);

/*
* Allocate memory for storing each line and keyword read from parameter file
* 3/96 te
*/
  if (parm->total > 0)
    ecalloc
    (
      (void **) &parm->line,
      parm->total,
      sizeof (STRING400),
      "parm line buffer",
      global.outfile
    );

  for
  (
    parm->total = 0;
    fgets (parm->line[parm->total], sizeof (STRING400), global.infile) != NULL;
    parm->total++
  );

  fseek (global.infile, 0, SEEK_END);
}


/* /////////////////////////////////////////////////////////// */

void get_parameter
(
  void *variable,
  PARM *parm,
  enum VARIABLE_TYPE variable_type,
  char variable_name[],
  char suggest[],
  int query
)
{
  int i, j, k;
  int find_value;
  int get_value;
  STRING400 parameter = "";
  STRING400 value = "";
  char **name = NULL;
  int name_total;
  char *name_ptr = NULL;

  find_value = FALSE;
  get_value = FALSE;

/*
* Count the number of backwards compatible alternatives in the variable name
* 12/97 te
*/
  for
  (
    name_ptr = variable_name, name_total = 1;
    name_ptr = strchr (name_ptr, '|');
    name_ptr++, name_total++
  );

/*
* Allocate space for each entry in variable name and copy it to name array
* 12/97 te
*/
  ecalloc
  (
    (void **) &name,
    name_total,
    sizeof (char *),
    "variable names",
    global.outfile
  );

  for (i = j = 0; i < name_total; i++)
  {
    ecalloc
    (
      (void **) &name[i],
      strlen (variable_name) + 1,
      sizeof (char),
      "variable names",
      global.outfile
    );

    for (k = 0; j < strlen (variable_name); j++, k++)
    {
      if (variable_name[j] == '|')
      {
        j++;
        break;
      }

      name[i][k] = variable_name[j];
    }
  }
  
/*
* Get a value from the user
* 3/96 te
*/
  if (query)
  {
/*
*   If we have buffered input lines, then scan them for a value
*   3/96 te
*/
    if ((global.infile != stdin) && (parm->total > 0))
    {
      for (i = 0; (i < name_total) && (find_value != TRUE); i++)
        for (j = 0; j < parm->total; j++)
          if
          (
            (sscanf (parm->line[j], "%200s", parameter) == 1) &&
            (!strcmp (parameter, name[i]))
          )
          {
            if (find_value)
              exit (fprintf (global.outfile,
                "ERROR get_parameter: "
                "Multiple entries of %s parameter in input file.\n", name[i]));

            find_value = TRUE;

            if (i > 0)
              fprintf (global.outfile,
                "\nWARNING get_parameter: "
                "%s not found in input file,\n"
                "  so %s value used instead.\n\n",
                name[0], name[i]);

            if (sscanf (parm->line[j], "%*s %200s", value) == 1)
              get_value = TRUE;

            else
              exit (fprintf (global.outfile,
                "ERROR get_parameter: "
                "Unable to get value for %s parameter.\n", name[i]));
          }
    }

/*
*   If we don't have buffered input lines, or the value wasn't found,
*   then try asking for the value interactively.
*   3/96 te
*/
    if (((global.infile == stdin) || !get_value) &&
      (global.outfile == stdout))
    {
      fprintf (global.outfile, "%-30s [%s] ", name[0], suggest);

      if (fgets (value, sizeof (STRING400), stdin) != NULL) {
        if (strcmp (value, "\n") == 0)
          sscanf (suggest, "%200s", value);
        get_value = TRUE;
      }
      else
        exit (fprintf (global.outfile,
          "ERROR get_parameter: "
          "Unable to get value for %s parameter.\n", name[0]));
    }
  }

/*
* If user not queried, the use the suggested value
* 3/96 te
*/
  else
  {
    sscanf (suggest, "%200s", value);
    find_value = TRUE;
    get_value = TRUE;
  }

/*
* Check if a value was read
* 3/96 te
*/
  if (!get_value)
    exit (fprintf (global.outfile,
      "ERROR get_parameter: Unable to get value for %s parameter.\n", name[0]));

/*
* Store the value in the variable
* 3/96 te
*/
  get_value = FALSE;

  switch (variable_type)
  {
    case Boolean:
      if (tolower (value[0]) == 'y')
      {
        *(int *) variable = TRUE;
        strcpy (value, "yes");
      }

      else
      {
        *(int *) variable = FALSE;
        strcpy (value, "no");
      }

      get_value = TRUE;

      break;

    case Integer:
      if (!strcmp (value, "<infinity>"))
      {
        *(int *) variable = INT_MAX;
        get_value = TRUE;
      }

      else if (sscanf (value, "%d", (int *) variable) == 1)
      {
        sprintf (value, "%d", *(int *) variable);
        get_value = TRUE;
      }

      break;

    case Real:
      if (!strcmp (value, "<infinity>"))
      {
        *(float *) variable = FLT_MAX;
        get_value = TRUE;
      }

      else if (sscanf (value, "%f", (float *) variable) == 1)
      {
        sprintf (value, "%g", *(float *) variable);
        get_value = TRUE;
      }

      break;

    case Character:
      if (*(int *) variable = (int) tolower (value[0]))
      {
        if (*(int *) variable == 'y')
          sprintf (value, "yes");

        else if (*(int *) variable == 'n')
          sprintf (value, "no");

        else
          sprintf (value, "%c", *(int *) variable);

        get_value = TRUE;
      }

      break;

    case String:
      if (sscanf (value, "%s", (char *) variable) == 1)
                               /* assume variable is big enough */
        get_value = TRUE;

      if ((!strncmp (value, "stdin", 5)) && (global.infile == stdin))
        exit (fprintf (global.outfile,
          "ERROR get_parameter: stdin stream already in use.\n"));

      else if ((!strncmp (value, "stdout", 6)) && (global.outfile == stdout))
        exit (fprintf (global.outfile,
          "ERROR get_parameter: stdout stream already in use.\n"));

      break;

    default:
      exit (fprintf (global.outfile,
        "ERROR get_parameter: Unknown data type for parameter.\n"));
  }

/*
* Check if value was stored in variable
* 3/96 te
*/
  if (!get_value)
    exit (fprintf (global.outfile,
      "ERROR get_parameter: Unable to store value for %s parameter.\n",
      name[0]));

/*
* Update the input file if:
*   An input file is being read AND either
*     The parameter was missing from the file, but was read interactively, OR
*     The parameter was NOT read from the file but running in verbose mode
* 3/96 te
*/
  if ((global.infile != stdin) && query && !find_value)
    fprintf (global.infile, "%-30s %s\n", name[0], value);

/*
* Report value stored in variable to output file, if;
*   The value was queried, but not entered interactively, or
*   The value was not queried and verbose mode
* 3/96 te
*/
  if ((query && find_value) || (!query && (global.output_volume == 'v')))
    fprintf (global.outfile, "%-30s %s\n", name[0], value);

/*
* Free memory for name array
* 12/97 te
*/
  for (i = 0; i < name_total; i++)
    efree ((void **) &name[i]);

  efree ((void **) &name);
}


/* //////////////////////////////////////////////////////////////////

Set limit on how much memory can be dynamically allocated for this run
11/95 te

////////////////////////////////////////////////////////////////// */

void set_memory_limit (void)
{
  time_t tp;
  struct rlimit total, resident;

#ifdef RLIMIT_RSS
/*
* Set the total memory limit to the resident amount (help prevent swapping)
* 11/95 te
*/

  getrlimit (RLIMIT_RSS, &resident);
  getrlimit (RLIMIT_DATA, &total);

  resident.rlim_cur = total.rlim_cur = MIN (resident.rlim_cur, total.rlim_cur);

  setrlimit (RLIMIT_RSS, &resident);
  setrlimit (RLIMIT_DATA, &total);

#else
/*
* Set the total memory limit
*/

  getrlimit (RLIMIT_DATA, &total);

#endif

/*
* Get the current time
* 11/95 te
*/
  time (&tp);

/*
* Output the information about this job
* 11/95 te
*/
  fprintf (global.outfile,
    "\n__________________Job_Information_________________\n");

  fprintf (global.outfile, "%-30s %s", "launch_time", ctime (&tp));
  fprintf (global.outfile, "%-30s %s\n", "host_name",
    getenv("HOST") ? getenv("HOST") : "unknown");
  fprintf (global.outfile, "%-30s %d\n", "memory_limit",
    total.rlim_cur);
  fprintf (global.outfile, "%-30s %s\n", "working_directory",
    getenv("PWD") ? getenv("PWD") : "unknown");
  fprintf (global.outfile, "%-30s %s\n", "user_name",
    getenv("USER") ? getenv("USER") : "unknown");
}

