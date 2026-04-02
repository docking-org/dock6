/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "define.h"
#include "global.h"

FILE *efopen
(
  const char *file_name,
  const char *mode,
  FILE *message_stream
)
{
  FILE	*file;

  if (!strncmp (file_name, "stdin", 5))
  {
    if (!strcmp (mode, "r"))
      file = stdin;

    else
      exit (fprintf (message_stream,
        "ERROR efopen: %s stream cannot be opened in %s mode.\n",
        file_name, mode));
  }
    
  else if (!strncmp (file_name, "stdout", 6))
  {
    if (!strcmp (mode, "w"))
      file = stdout;

    else
      exit (fprintf (message_stream,
        "ERROR efopen: %s stream cannot be opened in %s mode.\n",
        file_name, mode));
  }

  else
  {
    file = fopen (file_name, mode);

    if (file == NULL)
      exit (fprintf (message_stream,
        "ERROR efopen: Unable to open %s.  Terminating execution.\n",
        file_name));
  }

  return file;
}


/* ////////////////////////////////////////////////////////////////// */

FILE *rfopen
(
  const char *file_name,
  const char *mode,
  FILE *message_stream
)
{
  FILE	*file;

  if (!strncmp (file_name, "stdin", 5))
  {
    if (!strcmp (mode, "r"))
      file = stdin;

    else
      exit (fprintf (message_stream,
        "ERROR efopen: %s stream cannot be opened in %s mode.\n",
        file_name, mode));
  }
    
  else if (!strncmp (file_name, "stdout", 6))
  {
    if (!strcmp (mode, "w"))
      file = stdout;

    else
      exit (fprintf (message_stream,
        "ERROR efopen: %s stream cannot be opened in %s mode.\n",
        file_name, mode));
  }

  else
  {
    file = fopen (file_name, mode);

    if (file == NULL)
    {
      fprintf (message_stream,
        "WARNING rfopen: Unable to open %s, but will continue trying.\n",
        file_name);
      fflush (message_stream);
      while (file == NULL)
        file = fopen (file_name, mode);
    }
  }

  return file;
}


/* ////////////////////////////////////////////////////////////////// */

void efclose (FILE **file)
{
  if ((*file != stdin) && (*file != stdout) && (fclose (*file) != 0))
    exit (fprintf (global.outfile,
      "ERROR efclose: Unable to close file.\n"));

  *file = NULL;
}


/* ////////////////////////////////////////////////////////////////// */

void efread (void *ptr, size_t size, size_t nobj, FILE *stream)
{
  if (fread (ptr, size, nobj, stream) != nobj)
    exit (fprintf (stderr,
      "ERROR efread: Unable to read proper number of objects.\n"));
}


/* ////////////////////////////////////////////////////////////////// */

void rremove
(
  const char *file_name,
  FILE *message_stream
)
{
  if (remove (file_name))
  {
    fprintf (message_stream,
      "WARNING rremove: Unable to remove %s, but will continue trying.\n",
      file_name);
    fflush (message_stream);
    while (remove (file_name));
  }
}


/* ////////////////////////////////////////////////////////////////// */

void efwrite (void *ptr, size_t size, size_t nobj, FILE *stream)
{
  if (fwrite (ptr, size, nobj, stream) != nobj)
    exit (fprintf (stderr,
      "ERROR efwrite: Unable to write proper number of objects.\n"));
}


/* ////////////////////////////////////////////////////////////////// */

void emalloc
(
  void **ptr,
  const int size,
  const char *name,
  FILE *message_stream
)
{
  if (*ptr != NULL)
    exit (fprintf (message_stream,
      "ERROR emalloc: %s not freed prior to allocation.\n", name));

  if (size > 0)
  {
    *ptr = (void *) malloc (size);

    if (*ptr == NULL)
      exit (fprintf (message_stream,
        "ERROR emalloc: Insufficient memory for %s.\n", name));
  }
}


/* ////////////////////////////////////////////////////////////////// */

void ecalloc
(
  void **ptr,
  const int objectnum,
  const int objectsize,
  const char *name,
  FILE *message_stream
)
{
  if (*ptr != NULL)
    exit (fprintf (message_stream,
      "ERROR ecalloc: %s not freed prior to allocation.\n", name));

  if (objectnum * objectsize > 0)
  {
    *ptr = (void *) calloc (objectnum, objectsize);

    if (*ptr == NULL)
      exit (fprintf (message_stream,
        "ERROR ecalloc: Insufficient memory for %s.\n", name));
  }
}


/* ////////////////////////////////////////////////////////////////// */

void erealloc
(
  void **ptr,
  const int size,
  const char *name,
  FILE *message_stream
)
{
  void *new_ptr;

  new_ptr = (void *) realloc (*ptr, size);

  if ((size > 0) && (new_ptr == NULL))
    exit (fprintf (message_stream,
      "ERROR erealloc: Insufficient memory for %s.\n", name));

  else
    *ptr = new_ptr;
}


/* ////////////////////////////////////////////////////////////////// */

void efree (void **ptr)
{
  if (ptr == NULL)
    exit (fprintf (global.outfile,
      "ERROR efree: NULL address passed\n"));

  else if (*ptr != NULL)
  {
    free (*ptr);
    *ptr = NULL;
  }
}


/* ////////////////////////////////////////////////////////////////// */

char *subst_char (char *line, char old, char new)
{
  char *ptr;

  if (line)
    for (ptr = line; *ptr != '\0'; ptr++)
      if (*ptr == old)
        *ptr = new;

  return line;
}


/* ////////////////////////////////////////////////////////////////// */

char *strip_char (char *line, char character)
{
  if ((line) && (strchr (line, character)))
    strcpy (strchr (line, character), "");

  return line;
}


/* ////////////////////////////////////////////////////////////////// */

char *strip_newline (char *line)
{
  if (line[strlen (line) - 1] == '\n')
    line[strlen (line) - 1] = 0;

  return line;
}


/* ////////////////////////////////////////////////////////////////// */

char *white_line (char *line)
{
  int i;

  for (i = 0; i < strlen (line); i++)
    if (isspace (line[i]))
      line[i] = ' ';

  return line;
}


/* ////////////////////////////////////////////////////////////////// */

char *vstrcpy (char **copy_ptr, const char *original)
{
  if (original)
  {
    if (*copy_ptr)
    {
      if (strlen (original) > strlen (*copy_ptr))
      {
        efree ((void **) copy_ptr);
        ecalloc
        (
          (void **) copy_ptr,
          strlen (original) + 1,
          sizeof (char),
          original,
          stdout
        );
      }

      else
        memset (*copy_ptr, 0, strlen (*copy_ptr));
    }

    else
      ecalloc
      (
        (void **) copy_ptr,
        strlen (original) + 1,
        sizeof (char),
        original,
        stdout
      );

    strcpy (*copy_ptr, original);
  }

  else if (*copy_ptr)
  {
    efree ((void **) copy_ptr);
    *copy_ptr = NULL;
  }

  return *copy_ptr;
}


/* ////////////////////////////////////////////////////////////////// */

char *vstrcat (char **base_ptr, const char *addition)
{
  char *temp = NULL;

  if (addition)
  {
    if (*base_ptr)
    {
      temp = *base_ptr;
      *base_ptr = NULL;

      ecalloc
      (
        (void **) base_ptr,
        strlen (temp) + strlen (addition) + 1,
        sizeof (char),
        addition,
        stdout
      );

      strcpy (*base_ptr, temp);
      efree ((void **) &temp);
    }

    else
      ecalloc
      (
        (void **) base_ptr,
        strlen (addition) + 1,
        sizeof (char),
        addition,
        stdout
      );

    strcat (*base_ptr, addition);
  }

  return *base_ptr;
}


/* ////////////////////////////////////////////////////////////////// */

char *vchrcat (char **base_ptr, const char addition)
{
  char *temp = NULL;

  if (*base_ptr)
  {
    temp = *base_ptr;
    *base_ptr = NULL;

    ecalloc
    (
      (void **) base_ptr,
      strlen (temp) + 2,
      sizeof (char),
      "vchrcat intermediate",
      stdout
    );

    strcpy (*base_ptr, temp);
    efree ((void **) &temp);
  }

  else
    ecalloc
    (
      (void **) base_ptr,
      2,
      sizeof (char),
      "vchrcat intermediate",
      stdout
    );

  (*base_ptr)[strlen (*base_ptr)] = addition;

  return *base_ptr;
}


/* /////////////////////////////////////////////////////////////////// */

char *vfgets (char **line_ptr, FILE *file)
{
  STRING100 line = "";
  char *return_value = NULL;

  efree ((void **) line_ptr);

  while (return_value = fgets (line, 100, file))
  {
    vstrcat (line_ptr, line);

    if (line[strlen (line) - 1] == '\n')
      break;
  }

  if (return_value == NULL)
    return NULL;

  else
    return *line_ptr;
}


/* ///////////////////////////////////////////////////////////////////
Function to advance into a file until a particular record is found.

Return value:
  NULL: if end of file is reached before finding record
  char *: if the record is found, then return the record
/////////////////////////////////////////////////////////////////// */

char *find_record
(
  char **line,
  const char *record,
  FILE *file
)
{
  char *return_value = NULL;

  while (((return_value = vfgets (line, file)) != NULL) &&
    strncmp (*line, record, strlen (record)));

  return return_value;
}

