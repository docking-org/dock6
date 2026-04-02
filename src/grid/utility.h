/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/
FILE *efopen (const char *, const char *, FILE *);
FILE *rfopen (const char *, const char *, FILE *);
void efclose (FILE **);
void rremove (const char *, FILE *);

void efread (void *, size_t , size_t , FILE *);
void efwrite (void *, size_t , size_t , FILE *);

void emalloc (void **, const int, const char *, FILE *);
void ecalloc (void **, const int, const int, const char *, FILE *);
void erealloc (void **, const int, const char *, FILE *);
void efree (void **);

char *subst_char (char *line, char old, char new);
char *strip_char (char *line, char character);
char *strip_newline (char *line);
char *white_line (char *line);

char *vstrcpy (char **, const char *);
char *vstrcat (char **, const char *);
char *vchrcat (char **, const char);
char *vfgets (char **, FILE *);

char *find_record (char **, const char *, FILE *);
