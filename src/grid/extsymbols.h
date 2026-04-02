
/* Define the external symbols.  These must be defined so that */
/* the Fortran interface is portable.  This machine-dependency */
/* is localized in this file in order to simplify the C code.  */

/* After many years, three schemes of external name                    */
/* generation by Fortran compilers are still commonly used:            */
/* the name in lowercase with a trailing underscore, the name in       */
/* lowercase, and the name in uppercase                                */
/* in descending order of popularity.                                  */
/* These schemes are usually beyond the control of the compiler's      */
/* user. (Other schemes exist, notably two trailing underscores        */
/* by GNU compilers, but compiler options can control that behavior.)  */
/* Thus some mechanism to translate external names between Fortran     */
/* and C, which usually does no mangling, is still necessary.          */
/* We use a mechanism intended to decouple a platform name from the    */
/* scheme employed by that platform.  Two preprocessor names exist:    */
/* CLINK_CAPS  should be defined when external names are in uppercase; */
/* CLINK_PLAIN should be defined when external names are in lowercase. */
/* Defining one of these names occurs in the platform                  */
/* configuration files, if necessary.                                  */
/* When external names are in lowercase with a trailing underscore,    */
/* no special preprocessor name is necessary.  In addition, in         */
/* configuration files for GNU compilers the no second underscore      */
/* option, -fno-second-underscore, should be used.  Note that these    */
/* preprocessor names are the same as those used by the AMBER          */
/* software.                                                           */
/* Based on COLUMBUS, by Ron Shepard and Scott Brozell */


#if defined(CLINK_CAPS)

  /* Fortran symbols are converted to upper case with no trailing "_". */
#define EXTERNAL_GET_INDEX GET_INDEX
#define EXTERNAL_GET_VALUE GET_VALUE
#define EXTERNAL_ORIENT_GK ORIENT_GK
#define EXTERNAL_READGRID READGRID
#define EXTERNAL_TRANSFORM TRANSFORM
#define EXTERNAL_TRANSFORM_ATOM TRANSFORM_ATOM

#elif defined(CLINK_PLAIN)

  /* Fortran symbols are converted to lower case with no trailing "_". */
#define EXTERNAL_GET_INDEX get_index
#define EXTERNAL_GET_VALUE get_value
#define EXTERNAL_ORIENT_GK orient_gk
#define EXTERNAL_READGRID readgrid
#define EXTERNAL_TRANSFORM transform
#define EXTERNAL_TRANSFORM_ATOM transform_atom

#else

  /* default case: */
  /* Fortran symbols are converted to lower case, and a trailing "_" */
  /* is added automatically by the compiler.                         */
#define EXTERNAL_GET_INDEX get_index_
#define EXTERNAL_GET_VALUE get_value_
#define EXTERNAL_ORIENT_GK orient_gk_
#define EXTERNAL_READGRID readgrid_
#define EXTERNAL_TRANSFORM transform_
#define EXTERNAL_TRANSFORM_ATOM transform_atom_

#endif
