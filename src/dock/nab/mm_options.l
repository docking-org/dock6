%{
extern int mytaskid;

#ifdef flex

int  mmo_yyinput( char *, int );
#undef YY_INPUT
#define YY_INPUT(b, r, ms) (r = mmo_yyinput(b, ms))

#else

#undef input
#undef unput
#define input() (*gopts++)
#define unput(c) (*--gopts = c)

#endif

#undef ECHO
#define ECHO if(mytaskid==0)fprintf(yyout, "\tmm_options:  %s\n",yytext)

%}
%%

cut[=\ ][0-9eE\.\+\-]+		{ ECHO; cut  = atof( &yytext[4] ); }
scnb[=\ ][0-9eE\.\+\-]+		{ ECHO; scnb = atof( &yytext[5] ); }
scee[=\ ][0-9eE\.\+\-]+		{ ECHO; scee = atof( &yytext[5] ); }
ntpr[=\ ][1-9][0-9]*			{ ECHO; ntpr = atoi(&yytext[5]); }
nsnb[=\ ][0-9]+			{ ECHO; nsnb = atoi( &yytext[5] ); }
ntwx[=\ ][0-9]+			{ ECHO; ntwx = atoi( &yytext[5] ); }
nchk[=\ ][0-9]+			{ ECHO; nchk = atoi( &yytext[5] ); }
nchk2[=\ ][0-9]+		{ ECHO; nchk2 = atoi( &yytext[6] ); }

dt[=\ ][0-9eE\.\+\-]+		{ ECHO; dt  = atof( &yytext[3] ); }
t[=\ ][0-9eE\.\+\-]+		{ ECHO; t = atof( &yytext[2] ); }
tautp[=\ ][0-9eE\.\+\-]+	{ ECHO; tautp = atof( &yytext[6] ); }
temp0[=\ ][0-9eE\.\+\-]+	{ ECHO; temp0 = atof( &yytext[6] ); }
tempi[=\ ][0-9eE\.\+\-]+	{ ECHO; tempi = atof( &yytext[6] ); }
vlimit[=\ ][0-9eE\.\+\-]+	{ ECHO; vlimit = atof( &yytext[7] ); }
gamma_ln[=\ ][0-9eE\.\+\-]+	{ ECHO; gamma_ln = atof( &yytext[9] ); }
genmass[=\ ][0-9eE\.\+\-]+	{ ECHO; genmass = atof( &yytext[8] ); }
wcons[=\ ][0-9eE\.\+\-]+	{ ECHO; wcons = atof( &yytext[6] ); }
ntpr_md[=\ ][1-9][0-9]*		{ ECHO; ntpr_md = atoi(&yytext[8]); }
zerov[=\ ][0-9]+		{ ECHO; zerov = atoi( &yytext[6] ); }
k4d[=\ ][0-9eE\.\+\-]+		{ ECHO; k4d  = atof( &yytext[4] ); }

gb[=\ ][0125]			{ ECHO; gb = atoi( &yytext[3] ); }
gb_debug[=\ ][01]		{ ECHO; gb_debug = atoi( &yytext[9] ); }
nr_debug[=\ ][01]		{ ECHO; nr_debug = atoi( &yytext[9] ); }
gb2_debug[=\ ][01]		{ ECHO; gb2_debug = atoi( &yytext[10] ); }
e_debug[=\ ][01]		{ ECHO; e_debug = atoi( &yytext[8] ); }
epsext[=\ ][0-9eE\.\+\-]+	{ ECHO; epsext = atof( &yytext[7] ); }
kappa[=\ ][0-9eE\.\+\-]+	{ ECHO; kappa = atof( &yytext[6] ); }
rgbmax[=\ ][0-9eE\.\+\-]+	{ ECHO; rgbmax = atof( &yytext[7] ); }

surften[=\ ][0-9eE\.\+\-]+	{ ECHO; surften = atof( &yytext[8] ); }
gbsa[=\ ][01]			{ ECHO; gbsa = atoi( &yytext[5] ); }
dim[=\ ][34]			{ ECHO; dim = atoi( &yytext[4] ); }

diel[=\ ]Yammp				{ ECHO; dield = -5; }
diel[=\ ]R94				{ ECHO; dield = -4; }
diel[=\ ]RL94				{ ECHO; dield = -3; }
diel[=\ ]RL					{ ECHO; dield = -2; }
diel[=\ ]R					{ ECHO; dield =  0; }
diel[=\ ]C					{ ECHO; dield =  1; }
dielc[=\ ][0-9eE\.\+\-]+	{ ECHO; dielc = atof( &yytext[6] ); }

rattle[=\ ][01]			{ ECHO; irattle = atoi( &yytext[7] ); }

[\ \n\t,]+						;

.						{	
#ifdef flex
                if(mytaskid==0)fprintf( stderr, "Error: unknown mm_option:\n" );
		if(mytaskid==0)fprintf( stderr, "  first character of the option is '%s'.\n", yytext );
#else
		if(mytaskid==0)fprintf( stderr, "Error unknown mm_option: %s\n", --gopts );
#endif

		exit( 1 );
	}

%%

int yywrap(){
	return( 1 );
}

#ifdef flex


int  mmo_yyinput( buf, max_size )

char *buf;
int  max_size;
{
	int n;

	if( max_size > &gopts[mmoinputlim] - mmoinputptr ) 
		n =  &gopts[mmoinputlim] - mmoinputptr;
	else
		n = max_size;

	if( n>0 ){
		memcpy( buf, mmoinputptr, n );
		mmoinputptr += n;
	}
	return n;
}

#endif
