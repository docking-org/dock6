# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "ring.c"
# include "rotate.c"
# include "ac.c"
# include "charmm.c"
# include "mol2.c"
# include "mopcrt.c"
# include "divcrt.c"
# include "mopint.c"
# include "mopout.c"
# include "divout.c"
# include "gcrt.c"
# include "gzmat.c"
# include "gout.c"
# include "pdb.c"
# include "csd.c"
# include "mdl.c"
# include "alc.c"
# include "hin.c"
# include "prep.c"
# include "rst.c"
# include "jzmat.c"
# include "jcrt.c"
# include "jout.c"
# include "charge.c"

int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */
int ao_flag = 0;
/*For addtional file 
 1, only readin coordinates
 2, only readin charge
 3, only readin atom names
 4, only readin atom types
 5, only readin bond types 
*/
int atomtype_flag = 0;				/*judge atom type? */
int bondtype_flag = 0;				/*judge bond type? */
int default_flag = 0;				/*assign default information? */
int atomname_flag = 0;				/*assign atom name? */
int atomicnum_flag = 0;				/*judge atomic number according to atom name ? */
int adjustatomname_flag = 0;		/*adjust atom name? */
int duplicatedname_flag = 0;		/*check atom name duplication? */
int cartcoord_flag = 0;				/*generate coordinate from internal coordinate ? */
int connect_flag = 0;				/*judge atom connectivity and generate bond ? */
int divcon_flag = 1;  

int atomnum = 0;
int bondnum = 0;
int ringnum = 0;
ATOM *atom;
BOND *bond;
RING *ring;
AROM *arom;
MOLINFO minfo;
CONTROLINFO cinfo;

int atomnum_tmp = 0;
int bondnum_tmp = 0;
ATOM *atom_tmp;
BOND *bond_tmp;
char line[MAXCHAR];
char *system_env;
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char afilename[MAXCHAR];
char cfilename[MAXCHAR];

/*The following four functions, read_at(), write_at(), read_bt() and write_bt() are
used in amber sybyl interface development and are unreachable to the users
*/
int ra_flag = 0;
int rb_flag = 0;
int wa_flag = 0;
int wb_flag = 0;
char at_filename[MAXCHAR];
char bt_filename[MAXCHAR];

FILE *fpin;
FILE *fpout;

void usage()
{
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		printf("[31mUsage: antechamber -i  [0m input file name\n"
			   "[31m                   -fi [0m input file format\n"
			   "[31m                   -o  [0m output file name\n"
			   "[31m                   -fo [0m output file format\n"
			   "[31m                   -c  [0m charge method\n"
			   "[31m                   -cf [0m charge file name\n"
			   "[31m                   -nc [0m net molecular charge (int)\n"
			   "[31m                   -a  [0m additional file name\n"
			   "[31m                   -fa [0m additional file format\n"
			   "[31m                   -ao [0m additional file operation\n"
			   "[35m                        crd [0m: only read in coordinate\n"
			   "[35m                        crg[0m: only read in charge\n"
			   "[35m                        name  [0m: only read in atom name\n"
			   "[35m                        type  [0m: only read in atom type\n"
			   "[35m                        bond  [0m: only read in bond type \n"
			   "[31m                   -m  [0m mulitiplicity (2S+1), default is 1\n"
			   "[31m                   -rn [0m residue name, if not available in the input file, default is MOL\n"
			   "[31m                   -rf [0m residue toplogy file name in prep input file, default is molecule.res\n"
			   "[31m                   -ch [0m check file name in gaussian input file, default is molecule\n"
			   "[31m                   -mk [0m divcon or mopac keyword in a pair of quotation marks\n"
			   "[31m                   -gk [0m gaussian keyword in a pair of quotation marks\n"
			   "[31m                   -df [0m use divcon flag, 1 - always use divcon if $AMBERHOME is set (the default); 0 - not use divcon if $DOCK_HOME is set\n"
			   "[31m                   -at [0m atom type, can be gaff, amber, bcc and sybyl, default is gaff\n"
			   "[31m                   -du [0m check atom name duplications, can be yes(y) or no(n), default is yes\n"
			   "[31m                   -j  [0m atom type and bond type prediction index, default is 4 \n"
			   "[35m                        0    [0m: no assignment\n"
			   "[35m                        1    [0m: atom type \n"
			   "[35m                        2    [0m: full  bond types \n"
			   "[35m                        3    [0m: part  bond types \n"
			   "[35m                        4    [0m: atom and full bond type \n"
			   "[35m                        5    [0m: atom and part bond type \n"
			   "[31m                   -s  [0m status information, can be 0 (brief), 1 (the default) and 2 (verbose)\n"
			   "[31m                   -pf [0m remove the intermediate files: can be yes (y) and no (n), default is no\n"
			   "                   -i -o -fi and -fo must appear in command lines and the others are optional[0m");
		printf
			("\n\n	         	    [31m List of the File Formats [0m \n");
		printf
			("\n	 	file format type  abbre. index | file format type abbre. index");
		printf
			("\n		--------------------------------------------------------------- ");
		printf
			("\n		Antechamber        ac       1  | Sybyl Mol2         mol2    2 ");
		printf
			("\n		PDB                pdb      3  | Modifiled PDB      mpdb    4 ");
		printf
			("\n		AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 ");
		printf
			("\n		Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 ");
		printf
			("\n		Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 ");
		printf
			("\n		Gaussian Output    gout    11  | Mopac Output       mopout 12 ");
		printf
			("\n		Alchemy            alc     13  | CSD                csd    14 ");
		printf
			("\n		MDL                mdl     15  | Hyper              hin    16 ");
		printf
			("\n		AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 ");
		printf
			("\n		Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20");
		printf
			("\n		Divcon Input       divcrt  21  | Divcon Output      divout 22");
		printf
			("\n		--------------------------------------------------------------\n");
		printf
			("\n                AMBER restart file can only be read in as additional file\n");
		printf
			("\n	         	    [31m List of the Charge Methods [0m \n");
		printf
			("\n		charge method     abbre.  index | charge method      abbre. index");
		printf
			("\n		----------------------------------------------------------------  ");
		printf
			("\n		RESP               resp     1  |  AM1-BCC            bcc     2");
		printf
			("\n		CM1                cm1      3  |  CM2                cm2     4");
		printf
			("\n		ESP (Kollman)      esp      5  |  Mulliken           mul     6");
		printf
			("\n		Gasteiger          gas      7  |  Read in charge     rc      8");
		printf
			("\n		Write out charge   wc       9");
		printf
			("\n		----------------------------------------------------------------\n");
	} else {
		printf("Usage: antechamber -i   input file name\n");
		printf("                   -o   output file name\n");
		printf("                   -fi  input file format\n");
		printf("                   -fo  output file format\n");
		printf("                   -a   additional file name\n");
		printf("                   -fa  additional file format\n");
		printf("                   -ao  additional file options\n");
		printf("                        crd:  only read in coordinates\n");
		printf("                        crg:  only read in charges\n");
		printf("                        name: only read in atom names \n");
		printf("                        type: only read in atom type \n");
		printf("                        bond: only read in bond type \n");
		printf("                   -c   charge method\n");
		printf("                   -cf  charge filename\n");
		printf("                   -nc  net molecular charge (int)\n");
		printf
			("                   -m   mulitiplicity (2S+1), default is 1\n");
		printf
			("                   -rn  residue name, if not available in the input file, default is MOL\n");
		printf
			("                   -rf  residue topology file name in prep input file, default is molecule.res\n");
		printf
			("                   -ch  check file name in gaussian input file, default is molecule\n");
		printf
			("                   -mk  divcon or mopac keyword in a pair of quotation marks\n");
		printf
			("                   -gk  gaussian keyword in a pair of quotation marks\n");
		printf
			("                   -at  atom type: can be gaff, amber, bcc and sybyl, default is gaff\n");
		printf
		        ("	             -df  use divcon flag, 1 - always use divcon if $AMBERHOME is set (the default); 0 - not use divcon if $DOCK_HOME is set\n");
		printf
			("                   -du  check atom name duplications, can be yes(y) or no(n), default is no\n");
		printf
			("                   -j   atom type and bond type prediction index, default is 4\n");
		printf("                        0:  no assignment \n");
		printf("                        l  : atom type \n");
		printf("                        2  : full  bond types \n");
		printf("                        3  : part  bond types \n");
		printf("                        4  : atom and full bond type\n");
		printf("                        5  : atom and part bond type \n");
		printf
			("                   -s   status information, can be 0 (brief) ,1 (the default) and 2(verbose)\n");
		printf
			("                   -pf  remove the intermediate files: can be yes (y) or no (n), default is no\n");
		printf
			("                   -i -o -fi and -fo must appear in command lines and the others are optional");
		printf
			("\n\n	         	    List of the File Formats\n");
		printf
			("\n	 	file format type  abbre. index | file format type abbre. index");
		printf
                        ("\n            --------------------------------------------------------------- ");
                printf
                        ("\n            Antechamber        ac       1  | Sybyl Mol2         mol2    2 ");
                printf
                        ("\n            PDB                pdb      3  | Modifiled PDB      mpdb    4 ");
                printf
                        ("\n            AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 ");
                printf
                        ("\n            Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 ");
                printf
                        ("\n            Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 ");
                printf
                        ("\n            Gaussian Output    gout    11  | Mopac Output       mopout 12 ");
                printf
                        ("\n            Alchemy            alc     13  | CSD                csd    14 ");
                printf
                        ("\n            MDL                mdl     15  | Hyper              hin    16 ");
                printf
                        ("\n            AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 ");
                printf
                        ("\n            Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20");
                printf
                        ("\n            Divcon Input       divcrt  21  | Divcon Output      divout 22");
                printf
                        ("\n            --------------------------------------------------------------\n");
                printf
                        ("\n                AMBER restart file can only be read in as additional file\n");
                printf
                        ("\n                           List of the Charge Methods \n");
                printf
                        ("\n            charge method     abbre.  index | charge method      abbre. index");
                printf
                        ("\n            ----------------------------------------------------------------  ");
                printf
                        ("\n            RESP               resp     1  |  AM1-BCC            bcc     2");
                printf
                        ("\n            CM1                cm1      3  |  CM2                cm2     4");
                printf
                        ("\n            ESP (Kollman)      esp      5  |  Mulliken           mul     6");
                printf
                        ("\n            Gasteiger          gas      7  |  Read in charge     rc      8");
                printf
                        ("\n            Write out charge   wc       9");
                printf
                        ("\n            ----------------------------------------------------------------\n");
	}
}

void memory(int flag, int maxatom, int maxbond, int maxring)
{
	if (flag == 0) {
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stderr, "memory allocation error for *atom\n");
			exit(0);
		}
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stderr, "memory allocation error for *arom\n");
			exit(0);
		}
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stderr, "memory allocation error for *bond\n");
			exit(0);
		}
	}
/*flag = 1  <->atom
       = 2  <->bond
       = 3  <->arom
       = 4  <->atom + bond
       = 5  <->atom + arom 
       = 6  <->bond + arom
       = 7  <->atom + arom +bond
*/
	if (flag == 1 || flag == 4 || flag == 5 || flag == 7) {
		free(atom);
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stderr, "memory allocation error for *atom\n");
			exit(0);
		}
	}
	if (flag == 2 || flag == 4 || flag == 6 || flag == 7) {
		free(bond);
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stderr, "memory allocation error for *bond\n");
			exit(0);
		}
	}
	if (flag == 3 || flag == 5 || flag == 6 || flag == 7) {
		free(arom);
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stderr, "memory allocation error for *arom\n");
			exit(0);
		}
	}
	if (flag == 8) {
		free(ring);
		ring = (RING *) malloc(sizeof(RING) * maxring);
		if (ring == NULL) {
			fprintf(stderr, "memory allocation error for *ring\n");
			exit(0);
		}
	}
}

void judgebondtype(int atomnum, ATOM * atom, int bondnum, BOND * bond,
				   CONTROLINFO cinfo, MOLINFO minfo, int bondtype_flag)
{
	char tmpchar[MAXCHAR];
	char *system_env;
	int status = 0;

	wac("ANTECHAMBER_BOND_TYPE.AC0", atomnum, atom, bondnum, bond, cinfo,
		minfo);
	system_env = (char *) getenv("DOCK_HOME");
	if (system_env != NULL) {
		tmpchar[0] = '\0';
		strcpy(tmpchar, system_env);
		strcat(tmpchar, "/bin/bondtype");
	} else
		strcpy(tmpchar, "bondtype");

	strcat(tmpchar,
		   " -i ANTECHAMBER_BOND_TYPE.AC0 -o ANTECHAMBER_BOND_TYPE.AC -f ac -j ");
	if (bondtype_flag == 1)
		strcat(tmpchar, "part");
	else
		strcat(tmpchar, "full");

	if (cinfo.intstatus == 2)
		fprintf(stdout, "Running: %s\n", tmpchar);
	status = system(tmpchar);
        if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" in judgebondtype() of antechamber.c properly, exit\n", tmpchar);
                exit(0);
        }
	rac("ANTECHAMBER_BOND_TYPE.AC", &atomnum, atom, &bondnum, bond, &cinfo,
		&minfo);
}

int read_at(char *filename)  {
	FILE *fpin;
	int num = 0;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char line[MAXCHAR];

        if ((fpin = fopen(filename, "r")) == NULL) {
                fprintf(stderr, "Cannot open file %s, exit\n", filename);
                return 0;
        }
        for (;;) {
                if (fgets(line, LINELEN_MAX, fpin) == NULL)
                        break;
                        sscanf(line, "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
                        strcpy(atom[num++].ambername, tmpchar3);
        }
	fclose(fpin);
	return 0;
}

int read_bt(char *filename)  {
	FILE *fpin;
	int num = 0;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char line[MAXCHAR];
        if ((fpin = fopen(filename, "r")) == NULL) {
                fprintf(stderr, "Cannot open file %s, exit\n", filename);
                return 0;
        }
        for (;;) {
                if (fgets(line, LINELEN_MAX, fpin) == NULL)
                        break;
                        sscanf(line, "%s%s%s%ld", tmpchar1, tmpchar2, tmpchar3, &bond[num++].type);
        }
	fclose(fpin);
	return 0;

}

int write_at(char *filename)  {
	FILE *fpout;
	int i;
        if ((fpout = fopen(filename, "w")) == NULL) {
                fprintf(stderr, "Cannot open file %s, exit\n", filename);
                return 0;
        }
	for(i=0;i<atomnum;i++)
		fprintf(fpout, "%5d %5s %5s\n", i+1, atom[i].name, atom[i].ambername);
	fclose(fpout);
	return 0;

}

int write_bt(char *filename)  {
	FILE *fpout;
	int i;
        if ((fpout = fopen(filename, "w")) == NULL) {
                fprintf(stderr, "Cannot open file %s, exit\n", filename);
                return 0;
        }
	for(i=0;i<bondnum;i++)
		fprintf(fpout, "%5d %5d %5d %5d\n", i+1, bond[i].bondi+1, bond[i].bondj, bond[i].type);
	fclose(fpout);
	return 0;

}

int main(int argc, char *argv[])
{
	int i;
	int index;
	int status = 0;
	char tmpchar[MAXCHAR];

	if (argc == 2)
		if (strncmp(argv[1], "-h", 2) == 0
			|| strncmp(argv[1], "-H", 2) == 0) {
			usage();
			exit(0);
		}
	if (argc == 1) {
		usage();
		exit(0);
	}

/* 	set defaults information */
	default_minfo(&minfo);
	default_cinfo(&cinfo);
	atomtype_flag = 0;
	bondtype_flag = 0;
	default_flag = 0;
	atomname_flag = 0;
	atomicnum_flag = 0;
	adjustatomname_flag = 0;
	duplicatedname_flag = 1;
	cartcoord_flag = 0;
	connect_flag = 0;


	index = 0;
	for (i = 1; i < argc - 1; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			index++;
			strcpy(ifilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-o") == 0) {
			index++;
			strcpy(ofilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-fi") == 0) {
			index++;
			strcpy(cinfo.intype, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-j") == 0) {
			cinfo.prediction_index = atoi(argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-fo") == 0) {
			index++;
			strcpy(cinfo.outtype, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-a") == 0) {
			strcpy(afilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-fa") == 0) {
			strcpy(cinfo.atype, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-c") == 0) {
			strcpy(cinfo.chargetype, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-m") == 0) {
			minfo.multiplicity = atoi(argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-s") == 0) {
			cinfo.intstatus = atoi(argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-ra") == 0) {
			strcpy(at_filename, argv[i + 1]);
			ra_flag = 1;
			continue;
		} else if (strcmp(argv[i], "-rb") == 0) {
			strcpy(bt_filename, argv[i + 1]);
			rb_flag = 1;
			continue;
		} else if (strcmp(argv[i], "-wa") == 0) {
			strcpy(at_filename, argv[i + 1]);
			wa_flag = 1;
			continue;
		} else if (strcmp(argv[i], "-wb") == 0) {
			strcpy(bt_filename, argv[i + 1]);
			wb_flag = 1;
			continue;
		} else if (strcmp(argv[i], "-ao") == 0) {
			if (strcmp(argv[i + 1], "crd") == 0
				|| strcmp(argv[i + 1], "CRD") == 0) {
				ao_flag = 1;
				continue;
			}
			if (strcmp(argv[i + 1], "crg") == 0
				|| strcmp(argv[i + 1], "CRG") == 0) {
				ao_flag = 2;
				continue;
			}
			if (strcmp(argv[i + 1], "name") == 0
				|| strcmp(argv[i + 1], "NAME") == 0) {
				ao_flag = 3;
				continue;
			}
			if (strcmp(argv[i + 1], "type") == 0
				|| strcmp(argv[i + 1], "NAME") == 0) {
				ao_flag = 4;
				continue;
			}
			if (strcmp(argv[i + 1], "bond") == 0
				|| strcmp(argv[i + 1], "BOND") == 0) {
				ao_flag = 5;
				continue;
			}
		} else if (strcmp(argv[i], "-cf") == 0) {
			strcpy(cfilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-nc") == 0) {
			minfo.usercharge = atoi(argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-rn") == 0) {
			strcpy(minfo.longresname, argv[i + 1]);
			strncpy(minfo.resname, argv[i + 1], 3);
			cinfo.rnindex = 1;
			continue;
		} else if (strcmp(argv[i], "-ch") == 0) {
			strcpy(minfo.chkfile, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-at") == 0) {
			strcpy(minfo.atom_type_def, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-du") == 0) {
			if (strcmp(argv[i + 1], "Yes") == 0
				|| strcmp(argv[i + 1], "Y") == 0
				|| strcmp(argv[i + 1], "yes") == 0
				|| strcmp(argv[i + 1], "y") == 0) 
				duplicatedname_flag = 1;
			if (strcmp(argv[i + 1], "NO") == 0
				|| strcmp(argv[i + 1], "No") == 0
				|| strcmp(argv[i + 1], "N") == 0
				|| strcmp(argv[i + 1], "no") == 0
				|| strcmp(argv[i + 1], "n") == 0) 
				duplicatedname_flag = 0;
			continue;
		} else if (strcmp(argv[i], "-rf") == 0) {
			strcpy(minfo.resfilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-pf") == 0) {
			if (strcmp(argv[i + 1], "yes") == 0
				|| strcmp(argv[i + 1], "YES") == 0
				|| strcmp(argv[i + 1], "Y") == 0
				|| strcmp(argv[i + 1], "y") == 0
				|| strcmp(argv[i + 1], "Yes") == 0)
				cinfo.pfindex = 1;
			if (strcmp(argv[i + 1], "no") == 0
				|| strcmp(argv[i + 1], "NO") == 0
				|| strcmp(argv[i + 1], "N") == 0
				|| strcmp(argv[i + 1], "n") == 0
				|| strcmp(argv[i + 1], "No") == 0)
				cinfo.pfindex = 0;
			continue;
		} else if (strcmp(argv[i], "-mk") == 0) {
			strcpy(minfo.mkeyword, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-gk") == 0) {
			strcpy(minfo.gkeyword, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-df") == 0) {
			divcon_flag = atoi(argv[i + 1]);
			continue;
		} else {
			fprintf(stderr, "Flag not recognized: %s\n", argv[i]);
			fprintf(stderr,
					"Use antechamber -h for command-line syntax\n");
			exit(1);
		}
	}

	if (index != 4 && (strcmp(cinfo.chargetype, "wc") != 0)) {
		fprintf(stderr, "Need both input and output files & formats\n");
		fprintf(stderr, "Use antechamber -h for command-line syntax\n");
		exit(1);
	}

/* 	option of using the divcon program */

	if(divcon_flag == 0)
/*minfo.envtype == 3: both $AMBERHOME and $DOCK_HOME are set, use divcon*/
		if(minfo.envtype == 3) {
/*minfo.envtype == 4: both $AMBERHOME and $DOCK_HOME are set, use mopac.sh*/
			minfo.envtype = 4; 
			strcpy(minfo.mkeyword, "AM1 MMOK GEO-OK PRECISE");
		}

/* 	for connect.tpl and radius parameter files */
	system_env = (char *) getenv("DOCK_HOME");
	if (system_env != NULL) {
		minfo.connect_file[0] = '\0';
		strcpy(minfo.connect_file, system_env);
		strcat(minfo.connect_file, "/parameters/antechamber/CONNECT.TPL");
		minfo.radius_file[0] = '\0';
		strcpy(minfo.radius_file, system_env);
		strcat(minfo.radius_file, "/parameters/antechamber/RADIUS.DAT");
	}


/*      allocate memory using default parameters MAXATOM and MAXBOND */
	memory(0, MAXATOM, MAXBOND, MAXRING);

/******************************************/
/* 	The following codes readin input file */
/******************************************/

	if (strcmp("ac", cinfo.intype) == 0 || strcmp("1", cinfo.intype) == 0) {
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo,
					&minfo);
		}
		adjustatomname_flag = 1;
		atomicnum_flag = 1;
	}

	if (strcmp("mol2", cinfo.intype) == 0
		|| strcmp("2", cinfo.intype) == 0) {
		overflow_flag =
			rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo,
				  0);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo,
					  &minfo, 0);
		}
		default_flag = 1;
		atomicnum_flag = 1;
		adjustatomname_flag = 1;
	}
	if (strcmp("mopint", cinfo.intype) == 0
		|| strcmp("9", cinfo.intype) == 0) {
		overflow_flag = rmopint(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmopint(ifilename, &atomnum, atom, cinfo, minfo);
		}
		cartcoord_flag = 1;
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;

	}

	if (strcmp("mopcrt", cinfo.intype) == 0
		|| strcmp("10", cinfo.intype) == 0) {
		overflow_flag = rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("mopout", cinfo.intype) == 0
		|| strcmp("12", cinfo.intype) == 0) {
		overflow_flag = rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("gcrt", cinfo.intype) == 0
		|| strcmp("8", cinfo.intype) == 0) {
		overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
		}
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}
	if (strcmp("gzmat", cinfo.intype) == 0
		|| strcmp("7", cinfo.intype) == 0) {
		overflow_flag = rgzmat(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rgzmat(ifilename, &atomnum, atom, cinfo, minfo);
		}
		cartcoord_flag = 1;
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("gout", cinfo.intype) == 0
		|| strcmp("11", cinfo.intype) == 0) {
		overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
		}
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("jcrt", cinfo.intype) == 0
		|| strcmp("18", cinfo.intype) == 0) {
		overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
		}
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

        if (strcmp("jzmat", cinfo.intype) == 0
                || strcmp("19", cinfo.intype) == 0) {
                overflow_flag = rjzmat(ifilename, &atomnum, atom, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rjzmat(ifilename, &atomnum, atom, cinfo, minfo);
                }
                cartcoord_flag = 1;
                atomicnum_flag = 1;
                atomname_flag = 0;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }

        if (strcmp("jout", cinfo.intype) == 0
                || strcmp("20", cinfo.intype) == 0) {
                overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
                }
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }

	if (strcmp("pdb", cinfo.intype) == 0 || strcmp("3", cinfo.intype) == 0) {
		overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
		}
		adjustatomname_flag = 1;
		atomicnum_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("mpdb", cinfo.intype) == 0
		|| strcmp("4", cinfo.intype) == 0) {
		overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
		}
		adjustatomname_flag = 1;
		atomicnum_flag = 1;
		default_flag = 0;
		connect_flag = 1;
		bondtype_flag = 2;
	}
	if (strcmp("csd", cinfo.intype) == 0
		|| strcmp("14", cinfo.intype) == 0) {
		overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
		}
		atomicnum_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("mdl", cinfo.intype) == 0
		|| strcmp("15", cinfo.intype) == 0) {
		overflow_flag =
			rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo,
					 minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		bondtype_flag = 1;
	}

	if (strcmp("alc", cinfo.intype) == 0
		|| strcmp("13", cinfo.intype) == 0) {
		overflow_flag =
			ralc(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				ralc(ifilename, &atomnum, atom, &bondnum, bond, cinfo,
					 minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		bondtype_flag = 1;
		default_flag = 1;
	}

	if (strcmp("hin", cinfo.intype) == 0
		|| strcmp("16", cinfo.intype) == 0) {
		overflow_flag =
			rhin(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rhin(ifilename, &atomnum, atom, &bondnum, bond, cinfo,
					 minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		bondtype_flag = 1;
		default_flag = 1;
	}

	if (strcmp("prepi", cinfo.intype) == 0
		|| strcmp("5", cinfo.intype) == 0) {
		overflow_flag = rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		}
		atomicnum_flag = 1;
		bondtype_flag = 2;
		connect_flag = 0;
	}

	if (strcmp("prepc", cinfo.intype) == 0
		|| strcmp("5", cinfo.intype) == 0) {
		overflow_flag = rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		}
		atomicnum_flag = 1;
		bondtype_flag = 2;
		connect_flag = 1;
	}

	if (strcmp("rst", cinfo.intype) == 0
		|| strcmp("17", cinfo.intype) == 0) {
		fprintf(stderr,
				"RST (17) file format can only be additional file because it only has coordinate information\n");
		exit(0);
	}

        if (strcmp("divcrt", cinfo.intype) == 0
                || strcmp("21", cinfo.intype) == 0) {
                overflow_flag = rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
                }
                atomicnum_flag = 1;
                atomname_flag = 1;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }
                                                                                                                                                                                                           
        if (strcmp("divout", cinfo.intype) == 0
                || strcmp("22", cinfo.intype) == 0) {
                overflow_flag = rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
                }
                atomicnum_flag = 1;
                atomname_flag = 1;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }


/*****************************************************************************/
/*	assign atomtype_flag, bondtype_flag and charge_flag according to -j flag */
/*****************************************************************************/

	if (cinfo.prediction_index == 0) {
		atomtype_flag = 0;
		bondtype_flag = 0;
	}
	if (cinfo.prediction_index == 1)
		atomtype_flag = 1;
	if (cinfo.prediction_index == 2)
		bondtype_flag = 2;
	if (cinfo.prediction_index == 3)
		bondtype_flag = 1;
	if (cinfo.prediction_index == 4) {
		atomtype_flag = 1;
		bondtype_flag = 2;
	}
	if (cinfo.prediction_index == 5) {
		atomtype_flag = 1;
		bondtype_flag = 1;
	}

/*	reassign the connect_flag according to the output types 	*/
        if (strcmp("mopcrd", cinfo.outtype) == 0
                || strcmp("mopout", cinfo.outtype) == 0 
                || strcmp("gcrt", cinfo.outtype) == 0 
                || strcmp("gout", cinfo.outtype) == 0 
                || strcmp("jcrt", cinfo.outtype) == 0 
                || strcmp("jout", cinfo.outtype) == 0 
                || strcmp("pdb", cinfo.outtype) == 0 
                || strcmp("rst", cinfo.outtype) == 0 
                || strcmp("3", cinfo.outtype) == 0
                || strcmp("8", cinfo.outtype) == 0
                || strcmp("10", cinfo.outtype) == 0
                || strcmp("11", cinfo.outtype) == 0
                || strcmp("12", cinfo.outtype) == 0
                || strcmp("17", cinfo.outtype) == 0
                || strcmp("18", cinfo.outtype) == 0
                || strcmp("20", cinfo.outtype) == 0) {
		connect_flag = 0;
		bondtype_flag = 0;
		atomtype_flag = 0;
	}

        if (strcmp("mopout", cinfo.outtype) == 0
                || strcmp("gout", cinfo.outtype) == 0
                || strcmp("jout", cinfo.outtype) == 0
                || strcmp("rst", cinfo.outtype) == 0
                || strcmp("11", cinfo.outtype) == 0
                || strcmp("12", cinfo.outtype) == 0
                || strcmp("17", cinfo.outtype) == 0
                || strcmp("20", cinfo.outtype) == 0)
                duplicatedname_flag = 0;

        if (strcmp("mopint", cinfo.outtype) == 0
                || strcmp("gzmat", cinfo.outtype) == 0 
                || strcmp("jzmat", cinfo.outtype) == 0 
                || strcmp("7", cinfo.outtype) == 0 
                || strcmp("9", cinfo.outtype) == 0 
                || strcmp("19", cinfo.outtype) == 0 ) {
		bondtype_flag = 0;
		atomtype_flag = 0;
	}
/*     	the following code judge or assign atom name, atom type, bond type etc according to flags */
	if (adjustatomname_flag) {
		if(strcmp(cinfo.intype, "mol2")==0 || strcmp(cinfo.intype, "2")==0 ||
		   strcmp(cinfo.intype, "ac")==0 || strcmp(cinfo.intype, "1")==0)
			adjustatomname(atomnum, atom, 1);
		else
			adjustatomname(atomnum, atom, 0);
	}
	if (atomicnum_flag)
		atomicnum(atomnum, atom);
	if (atomname_flag)
		atomname(atomnum, atom);
	if (default_flag)
		default_inf(atomnum, atom, default_flag);
	if (cartcoord_flag)
		cartcoord(atomnum, atom);
	if (connect_flag) {
		overflow_flag =
			connect(minfo.connect_file, atomnum, atom, &bondnum, bond,
					cinfo.maxbond);
		if (overflow_flag) {
			cinfo.maxbond = bondnum + 10;
			memory(2, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				connect(minfo.connect_file, atomnum, atom, &bondnum, bond,
						cinfo.maxbond);
		}
	}
	if (bondtype_flag && bondnum > 0) {
		judgebondtype(atomnum, atom, bondnum, bond, cinfo, minfo,
					  bondtype_flag);
		if(cinfo.prediction_index ==2||cinfo.prediction_index ==3) cinfo.prediction_index = 0;
		if(cinfo.prediction_index ==4||cinfo.prediction_index ==5) cinfo.prediction_index = 1;
	}

	if (duplicatedname_flag)
		duplicatedname(atomnum, atom);
	if (atomtype_flag) {
		wac("ANTECHAMBER_AC.AC0", atomnum, atom, bondnum, bond, cinfo,
			minfo);
		system_env = (char *) getenv("DOCK_HOME");
		if (system_env != NULL) {
			tmpchar[0] = '\0';
			strcpy(tmpchar, system_env);
			strcat(tmpchar, "/bin/atomtype");
		} else
			strcpy(tmpchar, "atomtype");
		strcat(tmpchar, " -i ANTECHAMBER_AC.AC0 -o ANTECHAMBER_AC.AC -p ");
		strcat(tmpchar, minfo.atom_type_def);
		if (cinfo.intstatus == 2)
			fprintf(stdout, "\nRunning: %s\n", tmpchar);
		status = system(tmpchar);
	        if(status != 0) {
                	fprintf(stderr, "Error: cannot run \"%s\" in main() of antechamber.c properly, exit\n", tmpchar);
                	exit(0);
       		}
		rac("ANTECHAMBER_AC.AC", &atomnum, atom, &bondnum, bond, &cinfo,
			&minfo);
	}

/*     the following code readin or calculate charges */
/* 	usercharge info*/
	if (minfo.usercharge > -9999) {	/*read in charge with -nc flag */
		minfo.icharge = minfo.usercharge;
		minfo.dcharge = minfo.usercharge;
	}

	if (strcmp("resp", cinfo.chargetype) == 0
		|| strcmp("1", cinfo.chargetype) == 0)
		resp(ifilename, atomnum, atom, bondnum, bond, cinfo, minfo);
	if (strcmp("bcc", cinfo.chargetype) == 0
		|| strcmp("2", cinfo.chargetype) == 0)
		bcc(ifilename, atomnum, atom, bondnum, bond, arom, &cinfo, &minfo);
	if (strcmp("cm1", cinfo.chargetype) == 0
		|| strcmp("3", cinfo.chargetype) == 0)
		cm1(atomnum, atom, &cinfo, &minfo);
	if (strcmp("cm2", cinfo.chargetype) == 0
		|| strcmp("4", cinfo.chargetype) == 0)
		cm2(atomnum, atom, &cinfo, &minfo);
	if (strcmp("esp", cinfo.chargetype) == 0
		|| strcmp("5", cinfo.chargetype) == 0)
		esp(ifilename, atomnum, atom, cinfo, minfo);
	if (strcmp("mul", cinfo.chargetype) == 0
		|| strcmp("6", cinfo.chargetype) == 0)
		mul(ifilename, atomnum, atom, &cinfo, &minfo);
	if (strcmp("gas", cinfo.chargetype) == 0
		|| strcmp("7", cinfo.chargetype) == 0)
		gascharge(atomnum, atom, bondnum, bond, cinfo, &minfo);
	if (strcmp("rc", cinfo.chargetype) == 0
		|| strcmp("8", cinfo.chargetype) == 0)
		rcharge(cfilename, atomnum, atom, cinfo, &minfo);
	if (strcmp("wc", cinfo.chargetype) == 0
		|| strcmp("9", cinfo.chargetype) == 0)
		wcharge(cfilename, atomnum, atom, cinfo, minfo);

/*	judge the radii*/
	if (strcmp("mpdb", cinfo.intype) != 0
		&& strcmp("4", cinfo.intype) != 0)
		read_radius(minfo.radius_file, atomnum, atom);

/*    	read in additional files*/
/*	expand the array size a little bit; 
  	expand the array size a little bit 
  	since atom file type such as prepi may have additional atom records  
*/
	if (afilename && cinfo.atype) {
		cinfo.maxatom = atomnum + 10;
		cinfo.maxbond = bondnum + 10;
		atom_tmp = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom_tmp == NULL) {
			fprintf(stderr, "memory allocation error for *atom_tmp\n");
			exit(0);
		}
		bond_tmp = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
		if (bond_tmp == NULL) {
			fprintf(stderr, "memory allocation error for *bond_tmp\n");
			exit(0);
		}

		if (strcmp("ac", cinfo.atype) == 0
			|| strcmp("1", cinfo.atype) == 0) {
			overflow_flag =
				rac(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					bond_tmp, &cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}
		if (strcmp("mol2", cinfo.atype) == 0
			|| strcmp("2", cinfo.atype) == 0) {
			overflow_flag =
				rmol2(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					  bond_tmp, &cinfo, &minfo, 0);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}
		if (strcmp("mopint", cinfo.atype) == 0
			|| strcmp("9", cinfo.atype) == 0) {
			overflow_flag =
				rmopint(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}
		if (strcmp("mopcrt", cinfo.atype) == 0
			|| strcmp("10", cinfo.atype) == 0) {
			overflow_flag =
				rmopcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("mopout", cinfo.atype) == 0
			|| strcmp("12", cinfo.atype) == 0) {
			overflow_flag =
				rmopout(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}
		if (strcmp("gcrt", cinfo.atype) == 0
			|| strcmp("8", cinfo.atype) == 0) {
			overflow_flag =
				rgcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}
		if (strcmp("gzmat", cinfo.atype) == 0
			|| strcmp("7", cinfo.atype) == 0) {
			overflow_flag =
				rgzmat(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("gout", cinfo.atype) == 0
			|| strcmp("11", cinfo.atype) == 0) {
			overflow_flag =
				rgout(afilename, &atomnum_tmp, atom_tmp, cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("pdb", cinfo.atype) == 0
			|| strcmp("3", cinfo.atype) == 0) {
			overflow_flag =
				rpdb(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo, 0);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("mpdb", cinfo.atype) == 0
			|| strcmp("4", cinfo.atype) == 0) {
			overflow_flag =
				rpdb(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo, 1);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("csd", cinfo.atype) == 0
			|| strcmp("14", cinfo.atype) == 0) {
			overflow_flag =
				rcsd(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("mdl", cinfo.atype) == 0
			|| strcmp("15", cinfo.atype) == 0) {
			overflow_flag =
				rmdl(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					 bond_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("alc", cinfo.atype) == 0
			|| strcmp("13", cinfo.atype) == 0) {
			overflow_flag =
				ralc(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					 bond_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("hin", cinfo.atype) == 0
			|| strcmp("16", cinfo.atype) == 0) {
			overflow_flag =
				rhin(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp,
					 bond_tmp, cinfo, minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("prepi", cinfo.atype) == 0
			|| strcmp("5", cinfo.atype) == 0) {
			overflow_flag =
				rprepi(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, &cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}

		if (strcmp("prepc", cinfo.atype) == 0
			|| strcmp("6", cinfo.atype) == 0) {
			overflow_flag =
				rprepc(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}
		if (strcmp("rst", cinfo.atype) == 0
			|| strcmp("17", cinfo.atype) == 0) {
			overflow_flag = rrst(afilename, &atomnum_tmp, atom_tmp, cinfo);
			if (overflow_flag) {
				fprintf(stderr,
						"Overflow happens for additional files, exit");
				exit(0);
			}
		}
                if (strcmp("divcrt", cinfo.atype) == 0
                        || strcmp("10", cinfo.atype) == 0) {
                        overflow_flag =
                                rdivcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
                        if (overflow_flag) {
                                fprintf(stderr,
                                                "Overflow happens for additional files, exit");
                                exit(0);
                        }
                }
                                                                                                                                                                                                           
                if (strcmp("divout", cinfo.atype) == 0
                        || strcmp("12", cinfo.atype) == 0) {
                        overflow_flag =
                                rdivout(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
                        if (overflow_flag) {
                                fprintf(stderr,
                                                "Overflow happens for additional files, exit");
                                exit(0);
                        }
                }
		if (ao_flag == 1)
			for (i = 0; i < atomnum; i++) {
				atom[i].x = atom_tmp[i].x;
				atom[i].y = atom_tmp[i].y;
				atom[i].z = atom_tmp[i].z;
			}

		if (ao_flag == 2)
			for (i = 0; i < atomnum; i++)
				atom[i].charge = atom_tmp[i].charge;

		if (ao_flag == 3)
			for (i = 0; i < atomnum; i++)
				strcpy(atom[i].name, atom_tmp[i].name);
		if (ao_flag == 4)
			for (i = 0; i < atomnum; i++)
				strcpy(atom[i].ambername, atom_tmp[i].ambername);
		if (ao_flag == 5)
			for (i = 0; i < bondnum; i++)
				bond[i].type = bond_tmp[i].type;
		free(atom_tmp);
		free(bond_tmp);
	}

	if(atomtype_flag == 0 && (strcmp("charmm", cinfo.outtype) == 0||
	                         strcmp("ac", cinfo.outtype) == 0 || strcmp("1", cinfo.outtype) == 0 ||
	                         strcmp("mol2", cinfo.outtype) == 0 || strcmp("2", cinfo.outtype) == 0 ||
	                         strcmp("csd", cinfo.outtype) == 0 || strcmp("14", cinfo.outtype) == 0 ||
	                         strcmp("mdl", cinfo.outtype) == 0 || strcmp("15", cinfo.outtype) == 0 ||
	                         strcmp("alc", cinfo.outtype) == 0 || strcmp("13", cinfo.outtype) == 0 ||
	                         strcmp("hin", cinfo.outtype) == 0 || strcmp("16", cinfo.outtype) == 0)) {
        	memory(8, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
        	overflow_flag =
			ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
			  	 cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 0);
        	if (overflow_flag) {
			cinfo.maxring = ringnum + 10;
                	memory(8, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
        		overflow_flag =
				ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
			   		cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 0);
		}
	}
	if(ra_flag==1)	
		read_at(at_filename);
	if(rb_flag==1)	
		read_bt(bt_filename);

/*      write out files */
	if (strcmp("ac", cinfo.outtype) == 0
		|| strcmp("1", cinfo.outtype) == 0)
		wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
	if (strcmp("charmm", cinfo.outtype) == 0)
                wcharmm(ofilename, ifilename, atomnum, atom, bondnum, bond, cinfo, minfo);
	if (strcmp("mol2", cinfo.outtype) == 0
		|| strcmp("2", cinfo.outtype) == 0) 
		wmol2(ofilename, atomnum, atom, bondnum, bond, arom, cinfo, minfo);
	if (strcmp("mopint", cinfo.outtype) == 0
		|| strcmp("9", cinfo.outtype) == 0)
		wmopint(ofilename, atomnum, atom, minfo);
	if (strcmp("mopcrt", cinfo.outtype) == 0
		|| strcmp("10", cinfo.outtype) == 0)
		wmopcrt(ofilename, atomnum, atom, minfo);
	if (strcmp("mopout", cinfo.outtype) == 0
		|| strcmp("12", cinfo.outtype) == 0)
		wmopout(ofilename, atomnum, atom, cinfo, minfo);
	if (strcmp("gcrt", cinfo.outtype) == 0
		|| strcmp("8", cinfo.outtype) == 0)
		wgcrt(ofilename, atomnum, atom, minfo);
	if (strcmp("gzmat", cinfo.outtype) == 0
		|| strcmp("7", cinfo.outtype) == 0)
		wgzmat(ofilename, atomnum, atom, minfo);
	if (strcmp("gout", cinfo.outtype) == 0
		|| strcmp("11", cinfo.outtype) == 0)
		wgout();
	if (strcmp("jcrt", cinfo.outtype) == 0
		|| strcmp("18", cinfo.outtype) == 0)
		wjcrt(ofilename, atomnum, atom, minfo);
	if (strcmp("jzmat", cinfo.outtype) == 0
		|| strcmp("19", cinfo.outtype) == 0)
		wjzmat(ofilename, atomnum, atom, minfo);
	if (strcmp("jout", cinfo.outtype) == 0
		|| strcmp("20", cinfo.outtype) == 0)
		wjout();
	if (strcmp("pdb", cinfo.outtype) == 0
		|| strcmp("3", cinfo.outtype) == 0)
		wpdb(ofilename, atomnum, atom);
	if (strcmp("mpdb", cinfo.outtype) == 0
		|| strcmp("4", cinfo.outtype) == 0)
		wmpdb(ofilename, atomnum, atom);
	if (strcmp("csd", cinfo.outtype) == 0
		|| strcmp("14", cinfo.outtype) == 0)
		wcsd(ofilename, atomnum, atom);
	if (strcmp("mdl", cinfo.outtype) == 0
		|| strcmp("15", cinfo.outtype) == 0)
		wmdl(ofilename, atomnum, atom, bondnum, bond, cinfo);
	if (strcmp("alc", cinfo.outtype) == 0
		|| strcmp("13", cinfo.outtype) == 0)
		walc(ofilename, atomnum, atom, bondnum, bond);
	if (strcmp("hin", cinfo.outtype) == 0
		|| strcmp("16", cinfo.outtype) == 0)
		whin(ofilename, atomnum, atom, bondnum, bond);
	if (strcmp("prepi", cinfo.outtype) == 0
		|| strcmp("5", cinfo.outtype) == 0)
		wprep(ofilename, ifilename, atomnum, atom, bondnum, bond, cinfo,
			  &minfo, 1);
	if (strcmp("prepc", cinfo.outtype) == 0
		|| strcmp("6", cinfo.outtype) == 0)
		wprep(ofilename, ifilename, atomnum, atom, bondnum, bond, cinfo,
			  &minfo, 0);
	if (strcmp("rst", cinfo.outtype) == 0
		|| strcmp("17", cinfo.outtype) == 0)
		wrst(ofilename, atomnum, atom);
	if (strcmp("divcrt", cinfo.outtype) == 0
		|| strcmp("21", cinfo.outtype) == 0)
		wdivcrt(ofilename, atomnum, atom, minfo);
	if (strcmp("divout", cinfo.outtype) == 0
		|| strcmp("22", cinfo.outtype) == 0)
		wdivout(ofilename, atomnum, atom, cinfo, minfo);
/*	info(atomnum, atom, bondnum, bond, arom, cinfo, minfo); */
/*
        free(atom);
        free(bond);
        free(arom);
*/
	if(wa_flag==1)	
		write_at(at_filename);
	if(wb_flag==1)	
		write_bt(bt_filename);
	if (cinfo.pfindex == 1) {
		status = system("rm -f ANTECHAMBER* ATOMTYPE.INF BCCTYPE.INF NEWPDB.PDB PREP.INF");
	        if(status != 0) {
                	fprintf(stderr, "Error: cannot run \"%s\" in main() of antechamber.c properly, exit\n", "rm -f ANTECHAMBER* ATOMTYPE.INF BCCTYPE.INF NEWPDB.PDB PREP.INF");
                	exit(0);
        	}
	
	}
	printf("\n");
	return (0);
}
