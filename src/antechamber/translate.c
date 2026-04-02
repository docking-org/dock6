/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    translate                                                    *
*  Version: version 1.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "rotate.c"
# include "ac.c"
# include "pdb.c"
# include "mol2.c"
# include "prep.c"
# include "lsfit.c"

ATOM *atom;
AROM *arom;
BOND *bond;
ATOM *ref_atom;
AROM *ref_arom;
BOND *ref_bond;
ATOM tmpatom1;
ATOM tmpatom2;
RING ring[MAXRING];
MOLINFO minfo;
CONTROLINFO cinfo;
int atomnum = 0;
int bondnum = 0;
int ringnum;
char ifilename[MAXCHAR];
char rfilename[MAXCHAR];
char ofilename[MAXCHAR];
char line[MAXCHAR];
char *system_env;
int i, j, k;
int overflow_flag = 0;
int at1 = -99999;;
int at2 = -99999;
int at3 = -99999;
double vectx = -99999.0;
double vecty = -99999.0;
double vectz = -99999.0;
double coord_x1 = -99999.0;
double coord_y1 = -99999.0;
double coord_z1 = -99999.0;
double coord_x2 = -99999.0;
double coord_y2 = -99999.0;
double coord_z2 = -99999.0;
double degree = -99999.0;
double sum_coordx = 0.0;
double sum_coordy = 0.0;
double sum_coordz = 0.0;
double coordx;
double coordy;
double coordz;
double rmsd;
double w1, w2;
int iflag = 0;
int oflag = 0;
int rflag = 0;
int cflag = 0;
char command[MAXCHAR];
FILE *fp, *fpout;

int main(int argc, char *argv[]) {
	int i;
	int format;

	default_minfo(&minfo);
	default_cinfo(&cinfo);
	system_env = (char *) getenv("DOCK_HOME");
	if (system_env != NULL) {
		minfo.connect_file[0] = '\0';
		strcpy(minfo.connect_file, system_env);
		strcat(minfo.connect_file, "/parameters/antechamber/CONNECT.TPL");
	} else {
		strcpy(minfo.connect_file, "CONNECT.TPL");
	}

	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: translate -i [0m input file name (pdb, ac or mol2)\n"
				   "[31m                 -o [0m output file name\n"
				   "[31m                 -r [0m reference file name\n"
				   "[31m                 -f [0m file format\n"
				   "[31m                 -c [0m command (center, translate, rotate1, rotate2, match)\n"
				   "[31m                     center:     need -a1;\n"
				   "[31m                     translate:  need -vx, -vy and -vz;\n"
				   "[31m                     rotate1:    need -a1, -a2 and -d;\n"
				   "[31m                     rotate2:    need -x1, -y1, -z1, -x2, -y2, -z2 and -d;\n"
				   "[31m                     match:      need -r;\n"
				   "[31m                 -d [0m degree to be rotated\n"
				   "[31m                 -vx[0m x vector\n"
				   "[31m                 -vy[0m y vector\n"
				   "[31m                 -vz[0m z vector\n"
				   "[31m                 -a1[0m id of atom 1 (0 coordinate center)\n"
				   "[31m                 -a2[0m id of atom 2\n"
				   "[31m                 -x1[0m coord x for point 1\n"
				   "[31m                 -y1[0m coord y for point 1\n"
				   "[31m                 -z1[0m coord z for point 1\n"
				   "[31m                 -x2[0m coord x for point 2\n"
				   "[31m                 -y2[0m coord y for point 2\n"
				   "[31m                 -z2[0m coord z for point 2\n");
			exit(0);
		}
		if (argc < 5) {
			printf("[31mUsage: translate -i [0m input file name (pdb, ac or mol2)\n"
				   "[31m                 -o [0m output file name\n"
				   "[31m                 -r [0m reference file name\n"
				   "[31m                 -f [0m file format\n"
				   "[31m                 -c [0m command (center, translate, rotate1, rotate2, match)\n"
				   "[31m                     center:     need -a1;\n"
				   "[31m                     translate:  need -vx, -vy and -vz;\n"
				   "[31m                     rotate1:    need -a1, -a2 and -d;\n"
				   "[31m                     rotate2:    need -x1, -y1, -z1, -x2, -y2, -z2 and -d;\n"
				   "[31m                     match:      need -r;\n"
				   "[31m                 -d [0m degree to be rotated\n"
				   "[31m                 -vx[0m x vector\n"
				   "[31m                 -vy[0m y vector\n"
				   "[31m                 -vz[0m z vector\n"
				   "[31m                 -a1[0m id of atom 1 (0 coordinate center)\n"
				   "[31m                 -a2[0m id of atom 2\n"
				   "[31m                 -x1[0m coord x for point 1\n"
				   "[31m                 -y1[0m coord y for point 1\n"
				   "[31m                 -z1[0m coord z for point 1\n"
				   "[31m                 -x2[0m coord x for point 2\n"
				   "[31m                 -y2[0m coord y for point 2\n"
				   "[31m                 -z2[0m coord z for point 2\n");
			exit(0);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: translate -i  input file name (pdb, ac or mol2)\n"
				   "             -o  output file name\n"
				   "             -r  reference file name\n"
				   "             -f  file format\n"
				   "             -c  command (center, translate, rotate1, rotate2, rotate3, match)\n"
				   "		     center:     need -a1;\n"	
				   "		     translate: need -vx, -vy and -vz;\n"	
				   "		     rotate1:    need -a1, -a2 and -d;\n"	
				   "		     rotate2:    need -x1, -y1, -z1, -x2, -y2, -z2 and -d;\n"	
				   "		     match:      need -r;\n"	
				   "             -d  degree to be rotated\n"
				   "             -vx x vector\n"
				   "             -vy y vector\n"
				   "             -vz z vector\n"
				   "             -a1 id of atom 1 (0 coordinate center)\n"
				   "             -a2 id of atom 2\n"
				   "             -x1 coord x for point 1\n"
				   "             -y1 coord y for point 1\n"
				   "             -z1 coord z for point 1\n"
				   "             -x2 coord x for point 2\n"
				   "             -y2 coord y for point 2\n"
				   "             -z2 coord z for point 2\n");
			exit(0);
		}
		if (argc < 5) {
                        printf("Usage: translate -i  input file name (pdb, ac or mol2)\n"
                                   "             -o  output file name\n"
                                   "             -r  reference file name\n"
                                   "             -f  file format\n"
                                   "             -c  command (center, translate, rotate1, rotate2, rotate3, match)\n"
                                   "                 center:     need -a1;\n"         
                                   "                 translate: need -vx, -vy and -vz;\n"
                                   "                 rotate1:    need -a1, -a2 and -d;\n"
                                   "                 rotate2:    need -x1, -y1, -z1, -x2, -y2, -z2 and -d;\n"          
                                   "                 match:      need -r;\n"  
				   "             -d  degree to be rotated\n"
                                   "             -vx x vector\n"
                                   "             -vy y vector\n"
                                   "             -vz z vector\n"
                                   "             -a1 id of atom 1 (0 coordinate center)\n"
                                   "             -a2 id of atom 2\n"
                                   "             -x1 coord x for point 1\n"
                                   "             -y1 coord y for point 1\n"
                                   "             -z1 coord z for point 1\n"
                                   "             -x2 coord x for point 2\n"
                                   "             -y2 coord y for point 2\n"
                                   "             -z2 coord z for point 2\n");
                        exit(0);
		}
	}
	format = 0;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-c")== 0) { 
			strcpy(command, argv[i + 1]);
			cflag = 1;
		}
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			iflag = 1;
		}
		if (strcmp(argv[i], "-o") == 0) {
			strcpy(ofilename, argv[i + 1]);
			oflag = 1;
		}
		if (strcmp(argv[i], "-r") == 0) {
			strcpy(rfilename, argv[i + 1]);
			rflag = 1;
		}
		if (strcmp(argv[i], "-f") == 0) {
			if (strcmp(argv[i + 1], "ac") == 0)
				format = 0;
			if (strcmp(argv[i + 1], "pdb") == 0)
				format = 1;
			if (strcmp(argv[i + 1], "mol2") == 0)
				format = 2;
		}
		if (strcmp(argv[i], "-a1") == 0)  
			at1 = atoi(argv[i+1]);
		if (strcmp(argv[i], "-a2") == 0)  
			at2 = atoi(argv[i+1]);
		if (strcmp(argv[i], "-a3") == 0)  
			at3 = atoi(argv[i+1]);
		if (strcmp(argv[i], "-d") == 0)  
			degree = atof(argv[i+1]);
		if (strcmp(argv[i], "-vx") == 0)  
			vectx = atof(argv[i+1]);
		if (strcmp(argv[i], "-vy") == 0)  
			vecty = atof(argv[i+1]);
		if (strcmp(argv[i], "-vz") == 0)  
			vectz = atof(argv[i+1]);
		if (strcmp(argv[i], "-x1") == 0)  
			coord_x1 = atof(argv[i+1]);
		if (strcmp(argv[i], "-y1") == 0)  
			coord_y1 = atof(argv[i+1]);
		if (strcmp(argv[i], "-z1") == 0)  
			coord_z1 = atof(argv[i+1]);
		if (strcmp(argv[i], "-x2") == 0)  
			coord_x2 = atof(argv[i+1]);
		if (strcmp(argv[i], "-y2") == 0)  
			coord_y2 = atof(argv[i+1]);
		if (strcmp(argv[i], "-z2") == 0)  
			coord_z2 = atof(argv[i+1]);
	}

	if(cflag == 0) {
		printf("\nNo command specified, exit\n");
		exit(0);
	}
	if(iflag == 0) {
		printf("\nNo input file specified, exit\n");
		exit(0);
	}
	/*read in prep or ac file */
	atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom == NULL) {
		fprintf(stderr, "memory allocation error for *atom\n");
		exit(0);
	}
	bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond == NULL) {
		fprintf(stderr, "memory allocation error for *bond\n");
		exit(0);
	}

	if (format == 0)
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	if (format == 1)
		overflow_flag =
			 rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0); 	
	if (format == 2)
		overflow_flag =
			rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 1);
	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
		cinfo.maxbond = bondnum + 10;
		free(atom);
		free(bond);
		atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom == NULL) {
			fprintf(stderr, "memory allocation error for *atom\n");
			exit(0);
		}
		bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
		if (bond == NULL) {
			fprintf(stderr, "memory allocation error for *bond\n");
			exit(0);
		}

        	if (format == 0)
                	overflow_flag =
                        	rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        	if (format == 1)
                	overflow_flag =
                         	rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0); 
        	if (format == 2)
                	overflow_flag =
                        	rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 1);
	}

	if(strcmp(command, "center") == 0) {
		if(at1 <-99990) {
			printf("\nNo -a1 specified\n");
			exit(0);
		}
		if(at1 == 0) {
			for(i=0;i<atomnum;i++) {
				sum_coordx += atom[i].x;
				sum_coordy += atom[i].y;
				sum_coordz += atom[i].z;
			}
			coordx = sum_coordx/atomnum;
			coordy = sum_coordy/atomnum;
			coordz = sum_coordz/atomnum;
		}
		else {
			coordx = atom[at1-1].x;
			coordy = atom[at1-1].y;
			coordz = atom[at1-1].z;
		}
		for(i=0;i<atomnum;i++) {
			 atom[i].x-=coordx;
			 atom[i].y-=coordy;
			 atom[i].z-=coordz;
		}
		printf("\nThe molecule is translated: Vector X= %9.4lf, Vector Y= %9.4lf, Vector Z = %9.4lf\n", coordx, coordy, coordz);
	}	
	if(strcmp(command, "translate") == 0) {
		if(vectx <-99990) {
			printf("\nNo -vx specified\n");
			exit(0);
		}
		if(vecty <-99990) {
			printf("\nNo -vy specified\n");
			exit(0);
		}
		if(vectz <-99990) {
			printf("\nNo -vz specified\n");
			exit(0);
		}
		coordx =vectx; 
		coordy = vecty;
		coordz = vectz;
		for(i=0;i<atomnum;i++) {
			 atom[i].x-=coordx;
			 atom[i].y-=coordy;
			 atom[i].z-=coordz;
		}
		printf("\nThe molecule is translated: Vector X= %9.4lf, Vector Y= %9.4lf, Vector Z = %9.4lf\n", coordx, coordy, coordz);
	}	
	if(strcmp(command, "rotate1") == 0) {
		if(at1 <-99990) {
			printf("\nNo -a1 specified\n");
			exit(0);
		}
		if(at2 <-99990) {
			printf("\nNo -a2 specified\n");
			exit(0);
		}
		if(degree <-99990) {
			printf("\nNo -d specified\n");
			exit(0);
		}
		omegarotate(atom[at1-1], atom[at2-1], &w1, &w2);
		omegarotate2(atom, atomnum, atom[at2-1], degree, w1, w2);
	}	
	if(strcmp(command, "rotate2") == 0) {
		if(coord_x1 <-99990) {
			printf("\nNo -x1 specified\n");
			exit(0);
		}
		if(coord_y1 <-99990) {
			printf("\nNo -y1 specified\n");
			exit(0);
		}
		if(coord_z1 <-99990) {
			printf("\nNo -z1 specified\n");
			exit(0);
		}
		if(coord_x2 <-99990) {
			printf("\nNo -x2 specified\n");
			exit(0);
		}
		if(coord_y2 <-99990) {
			printf("\nNo -y2 specified\n");
			exit(0);
		}
		if(coord_z2 <-99990) {
			printf("\nNo -z2 specified\n");
			exit(0);
		}
		if(degree <-99990) {
			printf("\nNo -d specified\n");
			exit(0);
		}
		tmpatom1.x = coord_x1;
		tmpatom1.y = coord_y1;
		tmpatom1.z = coord_z1;
		tmpatom2.x = coord_x2;
		tmpatom2.y = coord_y2;
		tmpatom2.z = coord_z2;
		omegarotate(tmpatom1, tmpatom2, &w1, &w2);
		omegarotate2(atom, atomnum, tmpatom2, degree, w1, w2);
	}	
	if(strcmp(command, "match") == 0) {
		if(rflag == 0) {
			printf("\nNo -r specified\n");
			exit(0);
		}	
                ref_atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
                if (ref_atom == NULL) {
                        fprintf(stderr, "memory allocation error for *refatom\n");
                        exit(0);
                }
                ref_bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
                if (ref_bond == NULL) {
                        fprintf(stderr, "memory allocation error for *refbond\n");
                        exit(0);
                }
                if (format == 0)
                        overflow_flag =
                                rac(rfilename, &atomnum, ref_atom, &bondnum, ref_bond, &cinfo, &minfo);
                if (format == 1)
                        overflow_flag =
                                rpdb(rfilename, &atomnum, ref_atom, cinfo, minfo, 0);
                if (format == 2)
                        overflow_flag =
                                rmol2(rfilename, &atomnum, ref_atom, &bondnum, ref_bond, &cinfo, &minfo, 1);
		rmsd = lsfit(atom, ref_atom, atom, atomnum); 
		printf("\nThe rmsd is %9.4lf\n", rmsd);	
	}
	if(oflag == 1) {
		if(format == 0) wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
 		if(format == 1) wpdb(ofilename, atomnum, atom);
               	if(format == 2) wmol2(ofilename, atomnum, atom, bondnum, bond, arom, cinfo, minfo); 
	}
	return (0);
}
