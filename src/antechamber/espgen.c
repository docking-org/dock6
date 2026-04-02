/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    espgen                                                     * 
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
# define Bohr 0.529177249
# define MAXESP 10000
# define MAX_ATOM_CENTER 1000
int i, j, k;
int fail = 0;
int Found_Stationary = 0;
int espindex = 0;
int method = 0;
int tmpint1, tmpint2, tmpint3;
int index1 = 0;
int index2 = 0;
int index3 = 0;
int index0 = 0;
int maxesp = 0;
int status = 0;
char line[MAXCHAR];
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR],
	tmpchar4[MAXCHAR], tmpchar5[MAXCHAR];
char tmpchar[MAXCHAR] = "rm -f ";
double tmpfloat1, tmpfloat2, tmpfloat3, tmpfloat4;
double cord[MAX_ATOM_CENTER][3];
ESP *esp;
DM dipole;
QM quadrupole;
FILE *fpin, *fpout;

int main(int argc, char *argv[])
{
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: espgen -i  [0m input file name \n"
				   "[31m              -o  [0m output file name \n"
				   /*         "[31m              -dq [0m yes or no (optional, printing out dipole and\n"
				      "                   quadrupole moments or not, default is no\n" */
				);
			exit(0);
		}
		if (argc != 7 && argc != 5) {
			printf("[31mUsage: espgen -i  [0m input file name \n"
				   "[31m              -o  [0m output file name \n"
				   /*         "[31m              -dq [0m yes or no (optional, printing out dipole and\n"
				      "                   quadrupole moments or not, default is no\n" */
				);
			exit(0);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: espgen -i input file name\n");
			printf("              -o output file name \n");
/*
    printf("Usage: espgen -i  input_file_name  -o output file name \n");
    printf("              -dq yes or no (optional, printing out dipole and\n"); 
    printf("                  quadrupole moments or not, default is no\n");
*/
			exit(0);
		}
		if (argc != 7 && argc != 5) {
/*
    printf("Usage: espgen -i  input_file_name  -o output_file_name \n");
    printf("              -dq yes or no (optional, printing out dipole and\n"); 
    printf("                  quadrupole moments or not, default is no\n");
*/
			printf("Usage: espgen -i input file name\n");
			printf("              -o output file name \n");
			exit(0);
		}
	}

/* allocate memory for *esp */
	maxesp = MAXESP;
	esp = (ESP *) calloc(maxesp, sizeof(ESP));
	if (esp == NULL) {
		fprintf(stderr, "memory allocation error for *esp\n");
		exit(0);
	}

	method = 0;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-dq") == 0) {
			if (strcmp("YES", argv[i + 1]) == 0
				|| strcmp("yes", argv[i + 1]) == 0)
				method = 1;
			if (strcmp("Y", argv[i + 1]) == 0
				|| strcmp("y", argv[i + 1]) == 0)
				method = 1;
			if (strcmp("NO", argv[i + 1]) == 0
				|| strcmp("no", argv[i + 1]) == 0)
				method = 0;
			if (strcmp("N", argv[i + 1]) == 0
				|| strcmp("n", argv[i + 1]) == 0)
				method = 0;
		}
	}
	if ((fpin = fopen(ifilename, "r")) == NULL) {
		printf("\n Cannot open file %s, exit", ifilename);
		return (0);
	}
	if ((fpout = fopen(ofilename, "w")) == NULL) {
		printf("\n Cannot open file %s, exit", ofilename);
		return (0);
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL)
			break;
		sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5);
		if (strcmp("--", tmpchar1) == 0
			&& strcmp("Stationary", tmpchar2) == 0
			&& strcmp("point", tmpchar3) == 0
			&& strcmp("found.", tmpchar4) == 0)
			Found_Stationary = 1;
		if (Found_Stationary == 1)
			if (strcmp("ESP", tmpchar1) == 0
				&& strcmp("Fit", tmpchar2) == 0
				&& strcmp("Center", tmpchar3) == 0)
				espindex = 1;
	}
	if (espindex == 1)
		Found_Stationary = 2;
	if (Found_Stationary == 0)
		Found_Stationary = -1;
	tmpint1 = 0;
	tmpint2 = 0;
	tmpint3 = 0;
	rewind(fpin);

	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL)
			break;
		sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5);
		if (strcmp("--", tmpchar1) == 0
			&& strcmp("Stationary", tmpchar2) == 0
			&& strcmp("point", tmpchar3) == 0
			&& strcmp("found.", tmpchar4) == 0)
			Found_Stationary = 1;
		if (Found_Stationary == 2)
			continue;
		if (strcmp("Atomic", tmpchar1) == 0
			&& strcmp("Center", tmpchar2) == 0) {
			sscanf(&line[32], "%lf%lf%lf", &cord[tmpint1][0],
				   &cord[tmpint1][1], &cord[tmpint1][2]);
			tmpint1++;
			if(tmpint1 > MAX_ATOM_CENTER) {
				fprintf(stderr,"Error, the number of atomic center exceeds MAX_ATOM_CENTER defined in espgen.c, extend MAX_ATOM_CENTER and recompile the program, tmpint1, MAX_ATOM_CENTER\n");
				exit(0);
			}
		}
		if (strcmp("ESP", tmpchar1) == 0 && strcmp("Fit", tmpchar2) == 0 &&
			strcmp("Center", tmpchar3) == 0) {
			sscanf(&line[32], "%lf%lf%lf", &esp[tmpint2].x,
				   &esp[tmpint2].y, &esp[tmpint2].z);
			esp[tmpint2].x /= Bohr;
			esp[tmpint2].y /= Bohr;
			esp[tmpint2].z /= Bohr;
			tmpint2++;
			if (tmpint2 >= maxesp) {
				maxesp += MAXESP;
				printf
					("\nInfo: the number of the ESP exceeds the MAXESP(%d),extend the size and reallocate the memory automatically",
					 maxesp);
				esp = (ESP *) realloc(esp, maxesp * sizeof(ESP));
				if (esp == NULL) {
					fprintf(stderr,
							"Info: number of EPS exceeds MAXESP, reallocate memory for *esp automatically\n");
					exit(0);
				}
			}
		}
		if (strncmp("Fit", &line[6], 3) == 0) {
			sscanf(&line[12], "%lf", &esp[tmpint3++].esp);
/*     esp[tmpint3-1].esp/=Bohr; 
       esp[tmpint3-1].esp*=10.0; */
		}
	}
/*   printf("\n %10d%10d%10d", tmpint1,tmpint2,tmpint3);*/
	if (tmpint3 < 10)
		fail = 1;
	fprintf(fpout, "%5d%5d%5d\n", tmpint1, tmpint2, 0);
	for (i = 0; i < tmpint1; i++)
		fprintf(fpout, "%32.7E%16.7E%16.7E\n", cord[i][0] / Bohr,
				cord[i][1] / Bohr, cord[i][2] / Bohr);

	for (j = 0; j < tmpint2; j++)
		fprintf(fpout, "%16.7E%16.7E%16.7E%16.7E\n", esp[j].esp, esp[j].x,
				esp[j].y, esp[j].z);

	if (method == 1) {
		rewind(fpin);
		Found_Stationary = 0;
		for (;;) {
			if (fgets(line, MAXCHAR, fpin) == NULL) {
				if (Found_Stationary == 0) {
					Found_Stationary = 1;
					rewind(fpin);
				} else
					break;
			}
			sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
				   tmpchar4, tmpchar5);
			if (strcmp("--", tmpchar1) == 0
				&& strcmp("Stationary", tmpchar2) == 0
				&& strcmp("point", tmpchar3) == 0
				&& strcmp("found.", tmpchar4) == 0) {
				Found_Stationary = 1;
				continue;
			}
			if (Found_Stationary == 1)
				if (strcmp("Quadrupole", tmpchar1) == 0
					&& strcmp("moment", tmpchar2) == 0
					&& strcmp("(Debye-Ang):", tmpchar3) == 0) {
					index1 = 1;
					continue;
				}
			if (Found_Stationary == 1)
				if (strcmp("Octapole", tmpchar1) == 0
					&& strcmp("moment", tmpchar2) == 0
					&& strcmp("(Debye-Ang**2):", tmpchar3) == 0) {
					index2 = 1;
					continue;
				}
			if (Found_Stationary == 1 && index1 == 1 && index2 == 0
				&& index3 == 0) {
				sscanf(&line[7], "%lf", &quadrupole.xx);
				sscanf(&line[24], "%lf", &quadrupole.yy);
				sscanf(&line[41], "%lf", &quadrupole.zz);
				index3 = 1;
				continue;
			}
			if (Found_Stationary == 1 && index1 == 1 && index2 == 0
				&& index3 == 1) {
				sscanf(&line[8], "%lf", &quadrupole.xy);
				sscanf(&line[25], "%lf", &quadrupole.xz);
				sscanf(&line[43], "%lf", &quadrupole.yz);
				break;
			}
		}
		rewind(fpin);
		Found_Stationary = 0;
		for (;;) {
			if (fgets(line, LINELEN_MAX, fpin) == NULL) {
				if (Found_Stationary == 0) {
					Found_Stationary = 1;
					rewind(fpin);
				} else
					break;
			}
			sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
				   tmpchar4, tmpchar5);
			if (strcmp("--", tmpchar1) == 0
				&& strcmp("Stationary", tmpchar2) == 0
				&& strcmp("point", tmpchar3) == 0
				&& strcmp("found.", tmpchar4) == 0) {
				Found_Stationary = 1;
				continue;
			}
			if (Found_Stationary == 1)
				if (strcmp("Dipole", tmpchar1) == 0
					&& strcmp("moment", tmpchar2) == 0
					&& strcmp("(Debye):", tmpchar3) == 0) {
					index0 = 1;
					continue;
				}
			if (Found_Stationary == 1 && index0 == 1 && line[12] == '.'
				&& line[29] == '.' && line[46] == '.' && line[63] == '.') {
				sscanf(&line[8], "%lf", &dipole.x);
				sscanf(&line[25], "%lf", &dipole.y);
				sscanf(&line[43], "%lf", &dipole.z);
				sscanf(&line[60], "%lf", &dipole.total);
				break;
			}
		}
		fprintf(fpout, "%16.7E%16.7E%16.7E\n", dipole.x, dipole.y,
				dipole.z);
		fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xx,
				quadrupole.yy, quadrupole.zz);
		fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xy,
				quadrupole.xz, quadrupole.yz);
/*   fprintf(fpout,"\n\n"); */
	}
/*   fprintf(fpout,"\n\n");*/
	fclose(fpin);
	fclose(fpout);
	if (fail == 1) {
		strcat(tmpchar, ofilename);
		status = system(tmpchar);
		if(status != 0) {
                	fprintf(stderr, "Error: cannot run \"%s\" in main() of espgen.c properly, exit\n", tmpchar);
                	exit(0);
        	}
	}
	printf("\n");
/*
	free(esp);
*/
	return (0);
}
