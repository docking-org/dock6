/* GZMAT */
int rgzmat(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		   MOLINFO minfo)
{
        typedef struct {
                char name[30];
        } STRNAME;

	FILE *fpin;
	int i, j, index, index0;
	int overflow_flag = 0;
	int findindex;
	int numatom;
	STRNAME *bondstr;	
	STRNAME *anglestr;	
	STRNAME *twiststr;	
	char tmpchar1[10];
	char tmpchar2[10];
	char line[MAXCHAR];


	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
        bondstr = (STRNAME *) malloc(sizeof(STRNAME) * (cinfo.maxatom +10));
        if (bondstr == NULL) {
                fprintf(stderr, "memory allocation error for *bondstr in rgzmat()\n");
                exit(0);
        }
        anglestr = (STRNAME *) malloc(sizeof(STRNAME) * (cinfo.maxatom +10));
        if (anglestr == NULL) {
                fprintf(stderr, "memory allocation error for *anglestr in rgzmat()\n");
                exit(0);
        }
        twiststr = (STRNAME *) malloc(sizeof(STRNAME) * (cinfo.maxatom +10));
        if (twiststr == NULL) {
                fprintf(stderr, "memory allocation error for *twiststr in rgzmat()\n");
                exit(0);
        }

	initial(cinfo.maxatom, atom, minfo.resname);
	index = 0;
	index0 = 1;
	numatom = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s file.", cinfo.ifilename); */
			break;
		}
		if (spaceline(line) == 1) {
			index++;
			continue;
		}
		if (index >= 2)
			index++;
		if (index <= 3)
			continue;
		if (index >= 4) {
			if (spaceline(line) == 1 || strncmp(line, "Vari", 4) == 0
				|| strncmp(line, "vari", 4) == 0)
				index0 = -1;
			if (strncmp(line, "Const", 5) == 0
				|| strncmp(line, "const", 5) == 0)
				index0 = -1;
		}
		if (index == 4) {
			if (overflow_flag == 0)
				sscanf(line, "%s", atom[numatom].name);
			numatom++;
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
			continue;
		}
		if (index == 5) {
			if (overflow_flag == 0)
				sscanf(line, "%s%d%s", atom[numatom].name,
					   &atom[numatom].bondatom, bondstr[numatom].name);
			numatom++;
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
			continue;
		}
		if (index == 6) {
			if (overflow_flag == 0)
				sscanf(line, "%s%d%s%d%s", atom[numatom].name,
					   &atom[numatom].bondatom, bondstr[numatom].name,
					   &atom[numatom].angleatom, anglestr[numatom].name);
			numatom++;
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
			continue;
		}
		if (index0 != -1) {
			if (overflow_flag == 0)
				sscanf(line, "%s%d%s%d%s%d%s", atom[numatom].name,
					   &atom[numatom].bondatom, bondstr[numatom].name,
					   &atom[numatom].angleatom, anglestr[numatom].name,
					   &atom[numatom].twistatom, twiststr[numatom].name);
			numatom++;
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}

			continue;
		}
		if (index0 == -1) {
			sscanf(line, "%s%s", tmpchar1, tmpchar2);
			for (i = 1; i < numatom; i++) {
				findindex = 1;
				for (j = 0; j < strlen(bondstr[i].name); j++)
					if (bondstr[i].name[j] != tmpchar1[j]) {
						findindex = 0;
						break;
					}
				if (findindex == 1) {
					strcpy(bondstr[i].name, tmpchar2);
					break;
				}

			}
			for (i = 2; i < numatom; i++) {
				findindex = 1;
				for (j = 0; j < strlen(anglestr[i].name); j++)
					if (anglestr[i].name[j] != tmpchar1[j]) {
						findindex = 0;
						break;
					}
				if (findindex == 1) {
					strcpy(anglestr[i].name, tmpchar2);
					break;
				}
			}
			for (i = 3; i < numatom; i++) {
				findindex = 1;
				for (j = 0; j < strlen(twiststr[i].name); j++)
					if (twiststr[i].name[j] != tmpchar1[j]) {
						findindex = 0;
						break;
					}
				if (findindex == 1) {
					strcpy(twiststr[i].name, tmpchar2);
					break;
				}
			}
		}
	}
	atom[1].bondatom--;
	atom[2].bondatom--;
	atom[2].angleatom--;
	for (i = 3; i < numatom; i++) {
		atom[i].bondatom--;
		atom[i].angleatom--;
		atom[i].twistatom--;
	}
	for (i = 1; i < numatom; i++)
		atom[i].bond = atof(bondstr[i].name);
	for (i = 2; i < numatom; i++)
		atom[i].angle = atof(anglestr[i].name);
	for (i = 3; i < numatom; i++)
		atom[i].twist = atof(twiststr[i].name);
	*atomnum = numatom;
/* printf("\n atom number is  %5d", *atomnum); */
	fclose(fpin);
	free(bondstr);
	free(anglestr);
	free(twiststr);
	return overflow_flag;

}
void wgzmat(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
	FILE *fpout;
	int i;
	/* int index; */
	char tmpchar0[10];
	char tmpchar1[10];
	char tmpchar2[10];
	char tmpchar3[10];

	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
	intercoord(atomnum, atom);
        fprintf(fpout, "%s\n", "--Link1--");
        fprintf(fpout, "%s%s\n", "%chk=", minfo.chkfile);
	fprintf(fpout, "%s\n\n", minfo.gkeyword);
	fprintf(fpout, "%s\n\n", "remark line goes here");
	fprintf(fpout, "%d%4d\n", minfo.icharge, minfo.multiplicity);
	element(atomnum, atom);
	for (i = 0; i < atomnum; i++) {
		newitoa(i + 1, tmpchar0);
		if (i == 0) {
			fprintf(fpout, "%5s\n", atom[i].element);
			continue;
		}
		if (i == 1) {
			strcpy(tmpchar1, "b");
			strcat(tmpchar1, tmpchar0);
			fprintf(fpout, "%5s%5d%8s\n", atom[i].element,
					atom[i].bondatom + 1, tmpchar1);
			continue;
		}
		if (i == 2) {
			strcpy(tmpchar1, "b");
			strcat(tmpchar1, tmpchar0);
			fprintf(fpout, "%5s%5d%8s", atom[i].element,
					atom[i].bondatom + 1, tmpchar1);
			strcpy(tmpchar2, "a");
			strcat(tmpchar2, tmpchar0);
			fprintf(fpout, "%5d%8s\n", atom[i].angleatom + 1, tmpchar2);
			continue;
		}
		strcpy(tmpchar1, "b");
		strcat(tmpchar1, tmpchar0);
		fprintf(fpout, "%5s%5d%8s", atom[i].element, atom[i].bondatom + 1,
				tmpchar1);
		strcpy(tmpchar2, "a");
		strcat(tmpchar2, tmpchar0);
		fprintf(fpout, "%5d%8s", atom[i].angleatom + 1, tmpchar2);
		strcpy(tmpchar3, "t");
		strcat(tmpchar3, tmpchar0);
		fprintf(fpout, "%5d%8s\n", atom[i].twistatom + 1, tmpchar3);
	}

	fprintf(fpout, "Variables:\n");
	fprintf(fpout, "b2= %8.4lf\n", atom[1].bond);
	fprintf(fpout, "b3= %8.4lf\n", atom[2].bond);
	fprintf(fpout, "a3= %8.4lf\n", atom[2].angle);
	for (i = 3; i < atomnum; i++) {
		newitoa(i + 1, tmpchar0);
		strcpy(tmpchar1, "b");
		strcat(tmpchar1, tmpchar0);
		strcpy(tmpchar2, "a");
		strcat(tmpchar2, tmpchar0);
		strcpy(tmpchar3, "t");
		strcat(tmpchar3, tmpchar0);
		fprintf(fpout, "%s= %8.4lf\n", tmpchar1, atom[i].bond);
		fprintf(fpout, "%s= %8.4lf\n", tmpchar2, atom[i].angle);
		fprintf(fpout, "%s= %8.4lf\n", tmpchar3, atom[i].twist);
	}
	fprintf(fpout, "\n");
	fclose(fpout);
}
