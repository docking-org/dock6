/* Gaussian OUT */
int rgout(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		  MOLINFO *minfo)
{
	int i, j;
	int Found_Stationary = 0;
	int Standard = 0;
	int index;
	int overflow_flag = 0;
	int rwindex1 = 0;
	/* int rwindex2 = 0; */
	int rwindex3 = 0;
	/* int inputindex = 0; */
	int readinputindex = 0;
/* int index1, index2;*/
	long lindex;
	long lineindex;
	int chargeindex = 1;
/* double cord[MAXATOM][3];*/
/* char tmpchar[100];*/
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char line[MAXCHAR];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open %s, exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, (*minfo).resname);
	index = 0;
/*-- Stationary point found. */

	i = 0;
	lindex = 0;
	for (;;) {
		lindex++;
		sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5);
		if (chargeindex == 1 && strcmp(tmpchar1, "Charge") == 0) {
			sscanf(&line[9], "%lf", &(*minfo).dcharge);
			sscanf(&line[27], "%d", &(*minfo).multiplicity);
			(*minfo).icharge = (int) (*minfo).dcharge;
			chargeindex = 0;
		}
		if (fgets(line, MAXCHAR, fpin) == NULL) {
			if (Found_Stationary == 0) {
				Found_Stationary = 1;
				readinputindex = 1;
				lindex = 0;
				rwindex1 = 1;
				rewind(fpin);
			} else if (rwindex3 == 0) {
				rwindex3 = 1;
				rewind(fpin);
			} else if (rwindex1 == 1 && index == 0) {
				/*  rwindex2 = 1;  */
				readinputindex = 2;
				lindex = 0;
				rewind(fpin);
			} else
				break;
		}
		if (readinputindex == 0 && strcmp("Input", tmpchar1) == 0
			&& strcmp("orientation:", tmpchar2) == 0)
			lineindex = lindex;
		if (strcmp("--", tmpchar1) == 0
			&& strcmp("Stationary", tmpchar2) == 0
			&& strcmp("point", tmpchar3) == 0
			&& strcmp("found.", tmpchar4) == 0)
			Found_Stationary = 1;
		if (Found_Stationary == 1 && strcmp("Standard", tmpchar1) == 0
			&& strcmp("orientation:", tmpchar2) == 0)
			Standard = 1;
		if (Found_Stationary == 1 && Standard == 1
			&& strncmp(".", &line[28], 1) == 0
			&& strncmp(".", &line[40], 1) == 0
			&& strncmp(".", &line[52], 1) == 0) {
			/* if(i==0) printf("\n readin from standard coordinates"); */
			if (overflow_flag == 0) {
				sscanf(&line[10], "%d%lf%lf%lf", &atom[i].atomicnum,
					   &atom[i].x, &atom[i].y, &atom[i].z);
				atom[i].charge = 0.0;
			}
			index = 1;
			i++;
			if (i >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
		if (Found_Stationary == 1 && Standard == 1
			&& strncmp(".", &line[39], 1) == 0
			&& strncmp(".", &line[51], 1) == 0
			&& strncmp(".", &line[63], 1) == 0) {
			/* if(i==0) printf("\n readin from standard coordinates"); */
			if (overflow_flag == 0) {
				sscanf(&line[10], "%d", &atom[i].atomicnum);
				sscanf(&line[31], "%lf%lf%lf", &atom[i].x, &atom[i].y,
					   &atom[i].z);
				atom[i].charge = 0.0;
			}
			index = 1;
			i++;
			if (i >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
		if (readinputindex == 2 && (lindex - lineindex) > -3) {
			if (strncmp(".", &line[40], 1) == 0
				&& strncmp(".", &line[52], 1) == 0) {
				/*   if(i==0) printf("\n readin from input coordinates");  */
				if (overflow_flag == 0) {
					sscanf(&line[10], "%d%lf%lf%lf", &atom[i].atomicnum,
						   &atom[i].x, &atom[i].y, &atom[i].z);
					atom[i].charge = 0.0;
				}
				index = 1;
				i++;
				if (i >= cinfo.maxatom && overflow_flag == 0) {
					printf
						("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
					overflow_flag = 1;
				}
			}
			if (strncmp(".", &line[39], 1) == 0
				&& strncmp(".", &line[51], 1) == 0
				&& strncmp(".", &line[63], 1) == 0) {
				/*   if(i==0) printf("\n readin from input coordinates"); */
				if (overflow_flag == 0) {
					sscanf(&line[10], "%d", &atom[i].atomicnum);
					sscanf(&line[31], "%lf%lf%lf", &atom[i].x, &atom[i].y,
						   &atom[i].z);
					atom[i].charge = 0.0;
				}
				index = 1;
				i++;
				if (i >= cinfo.maxatom && overflow_flag == 0) {
					printf
						("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
					overflow_flag = 1;
				}
			}
		}
		if (index == 1 && strncmp(".", &line[40], 1) != 0
			&& strncmp(".", &line[52], 1) != 0
			&& strncmp(".", &line[51], 1) != 0
			&& strncmp(".", &line[62], 1) != 0)
			break;

	}
	*atomnum = i;
	fclose(fpin);
	element(*atomnum, atom);
	for (j = 0; j < *atomnum; j++)
		strcpy(atom[j].name, atom[j].element);
	return overflow_flag;
}


void wgout()
{
	printf
		("\n sorry, you may get the gaussian out by execute the Gaussian program");
}
