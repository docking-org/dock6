/* Gaussian OUT */
int rjout(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		  MOLINFO *minfo)
{
	int i = 0;
	int index =0;
	int overflow_flag = 0;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char line[MAXCHAR];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open %s, exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, (*minfo).resname);
	index = 0;

	for (;;) {
		strcpy(tmpchar1, "");
		strcpy(tmpchar2, "");
		strcpy(tmpchar3, "");
		strcpy(tmpchar4, "");
		tmpchar1[0]='\0';
		tmpchar2[0]='\0';
		tmpchar3[0]='\0';
		tmpchar4[0]='\0';
                if (fgets(line, MAXCHAR, fpin) == NULL) {
/*     printf("\nFinished reading %s file.", cinfo.ifilename); */
                        break;
                }
		sscanf(&line[0], "%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4);
		if(strcmp(tmpchar1, "net") ==0 && strcmp(tmpchar2, "molecular") ==0 && strcmp(tmpchar3, "charge:") ==0) {
		
		(*minfo).dcharge = atof(tmpchar4);
		continue;
	}
		if(strcmp(tmpchar1, "multiplicity:") ==0) {
			(*minfo).dcharge = atof(tmpchar2);
			continue;
		}
		if(strcmp(tmpchar1, "final") ==0 && strcmp(tmpchar2, "geometry:") ==0){
			index = 1;
			continue;
		}	
		if(index == 1 && strcmp(tmpchar1, "atom") == 0) {
			index = 2;
			continue;
		}
		if(index == 2 && strlen(line) <=4) break;
		if (overflow_flag == 0 && index == 2) {
			sscanf(line, "%s%lf%lf%lf", atom[i].name, &atom[i].x, &atom[i].y, &atom[i].z);
			i++;
		}
		if (i >= cinfo.maxatom && overflow_flag == 0) {
			printf
				("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
			overflow_flag = 1;
		}
	}
	*atomnum = i;
	fclose(fpin);
	element(*atomnum, atom);
	atomicnum(*atomnum, atom);
	return overflow_flag;
}


void wjout()
{
	printf
		("\n sorry, you may get the jaguar out by execute the jaguar program");
}
