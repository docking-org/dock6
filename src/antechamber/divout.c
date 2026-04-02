/* MOPAC OUTPUT */
int rdivout(char *filename, int *atomnum, ATOM * atom, CONTROLINFO *cinfo,
			MOLINFO *minfo)
{
	/* int i,j; */
	int flag;
	int number = 0;
	int tmpint;
	int overflow_flag = 0;
	FILE *fpout;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char tmpchar6[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1;
	double tmpfloat2;
	double tmpfloat3;

	if ((fpout = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open %s , exit\n", filename);
		exit(1);
	}
	initial((*cinfo).maxatom, atom, (*minfo).resname);
	flag = 0;

	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6);
		if (strcmp(tmpchar1, "TOTAL") == 0 
			&& strcmp(tmpchar2, "MULLIKEN") == 0 
			&& strcmp(tmpchar3, "CHARGE") == 0 
			&&  strcmp(tmpchar4, "=") == 0) {
			sscanf(line, "%s%s%s%s%lf", tmpchar1, tmpchar2, tmpchar3,
				   tmpchar4, &(*minfo).dcharge);
			(*minfo).icharge = (int) (*minfo).dcharge;
			continue;
		}
		if (strcmp(tmpchar1, "FINAL") == 0 && strcmp(tmpchar2, "QUANTITIES:") == 0) {
			flag = 1;	
			continue;
		}
		if (flag == 1 && strcmp(tmpchar1, "NUMBER") == 0
			&& strcmp("SYMBOL", tmpchar2) == 0
			&& strcmp("X", tmpchar3) == 0
			&& strcmp("Y", tmpchar4) == 0
			&& strcmp("Z", tmpchar5) == 0) {
			flag = 2;
			continue;
		}
		if(flag == 2) {
			flag = 3;
			continue;
		}
		if (flag == 3 && spaceline(line) == 1) {
			flag = 4;
			*atomnum = number;
			break;
		}
		if (flag == 3) {
			if (overflow_flag == 0) {
				sscanf(line, "%d%s%lf%lf%lf", &tmpint, tmpchar, &tmpfloat1,
					   &tmpfloat2, &tmpfloat3);
				strcpy(atom[number].name, tmpchar);
				atom[number].x = tmpfloat1;
				atom[number].y = tmpfloat2;
				atom[number].z = tmpfloat3;
				strcpy(atom[number].aa, (*minfo).resname);
			}
			number++;
			if (number >= (*cinfo).maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
	}
	fclose(fpout);
	return overflow_flag;
}

void wdivout(char *filename, int atomnum, ATOM * atom, CONTROLINFO cinfo,
			 MOLINFO minfo)
{
	char *system_env;
	char tmpchar[MAXCHAR];
	int status1 = 0;
	int status2 = 0;
	strcpy(minfo.mkeyword, "CARTESIAN AM1 STANDARD DIRECT OPT=BFGS XTEST=0.0001");	
	wdivcrt("divcon.in", atomnum, atom, minfo);
        system_env = (char *) getenv("AMBERHOME");
	if (system_env != NULL) {
		strcpy(tmpchar, system_env);
		strcat(tmpchar, "/exe/");
		strcat(tmpchar, "divcon");
	} else
		strcpy(tmpchar, "divcon");
	if (cinfo.intstatus == 1 || cinfo.intstatus == 2)
		fprintf(stderr, "\nRunning: %s\n", tmpchar);
	status1 = system(tmpchar);
	if(status1 == 0) {
		strcpy(tmpchar, "mv -f divcon.out ");
		strcat(tmpchar, filename);
		status2 = system(tmpchar);
		if(status2 != 0) {
			fprintf(stderr, "Error: cannot run \"%s\" in wdivout() of divout.c properly, exit\n", tmpchar);
			exit(0);
		}
	}
	if(status1 != 0) {
		fprintf(stderr, "Error: cannot run \"%s\" in wdivout() of divout.c properly, exit\n", tmpchar);
		exit(0);
	}
}
