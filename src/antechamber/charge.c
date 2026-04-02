# define CONVERG 0.00001
# define GASMAXITER 500
# define DAMPFACTOR 0.5
GASTEIGER gas[MAXGAS];

void rmopcharge(char *filename, int atomnum, ATOM atom[], MOLINFO *minfo)
{

	/* now modified for divcon.out output */

	int index;
	int number = 0;
	int tmpint;
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
	index = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6);
		if (strcmp("NET", tmpchar1) == 0 && strcmp("ATOMIC", tmpchar2) == 0
			&& strcmp("CHARGES", tmpchar3) == 0
			&& strcmp("AND", tmpchar4) == 0
			&& strcmp("DIPOLE", tmpchar5) == 0
			&& strcmp("CONTRIBUTIONS", tmpchar6) == 0) {
			index = 1;
			continue;
		}
		if (index == 1 && strcmp("ATOM", tmpchar1) == 0
			&& strcmp("NO.", tmpchar2) == 0) {
			index = 2;
			continue;
		}
		if (index == 2 && strcmp("DIPOLE", tmpchar1) == 0
			&& strcmp("X", tmpchar2) == 0 && strcmp("Y", tmpchar3) == 0
			&& strcmp("Z", tmpchar4) == 0
			&& strcmp("TOTAL", tmpchar5) == 0) {
			index = 3;
			break;
		}
		if (index == 2) {
			sscanf(line, "%d%s%lf%lf", &tmpint, tmpchar, &tmpfloat1,
			   	&tmpfloat2);
		 	atom[number].charge = tmpfloat1;
			number++;
		}
	}
	fclose(fpout);
	if (number == 0) {
		fprintf(stderr, "Error: unable to find mopac charges in %s\n", filename);
		exit(1);
	}
	if ((*minfo).usercharge < -9990)
		(*minfo).icharge = intcharge(atomnum, atom);
}


void rdivcharge(char *filename, int atomnum, ATOM atom[], MOLINFO *minfo, int flag)
{

	/* now modified for divcon.out output */

	int index;
	int number = 0;
	int tmpint;
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
	index = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6);
		/* parse DIVCON output  */
		if (strcmp("NO.", tmpchar1) == 0 && strcmp("SYMBOL", tmpchar2) == 0
			&& strcmp("MULLIKEN", tmpchar3) == 0
			&& strcmp("CM1", tmpchar4) == 0
			&& strcmp("CM2", tmpchar5) == 0 ){
			index = 1;
			continue;
		}
		if (index == 1 && strcmp("CHARGE", tmpchar1) == 0
			&& strcmp("CHARGE", tmpchar2) == 0) {
			index = 2;
			fgets(line, MAXCHAR, fpout);
			continue;
		}
		if (index == 2 && strcmp("\n", line) == 0 ) {
			index = 3;
			break;
		}
		if (index == 2) {
			sscanf(line, "%d%s%lf%lf%lf", &tmpint, tmpchar, &tmpfloat1,
			   	&tmpfloat2, &tmpfloat3);
			if(flag == 1) 
				atom[number].charge = tmpfloat1;
			if(flag == 2) 
				atom[number].charge = tmpfloat2;
			if(flag == 3) 
				atom[number].charge = tmpfloat3;
			number++;
		}
	}
	fclose(fpout);
	if (number == 0) {
		fprintf(stderr, "Error: unable to find divcon charges in %s\n", filename);
		exit(1);
	}
	if ((*minfo).usercharge < -9990)
		(*minfo).icharge = intcharge(atomnum, atom);
}

void rgaucharge(char *filename, char *chargemethod, int atomnum,
				ATOM atom[], MOLINFO *minfo)
{
	int chargeindex;
	int num;
	int Found_Stationary = 0;
	int index = 0;
	char line[MAXCHAR];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open %s, exit\n", filename);
		return;
	}
	if (strcmp(chargemethod, "mul") == 0)
		chargeindex = 1;
	else
		chargeindex = 2;
	num = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL)
			break;
		if (strncmp("-- Stationary point found", &line[4], 25) == 0) {
			Found_Stationary = 1; 
			continue; 
		}
		if (chargeindex == 1 && Found_Stationary == 1 && index == 0) {
			if(strncmp("Total atomic charges:", &line[1], 21) == 0) index = 1;
			if(strncmp("Mulliken atomic charges:", &line[1], 24) == 0) index = 1;
			continue;
		}
		if (chargeindex == 2 && Found_Stationary == 1 && index == 0 && 
			strncmp("Fitting point charges", &line[1], 21) == 0)
			index = 1;
		if (index ==1 && strncmp("Atomic charges with hydrogens summed into heavy atoms", 
			&line[1], 52) == 0)
			index = 2;
		if (Found_Stationary == 1 && index == 1) {
			if(line[11]=='.' ) {
				sscanf(&line[8], "%lf", &atom[num].charge);
				num = num + 1;
			}
			if(line[14]=='.') {
				sscanf(&line[11], "%lf", &atom[num].charge);
				num = num + 1;
			}
		}
		if (num > atomnum)
			break;
	}
	fclose(fpin);
	if ((*minfo).usercharge < -9990)
		(*minfo).icharge = intcharge(atomnum, atom);
}



/*CHARGE METHOD : READ CHARGE */
void rcharge(char *filename, int atomnum, ATOM atom[], CONTROLINFO cinfo,
			 MOLINFO *minfo)
{
	FILE *fpcharge;
	char line[MAXCHAR];
	int i;

	if ((fpcharge = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open charge file %s , exit\n", filename);
		exit(1);
	}
	i = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpcharge) == NULL) {
/*       printf("\nFinished reading file %s", filename); */
			break;
		}
		sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf", &atom[i].charge,
			   &atom[i + 1].charge, &atom[i + 2].charge,
			   &atom[i + 3].charge, &atom[i + 4].charge,
			   &atom[i + 5].charge, &atom[i + 6].charge,
			   &atom[i + 7].charge);
		i = i + 8;
	}
	fclose(fpcharge);
	if ((*minfo).usercharge < -9990.)
		(*minfo).icharge = intcharge(atomnum, atom);
}

/*CHARGE METHOD : WRITE CHARGE */

void wcharge(char *filename, int atomnum, ATOM atom[], CONTROLINFO cinfo,
			 MOLINFO minfo)
{
	FILE *fpcharge;
	int i;
	int num = 0;

	if ((fpcharge = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Cannot open charge file%s , exit\n", filename);
		exit(1);
	}
	for (i = 0; i < atomnum; i++) {
		fprintf(fpcharge, "%10.6lf", atom[i].charge);
		num++;
		if (num == 8) {
			num = 0;
			fprintf(fpcharge, "\n");
		}
	}
	fclose(fpcharge);
}



/*CHARGE METHOD : RESP CHARGE */

void resp(char *filename, int atomnum, ATOM * atom, int bondnum,
		  BOND * bond, CONTROLINFO cinfo, MOLINFO minfo)
{
	/* char tmpchar[200]; */
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char *system_env;
	int status = 0;

	if (strcmp(cinfo.intype, "gout") != 0
		&& strcmp(cinfo.intype, "11") != 0) {
		printf("\n Sorry, RESP charge needs a Gaussian output file (gout)");
		return;
	}
	wac("ANTECHAMBER_RESP.AC", atomnum, atom, bondnum, bond, cinfo, minfo);
	system_env = (char *) getenv("DOCK_HOME");
	if (system_env != NULL) {
		tmpchar1[0] = '\0';
		tmpchar2[0] = '\0';
		tmpchar3[0] = '\0';
		strcpy(tmpchar1, system_env);
		strcat(tmpchar1, "/bin/espgen -o ANTECHAMBER.ESP -i ");
		strcpy(tmpchar2, system_env);
		strcat(tmpchar2, "/bin/respgen");
                strcpy(tmpchar3, system_env);
                strcat(tmpchar3, "/bin/resp");
	} else {
		strcpy(tmpchar1, "espgen -o ANTECHAMBER.ESP -i ");
		strcpy(tmpchar2, "respgen");
		strcpy(tmpchar3, "resp");
	}
	strcat(tmpchar1, filename);
	if (cinfo.intstatus == 2)
		fprintf(stderr, "\nRunning: %s\n\n", tmpchar1);
	status = system(tmpchar1);
        if(status != 0) {
        	fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", tmpchar1);
                exit(0);
        }
	strcpy(tmpchar4, tmpchar2);
	strcat(tmpchar2,
		   " -i ANTECHAMBER_RESP.AC -o ANTECHAMBER_RESP1.IN -f resp1");
	if (cinfo.intstatus == 2)
		fprintf(stderr, "\nRunning: %s\n\n", tmpchar2);
	status = system(tmpchar2);
        if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", tmpchar2);
                exit(0);
        }
	strcat(tmpchar4,
		   " -i ANTECHAMBER_RESP.AC -o ANTECHAMBER_RESP2.IN -f resp2");
	if (cinfo.intstatus == 2)
		fprintf(stderr, "\nRunning: %s\n\n", tmpchar3);
	status = system(tmpchar4);
        if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", tmpchar3);
                exit(0);
        }
	if (cinfo.intstatus == 2)
		fprintf(stderr,
				" \nRunning: resp -O -i ANTECHAMBER_RESP1.IN -o ANTECHAMBER_RESP1.OUT -e ANTECHAMBER.ESP -t qout\n\n");
        strcpy(tmpchar4,tmpchar3);
        strcat(tmpchar4,
                   " -O -i ANTECHAMBER_RESP1.IN -o ANTECHAMBER_RESP1.OUT -e ANTECHAMBER.ESP -t qout");
        status = system(tmpchar4);
        if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", "resp -O -i ANTECHAMBER_RESP1.IN -o ANTECHAMBER_RESP1.OUT -e ANTECHAMBER.ESP -t qout");
                exit(0);
        }
	if (cinfo.intstatus == 2)
		fprintf(stderr,
				"resp -O -i ANTECHAMBER_RESP2.IN -o ANTECHAMBER_RESP2.OUT -e ANTECHAMBER.ESP -q qout -t QOUT\n\n");
        strcpy(tmpchar4,tmpchar3);
        strcat(tmpchar4,
		   " -O -i ANTECHAMBER_RESP2.IN -o ANTECHAMBER_RESP2.OUT -e ANTECHAMBER.ESP -q qout -t QOUT");
        status = system(tmpchar4);
        if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" in resp() of charge.c properly, exit\n", "resp -O -i ANTECHAMBER_RESP2.IN -o ANTECHAMBER_RESP2.OUT -e ANTECHAMBER.ESP -q qout -t QOUT");
                exit(0);
        }
	rcharge("QOUT", atomnum, atom, cinfo, &minfo);
}

/* CHARGE METHOD : BCC-AM1 */
void bccharge(int atomnum, ATOM atom[], int bondnum, BOND bond[],
			  AROM arom[], CONTROLINFO *cinfo, MOLINFO *minfo)
{
	int i;
	int atomnum_tmp;
	int bondnum_tmp;
	int status = 0;
	ATOM *atom_tmp;
	BOND *bond_tmp;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char *system_env;

	wac("ANTECHAMBER_AM1BCC_PRE.AC", atomnum, atom, bondnum, bond, *cinfo,
		*minfo);
	system_env = (char *) getenv("DOCK_HOME");
	if (system_env != NULL) {
		tmpchar1[0] = '\0';
		tmpchar2[0] = '\0';
		strcpy(tmpchar1, system_env);
		strcat(tmpchar1, "/bin/am1bcc");
		strcpy(tmpchar2, system_env);
		strcat(tmpchar2, "/parameters/antechamber/BCCPARM.DAT");
	} else {
		strcpy(tmpchar1, "am1bcc");
		strcpy(tmpchar2, "BCCPARM.DAT");
	}
	strcat(tmpchar1,
		   " -i ANTECHAMBER_AM1BCC_PRE.AC -o ANTECHAMBER_AM1BCC.AC -f ac -p ");
	strcat(tmpchar1, tmpchar2);
	strcat(tmpchar1, " -s ");
	newitoa((*cinfo).intstatus, tmpchar2);
	strcat(tmpchar1, tmpchar2);
/*
	if (strcmp((*cinfo).intype, "ac") == 0 || strcmp((*cinfo).intype, "1") == 0
		|| strcmp((*cinfo).intype, "mol2") == 0
		|| strcmp((*cinfo).intype, "2") == 0) {
		if ((*cinfo).prediction_index == -1)
			(*cinfo).prediction_index = 0;
		newitoa((*cinfo).prediction_index, tmpchar);
		strcat(tmpchar1, " -j ");
		strcat(tmpchar1, tmpchar);
	} else
		strcat(tmpchar1, " -j 4");
*/
	if ((*cinfo).prediction_index == -1)
		(*cinfo).prediction_index = 0;
	newitoa((*cinfo).prediction_index, tmpchar);
	strcat(tmpchar1, " -j ");
	strcat(tmpchar1, tmpchar);

	if ((*cinfo).intstatus == 2)
		fprintf(stderr, "\nRunning: %s\n", tmpchar1);
	status = system(tmpchar1);
        if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" in bccharge() of charge.c properly, exit\n", tmpchar1);
                exit(0);
        }
	(*cinfo).maxatom = atomnum + 10;
	(*cinfo).maxbond = bondnum + 10;
	atom_tmp = (ATOM *) malloc(sizeof(ATOM) * (*cinfo).maxatom);
	if (atom_tmp == NULL) {
		fprintf(stderr, "memory allocation error for *atom_tmp\n");
		exit(0);
	}
	bond_tmp = (BOND *) malloc(sizeof(BOND) * (*cinfo).maxbond);
	if (bond_tmp == NULL) {
		fprintf(stderr, "memory allocation error for *bond_tmp\n");
		exit(0);
	}
	rac("ANTECHAMBER_AM1BCC.AC", &atomnum_tmp, atom_tmp, &bondnum_tmp,
		bond_tmp, cinfo, minfo);

	for (i = 0; i < atomnum; i++)
		atom[i].charge = atom_tmp[i].charge;
	free(atom_tmp);
	free(bond_tmp);

}

/* CHARGE METHOD : BCC-AM1 */
void bcc(char *filename, int atomnum, ATOM atom[], int bondnum,
		 BOND bond[], AROM arom[], CONTROLINFO *cinfo, MOLINFO *minfo)
{
	char tmpchar[MAXCHAR];
	int status = 0;	

	if (strcmp((*cinfo).intype, "mopout") == 0
		|| strcmp((*cinfo).intype, "12") == 0)
		rmopcharge(filename, atomnum, atom, minfo);
	else if (strcmp((*cinfo).intype, "divout") == 0
		|| strcmp((*cinfo).intype, "22") == 0)
		rdivcharge(filename, atomnum, atom, minfo,1);
	else {
		if((*minfo).envtype == 1 || (*minfo).envtype == 4) {
			wmopcrt("mopac.in", atomnum, atom, *minfo);
			strcpy(tmpchar, "$DOCK_HOME");
			strcat(tmpchar, "/bin/mopac.sh");
		}
		if((*minfo).envtype == 2 || (*minfo).envtype == 3) {
			wdivcrt("divcon.in", atomnum, atom, *minfo);
			strcpy(tmpchar, "$AMBERHOME");
			strcat(tmpchar, "/bin/divcon");
		}
		if ((*minfo).envtype == 0) {
			wmopcrt("mopac.in", atomnum, atom, *minfo);
			strcpy(tmpchar, "mopac.sh");
		}
		if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
			fprintf(stderr, "\nRunning: %s\n", tmpchar);
		status = system(tmpchar);
		if(status != 0) {
			fprintf(stderr, "Error: cannot run \"%s\" of bcc() in charge.c properly, exit\n", tmpchar);
			exit(0);
		}
		if((*minfo).envtype == 2 || (*minfo).envtype == 3) 
			rdivcharge("divcon.out", atomnum, atom, minfo, 1);
		if((*minfo).envtype == 0 ||(*minfo).envtype == 1 ||(*minfo).envtype == 4) 
			rmopcharge("mopac.out", atomnum, atom, minfo);
	}
	bccharge(atomnum, atom, bondnum, bond, arom, cinfo, minfo);
}

/* CHARGE METHOD : CM1 */
void cm1(int atomnum, ATOM atom[], CONTROLINFO *cinfo, MOLINFO *minfo)
{
        char tmpchar[200];
	int status = 0;

	if((*minfo).envtype == 1 || (*minfo).envtype == 4) {
        	wmopcrt("mopac.in", atomnum, atom, *minfo);
                strcpy(tmpchar, "$DOCK_HOME");
                strcat(tmpchar, "/bin/mopac.sh");
	}
	if((*minfo).envtype == 2 || (*minfo).envtype == 3) {
        	wdivcrt("divcon.in", atomnum, atom, *minfo);
                strcpy(tmpchar, "$AMBERHOME");
                strcat(tmpchar, "/bin/divcon");
	}
	if((*minfo).envtype == 0) {
        	wmopcrt("mopac.in", atomnum, atom, *minfo);
                strcpy(tmpchar, "mopac.sh");
	}
         if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
         	fprintf(stderr, "\nRunning: %s\n", tmpchar);
        status = system(tmpchar);
	if(status != 0) {
		fprintf(stderr, "Error: cannot run \"%s\" of cm1() in charge.c properly, exit\n", tmpchar);
		exit(0);
	}
	if((*minfo).envtype == 2 || (*minfo).envtype == 3) 
		rdivcharge("divcon.out", atomnum, atom, minfo, 2);
	if((*minfo).envtype == 0 ||(*minfo).envtype == 1 ||(*minfo).envtype == 4) {
		rmopcharge("mopac.out", atomnum, atom, minfo);
		fprintf(stderr, "Warning: if the mopac program cannot generate the CM1 charges, Mulliken charges are generated\n");
	}
}

/* CHARGE METHOD : CM2 */
void cm2(int atomnum, ATOM atom[], CONTROLINFO *cinfo, MOLINFO *minfo)
{
        char tmpchar[200];
	int status = 0;

	if((*minfo).envtype == 1 || (*minfo).envtype == 4) {
        	wmopcrt("mopac.in", atomnum, atom, *minfo);
                strcpy(tmpchar, "$DOCK_HOME");
                strcat(tmpchar, "/bin/mopac.sh");
	}
	if((*minfo).envtype == 2 || (*minfo).envtype == 3) {
        	wdivcrt("divcon.in", atomnum, atom, *minfo);
                strcpy(tmpchar, "$AMBERHOME");
                strcat(tmpchar, "/bin/divcon");
	}
	if((*minfo).envtype == 0) {
        	wmopcrt("mopac.in", atomnum, atom, *minfo);
                strcpy(tmpchar, "mopac.sh");
	}
         if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
         	fprintf(stderr, "\nRunning: %s\n", tmpchar);
        status = system(tmpchar);
	if(status != 0) {
		fprintf(stderr, "Error: cannot run \"%s\" of cm2() in charge.c properly, exit\n", tmpchar);
		exit(0);
	}
	if((*minfo).envtype == 2 || (*minfo).envtype == 3) 
		rdivcharge("divcon.out", atomnum, atom, minfo, 3);
	if((*minfo).envtype == 0 ||(*minfo).envtype == 1 ||(*minfo).envtype == 4) {
		rmopcharge("mopac.out", atomnum, atom, minfo);
		fprintf(stderr, "Warning: if the mopac program cannot generate the CM2 charges, Mulliken charges are generated\n");
	}
}

/* CHARGE METHOD : Mulliken  */
void mul(char *filename, int atomnum, ATOM atom[], CONTROLINFO *cinfo,
		 MOLINFO *minfo)
{
	char tmpchar[200];
	int status = 0;

	if (strcmp((*cinfo).intype, "gout") == 0
		|| strcmp((*cinfo).intype, "11") == 0)
		rgaucharge(filename, "mul", atomnum, atom, minfo);
	else if (strcmp((*cinfo).intype, "mopout") == 0
			 || strcmp((*cinfo).intype, "12") == 0)
		rmopcharge(filename, atomnum, atom, minfo);
	else if (strcmp((*cinfo).intype, "divout") == 0
			 || strcmp((*cinfo).intype, "22") == 0)
		rdivcharge(filename, atomnum, atom, minfo,1);
	else {
               if((*minfo).envtype == 1 || (*minfo).envtype == 4) {
                        wmopcrt("mopac.in", atomnum, atom, *minfo);
                        strcpy(tmpchar, "$DOCK_HOME");
                        strcat(tmpchar, "/bin/mopac.sh");
                }
                if((*minfo).envtype == 2 || (*minfo).envtype == 3) {
                        wdivcrt("divcon.in", atomnum, atom, *minfo);
                        strcpy(tmpchar, "$AMBERHOME");
                        strcat(tmpchar, "/bin/divcon");
                }
                if((*minfo).envtype == 0) {
                        wmopcrt("mopac.in", atomnum, atom, *minfo);
                        strcpy(tmpchar, "mopac.sh");
                }
                if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
                        fprintf(stderr, "\nRunning: %s\n", tmpchar);
                status = system(tmpchar);
		if(status != 0) {
			fprintf(stderr, "Error: cannot run \"%s\" of mul() in charge.c properly, exit\n", tmpchar);
			exit(0);
		}
                if((*minfo).envtype == 2 || (*minfo).envtype == 3) 
                        rdivcharge("divcon.out", atomnum, atom, minfo, 1);
		if((*minfo).envtype == 0 ||(*minfo).envtype == 1 ||(*minfo).envtype == 4) 
                        rmopcharge("mopac.out", atomnum, atom, minfo);
        }
}

/* CHARGE METHOD : ESP(Kollman) */
void esp(char *filename, int atomnum, ATOM atom[], CONTROLINFO cinfo,
		 MOLINFO minfo)
{
	if (strcmp(cinfo.intype, "gout") == 0)
		rgaucharge(filename, "esp", atomnum, atom, &minfo);
	else
		printf("\n Kollman ESP charges IS Only Used for Gaussian Output");
}

/* CHARGE METHOD : Gasteiger charge */
void rgasparm(char *filename, int *gasparmnum, GASTEIGER gas[])
{
	FILE *fp;
	int num = 0;
	char line[MAXCHAR];

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		return;
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp(line, "GASPARM", 7) == 0) {
			gas[num].name[0] = '\0';
			sscanf(&line[8], "%s%lf%lf%lf%lf%lf", gas[num].name, &gas[num].a,
				   &gas[num].b, &gas[num].c, &gas[num].d, &gas[num].charge);
			num++;
		}
	}
	*gasparmnum = num;
	fclose(fp);
}
void assign(int atomnum, ATOM atom[], int gasparmnum, GASTEIGER gas[],
			int gasparmindex[])
{
	int i, j;
	int flag;
	for (i = 0; i < atomnum; i++) {
		flag = 0;
		for (j = 0; j < gasparmnum; j++)
			if (strcmp(gas[j].name, atom[i].ambername) == 0) {
				gasparmindex[i] = j;
				flag = 1;
				break;
			}
		if (flag == 0) {
			fprintf(stderr,
					"\nNo Gasteiger parameter for atom[%d]:%s:%s, exit(0)\n",
					i, atom[i].name, atom[i].ambername);
			exit(0);
		}
	}
}

double rmscal(int atomnum, double gaschargep[], double gaschargea[])
{
	double total = 0.0;
	double rms;
	int i;

	for (i = 0; i < atomnum; i++) {
		total +=
			(gaschargep[i] - gaschargea[i]) * (gaschargep[i] -
											   gaschargea[i]);
		gaschargep[i] = gaschargea[i];
	}
	rms = total / atomnum;
	rms = sqrt(rms);
	return rms;
}

void gasiter(int atomnum, ATOM atom[], int gasparmnum, GASTEIGER gas[])
{
	int i, j;
	int iteration = 0;
	double q;
	double xx;
	double rmsd;
	int *gasparmindex;
	double *x;
	double *gaschargep;
	double *gaschargea;

	gasparmindex = (int *) malloc(sizeof(int) * atomnum);
	if (gasparmindex == NULL) {
		fprintf(stderr,
				"memory allocation error for gasparmindex in gasiter(), exit\n");
		exit(0);
	}
	x = (double *) malloc(sizeof(double) * atomnum);
	if (x == NULL) {
		fprintf(stderr,
				"memory allocation error for x in gasiter(), exit\n");
		exit(0);
	}
	gaschargep = (double *) malloc(sizeof(double) * atomnum);
	if (gaschargep == NULL) {
		fprintf(stderr,
				"memory allocation error for gaschargep in gasiter(), exit\n");
		exit(0);
	}
	gaschargea = (double *) malloc(sizeof(double) * atomnum);
	if (gaschargea == NULL) {
		fprintf(stderr,
				"memory allocation error for gaschargea in gasiter(), exit\n");
		exit(0);
	}
	assign(atomnum, atom, gasparmnum, gas, gasparmindex);
	for (i = 0; i < atomnum; i++) {
		gaschargep[i] = gas[gasparmindex[i]].charge; 
		gaschargea[i] = gas[gasparmindex[i]].charge;
	}
	do {
		for (i = 0; i < atomnum; i++) {
			x[i] =
				gas[gasparmindex[i]].a +
				gas[gasparmindex[i]].b * gaschargep[i];
			x[i] += gas[gasparmindex[i]].c * gaschargep[i] * gaschargep[i];
			if (x[i] == 0.0)
				x[i] = 0.0000000001;
		}
		for (i = 0; i < atomnum; i++)
			for (j = i + 1; j < atomnum; j++)
				if (atom[i].con[0] == j || atom[i].con[1] == j
					|| atom[i].con[2] == j || atom[i].con[3] == j
					|| atom[i].con[4] == j || atom[i].con[5] == j) {
					if (x[i] <= x[j]) {
						xx = gas[gasparmindex[i]].d ;
						q = (x[j] - x[i]) / xx * pow(DAMPFACTOR,
													 iteration + 1);
						gaschargea[i] += q;
						gaschargea[j] -= q;
					}
					if (x[i] > x[j]) {
						xx = gas[gasparmindex[j]].d; 
						q = (x[i] - x[j]) / xx * pow(DAMPFACTOR,
													 iteration + 1);
						gaschargea[i] -= q;
						gaschargea[j] += q;
					}
				}
		iteration++;
/*
		printf("\nIteration %5d", iteration);
		for (i = 0; i < atomnum; i++) 
			printf("\n%5d %5s %8.4lf %8.4lf", i+1, atom[i].name, gaschargep[i], gaschargea[i]);
*/
		rmsd = rmscal(atomnum, gaschargep, gaschargea);
	} while (rmsd > CONVERG && iteration < GASMAXITER);

	for (i = 0; i < atomnum; i++)
		atom[i].charge = gaschargea[i];
	free(gasparmindex);
	free(x);
	free(gaschargep);
	free(gaschargea);
}

void gascharge(int atomnum, ATOM atom[], int bondnum, BOND bond[],
			   CONTROLINFO cinfo, MOLINFO *minfo)
{
	int i;
	int atomnum_tmp;
	int bondnum_tmp;
	int status = 0;
	ATOM *atom_tmp;
	BOND *bond_tmp;

	int gasparmnum = 0;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char *system_env;
	GASTEIGER *gas;

	gas = (GASTEIGER *) malloc(sizeof(GASTEIGER) * MAXGAS);
	if (gas == NULL) {
		fprintf(stderr,
				"memory allocation error for gas in gascharge(),increase GASMAX in define.h exit\n");
		exit(0);
	}
	system_env = (char *) getenv("DOCK_HOME");
	if (system_env != NULL) {
		(*minfo).gfilename[0] = '\0';
		strcpy((*minfo).gfilename, system_env);
		strcat((*minfo).gfilename, "/parameters/antechamber/GASPARM.DAT");
	} else {
		/*    printf("DOCK_HOME environment not set, using GASPARM.DAT in the current directory\n"); */
		strcpy((*minfo).gfilename, "GASPARM.DAT");
	}
	rgasparm((*minfo).gfilename, &gasparmnum, gas);	/*principle parameter file */
	wac("ANTECHAMBER_GAS.AC", atomnum, atom, bondnum, bond, cinfo, (*minfo));
	system_env = (char *) getenv("DOCK_HOME");
	if (system_env != NULL) {
		tmpchar1[0] = '\0';
		tmpchar2[0] = '\0';
		strcpy(tmpchar1, system_env);
		strcpy(tmpchar2, system_env);
		strcat(tmpchar1, "/bin/atomtype");
		strcat(tmpchar2, "/parameters/antechamber/ATOMTYPE_GAS.DEF");
	} else {
		strcpy(tmpchar1, "atomtype");
		strcpy(tmpchar2, "ATOMTYPE_GAS.DEF");
	}
	strcat(tmpchar1,
		   " -i ANTECHAMBER_GAS.AC -o ANTECHAMBER_GAS_AT.AC -d ");
	strcat(tmpchar1, tmpchar2);
	if (cinfo.intstatus == 2)
		fprintf(stderr, "Running: %s\n", tmpchar1);
	status = system(tmpchar1);
	if(status != 0) {
        	fprintf(stderr, "Error: cannot run \"%s\" of gascharge() in charge.c properly, exit\n", tmpchar1);
                exit(0);
	}
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

	rac("ANTECHAMBER_GAS_AT.AC", &atomnum_tmp, atom_tmp, &bondnum_tmp,
		bond_tmp, &cinfo, minfo);
	gasiter(atomnum_tmp, atom_tmp, gasparmnum, gas);
	for (i = 0; i < atomnum; i++)
		atom[i].charge = atom_tmp[i].charge;
	free(atom_tmp);
	free(bond_tmp);
	free(gas);
}



void write_sybyl_bat(char *str)
{
	FILE *fpout;
	char tmpchar[MAXCHAR];
	char *system_env;
	int status = 0;

	if ((fpout = fopen("antechamber_sybyl.bat", "w")) == NULL) {
		fprintf(stderr, "Cannot open antechamber_sybyl.bat , exit\n");
		exit(1);
	}
	system_env = (char *) getenv("DOCK_HOME");
	if (system_env != NULL) {
		tmpchar[0] = '\0';
		strcpy(tmpchar, system_env);
		strcat(tmpchar, "/parameters/antechamber/charge.spl");
	} else
		strcpy(tmpchar, "charge.spl");

	fprintf(fpout, "%s", "#!/bin/csh");
	fprintf(fpout, "\n%s", "sybyl << @");
	fprintf(fpout, "\n%s", "take ");
	fprintf(fpout, "%s", tmpchar);
	if (strcmp(str, "gas1") == 0)
		fprintf(fpout, "\n%s", "CHARGE GASTEIGER");
	if (strcmp(str, "gas2") == 0)
		fprintf(fpout, "\n%s", "CHARGE GAST_HUCK");
	if (strcmp(str, "del") == 0)
		fprintf(fpout, "\n%s", "CHARGE DELRE");
	if (strcmp(str, "pull") == 0)
		fprintf(fpout, "\n%s", "CHARGE PULLMAN ");
	if (strcmp(str, "huc") == 0)
		fprintf(fpout, "\n%s", "CHARGE HUCKEL");
	if (strcmp(str, "mmff") == 0)
		fprintf(fpout, "\n%s", "CHARGE MMFF94");
	fprintf(fpout, "\n%s", "QUIT YES");
	fprintf(fpout, "\n%s\n\n", "@");
	fclose(fpout);
	status = system("chmod +x antechamber_sybyl.bat");
	if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" of write_sybyl_bat() in charge.c properly, exit\n", "chmod +x antechamber_sybyl.bat");
                exit(0);
        }
	status = system("antechamber_sybyl.bat");
	if(status != 0) {
                fprintf(stderr, "Error: cannot run \"%s\" of write_sybyl_bat() in charge.c properly, exit\n", "antechamber_sybyl.bat");
                exit(0);
        }

}


/* CHARGE METHOD : Gasteiger, using sybyl */
void gas1(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("gas1");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD : Del Re, using sybyl*/
void del(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		 CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("del");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD : Pullman, using sybyl */
void pull(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("pull");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD :Gasteiger-Huckel, using sybyl*/
void gas2(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("gas2");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD : Huckel, using sybyl */
void huc(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		 CONTROLINFO cinfo, MOLINFO minfo)
{

	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("huc");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}

/* CHARGE METHOD :MMFF94, using sybyl*/
void mmff(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom,
		  cinfo, minfo);
	write_sybyl_bat("mmff");
	rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo,
		  &minfo, 0);
}
