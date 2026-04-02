
/* CHARMM */
void wcharmm(char *filename, char *ifilename, int atomnum, ATOM * atom,
                  int bondnum, BOND * bond, CONTROLINFO cinfo, MOLINFO minfo)
{
	FILE *fpin;
	FILE *fpout;
	char tmpchar[200];
	int chargeindex;
	int typeindex;
	int changelineindex;
	int num;
	int tmpint1, tmpint2, tmpint3, tmpint4;
	int i;
	char tmpchar1[10], tmpchar2[10], tmpchar3[10];
	double charge;
	double tmpfloat1, tmpfloat2, tmpfloat3;
	char *system_env;
/*        char *system_env0; */
	char line[MAXCHAR];

	wac("ANTECHAMBER_PREP.AC0", atomnum, atom, bondnum, bond, cinfo,
		minfo);
	system("cp -rf ANTECHAMBER_PREP.AC0 ANTECHAMBER_PREP.AC");

/*part1: if intype is not prepi, prepc and ac, judge atom type*/
	if ((strcmp(cinfo.intype, "prepi") != 0 &&
		 strcmp(cinfo.intype, "prepc") != 0 &&
		 strcmp(cinfo.intype, "5") != 0 &&
		 strcmp(cinfo.intype, "6") != 0 &&
		 strcmp(cinfo.intype, "ac") != 0 &&
		 strcmp(cinfo.intype, "1") != 0) || cinfo.prediction_index == 1
		|| cinfo.prediction_index == 3) {
	        system_env = (char *) getenv("DOCK_HOME");
                if (system_env == NULL)
		  system_env = (char *) getenv("AMBERHOME");
		if (system_env != NULL) {
			tmpchar[0] = '\0';
			strcpy(tmpchar, system_env);
			strcat(tmpchar, "/bin/atomtype");
		} else
			strcpy(tmpchar, "atomtype");
		strcat(tmpchar,
			   " -i ANTECHAMBER_PREP.AC0 -o ANTECHAMBER_PREP.AC -p ");
		strcat(tmpchar, minfo.atom_type_def);
                if(cinfo.intstatus == 2)
			fprintf(stderr, "\nRunning: %s\n", tmpchar);
		system(tmpchar);
	}

	system_env = (char *) getenv("DOCK_HOME");
        if (system_env == NULL)
	  system_env = (char *) getenv("AMBERHOME");
	if (system_env != NULL) {
		tmpchar[0] = '\0';
		strcpy(tmpchar, system_env);
		strcat(tmpchar, "/bin/charmmgen");
	} else
		strcpy(tmpchar, "charmmgen");
	strcat(tmpchar, " -i ANTECHAMBER_PREP.AC -f ac -o ");
	strcat(tmpchar, filename);
	strcat(tmpchar, " -rn \"");
	strcat(tmpchar, atom[0].aa);
	strcat(tmpchar, "\" -rf ");
	strcat(tmpchar, minfo.resfilename);
	if(cinfo.intstatus == 2)
       	     	fprintf(stderr, "\nRunning: %s\n", tmpchar);
       system(tmpchar);
/*  system("rm -f ANTECHAMBER_PREP.AC");*/
}

