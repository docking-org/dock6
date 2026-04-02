void default_minfo(MOLINFO * minfo)
{
        char *system_env;
	(*minfo).envtype = 0;
        system_env = (char *) getenv("DOCK_HOME");
        if (system_env != NULL) 
		(*minfo).envtype = 1;

/* In default, divcon is used if $AMBERHOME is set*/
        system_env = (char *) getenv("AMBERHOME");
        if (system_env != NULL) {
		if ((*minfo).envtype == 1) 
			(*minfo).envtype = 3;
		else
			(*minfo).envtype = 2;
	}
	if ((*minfo).envtype == 0) 
		fprintf(stderr, "Warning: $AMBERHOME and $DOCK_HOME enviornment strings are not set, use \"mopac.sh\" in the work directory\n");
	if((*minfo).envtype == 2 || (*minfo).envtype == 3)
        	strcpy((*minfo).mkeyword, "CARTESIAN AM1 STANDARD DIRECT OPT=BFGS XTEST=0.0001");
	else
        	strcpy((*minfo).mkeyword, "AM1 MMOK GEO-OK PRECISE");
	strcpy((*minfo).gkeyword,
		   "#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt");
	strcpy((*minfo).resname, "MOL");
	strcpy((*minfo).atom_type_def, "gaff");
	strcpy((*minfo).resfilename, "molecule.res");
	strcpy((*minfo).chkfile, "molecule");
	strcpy((*minfo).gfilename, "GASPARM.DAT");
	strcpy((*minfo).connect_file, "CONNECT.TPL");
	strcpy((*minfo).radius_file, "RADIUS.DAT");
	strcpy((*minfo).inf_filename, "ATOMTYPE.INF");
	(*minfo).longresname[0] = '\0';
	(*minfo).multiplicity = 1;
	(*minfo).icharge = 0;
	(*minfo).usercharge = -9999;
	(*minfo).dcharge = -9999.0;

}

void default_cinfo(CONTROLINFO * cinfo)
{
	(*cinfo).intype[0] = '\0';
	(*cinfo).outtype[0] = '\0';
	(*cinfo).atype[0] = '\0';
	(*cinfo).chargetype[0] = '\0';
	(*cinfo).rnindex = 0;
	(*cinfo).intstatus = 1;
	(*cinfo).pfindex = 0;
	(*cinfo).prediction_index = 4;
	(*cinfo).bpindex = 1;
	(*cinfo).maxatom = MAXATOM;
	(*cinfo).maxbond = MAXBOND;
	(*cinfo).maxring = MAXRING;
}



void default_inf(int atomnum, ATOM atom[], int index)
{
	int i;
	for (i = 0; i < atomnum; i++) {
		strcpy(atom[i].chain, " ");
		atom[i].ter = -1;
		if (index == 1)
			strcpy(atom[i].ambername, atom[i].name);
	}
}



int intcharge(int atomnum, ATOM atom[])
{
	int i;
	int icharge = 0;
	double fraction = 0.0;		/*decimal fraction */
	double dcharge = 0.0;		/*double precision charge */
	double tmpf;


	for (i = 0; i < atomnum; i++)
		dcharge += atom[i].charge;
	fraction = modf(dcharge, &tmpf);
	icharge = (int) tmpf;
	if (fabs(fraction) >= 0.50) {
		if (dcharge < 0)
			icharge--;
		if (dcharge > 0)
			icharge++;
	}
	return icharge;
}

void formula(int atomnum, ATOM atom[], char *form)
{
	int i, j;
	int countatom[150];
	char tmpchar[10];
	char elemname[150][5];

	for (i = 0; i < 150; i++) {
		countatom[i] = 0;
		strcpy(elemname[i], "X");
	}
        strcpy(elemname[1], "H");
        strcpy(elemname[3], "Li");
        strcpy(elemname[4], "Be");
        strcpy(elemname[5], "B");
        strcpy(elemname[6], "C");
        strcpy(elemname[7], "N");
        strcpy(elemname[8], "O");
        strcpy(elemname[9], "F");
        strcpy(elemname[11], "Na");
        strcpy(elemname[12], "Mg");
        strcpy(elemname[13], "Al");
        strcpy(elemname[14], "Si");
        strcpy(elemname[15], "P");
        strcpy(elemname[16], "S");
        strcpy(elemname[17], "Cl");
        strcpy(elemname[19], "K");
        strcpy(elemname[20], "Ca");
        strcpy(elemname[21], "Sc");
        strcpy(elemname[22], "Ti");
        strcpy(elemname[23], "V");
        strcpy(elemname[24], "Cr");
        strcpy(elemname[25], "Mn");
        strcpy(elemname[26], "Fe");
        strcpy(elemname[27], "Co");
        strcpy(elemname[28], "Ni");
        strcpy(elemname[29], "Cu");
        strcpy(elemname[30], "Zn");
        strcpy(elemname[31], "Ga");
        strcpy(elemname[32], "Ge");
        strcpy(elemname[33], "As");
        strcpy(elemname[34], "Se");
        strcpy(elemname[35], "Br");
        strcpy(elemname[38], "Sr");
        strcpy(elemname[44], "Ru");
        strcpy(elemname[45], "Rh");
        strcpy(elemname[46], "Pd");
        strcpy(elemname[47], "Ag");
        strcpy(elemname[48], "Cd");
        strcpy(elemname[53], "I");
        strcpy(elemname[56], "Ba");
        strcpy(elemname[78], "Pt");
        strcpy(elemname[79], "Au");
        strcpy(elemname[80], "Hg");
        strcpy(elemname[80], "Tl");
        strcpy(elemname[82], "Pb");
	strcpy(form, "");
	for (j = 0; j < atomnum; j++)
		countatom[atom[j].atomicnum]++;
	for (i = 0; i < 150; i++)
		if (countatom[i] >= 1) {
			strcat(form, elemname[i]);
			newitoa(countatom[i], tmpchar);
			strcat(form, tmpchar);
			strcat(form, " ");
		}
}

void initial(int num, ATOM * atom, char *resname)
{
	int i;
	for (i = 0; i < num; i++) {
		atom[i].connum = 0;
		atom[i].resno = 1;
		atom[i].con[0] = -1;
		atom[i].con[1] = -1;
		atom[i].con[2] = -1;
		atom[i].con[3] = -1;
		atom[i].con[4] = -1;
		atom[i].con[5] = -1;
		atom[i].charge = 0.0;
		atom[i].improper = 0;
		strcpy(atom[i].aa, resname);
	}
}


void read_radius(char *radius_parm, int atomnum, ATOM atom[])
{
	typedef struct {
		char name[5];
		double radius;
	} RADII;
	RADII radii[200];
	int i, j, k;
	int index;
	int num = 0;
	int parmnum;
	double tmpfloat;
	char line[MAXCHAR];
	char tmpchar[MAXCHAR];
	FILE *fpin;


	if ((fpin = fopen(radius_parm, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", radius_parm);
		return;
	}
	for (;;) {
		if (fgets(line, LINELEN_MAX, fpin) == NULL)
			break;
		if (strncmp(line, "RADIUS", 6) == 0) {
			sscanf(&line[7], "%s%lf", tmpchar, &tmpfloat);
			strcpy(radii[num].name, tmpchar);
			radii[num++].radius = tmpfloat;
		}
	}
	fclose(fpin);
	parmnum = num;

	for (i = 0; i < atomnum; i++) {
		index = 0;
		for (j = 0; j < parmnum; j++)
			if (strcmp(radii[j].name, atom[i].ambername) == 0) {
				atom[i].radius = radii[j].radius;
				index = 1;
				break;
			}
		if (index == 0) {
			tmpchar[0] = '\0';
			num = 0;
			for (k = 0; k <= strlen(atom[i].ambername); k++)
				if ((atom[i].ambername[k] <= 'Z'
					 && atom[i].ambername[k] >= 'A')
					|| (atom[i].ambername[k] <= 'z'
						&& atom[i].ambername[k] >= 'a'))
					tmpchar[num++] = atom[i].ambername[k];
			tmpchar[num] = '\0';
			for (j = 0; j < parmnum; j++)
				if (strcmp(tmpchar, radii[j].name) == 0) {
					atom[i].radius = radii[j].radius;
					break;
				}
		}
	}
}



void adjustatomname(int atomnum, ATOM atom[], int check_atom_type)
{								/* just strips blanks out of the atomname */
	int i, j;
	int num, num1, num2;
	int tmpint;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char filename[MAXCHAR];
	char line[MAXCHAR];
	char *system_env;
	FILE *fp;
	int index;

	for (i = 0; i < atomnum; i++) {
		num = 0;
		tmpchar[0] = '\0';
		for (j = 0; j < strlen(atom[i].name); j++) {
			if (atom[i].name[j] == ' ')
				continue;
			tmpchar[num] = atom[i].name[j];
			if (tmpchar[num] == '\0')
				break;
			num++;
		}
		tmpchar[num] = '\0';
		strcpy(atom[i].name, tmpchar);
	}
	if(check_atom_type == 1) {
	        system_env = (char *) getenv("DOCK_HOME");
        	if (system_env != NULL) {
                	strcpy(filename, "");
                	strcpy(filename, system_env);
                	strcat(filename, "/parameters/antechamber/CORR_NAME_TYPE.DAT");
        	} else
                	strcpy(filename, "CORR_NAME_TYPE.DAT");
		if ((fp = fopen(filename, "r")) == NULL) {
                	fprintf(stderr, "Cannot open file %s, exit\n", filename);
                	exit(0);
        	}
		for (i = 0; i < atomnum; i++) {
			rewind(fp);	
			for(;;) {
			  	if (fgets(line, MAXCHAR, fp) == NULL)  break;
				sscanf(line, "%s%s", tmpchar1, tmpchar2);
				tmpint = strlen(tmpchar2);
				if(strcmp(tmpchar1, atom[i].ambername) == 0) {
					if(strncmp(atom[i].name, tmpchar2, tmpint)!=0)
						for(j=0;j<tmpint;j++)	
							atom[i].name[j] = tmpchar2[j];
					break;
				}
			}
		}
		fclose(fp);
	}
	for (i = 0; i < atomnum; i++) {
		num1 = 0;
		num2 = 0;
		index = 0;
		strcpy(tmpchar1, "");
		strcpy(tmpchar2, "");

		for (j = 0; j < strlen(atom[i].name); j++) {
			if (index == 0 && atom[i].name[0] >= '0'
				&& atom[i].name[0] <= '9')
				index = 1;
			if (index == 1
				&& (atom[i].name[j] >= '0' && atom[i].name[j] <= '9'))
				tmpchar1[num1++] = atom[i].name[j];
			if (index == 1
				&& (atom[i].name[j] < '0' || atom[i].name[j] > '9'))
				index = 2;
			if (index == 0 || index == 2)
				tmpchar2[num2++] = atom[i].name[j];
		}
		tmpchar2[num2] = '\0';
		if (index == 2) {
			tmpchar1[num1] = '\0';
			strcat(tmpchar2, tmpchar1);
		}
		strcpy(atom[i].name, tmpchar2);
	}
}


void atomname(int atomnum, ATOM * atom)
{
	int i, j;
	int countatom[150];
	char tmpchar[MAXCHAR];

	for (i = 0; i < 150; i++)
		countatom[i] = 0;
	for (j = 0; j < atomnum; j++) {
		countatom[atom[j].atomicnum]++;
		newitoa(countatom[atom[j].atomicnum], tmpchar);
		strcat(atom[j].name, tmpchar);
	}
}

void atomicnum(int atomnum, ATOM * atom)
{
	int i;
	for (i = 0; i < atomnum; i++) {
		switch (atom[i].name[0]) {
		case 'H':
			if (atom[i].name[1] == 'g') 
				atom[i].atomicnum = 80;
			else
				atom[i].atomicnum = 1;
			break;
		case 'C':
			if (atom[i].name[1] == 'l' || atom[i].name[1] == 'L') 
				atom[i].atomicnum = 17;
			else if (atom[i].name[1] == 'a') 
				atom[i].atomicnum = 20;
			else if (atom[i].name[1] == 'o') 
				atom[i].atomicnum = 27;
			else if (atom[i].name[1] == 'u')  
				atom[i].atomicnum = 29;
			else if (atom[i].name[1] == 'r')  
				atom[i].atomicnum = 24;
			else if (atom[i].name[1] == 'd')  
				atom[i].atomicnum = 48;
			else 
				atom[i].atomicnum = 6;
			break;
		case 'M':
			if (atom[i].name[1] == 'n')
				atom[i].atomicnum = 25;
			else if (atom[i].name[1] == 'g')
				atom[i].atomicnum = 12;
			break;
		case 'N':
			if (atom[i].name[1] == 'i') 
				atom[i].atomicnum = 28;
			else if (atom[i].name[1] == 'a') 
				atom[i].atomicnum = 11;
			else 
				atom[i].atomicnum = 7;
			break;
		case 'O':
			atom[i].atomicnum = 8;
			break;
		case 'F':
			if (atom[i].name[1] == 'e')
				atom[i].atomicnum = 26;
			else
				atom[i].atomicnum = 9;
			break;
		case 'B':
			if (atom[i].name[1] == 'r' || atom[i].name[1] == 'R') 
				atom[i].atomicnum = 35;
			else if (atom[i].name[1] == 'e') 
				atom[i].atomicnum = 4;
			else if (atom[i].name[1] == 'a') 
				atom[i].atomicnum = 56;
			else 
				atom[i].atomicnum = 5;
			break;
		case 'S':
			if (atom[i].name[1] == 'i' || atom[i].name[1] == 'I') 
				atom[i].atomicnum = 14;
			else if (atom[i].name[1] == 'c') 
				atom[i].atomicnum = 21;
			else if (atom[i].name[1] == 'e') 
				atom[i].atomicnum = 34;
			else if (atom[i].name[1] == 'r') 
				atom[i].atomicnum = 38;
			else 
				atom[i].atomicnum = 16;
			break;
		case 'P':
			if (atom[i].name[1] == 'd') 
				atom[i].atomicnum = 46;
			else if (atom[i].name[1] == 't') 
				atom[i].atomicnum = 78;
			else if (atom[i].name[1] == 'b') 
				atom[i].atomicnum = 82;
			else
				atom[i].atomicnum = 15;
			break;
		case 'Z':
			if (atom[i].name[1] == 'n') 
				atom[i].atomicnum = 30;
			break;
		case 'K':
			atom[i].atomicnum = 19;
			break;
		case 'I':
			atom[i].atomicnum = 53;
			break;
		case 'V':
			atom[i].atomicnum = 23;
			break;
		case 'R':
			if (atom[i].name[1] == 'u')
				atom[i].atomicnum = 44;
			if (atom[i].name[1] == 'h')
				atom[i].atomicnum = 45;
			break;
		case 'T':
			if (atom[i].name[1] == 'i')
				atom[i].atomicnum = 22;
			else if (atom[i].name[1] == 'l')
				atom[i].atomicnum = 81;
			break;
		case 'G':
			if (atom[i].name[1] == 'a')
				atom[i].atomicnum = 31;
			else if (atom[i].name[1] == 'e')
				atom[i].atomicnum = 32;
			break;
		case 'A':
			if (atom[i].name[1] == 's')
				atom[i].atomicnum = 33;
			else if (atom[i].name[1] == 'l')
				atom[i].atomicnum = 13;
			else if (atom[i].name[1] == 'g')
				atom[i].atomicnum = 47;
			else if (atom[i].name[1] == 'u')
				atom[i].atomicnum = 79;
			break;
		case 'l':
			if (atom[i].name[1] == 'p')
				atom[i].atomicnum = 0;
			break;
		case 'L':
			if (atom[i].name[1] == 'i') 
				atom[i].atomicnum = 3;
			else if (atom[i].name[1] == 'P') 
				atom[i].atomicnum = 0;
			break;
		case 'E':
			if (atom[i].name[1] == 'P')
				atom[i].atomicnum = 0;
			break;
		default:
			printf("\n Unrecognized atomic name %5s, exit", atom[i].name);
		}
	}
}

void element(int atomnum, ATOM * atom)
{
	int i;
	for (i = 0; i < atomnum; i++)
		switch (atom[i].atomicnum) {
		case 1:
			strcpy(atom[i].element, "H");
			break;
		case 3:
			strcpy(atom[i].element, "Li");
			break;
		case 4:
			strcpy(atom[i].element, "Be");
			break;
		case 5:
			strcpy(atom[i].element, "B");
			break;
		case 6:
			strcpy(atom[i].element, "C");
			break;
		case 7:
			strcpy(atom[i].element, "N");
			break;
		case 8:
			strcpy(atom[i].element, "O");
			break;
		case 9:
			strcpy(atom[i].element, "F");
			break;
		case 11:
			strcpy(atom[i].element, "Na");
			break;
		case 12:
			strcpy(atom[i].element, "Mg");
			break;
		case 13:
			strcpy(atom[i].element, "Al");
			break;
		case 14:
			strcpy(atom[i].element, "Si");
			break;
		case 15:
			strcpy(atom[i].element, "P");
			break;
		case 16:
			strcpy(atom[i].element, "S");
			break;
		case 17:
			strcpy(atom[i].element, "Cl");
			break;
		case 19:
			strcpy(atom[i].element, "K");
			break;
		case 20:
			strcpy(atom[i].element, "Ca");
			break;
		case 21:
			strcpy(atom[i].element, "Sc");
			break;
		case 22:
			strcpy(atom[i].element, "Ti");
			break;
		case 23:
			strcpy(atom[i].element, "V");
			break;
		case 24:
			strcpy(atom[i].element, "Cr");
			break;
		case 25:
			strcpy(atom[i].element, "Mn");
			break;
		case 26:
			strcpy(atom[i].element, "Fe");
			break;
		case 27:
			strcpy(atom[i].element, "Co");
			break;
		case 28:
			strcpy(atom[i].element, "Ni");
			break;
		case 29:
			strcpy(atom[i].element, "Cu");
			break;
		case 30:
			strcpy(atom[i].element, "Zn");
			break;
		case 31:
			strcpy(atom[i].element, "Ga");
			break;
		case 32:
			strcpy(atom[i].element, "Ge");
			break;
		case 33:
			strcpy(atom[i].element, "As");
			break;
		case 34:
			strcpy(atom[i].element, "Se");
			break;
		case 35:
			strcpy(atom[i].element, "Br");
			break;
		case 45:
			strcpy(atom[i].element, "Rh");
			break;
		case 46:
			strcpy(atom[i].element, "Pd");
			break;
		case 47:
			strcpy(atom[i].element, "Ag");
			break;
		case 48:
			strcpy(atom[i].element, "Cd");
			break;
		case 53:
			strcpy(atom[i].element, "I");
			break;
		case 56:
			strcpy(atom[i].element, "Ba");
			break;
		case 78:
			strcpy(atom[i].element, "Pt");
			break;
		case 79:
			strcpy(atom[i].element, "Au");
			break;
		case 80:
			strcpy(atom[i].element, "Hg");
			break;
		case 81:
			strcpy(atom[i].element, "Tl");
			break;
		case 82:
			strcpy(atom[i].element, "Pb");
			break;
		default:
			strcpy(atom[i].element, "du");
			break;
		}
}


void duplicatedname(int atomnum, ATOM * atom)
{
	int i, j, k;
	int index;
	int id;
	char tmpchar[MAXCHAR];
	char name[MAXCHAR];
	char elemname[150][5];
	for (i = 0; i < 150; i++)
		strcpy(elemname[i], "X");
        strcpy(elemname[1], "H");
        strcpy(elemname[3], "Li");
        strcpy(elemname[4], "Be");
        strcpy(elemname[5], "B");
        strcpy(elemname[6], "C");
        strcpy(elemname[7], "N");
        strcpy(elemname[8], "O");
        strcpy(elemname[9], "F");
        strcpy(elemname[11], "Na");
        strcpy(elemname[12], "Mg");
        strcpy(elemname[13], "Al");
        strcpy(elemname[14], "Si");
        strcpy(elemname[15], "P");
        strcpy(elemname[16], "S");
        strcpy(elemname[17], "Cl");
        strcpy(elemname[19], "K");
        strcpy(elemname[20], "Ca");
        strcpy(elemname[21], "Sc");
        strcpy(elemname[22], "Ti");
        strcpy(elemname[23], "V");
        strcpy(elemname[24], "Cr");
        strcpy(elemname[25], "Mn");
        strcpy(elemname[26], "Fe");
        strcpy(elemname[27], "Co");
        strcpy(elemname[28], "Ni");
        strcpy(elemname[29], "Cu");
        strcpy(elemname[30], "Zn");
        strcpy(elemname[31], "Ga");
        strcpy(elemname[32], "Ge");
        strcpy(elemname[33], "As");
        strcpy(elemname[34], "Se");
        strcpy(elemname[35], "Br");
        strcpy(elemname[45], "Rh");
        strcpy(elemname[46], "Pd");
        strcpy(elemname[47], "Ag");
        strcpy(elemname[48], "Cd");
        strcpy(elemname[53], "I");
        strcpy(elemname[56], "Ba");
        strcpy(elemname[78], "Pt");
        strcpy(elemname[79], "Au");
        strcpy(elemname[80], "Hg");
        strcpy(elemname[80], "Tl");
        strcpy(elemname[82], "Pb");

	for (i = 0; i < atomnum; i++)
		for (j = i + 1; j < atomnum; j++) 
			if (strcmp(atom[i].name, atom[j].name) == 0) {
				id = 1;
				index = 1;
				while (index) {
					newitoa(id, tmpchar);
					strcpy(name, elemname[atom[j].atomicnum]);
					strcat(name, tmpchar);
					for (k = 0; k < atomnum; k++)
						if (strncmp(name, atom[k].name,strlen(name)) == 0) {
							index = 0;
							break;
						}
					if (index == 1) 
						break;	
					if (index == 0) {
						id++;
					 	index = 1;	
					}
				}
				strcpy(atom[j].name, name);
			}
}

void info(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	int i;
	printf("\natomnum: %d", atomnum);
	printf("\nbondnum: %d", bondnum);

	printf("\n-----------------------------------------\n");
	for (i = 0; i < atomnum; i++)
		printf("ATOM%7d  %-4s%-4s%5d%5s%5s%12.3f%8.3f%8.3f%8d\n",
			   i + 1, atom[i].name, atom[i].aa, atom[i].resno,
			   atom[i].ambername, atom[i].element, atom[i].x, atom[i].y,
			   atom[i].z, atom[i].atomicnum);

	printf("\n-----------------------------------------\n");

	for (i = 0; i < bondnum; i++)
		printf("BOND  %5d  %5d  %5d  %5s  %5s  %5d\n", i + 1,
			   bond[i].bondi, bond[i].bondj, atom[bond[i].bondi].name,
			   atom[bond[i].bondj].name, bond[i].type);

	printf("\n-----------------------------------------\n");

	for (i = 0; i < atomnum; i++) {
		printf("\nFor Atom %d %s\n\n", i + 1, atom[i].name);
		printf("RING nr: %5d\n", arom[i].nr);
		printf("RING ar1=%d, ar2=%d, ar3=%d, ar4=%d, ar5=%d \n",
			   arom[i].ar1, arom[i].ar2, arom[i].ar3, arom[i].ar4,
			   arom[i].ar5);
		printf("RING rg3=%d, rg4=%d, rg5=%d, rg6=%d, rg7=%d \n",
			   arom[i].rg[3], arom[i].rg[4], arom[i].rg[5], arom[i].rg[6],
			   arom[i].rg[7]);
	}
	printf("\n-----------------------------------------\n");

	printf("\ngkeyword: %s", minfo.gkeyword);
	printf("\nmkeyword: %s", minfo.mkeyword);
	printf("\nresname: %s", minfo.resname);
	printf("\natom_type_def: %s", minfo.atom_type_def);
	printf("\nresfilename: %s", minfo.resfilename);
	printf("\ngfilename: %s", minfo.gfilename);
	printf("\nconnect_file: %s", minfo.connect_file);
	printf("\nradius_file: %s", minfo.radius_file);
	printf("\ninf_filename: %s", minfo.inf_filename);
	printf("\nlongresname: %s", minfo.longresname);
	printf("\nmultiplicity: %d", minfo.multiplicity);
	printf("\nicharge: %d", minfo.icharge);
	printf("\nusercharge: %d", minfo.usercharge);
	printf("\ndcharge: %lf", minfo.dcharge);

	printf("\n-----------------------------------------\n");
	printf("\nintype: %s", cinfo.intype);
	printf("\nouttype: %s", cinfo.outtype);
	printf("\natype: %s", cinfo.atype);
	printf("\nchargetype: %s", cinfo.chargetype);
	printf("\nmaxatom: %d", cinfo.maxatom);
	printf("\nmaxbond: %d", cinfo.maxbond);
	printf("\nmaxring: %d", cinfo.maxring);
	printf("\nrnindex: %d", cinfo.rnindex);
	printf("\nintstatus: %d", cinfo.intstatus);
	printf("\npfindex: %d", cinfo.pfindex);
	printf("\nprediction_index %d", cinfo.prediction_index);
	printf("\nbpindex: %d", cinfo.bpindex);
	printf("\n-----------------------------------------\n\n");


}
